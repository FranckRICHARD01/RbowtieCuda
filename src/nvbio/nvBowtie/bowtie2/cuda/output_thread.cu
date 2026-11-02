/*
 * nvbio
 * Copyright (c) 2011-2014, NVIDIA CORPORATION. All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of the NVIDIA CORPORATION nor the
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#include <nvBowtie/bowtie2/cuda/defs.h>
#include <nvBowtie/bowtie2/cuda/params.h>
#include <nvBowtie/bowtie2/cuda/stats.h>
#include <nvBowtie/bowtie2/cuda/output_thread.h>
#include <nvbio/io/output/output_utils.h>
#include <nvbio/basic/threads.h>
#include <nvbio/basic/atomics.h>
#include <nvbio/basic/timer.h>
#include <nvbio/basic/exceptions.h>
#include <chrono>
#include <thread>

namespace nvbio {
namespace bowtie2 {
namespace cuda {


//
// Enqueue a fully read-back host batch to the writer-thread queue.
// The batch is deep-copied into a ring slot: the writer thread formats and
// writes it asynchronously, while the compute thread recycles its scratch
// buffers and the input thread recycles the read storage for this batch.
//
void OutputThreadSE::submit( io::HostOutputBatchSE& cpu_batch)
{
    std::unique_lock<std::mutex> lock( m_queue_lock );

    m_cv_free.wait( lock, [this](){ return !m_free_idx.empty(); });

    const uint32 idx = m_free_idx.front();
    m_free_idx.pop_front();

    QueueSlot& slot = m_slots[idx];

    slot.batch  = cpu_batch;             // deep copy of the host output-batch containers
    slot.reads  = *cpu_batch.read_data;  // deep copy of the host read storage for this batch
    slot.batch.read_data = &slot.reads;  // bind the batch views to the slot-owned storage

    m_ready_idx.push_back( idx );

    lock.unlock();
    m_cv_ready.notify_one();
}

//
// Mark the queue closed: wakes the writer thread, which drains any remaining
// published batches before exiting; a subsequent join() then returns.
//
void OutputThreadSE::close_and_drain()
{
    {
        std::lock_guard<std::mutex> lock( m_queue_lock );
        m_closed = true;
    }
    m_cv_ready.notify_all();
    m_cv_free.notify_all();
}

void OutputThreadSE::run()
{
    log_verbose( stderr, "starting background output thread\n" );

    try
    {
        while (true)
        {
            uint32 idx = 0;
            {
                std::unique_lock<std::mutex> lock( m_queue_lock );

                m_cv_ready.wait( lock, [this](){ return !m_ready_idx.empty() || m_closed; });

                if ( m_ready_idx.empty() && m_closed )
                    break;

                idx = m_ready_idx.front();
                m_ready_idx.pop_front();
            }

            QueueSlot& slot = m_slots[idx];

            //
            // The writer thread is the only one touching the output file: it
            // formats the SAM records and writes them in batch-submission
            // order, preserving the exact ordering of the synchronous path.
            //
            if ( output_file )
                output_file->process( slot.batch );

            {
                std::lock_guard<std::mutex> lock( m_queue_lock );
                m_free_idx.push_back( idx );
            }
            m_cv_free.notify_one();
        }

        log_verbose( stderr, "stop background output thread\n" );
    }
    catch (nvbio::bad_alloc &e)
    {
        log_error(stderr, "caught a nvbio::bad_alloc exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (nvbio::logic_error &e)
    {
        log_error(stderr, "caught a nvbio::logic_error exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (nvbio::runtime_error &e)
    {
        log_error(stderr, "caught a nvbio::runtime_error exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (std::bad_alloc &e)
    {
        log_error(stderr, "caught a std::bad_alloc exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (std::logic_error &e)
    {
        log_error(stderr, "caught a std::logic_error exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (std::runtime_error &e)
    {
        log_error(stderr, "caught a std::runtime_error exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (...)
    {
        log_error(stderr, "caught an unknown exception!\n");
        exit(1);
    }
}


void OutputThreadPE::run()
{
    log_verbose( stderr, "starting background paired-end output thread\n" );
    
    while ( !write_thread_data )
    {
        std::this_thread::sleep_for(std::chrono::microseconds(1));
    }

    try
    {
        while (!stop_thread)
        {
            if ( write_thread_data && write_thread_data[0] != 0 )
            {
                log_debug(stderr, "      writing output file.\n");

                fwrite(write_thread_data, 1, write_thread_data_size, fp);
            
                write_thread_data[0]=0;
            }       
        }

        if (write_thread_data)
        {
            delete [] write_thread_data;
            write_thread_data = nullptr;
        }     

        log_verbose( stderr, "stop background paired-end output thread\n" );
    }
    catch (nvbio::bad_alloc &e)
    {
        log_error(stderr, "caught a nvbio::bad_alloc exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (nvbio::logic_error &e)
    {
        log_error(stderr, "caught a nvbio::logic_error exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (nvbio::runtime_error &e)
    {
        log_error(stderr, "caught a nvbio::runtime_error exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (std::bad_alloc &e)
    {
        log_error(stderr, "caught a std::bad_alloc exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (std::logic_error &e)
    {
        log_error(stderr, "caught a std::logic_error exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (std::runtime_error &e)
    {
        log_error(stderr, "caught a std::runtime_error exception:\n");
        log_error(stderr, "  %s\n", e.what());
        exit(1);
    }
    catch (...)
    {
        log_error(stderr, "caught an unknown exception!\n");
        exit(1);
    }
}

// release a batch
//
void OutputThreadPE::release()
{
    // push back to the free pool
    ScopedLock lock( &m_free_pool_lock );
  //  m_free_pool1.push( write_data.first );
}

} // namespace io
} // namespace nvbio
}
