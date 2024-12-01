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


void OutputThreadSE::run()
{
   log_verbose( stderr, "starting background output thread\n" );

   while ( !write_thread_data )
   {
        std::this_thread::sleep_for(std::chrono::microseconds(1));
   }

   try
   {
        while (!stop_thread)
        {
            if ( write_thread_data &&  write_thread_data[0] != 0 )
            {
                log_debug(stderr, "      writing output file.\n");

                fwrite(write_thread_data, 1, write_thread_data_size, fp);
            
                write_thread_data[0]=0;
            }   

            if (write_thread_data)
            {
                delete [] write_thread_data;
                write_thread_data = nullptr;
            }       
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

// release a batch
//
void OutputThreadSE::release()
{
    // push back to the free pool
    ScopedLock lock( &m_free_pool_lock );
    //m_free_pool.push( write_data );
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
