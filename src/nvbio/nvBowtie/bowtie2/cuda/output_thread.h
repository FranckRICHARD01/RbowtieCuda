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

#pragma once

#include <nvBowtie/bowtie2/cuda/defs.h>
#include <nvBowtie/bowtie2/cuda/stats.h>
#include <nvbio/basic/threads.h>
#include <nvbio/basic/timer.h>
#include <nvbio/io/sequence/sequence.h>
#include <nvbio/basic/numbers.h>
#include <stack>
#include <deque>

namespace nvbio {
namespace bowtie2 {
namespace cuda {

//
// A class implementing a background input thread, providing
// a set of input read-streams which are read in parallel to the
// operations performed by the main thread.
//

struct OutputThreadSE : public Thread<OutputThreadSE> 
{
 static const uint32 BUFFERS = 2;

    OutputThreadSE( FILE * file_pointer) :
    write_thread_data(nullptr),
    stop_thread(false),
    fp( file_pointer )
    { } 
    
    void run();

    // release a batchs
    //
    void release();

    uint32 mystrcat( char* dest, const char* src, bool init = false )
    {
        static uint32 p = 0;

        if (init) { p = 0; dest[0] = 0; return 0; }

        //if (src[0] == 0) return p;

        dest += p;

        while (*src) { *dest++ = *src++; p++; }

        --dest;

        return p;
    }

    void createCache()
    {       
        if (!write_thread_data)
        {
            log_debug(stderr, "write thread: create cache memory\n");
            write_thread_data = new char[MAX_CHAR];

            #pragma parallel for
            for(uint32 i = 0; i < MAX_CHAR; i++) write_thread_data[i]=0;
        }  
    }

    
    char * write_thread_data;

    uint32 write_thread_data_size;

    bool stop_thread;


private:

    Mutex                                m_free_pool_lock;
    
    FILE * fp;   
};

//
// A class implementing a background input thread, providing
// a set of input read-streams which are read in parallel to the
// operations performed by the main thread.
//

struct OutputThreadPE : public Thread<OutputThreadPE>
{
    static const uint32 BUFFERS = 2;

    OutputThreadPE( FILE * file_pointer) :
    write_thread_data(nullptr),
    stop_thread(false),
    fp( file_pointer )
    {} 
    
    void run();

    // release a batchs
    //
    void release();

    uint32 mystrcat( char* dest, const char* src, bool init = false )
    {
        static uint32 p = 0;

        if (init) { p = 0; dest[0] = 0; return 0; }

        if (src[0] == 0) return p;

        dest += p;

        while (*src) { *dest++ = *src++; p++; }

        --dest;

        return p;
    }

    void createCache()
    {       
        if (!write_thread_data)
        {
            log_debug(stderr, "write thread: create cache memory\n");
            write_thread_data = new char[MAX_CHAR];

            #pragma parallel for
            for(uint32 i = 0; i < MAX_CHAR; i++) write_thread_data[i]=0;
        }  
    }

    char * write_thread_data;

    uint32 write_thread_data_size;

    bool stop_thread;

private:

    Mutex                                m_free_pool_lock;
    
    FILE * fp;   
};

} // namespace io
} // namespace nvbio
}



