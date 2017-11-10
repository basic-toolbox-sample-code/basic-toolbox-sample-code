/* Copyright (c) 2014-2017, Intel Corporation
   All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "stdafx.h"
#define _VARIADIC_MAX 10
#include <iostream>
#include <thread>
#include <random>
#include <functional>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <mutex>
#include <memory>
#include <atomic>
#include <numeric>

#ifdef TBB
#include "tbb/tbb.h"
#endif

#ifdef NUMA
#include <pthread.h>
#endif

#ifdef STD_PARALLEL_MODE
#include <parallel/algorithm>
#endif

#include "types.h"

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cerr << "Parameters: num_threads [input_size]" << endl;
        return -1;
    }
    const unsigned p = atoi(argv[1]);
    size_t n = p*8ULL*1024ULL*1024ULL;
    if(argc >= 3) {
        n = ((size_t)atoi(argv[2]))*((size_t)p);
    }
    mt19937 rndEngine;

    cout << "Number of threads: "<< p << endl;
    cout << "Number of elements: "<< n << endl;

    #ifdef  STD_PARALLEL_MODE
    omp_set_num_threads(p);
    #endif 

    #ifdef TBB
    tbb::task_scheduler_init init(p);
    #endif

    try {
        cout <<"Fill the input with random data" << endl;
        uniform_int_distribution<KeyType> KeyDistribution(0, (numeric_limits<KeyType>::max)() - 1);

#ifndef NUMA
        // naive (non-NUMA) optimized inizialization

        vector<Element> s(n);
        auto sBegin = s.begin(), sEnd = s.end();
	#ifdef  STD_PARALLEL_MODE
        __gnu_parallel::generate(sBegin, sEnd, [&]{ return 0ULL; } ); // only touch the array, can't use the thread-unsafe rnd number generator
        std::generate(sBegin, sEnd, bind(KeyDistribution, rndEngine)); // single-threaded
	#else 
           #ifdef TBB
           tbb::parallel_for(size_t(0), n, size_t(1) , [=](size_t i) { *(sBegin + i) = 0ULL;}); // only touch the array, can't use the thread-unsafe rnd number generator
           #endif
           std::generate(sBegin, sEnd, bind(KeyDistribution, rndEngine)); // single-threaded
	#endif
#else
        // NUMA optimized initialization (Linux only)

        // reserve large memory space from Linux without touching/initializing it
        shared_ptr<Element> s(new Element[n]);
        Element * sBegin = s.get();
        Element * sEnd = sBegin + n;

        // each thread initializes its own area in the memory causing lazy physical memory allocation
        // on the socket it runs (according to the Linux first-touch NUMA allocation policy)
        #ifdef  STD_PARALLEL_MODE
        __gnu_parallel::generate(sBegin, sEnd, [&]{ return 0ULL; } ); // only touch the array, can't use the thread-unsafe rnd number generator
        std::generate(sBegin, sEnd, bind(KeyDistribution, rndEngine)); // single-threaded
        #else
            #ifdef TBB
            tbb::parallel_for(size_t(0), n, size_t(1) , [=](size_t i) { *(sBegin + i) = 0ULL;}); // only touch the array, can't use the thread-unsafe rnd number generator
            #endif
            std::generate(sBegin, sEnd, bind(KeyDistribution, rndEngine)); // single-threaded
        #endif
#endif
        Element checksum = accumulate(sBegin, sEnd, Element(0));

        cout << "Sorting started." << endl;
        chrono::steady_clock::time_point start = chrono::steady_clock::now();
#ifdef  STD_SORT
        sort(sBegin,sEnd);
#else
	#ifdef STD_PARALLEL_MODE
	__gnu_parallel::sort(sBegin,sEnd);
	#else
            #ifdef TBB
                tbb::parallel_sort(sBegin,sEnd);
            #endif
	#endif
#endif
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        cout << "Sorting completed." << endl;

        // validate the result
        bool correct;
        cout << "The data is "<< ((correct = is_sorted(sBegin, sEnd))?"sorted":"NOT sorted" ) << endl;
        Element newChecksum = accumulate(sBegin, sEnd, Element(0));
        correct &= (newChecksum.key == checksum.key);
        cout << "The checksum is "<< ((newChecksum.key == checksum.key)?"correct":"WRONG" ) << endl;

        const long long timeMs = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        if(timeMs == 0) {
            cout << "zero Ms" << endl;
            exit(-1);
        }
        const long long tput = 1000ULL*n/timeMs;
        cout << "Sorting took " << timeMs << " ms -> " << tput << " elements per second." << endl;

        if(!correct) {
            for(size_t i=0; i < n; ++i) {
                cout << i << " " << (sBegin + i)->key << endl;
            }
        }
        else {
           cout << ">>>" << n << ";" << p << ";" << timeMs << ";" <<tput << endl;
        }

    } catch (bad_alloc & )
    {
        cerr << "ERROR: An out of memory exception thrown" << endl;
        return -1;
    }

    return 0;
}


