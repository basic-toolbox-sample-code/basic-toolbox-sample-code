/* Copyright (c) 2015-2017, Intel Corporation
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
#include <stdlib.h>
#include <sys/mman.h>

#include "barrier/barrier.h"

#ifdef USE_POOL_MEMORY
#define TBB_PREVIEW_MEMORY_POOL 1
#include "tbb/memory_pool.h"
#endif

#ifdef TBB
#include "tbb/tbb.h"
#endif

#define MAX_THREADS (1024*1024)

#ifdef NUMA
#include <pthread.h>
#endif


#include "types.h"

void pinToCore(size_t core)
{
#ifdef NUMA
   // pin the current thread to the specified core to make NUMA memory access mapping deterministic
   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);

#if 1 
   int socket = core % 4;
   int ht = core / (4*18);
   core = (core / 4) % 18;
   core = ht*72 + socket * 18 + core;
#endif
   
   CPU_SET(core, &cpuset);

   int res = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
   if(res) {
      cerr << "could not set thread affinity to "<< core << endl;
      throw exception();
   }
#endif
}

char pad0[64];
#ifdef USE_POOL_MEMORY
    #ifdef NUMA
    std::vector<tbb::fixed_pool *> memPool(MAX_THREADS,NULL);
    #else
    tbb::fixed_pool * memPool = NULL;
    #endif
typedef tbb::memory_pool_allocator<Element> PoolAllocT;
// buckets array indexed as [iPE][bucket]
typedef vector<Element,PoolAllocT> BucketT;
typedef vector<BucketT> PEBucketT;
#else
#define  PoolAllocT(x) 
// buckets array indexed as [iPE][bucket]
typedef vector<Element> BucketT;
typedef vector<BucketT> PEBucketT;
#endif


typedef BinomialBarrier BarrierType;

template <class Iterator>
void pSampleSort(const Iterator & s, const size_t n, const unsigned p)
{
    mt19937 rndEngine;
    BarrierType barrier(p); // for barrier synchronization
    /////// Choose random samples
    vector<KeyType> S; // random sample of a elements from s
    uniform_int_distribution<size_t> SampleDistribution(0, n-1);
    const int a = (int)(16*log(p)/log(2.)); // oversampling ratio
    for(size_t i=0; i < (size_t)(a+1)*p - 1; ++i)
        S.push_back((s + SampleDistribution(rndEngine))->key);
    sort(S.begin(),S.end()); // sort samples
    for(size_t i=0; i < p-1 ; ++i) // select splitters
        S[i] = S[(a+1)*(i+1)];
    S.resize(p-1);
    // buckets array indexed as [iPE][bucket]
    vector<PEBucketT> buckets;
    for(size_t i=0; i < p ; ++i) {
#ifdef NUMA
        buckets.push_back(PEBucketT(p,BucketT(PoolAllocT(*memPool[i]))) );
#else
        buckets.push_back(PEBucketT(p,BucketT(PoolAllocT(*memPool))) );
#endif
    }
    vector<size_t> bucketSize(p, 0);
    vector<thread> threads(p);
    for (size_t i = 0; i < p; ++i) { /// go parallel
        threads[i] = thread( [&] (const unsigned iPE) {
                pinToCore(iPE); // no-op if not NUMA
                auto current = s + iPE*n/p, end = s + (iPE+1)*n/p;
                auto & myBuckets = buckets[iPE]; // buckets array indexed as [iPE][bucket]
                while(current != end) {
                   // do binary search
                   const size_t i = upper_bound(S.begin(), S.end(), current->key)
                         - S.begin();
                   myBuckets[i].push_back(*current++);
                }
                barrier.wait(iPE, p);
                // now each thread works on "iPE" bucket. First compute the total size of bucket iPE:
                size_t bSize = 0; // use local variable to avoid false-sharing
                for(const auto & PEbuckets: buckets) bSize += PEbuckets[iPE].size();
                bucketSize[iPE] = bSize;
                barrier.sync(iPE);
                // find the bucket begin position in the s by summing up the
                // sizes of the previous buckets (<iPE)
                auto bucketBegin = s;
                for(size_t b = 0; b < iPE ; ++b) bucketBegin += bucketSize[b];
                // copy the bucket 'iPE' from all PEs into s
                auto currentOut = bucketBegin;
                for(const auto & PEbuckets: buckets)
                    currentOut = copy(PEbuckets[iPE].cbegin(), PEbuckets[iPE].cend(), currentOut);
                // sort the bucket
                sort(bucketBegin, currentOut);
            }, i);
    }
    for (auto & t : threads) t.join();
}

#ifdef NUMA

//#define USE_HUGE_PAGES

// used in NUMA version only
// each thread initializes its own area in the memory and causes lazy physical memory allocation
// on the socket it runs (according to the Linux first-touch NUMA allocation policy)
void initialize(const size_t iPE, Element * begin, Element * end, const size_t capacity)
{
    pinToCore(iPE); // no-op if not NUMA
    mt19937 rndEngine((unsigned long long)begin);
    uniform_int_distribution<KeyType> KeyDistribution(0, (numeric_limits<KeyType>::max)() - 1);
    generate(begin, end, bind(KeyDistribution, rndEngine));
#ifdef USE_POOL_MEMORY
    if(iPE == 0) cout << "Allocate Pool Memory (NUMA)" << std::endl;
#ifdef USE_HUGE_PAGES
    char * buffer = (char*) mmap(NULL, capacity, PROT_READ | PROT_WRITE,
                       MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, 0, 0);
    if(buffer == MAP_FAILED) { cout << "mmap failed for capacity " << capacity << std::endl; exit(-1); }
#else
    char * buffer = new char[capacity];
#endif
    std::fill(buffer, buffer + capacity, 0); // allocate physical memory by touching the array
    memPool[iPE] = new tbb::fixed_pool(buffer, capacity); 
#endif
}

#endif

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

    try {
        cout <<"Fill the input with random data" << endl;
        uniform_int_distribution<KeyType> KeyDistribution(0, (numeric_limits<KeyType>::max)() - 1);

        const size_t poolCapacity = 3ULL* std::max(n,(size_t)(1024ULL*1024ULL)) *sizeof(Element);
#ifndef NUMA
        // naive (non-NUMA) optimized inizialization

        vector<Element> s(n);
        auto sBegin = s.begin(), sEnd = s.end();
        std::generate(sBegin, sEnd, bind(KeyDistribution, rndEngine));

        #ifdef USE_POOL_MEMORY
        cout << "Allocate Pool Memory." << endl;
        char * buffer = new char[poolCapacity];
        std::fill(buffer, buffer + poolCapacity, 0); // allocate physical memory by touching the array
        cout << "buffer: "<< (unsigned long long) buffer << " capacity: "<<  poolCapacity << std::endl; 
        memPool = new tbb::fixed_pool(buffer, poolCapacity);
        #endif
#else
        // NUMA optimized initialization (Linux only)

        // reserve large memory space from Linux without touching/initializing it
        shared_ptr<Element> s(new Element[n]);
        Element * sBegin = s.get();
        Element * sEnd = sBegin + n;

        vector<thread> threads(p);
        for (size_t i = 0; i < p; ++i) {
            threads[i] = thread( initialize, i, sBegin + i*n/p, sBegin + (i+1)*n/p, poolCapacity/p);
        }
        for (auto & t : threads) {
            t.join();
        }
#endif
        Element checksum = accumulate(sBegin, sEnd, Element(0));

#ifdef PROFILE
        system("/usr/bin/sh start_profiling.sh");
#endif
        cout << "Sorting started." << endl;
        chrono::steady_clock::time_point start = chrono::steady_clock::now();
        pSampleSort(sBegin, n, p);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        cout << "Sorting completed." << endl;
#ifdef PROFILE
        system("/usr/bin/sh stop_profiling.sh");
#endif

        // validate the result
        bool correct;
        cout << "The data is "<< ((correct = is_sorted(sBegin, sEnd))?"sorted":"NOT sorted" ) << endl;
        Element newChecksum = accumulate(sBegin, sEnd, Element(0));
        correct &= (newChecksum.key == checksum.key);
        cout << "The checksum is "<< ((newChecksum.key == checksum.key)?"correct":"WRONG" ) << endl;

        const long long timeMs = chrono::duration_cast<chrono::milliseconds>(end - start).count();
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

#ifdef USE_POOL_MEMORY
    #ifdef NUMA
      for (size_t i = 0; i < MAX_THREADS; ++i) {
         if(memPool[i]) delete memPool[i];
      }
    #else
      delete memPool;
    #endif
#endif

    return 0;
}


