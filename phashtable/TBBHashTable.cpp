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

#ifndef __linux__
#include "stdafx.h"
#endif

#include <atomic>
#include <vector>
#include <functional>
#include <limits>
#include <assert.h>
#include <iostream>
#include <xmmintrin.h>
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

#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_hash_map.h"

#ifdef __linux__
#include <malloc.h>
#define _aligned_malloc(size, al) memalign(al, size) 
#define _aligned_free free
#endif

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

#include "backward/hashtable.h"

using namespace std;

struct Overwrite {};
struct Increment {};
struct Decrement {};


using namespace std;
int insertPerc;
size_t opsPerThreads, initialSize;
typedef long long int Key;
typedef long long int Data;
#ifdef TBB_HASH_MAP
typedef tbb::concurrent_hash_map<Key,Data> HTableType;
#else
typedef tbb::concurrent_unordered_map<Key,Data> HTableType;
#endif
HTableType * htable = NULL;
typedef std::pair<Key,Data> MyElement;
atomic<bool> go(false);


void worker( const size_t iPE // the ID of PE or thread
           )
{
       mt19937 rndEngine(iPE);
       uniform_int_distribution<Key> KeyDistribution(0, (numeric_limits<Key>::max)() - 1);
       auto valueRnd = bind(KeyDistribution, rndEngine);
       uniform_int_distribution<int> OpDistribution(0, 100);
       auto opRnd = bind(OpDistribution, rndEngine);

       while(!go);

       for(size_t i =0 ; i< opsPerThreads; ++i) {
           if(insertPerc <= opRnd()) {
#ifdef TBB_HASH_MAP
              HTableType::const_accessor result;
              htable->find(result, valueRnd());
#else
              htable->find(valueRnd());
#endif
           } else {
              htable->insert(MyElement(valueRnd(),0));
           }
       }
}

#ifdef CONTENTION

void workerContention( const size_t iPE // the ID of PE or thread
           )
{
       mt19937 rndEngine(iPE);
       uniform_int_distribution<Key> KeyDistribution(0, (numeric_limits<Key>::max)() - 1);
       auto valueRnd = bind(KeyDistribution, rndEngine);
       uniform_int_distribution<int> OpDistribution(0, 100);
       auto opRnd = bind(OpDistribution, rndEngine);

       while(!go);

       for(size_t i =0 ; i< opsPerThreads; ++i) {
           if(insertPerc <= opRnd()) {
#ifdef TBB_HASH_MAP
              HTableType::accessor result; // write lock
              if(!htable->insert(result, MyElement(valueRnd(), 0)))
              		++(result->second);
#else
              std::pair<HTableType::iterator, bool> result =
                   htable->insert(MyElement(valueRnd(), 0));
              if(!result.second)
                    __sync_fetch_and_add(&(result.first->second),1);
#endif

           } else {
#ifdef TBB_HASH_MAP
              HTableType::accessor result; // write lock
              if(!htable->insert(result, MyElement(0xabcdabcd, 0)))
                        ++(result->second);
#else
              std::pair<HTableType::iterator, bool> result =
                   htable->insert(MyElement(0xabcdabcd, 0));
              if(!result.second)
                    __sync_fetch_and_add(&(result.first->second),1);
#endif

           }
       }
}

#endif


#ifdef READ_CONTENTION

volatile Data dummy = 0;

void workerReadContention( const size_t iPE // the ID of PE or thread
           )
{
       mt19937 rndEngine(iPE);
       uniform_int_distribution<Key> KeyDistribution(0, initialSize - 1);
       auto valueRnd = bind(KeyDistribution, rndEngine);
       uniform_int_distribution<int> OpDistribution(0, 100);
       auto opRnd = bind(OpDistribution, rndEngine);

       while(!go);

       Data sum = 0;

       for(size_t i =0 ; i< opsPerThreads; ++i) {
           if(insertPerc <= opRnd()) {
#ifdef TBB_HASH_MAP
              HTableType::const_accessor result;
              htable->find(result, valueRnd());
              sum += result->second;
#else
              sum += htable->find(valueRnd())->second;
#endif
           } else {
#ifdef TBB_HASH_MAP
              HTableType::const_accessor result;
              htable->find(result, 42);
              sum += result->second;
#else
              sum += htable->find(42)->second;
#endif

           }
       }

       dummy += sum;
}

#endif


int main(int argc, char* argv[])
{
        if(argc < 5)
        {
             cout << "Parameters: threads min_capacity initial_size ops inserts% " << endl;
             return -1;
        }

        size_t p = atoi(argv[1]);
        size_t capacity = atoi(argv[2]);
        initialSize = atoi(argv[3]);
        size_t n = atoi(argv[4]);
        insertPerc =atoi(argv[5]);
#ifndef CONTENTION
#ifndef READ_CONTENTION
        if(insertPerc > 0) // n is the number of inserts
             n = double(n)*100.0/double(insertPerc);
        else
             n *= 10ULL;
#endif
#endif
        opsPerThreads = n/p;

        htable = new HTableType(capacity);
        cout << "Using capacity of "<< capacity << " elements. "<< capacity*sizeof(MyElement)/(1024ULL*1024ULL)<< " Mbytes " << endl;
        mt19937 rndEngine(0xabcdef);
        uniform_int_distribution<Key> KeyDistribution(0, (numeric_limits<Key>::max)() - 1);
        auto rnd = bind(KeyDistribution, rndEngine);


        cout << "Filling up to "<< initialSize << " elements." << endl;

#ifdef READ_CONTENTION
       for(size_t i = 0; i< initialSize; ++i) {
            htable->insert(MyElement(i, i));
        }
#else
        for(size_t i = 0; i< initialSize; ++i) {
            htable->insert(MyElement(rnd(), 0));
        }
#endif

        vector<thread> threads(p);
        for (size_t i = 0; i < p; ++i) {
#ifdef CONTENTION
            threads[i] = thread( workerContention, i);
#else
        #ifdef READ_CONTENTION
            threads[i] = thread( workerReadContention, i);
        #else
            threads[i] = thread( worker, i);
        #endif
#endif
        }

#ifdef PROFILE
        system("/usr/bin/sh start_profiling.sh");
#endif
        cout << "Test started." << endl;
        chrono::steady_clock::time_point start = chrono::steady_clock::now();

        go = true;

        for (auto & t : threads) {
            t.join();
        }

        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        cout << "Test completed completed." << endl;
#ifdef PROFILE
        system("/usr/bin/sh stop_profiling.sh");
#endif

        const long long timeMs = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        const long long tput = 1000ULL*n/timeMs;
        cout << "Elapsed time " << timeMs << " ms -> " << tput << " operations per second." << endl;

        cout << ">>>" << p << ";" << capacity << ";" << initialSize <<";"<< n << ";"<< insertPerc  <<";"<< timeMs << ";" <<tput << endl;

	return 0;
}


