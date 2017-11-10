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

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <chrono>
#include <thread>
#include <atomic>
#include <assert.h>
#include <mutex>
#include <pthread.h>
#ifdef OMP
#include <omp.h>
#endif
#include "barrier.h"
#include <pthread.h>

//#define HT

void pinToCore(size_t core)
{
   // pin the current thread to the specified core to make NUMA memory access mapping deterministic
#ifndef HT 
   size_t s = core / 16;
   size_t c = core % 16;
   core = c + s*18;
#else
   size_t s = core / 32;
   size_t c = (core % 32)/2;
   size_t h = core % 2;
   core = c + s*18 + h*18*4;
#endif

   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(core, &cpuset);

   int res = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
   if(res) {
      cerr << "could not set thread affinity to "<< core << endl;
      throw exception();
   }
}

chrono::steady_clock::time_point start;

#ifdef PROFILE
//#define NOINLINE
#endif

char pad0[64];
Barrier barriers[63];
char pad1[64];
FolkBarrier2 folkBarrier2;
char pad2[64];
BinomialBarrier * binBarrier = NULL;
char padabcd[64];
pthread_barrier_t posixBarrier;
char paddef[64];

#ifdef FOLK
#define BARRIER(iPE) barriers[barrierCnt++].wait(p);
#endif

#ifdef BIN
#define BARRIER(iPE) binBarrier->sync(iPE); barrierCnt++;
#endif

#ifdef POSIX_BAR
#define BARRIER(iPE) pthread_barrier_wait(&posixBarrier);
#endif

#ifdef FOLK2
#define BARRIER(iPE) folkBarrier2.wait(p);
#endif

#ifdef OMP
#define BARRIER(iPE) _Pragma("omp barrier")
#endif

#define DEBUG(x) 

double * a;
double * b;
bool go = false;

void worker(const size_t iPE, // the ID of PE or thread
    const size_t n,
    const size_t p,
    const size_t k
    )
{
    pinToCore(iPE);
    while(!go) asm volatile("pause": : :"memory");
    for(int i = 0; i< k; i+=2)
    {
        int barrierCnt = 0;
        size_t j = 1 + iPE*n/p;
        const size_t end = 1 + (iPE+1)*n/p;
        // if(iPE == 10) std::cout << j << " " << end << endl;
        for(;j<end;++j) b[j] = (a[j-1] + a[j] + a[j+1])/3.0;
        BARRIER(iPE)
        j = 1 + iPE*n/p;
        for(;j<end;++j) a[j] = (b[j-1] + b[j] + b[j+1])/3.0;
        BARRIER(iPE)
    }
}

int main(int argc, char* argv[])
{
    if(argc < 4)
    {
         cerr << "Parameters: n threads k"<< endl;
         return -1;
    }
    size_t n = atoi(argv[1]);
    size_t p = atoi(argv[2]);
    size_t k = atoi(argv[3]);
    a = new double[n+2];
    b = new double[n+2];
    std::fill(a, a + n + 2, 0);
    std::fill(b, b + n + 2, 0);
    for(int i = 0; i <  n + 2; ++i)
    {
         a[i] = i*4.0;
         b[i] = i*100;
    }
    binBarrier = new BinomialBarrier(p);
    pthread_barrier_init(&posixBarrier, NULL, p);
#ifndef OMP
    vector<thread> threads(p);
    for (size_t i = 0; i < p; ++i) { // go parallel
        threads[i] = thread(worker, i, n, p, k);
    }
#endif
#ifdef PROFILE
    system("/usr/bin/sh start_profiling.sh");
#endif
    cout << "Test started." << endl;
    start = chrono::steady_clock::now();
    go = true;
#ifndef OMP
    for (auto & t : threads) t.join();
#else
    #pragma omp parallel
    {
        worker(omp_get_thread_num(), n, p, k);
    }
#endif
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Test completed." << endl;
#ifdef PROFILE
    system("/usr/bin/sh stop_profiling.sh");
#endif

    double sum = 0;
    for(int i = 0; i <  n + 2; ++i)
    {
         sum += a[i];
         sum += b[i];
    }
    std::cout << "checksum: "<< sum << std::endl;

    //delete binBarrier;
    const long long timeMs = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "Elapsed time " << timeMs << " ms" << endl;
    cout << ">>>"<< n <<";"<<p<<";"<<k<<";"<<timeMs<<endl;
    return 0;
}

