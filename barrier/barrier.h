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

#include <vector>
#include <algorithm>
#include <atomic>

using namespace std;

inline void relax() {
    asm volatile("pause": : :"memory");
}

class Barrier {
    atomic<int> a __attribute__((aligned(64)));
    char padding[64 - sizeof(atomic<int>)];
public:
    Barrier() : a(0) {}
    void
#ifdef NOINLINE
__attribute__ ((noinline)) 
#endif
            wait(int p) {
        int pp = p;
        a.compare_exchange_strong(pp, 0);
        ++a;
        while (a != p) {
		relax();
       }
    }
};

class FolkBarrier2 
{
public:
    FolkBarrier2(): c(0), epoch(0) {}

    void wait(int p)
    {
        register const int startEpoch = epoch;
        if (c.fetch_add(1)== p-1)
        {
            c = 0;
            ++epoch; 
        }
        else
        {
            while(epoch == startEpoch)
                relax();
        }
        //atomic_thread_fence(memory_order_seq_cst);
    }

protected:
    atomic<int> c __attribute__((aligned(64)));
    char pad[64 - sizeof(atomic<int>)];
    volatile int epoch __attribute__((aligned(64)));
};

class BinomialBarrier
{
	BinomialBarrier();
//        typedef  std::pair<atomic<int> , char [64-sizeof(atomic<int>)]> alignedAtomicInt;
        typedef  std::pair<volatile int , char [64-sizeof(int)]> alignedAtomicInt;
        struct EpochNumChildren {
           EpochNumChildren() : epoch(0), numChildren(0) {}
           int epoch, numChildren; // co-locate to avoid additional cache miss
           char padding[64 - 2*sizeof(int)];
        };
        template <class It>
	void initNumChildren(It begin, int size)
	{
		if(size >= 2)
		for(int i = 1; i < size; i *= 2)
		{
			++(begin->numChildren);
			int child = i;
			int childSize = i;
			if(size < child + childSize) childSize = size - child;	
			initNumChildren(begin + child, childSize); 
		}
	}
public:

#define USE_VECTOR 1 
#if USE_VECTOR
    BinomialBarrier(int p) : readyEpoch(p), peData(p) {
                initNumChildren(peData.begin(), p);
                for(auto && e : readyEpoch) {
                        e.first = 0;
                }
    }
#else
    BinomialBarrier(unsigned p) : peData(new EpochNumChildren[p]), readyEpoch(new alignedAtomicInt[p]) {
                std::fill(peData, peData + p, EpochNumChildren());
                initNumChildren(peData, p);
                for(int i = 0; i < p ; ++i)
                    readyEpoch[i].first = 0;
    }
    ~BinomialBarrier() {
        delete [] peData;
        delete [] readyEpoch;
    }
#endif

	void
#ifdef NOINLINE
 __attribute__ ((noinline)) 
#endif
            sync(const int iPE) {
		register const int myEpoch = ++(peData[iPE].epoch);
                register const int numC = peData[iPE].numChildren;
		for(int i=0; i < numC; ++i) {
                        auto & e = readyEpoch[iPE + (1<<i)].first;
			while(e != myEpoch)
				relax();
		}
		readyEpoch[iPE].first = myEpoch;
		atomic_thread_fence(memory_order_seq_cst);
                auto & e = readyEpoch[0].first;
		while(e != myEpoch)
			relax();
	}
        void wait(const int iPE, int) {
            sync(iPE);
        }
private:
#if USE_VECTOR
        vector<alignedAtomicInt> readyEpoch; 
	vector<EpochNumChildren> peData;
#else
        alignedAtomicInt * readyEpoch;
        EpochNumChildren * peData;
#endif
};


