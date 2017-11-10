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
#include "barrier/barrier.h"

typedef int NodeId;

#define SEQ_PREFIX_SUM

#define INVALID_NODE_ID (numeric_limits<NodeId>::min())

#define check_align8(x) assert((unsigned long long int)(&(x[0])) % 8 == 0);

using namespace std;

#include <pthread.h>

void pinToCore(size_t core)
{
   // pin the current thread to the specified core to make NUMA memory access mapping deterministic
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

class Graph {
    NodeId n;
    size_t m;
    std::vector<size_t> begin;
    std::vector<NodeId> adj;
    Graph();
    NodeId readNumber(ifstream & f) {
        NodeId result;
        if (!(f >> result))
        {
            cerr << "Error reading file" << endl;
            throw exception();
        }
        return result;
    }
public:
    Graph(const string & filename, int p) {
        ifstream f(filename);
        std::string formatName;
        f >> formatName;
        if (formatName != "AdjacencyGraph") {
            cerr << "Wrong format "<< formatName<< ", expecting AdjacencyGraph format." << endl;
            throw exception();
        }
        n = readNumber(f);
        cout << "n=" << n << endl;
        m = readNumber(f);
        cout << "m=" << m << endl;
        begin.reserve(n+1);
        size_t i;
        cout << "Reading nodes" << endl;
        for (i = 0; i < n; ++i) {
            begin.push_back(readNumber(f));
            if ((i % 200000) == 0)
                cout << "*"<< flush;
        }
        cout << endl;
        begin.push_back(m);
        adj.reserve(m);
        cout << "Reading edges" << endl;
        for (i = 0; i < m; ++i) {
            adj.push_back(readNumber(f));
            if ((i % 200000) == 0)
                cout << "*" << flush;
        }
        cout << endl;
        check_align8(begin); // must be aligned for optimal performance
        check_align8(adj); // must be aligned for optimal performance
    }
    virtual ~Graph() {}
    vector<NodeId>::const_iterator beginNeighbor(const NodeId v) const {
        return adj.cbegin() + begin[v];
    }
    vector<NodeId>::const_iterator endNeighbor(const NodeId v) const {
        return adj.cbegin() + begin[v + 1];
    }
    NodeId getN() const { return n; }
    size_t getM() const { return m; }
};

void seqBFS(const Graph & g, const NodeId s, vector<NodeId> & d, vector<NodeId> & parent) {
#ifdef WRITE_D
    d.resize(g.getN());
    check_align8(d); // must be aligned for optimal performance
    fill(d.begin(), d.end(), (numeric_limits<NodeId>::max)());
#endif
    parent.resize(g.getN());
    check_align8(parent); // must be aligned for optimal performance
    fill(parent.begin(), parent.end(), INVALID_NODE_ID);

#ifdef WRITE_D
    d[s] = 0;
#endif
    parent[s] = s;

    vector<NodeId> Q,Qp;
    Q.reserve(g.getN());
    Qp.reserve(g.getN());
    Q.push_back(s);
    
    for (NodeId l = 1; !Q.empty(); ++l) {
        cout << "computing level " << l << " Q's size is " << Q.size() << endl;
        for (const NodeId u: Q) {
            auto vIter = g.beginNeighbor(u);
            auto vEnd = g.endNeighbor(u);
            for (; vIter != vEnd; ++vIter) {
                const NodeId v = *vIter;
                if (parent[v] == INVALID_NODE_ID) {
                    parent[v] = u;
#ifdef WRITE_D
                    d[v] = l;
#endif
                    Qp.push_back(v);
                }
            }
        }
        swap(Q, Qp); // avoid copying
        Qp.clear();
    }
}

char pad0[64];
Barrier barriers[63];
char pad1[64];

BinomialBarrier * barrier = NULL;

#define BARRIER(iPE) barriers[barrierCnt++].wait(p);

//#define BARRIER(iPE) barrier->sync(iPE); barrierCnt++;

template <class Iterator, class F>
void prefixSum(Iterator outBegin, Iterator outEnd, int iPE, int p, Iterator tmp, F f) {
    const size_t begin = (outEnd-outBegin)*iPE/p;
    const size_t end = (outEnd-outBegin)*(iPE+1)/p;
    size_t sum = 0, i;
    for (i= begin; i != end ; ++i) {
        *(outBegin + i) = (sum += f(i)); 
    }
    *(tmp + iPE) = sum;
    int barrierCnt = 62;
    BARRIER(iPE)
    size_t offset = 0;
    for(i=0; i< iPE; ++i)
	offset += *(tmp + i);
    for(i=begin; i!=end; ++i)
        *(outBegin + i) += offset;
}

// pad to avoid false-sharing
typedef pair<vector<NodeId>, char [64-sizeof(vector<NodeId>)] > PaddedVector;

#define DEBUG(x) 

void worker(const unsigned iPE, // the ID of PE or thread
    const int p,
    const NodeId s,
    const Graph & g,
    vector<NodeId> & Q,
    vector<PaddedVector> & Qp,
    vector<size_t> & sigma,
    atomic<bool> & done,
    vector<NodeId> & d,
    vector<NodeId> & parent,
    vector<size_t> & tmp
    )
{
    chrono::steady_clock::time_point last = start;    
    const size_t beginI = iPE*g.getN()/p;
    const size_t endI = (iPE + 1)*g.getN()/p;
#ifdef WRITE_D
    fill(d.begin() + beginI, d.begin() + endI, (numeric_limits<NodeId>::max)());
#endif
    fill(parent.begin() + beginI, parent.begin() + endI, INVALID_NODE_ID);
    if(s >= beginI && s < endI) {
#ifdef WRITE_D
        d[s] = 0;
#endif
        parent[s] = s;
    }

    NodeId l = 1;

    for (; !done; ++l) {
        int barrierCnt = 1;
        DEBUG(if (iPE == 0) cout << "computing level " << l << " Q's size is " << Q.size() << endl;)

        if(Q.size()/p < 10) // too little work for a PE
        { // PE 0 does all the work
            if (iPE == 0) {
                sigma[0] = g.endNeighbor(Q[0]) - g.beginNeighbor(Q[0]);
                for (NodeId j = 1; j < Q.size(); ++j)
                    sigma[j] = sigma[j-1] + (g.endNeighbor(Q[j]) - g.beginNeighbor(Q[j]));
            }
        }
	else
        prefixSum(sigma.begin(), sigma.end(), iPE, p, tmp.begin(), 
                  [&] (int j) {
                    return g.endNeighbor(Q[j]) - g.beginNeighbor(Q[j]);
               }
           );
        BARRIER(iPE)
        DEBUG(if(iPE == 0) { cout << "Delay: " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - last).count()<<  " ms for barrier "<< (barrierCnt-1) << endl; last = chrono::steady_clock::now(); })

        Qp[iPE].first.clear();
        size_t ml = sigma.back();
        Qp[iPE].first.reserve(2*ml/p);
        size_t curQ = upper_bound(sigma.cbegin(), sigma.cend(), iPE*ml / p) - sigma.cbegin();
        const size_t endQ = upper_bound(sigma.cbegin(), sigma.cend(), (iPE+1)*ml / p) - sigma.cbegin();
        //cout << ml << " " << curQ << " " << endQ << endl;
        for (; curQ != endQ; ++curQ) {
            const NodeId u = Q[curQ];
            auto vIter = g.beginNeighbor(u);
            auto vEnd = g.endNeighbor(u);
            for (; vIter != vEnd; ++vIter) {
                const NodeId v = *vIter;
#ifdef USE_CAS

#ifdef __linux__
                if(parent[v] == INVALID_NODE_ID && __sync_bool_compare_and_swap(&parent[v], INVALID_NODE_ID, u)) {
#else
                if(parent[v] == INVALID_NODE_ID && _InterlockedCompareExchange((long *)&(parent[v]), u, INVALID_NODE_ID) == INVALID_NODE_ID) {
#endif

#else
                if(parent[v] == INVALID_NODE_ID) {
                    parent[v] = u;
#endif
#ifdef WRITE_D
                    d[v] = l;
#endif
                    Qp[iPE].first.push_back(v);
                }
            }
        }
        BARRIER(iPE)
        DEBUG(if(iPE == 0) { cout << "Delay: " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - last).count()<<  " ms for barrier "<< (barrierCnt-1) << endl; last = chrono::steady_clock::now(); })
        size_t outPos = 0;
        for(int j=0; j < iPE; ++j) outPos += Qp[j].first.size();
        if(iPE == p-1) {
            Q.resize(outPos + Qp[iPE].first.size());
            sigma.resize(Q.size());
            if (Q.empty()) done = true;
        }
        BARRIER(iPE)
        DEBUG(if(iPE == 0) { cout << "Delay: " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - last).count()<<  " ms for barrier "<< (barrierCnt-1) << endl; last = chrono::steady_clock::now(); } )

        copy(Qp[iPE].first.cbegin(), Qp[iPE].first.cend(), Q.begin() + outPos);

        BARRIER(iPE)
        DEBUG(if(iPE == 0) { cout << "Delay: " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - last).count()<<  " ms for barrier "<< (barrierCnt-1) << endl; last = chrono::steady_clock::now(); })
    }
}

void pBFS(unsigned p, const Graph & g, const NodeId s, vector<NodeId> & d, vector<NodeId> & parent)
{
#ifdef WRITE_D
    d.resize(g.getN());
    check_align8(d); // must be aligned for optimal performance
#endif
    parent.resize(g.getN());
    check_align8(parent); // must be aligned for optimal performance
    vector<NodeId> Q;
    Q.reserve(g.getN());    
    Q.push_back(s);
    vector<PaddedVector> Qp(p);

    vector<size_t> sigma(1);
    sigma.reserve(g.getN());
    vector<size_t> tmp(p);
    atomic<bool> done(false);
    vector<thread> threads(p);
    for (unsigned i = 0; i < p; ++i) { // go parallel
        threads[i] = thread(worker, i, p, s, ref(g), ref(Q), ref(Qp), ref(sigma), ref(done), ref(d), ref(parent), ref(tmp));
    }
    for (auto & t : threads) t.join();
}

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
         cerr << "Parameters: input threads"<< endl;
         return -1;
    }
    int p = atoi(argv[2]);
    Graph g(argv[1], p);
    barrier = new BinomialBarrier(p);
    vector<NodeId> d, parent;
#ifdef PROFILE
    system("/usr/bin/sh start_profiling.sh");
#endif
    cout << "Test started." << endl;
    start = chrono::steady_clock::now();
#ifdef SEQUENTIAL
    seqBFS(g, 0, d, parent);
#else
    pBFS(p, g, 0, d, parent);
#endif
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Test completed." << endl;
#ifdef PROFILE
    system("/usr/bin/sh stop_profiling.sh");
#endif

    delete barrier;
    const long long timeMs = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "Elapsed time " << timeMs << " ms" << endl;
    cout << ">>>"<< argv[1] <<";"<<g.getN()<<";"<<g.getM()<<";"<<p<<";"<<timeMs<<endl;
    return 0;
}

