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

class MyElement
{
public:
	typedef long long int Key;
	typedef long long int Data;
private:
	Key key;
	Data data;
	template <class T> T & asWord() { // a helper cast
		return *reinterpret_cast<T *> (this);
	}
	template <class T> const T & asWord() const {
		return *reinterpret_cast<const T *> (this);
	}
public:

    template <class F>
    void update(const MyElement & desired, F f) {
        data = f(*this,desired).data;
    }
    void update(const MyElement & desired, Overwrite f) {
        update(desired);
    }
    void update(const MyElement & desired, Increment f) {
        ++data;
    }
    void update(const MyElement & desired, Decrement f) {
        --data;
    }
	MyElement() {}
	MyElement(const Key & k, const Data & d) : key(k), data(d) {}
	Key getKey() const { return key; }
	Data getData() const { return data; } // not needed by HashTable
	static MyElement getEmptyValue() { return MyElement(numeric_limits<Key>::max(), 0); }
	bool isEmpty() const { return key == numeric_limits<Key>::max(); };
	bool isLockFree() { return true; }

    bool CAS(MyElement & expected, const MyElement & desired) { 
#ifdef __linux__
        return __sync_bool_compare_and_swap_16(&asWord<__int128>(),
			expected.asWord<__int128>(), desired.asWord<__int128>());
#else
        return _InterlockedCompareExchange128(&key, desired.data, desired.key, &expected.key) == 1;
#endif
    }
    template <class F>
    bool atomicUpdate(MyElement & expected, const MyElement & desired, F f) {
#ifdef __linux
        return  __sync_bool_compare_and_swap(&data, expected.data, f(expected,desired).data );
#else
        return _InterlockedCompareExchange64(&data,  f(expected,desired).data, expected.data) == 1;
#endif
    }
    bool atomicUpdate(MyElement & expected, const MyElement & desired, Overwrite f) {
        update(desired);
        return true;
    }
    bool atomicUpdate(MyElement & expected, const MyElement & desired, Increment f) {
#ifdef __linux__
        __sync_fetch_and_add(&data,1);
#else
        _InterlockedIncrement64(&data);
#endif
        return true;
    }
    bool atomicUpdate(MyElement & expected, const MyElement & desired, Decrement f) {
#ifdef __linux__
        __sync_fetch_and_sub(&data,1);
#else
        _InterlockedDecrement64(&data);
#endif
        return true;
    }
	MyElement(const MyElement & e) {
		/* using 'movups' instruction to implement reading using at most two atomic quad-word aligned loads (key,data)
		   If the hash table does not allow deletes this level of read atomicity is enough for this Element implementation
		   'e' is empty              :  1. empty 'key' copied 
                                        2. an insert writes a value into e
                                        3. the new 'data' copied -> copied element remains empty because the 'key' is empty
		   'e' contains a valid value:  1. valid 'key' copied 
                                        2. an update writes a new value into e 
                                        3. new 'data' copied -> returned the new element, is valid because the 'key' must remain the same
        */                                         
                asWord<__m128i>() = _mm_loadu_si128((__m128i*)&e);
	}
    MyElement & operator = (const MyElement & e) {
        asWord<__m128i>() = _mm_loadu_si128((__m128i*)&e);
        return *this;
    }
	void update(const MyElement & e) { data = e.data; }
};

typedef MyElement::Key Key;
typedef MyElement::Data Data;

struct SquareExisting
{
    MyElement operator () (const MyElement & existingValue, const MyElement & newValue) {
        return MyElement(existingValue.getKey(), existingValue.getData()*existingValue.getData());
    }
};

struct CASIncrement 
{
    MyElement operator () (const MyElement & existingValue, const MyElement & newValue) {
        return MyElement(0, existingValue.getData() + 1);
    }
};

template <class Element>
class HashTable
{
	typedef typename Element::Key Key;
	inline size_t h(const Key & k) const { return hash(k) & mask; }
        enum { maxDist = 100 };
public:
    HashTable(size_t logSize = 24) : mask((1ULL << logSize) - 1){
        	t = (Element *)_aligned_malloc((mask + 1)*sizeof(Element), 16);
        	if (t == NULL) std::bad_alloc();
		std::fill(t, t + mask + 1, Element::getEmptyValue());
	}
    void clean() {
        std::fill(t, t + mask + 1, Element::getEmptyValue());
    }
    virtual ~HashTable() { if (t) _aligned_free(t); }

#ifdef SEQUENTIAL
        bool insert(const Element & e) {
                const Key k = e.getKey();
                const size_t H = h(k), end = H + maxDist;
                for (size_t i = H; i < end; ++i) {
                        Element & current = t[i&mask];
                        if (current.getKey() == k) return false; // key already exists
                        if (current.isEmpty()) { // found free space
                                current = e;
                                return true;
                        }
                }
                throw bad_alloc(); // no space found for the element, use a table with a larger capacity
                return false;
        }
        Element find(const Key & k) {
                const size_t H = h(k), end = H + maxDist;
                for (size_t i = H; i < end; ++i) {
                        const Element & e = t[i&mask];
                        if (e.isEmpty() || (e.getKey() == k)) return e;
                }
                return Element::getEmptyValue();
        }
        bool update(const Element & e) {
                const Key k = e.getKey();
                const size_t H = h(k), end = H + maxDist;
                for (size_t i = H; i < end; ++i) {
                        Element & current = t[i&mask];
                        if (current.getKey() == k) {
                                current.update(e);
                                return true;
                        }
                        if (current.isEmpty()) return false;
                }
                return false;
        }
    template <class F = Overwrite>
    bool insertOrUpdate(const Element & e, F f = F()) {
        const Key k = e.getKey();
        const size_t H = h(k), end = H + maxDist;
        for (size_t i = H; i < end; ++i) {
            Element & current = t[i&mask];
            if (current.getKey() == k) { // key already exists
                current.update(e,f);
                return true;
            }
            if (current.isEmpty()) { // found free space
                current = e;
                return true;
            }
        }
        throw bad_alloc(); // no space found for the element, use a table with a larger capacity
        return false;
    }

#else
	bool insert(const Element & e) {
		const Key k = e.getKey();
                const size_t H = h(k), end = H + maxDist;
                for (size_t i = H; i < end; ++i) {
			// copy the element guaranteeing that a concurrent update of the source will not result in an inconsistent state
			Element current(t[i&mask]);
			if (current.getKey() == k) return false; // key already exists
			if (current.isEmpty()) { // found free space
				if (t[i&mask].CAS(current, e)) return true;
				// potentially collided with another insert
				--i; // need to reinspect position i
			}
		}
		throw bad_alloc(); // no space found for the element, use a table with a larger capacity
		return false;
	}
	Element find(const Key & k) {
                const size_t H = h(k), end = H + maxDist;
                for (size_t i = H; i < end; ++i) {
			const Element e(t[i&mask]);
			if (e.isEmpty() || (e.getKey() == k)) return e;
		}
		return Element::getEmptyValue();
	}
	bool update(const Element & e) {
		const Key k = e.getKey();
                const size_t H = h(k), end = H + maxDist;
                for (size_t i = H; i < end; ++i) {
			const Element current(t[i&mask]);
			if (current.getKey() == k) {
				t[i&mask].update(e);
				return true;
			}
			if (current.isEmpty()) return false;
		}
		return false;
	}
     template <class F = Overwrite>
     bool insertOrUpdate(const Element & e, F f = F()) {
#ifdef USE_TSX
	if(_xbegin() == _XBEGIN_STARTED) // successful transaction start
	{
		const Key k = e.getKey();
                const size_t H = h(k), end = H + maxDist;
                for (size_t i = H; i < end; ++i) {
			Element & current = t[i&mask];
                        if (current.getKey() == k) { //key already exists
                                current.update(e,f);
                                _xend();
                                return false;
                        }
                        if (current.isEmpty()) { //found free space
                                current = e;
                                _xend();
                                return true;
                        }
        	}
		_xend();
		//no space found for the element, 
		//use a table with a larger capacity
		throw bad_alloc();
    	}
#endif
        const Key k = e.getKey();
        const size_t H = h(k), end = H + maxDist;
        for (size_t i = H; i < end; ++i) {
            // copy the element guaranteeing that a concurrent update of the source will not result in an inconsistent state
            Element current(t[i&mask]);
            if (current.getKey() == k) { // key already exists
                while (!t[i&mask].atomicUpdate(current, e, f)) {
                    // potentially collided with another update
                    current = t[i&mask]; // need to reinspect position i
                }
                return false;
            }
            if (current.isEmpty()) { // found free space
                if (t[i&mask].CAS(current, e)) return true;
                // potentially collided with another insert
                --i; // need to reinspect position i
            }
        }
        throw bad_alloc(); // no space found for the element, use a table with a larger capacity
        return false;
    }

#endif

private:
	// the hash table array must be 16 byte aligned to work with the x86 16-byte compare & exchange instruction
	// then also key, data are 8-byte aligned to guarantee the atomicity of 8-byte loads
	Element * t;
	std::hash<Key> hash;
	const size_t mask;
};

template <class HashTableType>
void unitTest1()
{
	HashTableType * table = new HashTableType();
	assert(table->insert(MyElement(10, 20)));
	assert(!table->insert(MyElement(10, 20)));
	MyElement result = table->find(10);
	assert(!result.isEmpty());
	assert(result.getKey() == 10 && result.getData() == 20);
	result = table->find(20);
	assert(result.isEmpty());
    assert(table->insertOrUpdate(MyElement(20, 2), Overwrite()));
    result = table->find(20);
    assert(result.getKey() == 20 && result.getData() == 2);
    assert(!table->insertOrUpdate(MyElement(20, 3)));
    result = table->find(20);
    assert(result.getKey() == 20 && result.getData() == 3);
    assert(!table->insertOrUpdate(MyElement(20, 3), SquareExisting()));
    result = table->find(20);
    assert(result.getKey() == 20 && result.getData() == 9);
    assert(!table->insertOrUpdate(MyElement(20, 3), Increment()));
    result = table->find(20);
    assert(result.getKey() == 20 && result.getData() == 10);
    assert(!table->insertOrUpdate(MyElement(20, 3), Decrement()));
    result = table->find(20);
    assert(result.getKey() == 20 && result.getData() == 9);
    assert(!table->insertOrUpdate(MyElement(20, 3), 
            [] (const MyElement & existingValue, const MyElement & newValue) {
                return MyElement(existingValue.getKey(), 
                    existingValue.getData()*existingValue.getData()*existingValue.getData());
            } 
        ));
    result = table->find(20);
    assert(result.getKey() == 20 && result.getData() == 9*9*9);
	delete table;
}

template <class HashTableType>
void unitTest2()
{
	HashTableType * table = new HashTableType();
	for (int i = 0; i< 15485867 * 0.9; ++i)
	{
		assert(!table->update(MyElement(i, i)));
	}
	for (int i = 0; i< 15485867 * 0.9; ++i)
	{
		assert(table->insert(MyElement(i, i + 1)));
	}
	for (int i = 0; i< 15485867 * 0.9; ++i)
	{
		MyElement result = table->find(i);
		assert(!result.isEmpty());
		assert(result.getKey() == i && result.getData() == i + 1);
	}
	for (int i = 0; i< 15485867 * 0.9; ++i)
	{
		assert(!table->insert(MyElement(i, i + 1)));
	}
	for (int i = 0; i< 15485867 * 0.9; ++i)
	{
		assert(table->update(MyElement(i, i + 2)));
	}
	for (int i = 0; i< 15485867 * 0.9; ++i)
	{
		MyElement result = table->find(i);
		assert(!result.isEmpty());
		assert(result.getKey() == i && result.getData() == i + 2);
	}
	delete table;
}


using namespace std;
int insertPerc;
size_t opsPerThreads, initialSize;
HashTable<MyElement> * htable = NULL;
atomic<bool> go(false);


void worker( const size_t iPE // the ID of PE or thread
           )
{
       mt19937 rndEngine(iPE);
       uniform_int_distribution<MyElement::Key> KeyDistribution(0, (numeric_limits<MyElement::Key>::max)() - 1);
       auto valueRnd = bind(KeyDistribution, rndEngine);
       uniform_int_distribution<int> OpDistribution(0, 100);
       auto opRnd = bind(OpDistribution, rndEngine);

       while(!go);

       for(size_t i =0 ; i< opsPerThreads; ++i) {
           if(insertPerc <= opRnd()) {
              htable->find(valueRnd());
           } else {
              htable->insert(MyElement(valueRnd(),0));
           }
       }
}

#ifdef CASINCREMENT
typedef CASIncrement IncrementT;
#else
typedef Increment IncrementT;
#endif

#ifdef CONTENTION

void workerContention( const size_t iPE // the ID of PE or thread
           )
{
       mt19937 rndEngine(iPE);
       uniform_int_distribution<MyElement::Key> KeyDistribution(0, (numeric_limits<MyElement::Key>::max)() - 1);
       auto valueRnd = bind(KeyDistribution, rndEngine);
       uniform_int_distribution<int> OpDistribution(0, 100);
       auto opRnd = bind(OpDistribution, rndEngine);

       while(!go);

       for(size_t i =0 ; i< opsPerThreads; ++i) {
           if(insertPerc <= opRnd()) {
              htable->insertOrUpdate(MyElement(valueRnd(), 0), IncrementT());
           } else {
              htable->insertOrUpdate(MyElement(0xabcdabcd, 0), IncrementT());
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
              sum += htable->find(valueRnd()).getData();
           } else {
              sum += htable->find(42).getData();
           }
       }

       dummy += sum;
}

#endif

int main(int argc, char* argv[])
{
    //static_assert(sizeof(MyElement) == sizeof(atomic<long long int>), "MyElement does not fit into the single atomic 64-bit word");
    assert(MyElement().isLockFree());

    unitTest1<HashTable<MyElement> >();
    unitTest2<HashTable<MyElement> >();
    cout << "tests passed, enter x to exit" << endl;

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

    htable = new HashTable<MyElement>(capacity);
    capacity = 1ULL << capacity;
    cout << "Using capacity of "<< capacity << " elements. "<< capacity*sizeof(MyElement)/(1024ULL*1024ULL)<< " Mbytes " << endl;

    const size_t iterations = 1;

    for(int j=0;j<iterations;++j) {

        mt19937 rndEngine(0xabcdef);
        uniform_int_distribution<MyElement::Key> KeyDistribution(0, (numeric_limits<MyElement::Key>::max)() - 1);
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
        htable->clean();
        
   }

   return 0;
}


