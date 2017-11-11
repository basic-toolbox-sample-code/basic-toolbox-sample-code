all: psamplesort stdsort parallelmode parallelmode_numa tbb_numa psamplesort_pool psamplesort_numa_pool 

TBB_INCL=tbb/include/
TBB_LIB=tbb/build/linux_intel64_gcc_cc4.3_libc2.11.3_kernel3.13.5_release/

TBB_OPTIONS=-I${TBB_INCL} -L${TBB_LIB} -ltbb -ltbbmalloc

OPTIONS=-march=corei7-avx -O3 -g -pthread -std=c++11

psamplesort: psamplesort.cpp
	g++ ${OPTIONS} psamplesort.cpp -o psamplesort.x

psamplesort_pool: psamplesort.cpp
	g++ -DUSE_POOL_MEMORY  ${OPTIONS} ${TBB_OPTIONS} psamplesort.cpp -o psamplesort_pool.x

psamplesort_numa_pool: psamplesort.cpp
	g++ -DNOPROFILE -DUSE_POOL_MEMORY -DNUMA ${OPTIONS} ${TBB_OPTIONS} psamplesort.cpp -o psamplesort_numa_pool.x

stdsort: othersort.cpp
	g++ -DSTD_SORT ${OPTIONS} othersort.cpp -o stdsort.x

parallelmode: othersort.cpp
	g++ -DSTD_PARALLEL_MODE -fopenmp ${OPTIONS} othersort.cpp -o std_parallel_mode.x

parallelmode_numa: othersort.cpp
	g++ -DNOPROFILE -DNUMA -DSTD_PARALLEL_MODE -fopenmp ${OPTIONS} othersort.cpp -o std_parallel_mode_numa.x

tbb_numa: othersort.cpp
	g++ -DTBB -DNUMA ${OPTIONS} ${TBB_OPTIONS} othersort.cpp -o tbb_numa.x

clean:
	rm -f *.x *.o
