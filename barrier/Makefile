all: barrier_bin barrier_folk2 barrier_omp barrier_folk2_gcc barrier_bin_gcc

COMMONGCC= -DNOPROFILE -march=corei7-avx -mavx -O3 -g -pthread -std=c++11
COMMON= -DPROFILE -O3 -march=core-avx2 -xCORE-AVX2 -axCORE-AVX2 -g -pthread -std=c++11

barrier_folk: barrier.cpp 
	g++ -DFOLK $(COMMON) barrier.cpp -o barrier_folk.x 

barrier_folk2: barrier.cpp
	icpc -DFOLK2 $(COMMON) barrier.cpp -o barrier_folk2.x

barrier_bin: barrier.cpp
	icpc -DBIN $(COMMON) barrier.cpp -o barrier_bin.x

barrier_folk2_gcc: barrier.cpp
	g++ -DFOLK2 $(COMMONGCC) barrier.cpp -o barrier_folk2gcc.x

barrier_bin_gcc: barrier.cpp
	g++ -DBIN $(COMMONGCC) barrier.cpp -o barrier_bingcc.x

barrier_posix: barrier.cpp
	g++ -DPOSIX_BAR $(COMMON) barrier.cpp -o barrier_posix.x 

barrier_omp: barrier.cpp
	icpc -DOMP $(COMMON) barrier.cpp -o barrier_omp.x -qopenmp

clean:
	rm -f *.x *.o
