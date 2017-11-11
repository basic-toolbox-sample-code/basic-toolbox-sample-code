all: pbfs seqbfs pbfs_cas

pbfs: pbfs.cpp 
	g++ -DNOSTATS -DNOPROFILE -march=corei7-avx -O2 -g -pthread -std=c++11 pbfs.cpp -o pbfs.x 

pbfs_cas: pbfs.cpp
	g++ -DNOPROFILE -DUSE_CAS -march=corei7-avx -O2 -g -pthread -std=c++11 pbfs.cpp -o pbfs_cas.x

seqbfs: pbfs.cpp
	g++ -DNOPROFILE -march=corei7-avx -O2 -g -pthread -std=c++11 pbfs.cpp -o seqbfs.x -DSEQUENTIAL 

clean:
	rm -f *.x *.o
