all: phashtable_c_tsx phashtable phashtable_seq phashtable_c phashtable_c_seq phashtable_cas tbb_unordered_map tbb_hash_map phashtable_rc phashtable_rc_seq tbb_hash_map_rc tbb_unordered_map_rc tbb_hash_map_c tbb_unordered_map_c

TBB=-std=c++11 -Itbb/include/ -Ltbb/build/linux_intel64_gcc_cc4.3_libc2.11.3_kernel3.13.5_release/ -ltbb -ltbbmalloc 

GCC=g++
OPT=-DNOPROFILE -march=corei7-avx -O3 -g -pthread -std=c++11 -mrtm

phashtable: PHashTable.cpp 
	$(GCC) $(OPT) PHashTable.cpp -o phashtable.x 

phashtable_seq: PHashTable.cpp
	$(GCC) -DSEQUENTIAL $(OPT) -std=c++11 PHashTable.cpp -o phashtable_seq.x

phashtable_c: PHashTable.cpp
	$(GCC) -DCONTENTION $(OPT) PHashTable.cpp -o phashtable_c.x

phashtable_c_tsx: PHashTable.cpp
	$(GCC) -DCONTENTION -DUSE_TSX $(OPT) PHashTable.cpp -o phashtable_c_tsx.x

phashtable_rc: PHashTable.cpp
	$(GCC) -DREAD_CONTENTION $(OPT) PHashTable.cpp -o phashtable_rc.x

phashtable_rc_seq: PHashTable.cpp
	 $(GCC) -DREAD_CONTENTION -DSEQUENTIAL $(OPT) PHashTable.cpp -o phashtable_rc_seq.x

phashtable_c_seq: PHashTable.cpp
	$(GCC) -DCONTENTION -DSEQUENTIAL $(OPT) PHashTable.cpp -o phashtable_c_seq.x	

phashtable_cas: PHashTable.cpp
	$(GCC) -DCONTENTION -DCASINCREMENT $(OPT) PHashTable.cpp -o phashtable_cas.x

tbb_unordered_map: TBBHashTable.cpp
	 $(GCC) $(OPT) $(TBB) TBBHashTable.cpp -o tbb_unordered_map.x

tbb_hash_map: TBBHashTable.cpp
	 $(GCC) -DTBB_HASH_MAP $(OPT) $(TBB) TBBHashTable.cpp -o tbb_hash_map.x

tbb_unordered_map_rc: TBBHashTable.cpp
	$(GCC) -DREAD_CONTENTION $(OPT) $(TBB) TBBHashTable.cpp -o tbb_unordered_map_rc.x

tbb_hash_map_rc: TBBHashTable.cpp
	$(GCC) -DREAD_CONTENTION -DTBB_HASH_MAP $(OPT) $(TBB) TBBHashTable.cpp -o tbb_hash_map_rc.x

tbb_unordered_map_c: TBBHashTable.cpp
	$(GCC) -DCONTENTION $(OPT) $(TBB) TBBHashTable.cpp -o tbb_unordered_map_c.x

tbb_hash_map_c: TBBHashTable.cpp
	$(GCC) -DCONTENTION -DTBB_HASH_MAP $(OPT) $(TBB) TBBHashTable.cpp -o tbb_hash_map_c.x


clean:
	rm -f *.x *.o
