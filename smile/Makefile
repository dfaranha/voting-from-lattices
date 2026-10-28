CC ?= /usr/bin/cc
CPP = g++
CFLAGS += -Wall -Wextra  -Wredundant-decls \
  -Wshadow -Wpointer-arith -O3 -fomit-frame-pointer -march=native
NISTFLAGS += -Wno-unused-result -O3 -fomit-frame-pointer
RM = /bin/rm

SOURCES = polyvec.c poly.c ntt.c util.c
SOURCESKECCAK = $(SOURCES) fips202.c symmetric-shake.c
HEADERS = params.h polyvec.h poly.h ntt.h reduce.h util.h rej.h
HEADERSKECCAK = $(HEADERS) fips202.h

.PHONY: all speed shared clean

all: \
  test_rlwe \
  test_shuffle


%_short.o: %.c params.h Makefile
	$(CC) $(CFLAGS) -DSHORT=1 -c $< -o $@

%_shuffle.o: %.c params.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

test_mx: $(SOURCESKECCAK) rlwe.c comm_shuffle.o comm_short.o shortness_proof_short.o mxproof.c test_shuffle.c randombytes.c 
	$(CC) $(CFLAGS) $^ -o test_mx -lm

test_shuffle: $(SOURCESKECCAK) rlwe.c comm_shuffle.o comm_short.o shortness_proof_short.o mxproof.c test_shuffle.c randombytes.c 
	$(CC) $(CFLAGS) -DNO_SHORT=1 $^ -o test_shuffle -lm

test_shortness: $(SOURCESKECCAK) randombytes.c shortness_proof_short.o comm_short.o comm_shuffle.o rlwe.c test_shortness.c
	$(CC) $(CFLAGS) $^ -o test_shortness -lm

test_rlwe: $(SOURCESKECCAK) randombytes.c rlwe.c test_rlwe.c 
	$(CC) $(CFLAGS) $^ -o test_rlwe -lm

test_bgv: $(SOURCESKECCAK) randombytes.c bgv.c test_bgv.c 
	$(CC) $(CFLAGS) $^ -o test_bgv -lm

test_smile: $(SOURCESKECCAK) randombytes.c bgv.c comm.c gaussian_ct.cpp rej.c smile.c test_smile.c 
	$(CPP) $(CFLAGS) -c gaussian_ct.cpp -o gaussian.o 
	$(CPP) $(CFLAGS) -DGAUSSIAN test_smile.c gaussian.o $(SOURCESKECCAK) randombytes.c bgv.c comm.c rej.c smile.c -o test_smile -lm -lgmp -lmpfr

test_gaussian: randombytes.c cpucycles.c params.h gaussian_ct.cpp
	$(CPP) $(CFLAGS) -DMAIN=1 $^ -o test_gaussian.o -lm -march=native 

clean:
	-$(RM) -rf *.gcno *.gcda *.gch *.lcov *.o *.so
	-$(RM) -rf test_rlwe
	-$(RM) -rf test_bgv
	-$(RM) -rf test_shuffle
	-$(RM) -rf test_shortness
	-$(RM) -rf test_mx
	-$(RM) -rf test_smile
	-$(RM) -rf test_gaussian
