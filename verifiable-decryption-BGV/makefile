
CPP=g++

FLAGS=-g -O2

AUTOFLAGS=-std=c++11 -pthread -march=native

LDLIBS=-lntl -lgmp -lm

FILES= main.cpp encryption.cpp commitment.cpp helper.cpp linearityproof.cpp shortnessproof.cpp fiatshamir.cpp

run:
	clear
	gcc $(FLAGS) -c fips202.c -o fips202.o
	$(CPP) $(FLAGS) $(AUTOFLAGS) $(FILES) fips202.o -o main $(LDLIBS)
	./main

clean:
	rm -rf main fips202.o