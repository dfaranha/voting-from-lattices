# Cumulative voting with homomorphic tallying. The ring arithmetic is
# include/zq.hpp; GMP is the only library the voter's proof needs.
CXX      = g++
CC       = gcc
OPT      = -O3 -march=native -mtune=native
CXXFLAGS = -std=c++17 $(OPT) -Wall -I include -I lib $(CONFIG)
CFLAGS   = $(OPT) -Wall -I lib
# Extra defines for a build, e.g. make CONFIG="-DNV=10 -DNC=8 -DLC=4 -DLS=8".
CONFIG   =
LIBS     = -lgmpxx -lgmp


OBJ      = obj
LIBOBJ   = $(OBJ)/fips202.o $(OBJ)/sp800-185.o $(OBJ)/cpucycles.o

.PHONY: all clean
all: vericrypt

$(OBJ):
	mkdir -p $(OBJ)

$(OBJ)/%.o: lib/%.c | $(OBJ)
	$(CC) $(CFLAGS) -c $< -o $@

vericrypt: main.cpp src/vericrypt.cpp include/*.hpp $(LIBOBJ)
	$(CXX) $(CXXFLAGS) main.cpp src/vericrypt.cpp $(LIBOBJ) -o $@ $(LIBS)

clean:
	rm -rf $(OBJ) vericrypt
