# Cumulative voting with homomorphic tallying. The ring arithmetic is
# include/zq.hpp; GMP is the only library the voter's proof needs. election
# also links the proof of correct decryption of ../verifiable-decryption-BGV,
# which needs NTL.
CXX      = g++
CC       = gcc
OPT      = -O3 -march=native -mtune=native
CXXFLAGS = -std=c++17 $(OPT) -Wall -I include -I lib $(CONFIG)
CFLAGS   = $(OPT) -Wall -I lib
# Extra defines for a build, e.g. make CONFIG="-DNV=10 -DNC=8 -DLC=4 -DLS=8".
CONFIG   =
LIBS     = -lgmpxx -lgmp

VERDEC   = ../verifiable-decryption-BGV
VERDEC_SRC = encryption.cpp commitment.cpp helper.cpp linearityproof.cpp \
             shortnessproof.cpp
VERDEC_OBJ = $(addprefix $(OBJ)/verdec_,$(VERDEC_SRC:.cpp=.o))

OBJ      = obj
LIBOBJ   = $(OBJ)/fips202.o $(OBJ)/sp800-185.o $(OBJ)/cpucycles.o

.PHONY: all clean
all: vericrypt election

$(OBJ):
	mkdir -p $(OBJ)

$(OBJ)/%.o: lib/%.c | $(OBJ)
	$(CC) $(CFLAGS) -c $< -o $@

vericrypt: main.cpp src/vericrypt.cpp include/*.hpp $(LIBOBJ)
	$(CXX) $(CXXFLAGS) main.cpp src/vericrypt.cpp $(LIBOBJ) -o $@ $(LIBS)

# verifiable-decryption-BGV is compiled as its own makefile does, -O2, and
# never in the same translation unit as this code: see include/verdec_bridge.h.
$(OBJ)/verdec_%.o: $(VERDEC)/%.cpp $(VERDEC)/*.h | $(OBJ)
	$(CXX) -std=c++17 -O2 -march=native -pthread -I $(VERDEC) -c $< -o $@

$(OBJ)/verdec_bridge.o: src/verdec_bridge.cpp include/verdec_bridge.h $(VERDEC)/*.h | $(OBJ)
	$(CXX) -std=c++17 -O2 -march=native -pthread -I include -I $(VERDEC) -c $< -o $@

election: election.cpp src/vericrypt.cpp include/*.hpp include/verdec_bridge.h \
		$(LIBOBJ) $(OBJ)/verdec_bridge.o $(VERDEC_OBJ)
	$(CXX) $(CXXFLAGS) -pthread election.cpp src/vericrypt.cpp $(LIBOBJ) \
		$(OBJ)/verdec_bridge.o $(VERDEC_OBJ) -o $@ -lntl $(LIBS) -lm

clean:
	rm -rf $(OBJ) vericrypt election
