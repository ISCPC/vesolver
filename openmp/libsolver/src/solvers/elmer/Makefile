WITH_MKL=true
WITH_MKL_SBLAS=true
WITH_ONEAPI=false

COMMON_CXX_SRC = main.cpp SpMatrix.cpp LinearSolverC.cpp
COMMON_CXX_SRCDIR = ../standalone
COMMON_C_SRC = Matrix.c
COMMON_C_SRCDIR = ../../common

CXXFLAGS = -g -D_STANDALONE
CFLAGS = -g -D_STANDALONE

include ../standalone/Make.icc
#include ../standalone/Make.gcc
#include ../standalone/Make.SXAT
#include ../standalone/Make.fcc

#######
SOLVERDIR = .

CMDS = bicgstab2 elmer_cg
SOLVER_SRC = bicgstab2.c elmer_cg.c
SOLVER_OBJ = bicgstab2.o elmer_cg.o

all: $(CMDS)

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

bicgstab2: $(COMMON_OBJS) bicgstab2.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

elmer_cg: $(COMMON_OBJS) elmer_cg.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
