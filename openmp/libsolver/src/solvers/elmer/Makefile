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

CMDS = bicgstab2 cg
SOLVER_SRC = bicgstab2.c cg.c
SOLVER_OBJ = bicgstab2.o cg.o

all: $(CMDS)

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

bicgstab2: $(COMMON_OBJS) bicgstab2.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

cg: $(COMMON_OBJS) cg.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
