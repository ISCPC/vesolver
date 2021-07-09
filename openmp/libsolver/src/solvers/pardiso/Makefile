WITH_MKL=true
WITH_MKL_SBLAS=true
WITH_ONEAPI=true

COMMON_CXX_SRC = main.cpp SpMatrix.cpp LinearSolverC.cpp
COMMON_CXX_SRCDIR = ../standalone
COMMON_C_SRC = Matrix.c
COMMON_C_SRCDIR = ../../common

CXXFLAGS = -g -D_STANDALONE
CFLAGS = -g -D_STANDALONE

include ../standalone/Make.icc
#include ../standalone/Make.gcc

#######
SOLVERDIR = .

CMDS = pardiso
SOLVER_SRC = pardiso.c
SOLVER_OBJ = pardiso.o

all: $(CMDS)

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

pardiso: $(COMMON_OBJS) pardiso.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
