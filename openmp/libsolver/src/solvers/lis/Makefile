WITH_MKL=true
WITH_MKL_SBLAS=true
WITH_ONEAPI=true

COMMON_CXX_SRC = main.cpp SpMatrix.cpp LinearSolverC.cpp
COMMON_CXX_SRCDIR = ../standalone
COMMON_C_SRC = Matrix.c
COMMON_C_SRCDIR = ../../common

#LIS_ROOT = /home/u/u906/local
LIS_ROOT = /home/u/u906/work/LS-Particle/ruri/2DRC-17_IMPS-SPP-PGC-KONDO_L50mic/Lis-2.0.20

CXXFLAGS = -g -D_STANDALONE -I$(LIS_ROOT)/include
#CXXFLAGS += -DSORT_CSR
CFLAGS = -g -D_STANDALONE -I$(LIS_ROOT)/include
LDFLAGS = -L$(LIS_ROOT)/lib -llis -nofor-main -lstdc++

include ../standalone/Make.icc
#include ../standalone/Make.gcc

#######
SOLVERDIR = .

CMDS = lis
SOLVER_SRC = lis.c
SOLVER_OBJ = lis.o

all: $(CMDS)

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

lis: $(COMMON_OBJS) lis.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
