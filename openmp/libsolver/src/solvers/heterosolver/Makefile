WITH_MKL=true
WITH_MKL_SBLAS=true

COMMON_CXX_SRC = main.cpp SpMatrix.cpp LinearSolverC.cpp
COMMON_CXX_SRCDIR = ../standalone
COMMON_C_SRC = Matrix.c
COMMON_C_SRCDIR = ../../common

CXXFLAGS = -g -D_STANDALONE
CFLAGS = -g -D_STANDALONE
LDFLAGS = -lheterosolver_openmp -lblas_sequential -fopenmp

include ../standalone/Make.SXAT

#######
SOLVERDIR = .

CMDS = hs_c
SOLVER_SRC = hs.c
SOLVER_OBJ = hs.o

all: $(CMDS)

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

hs_c: $(COMMON_OBJS) hs.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
