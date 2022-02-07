WITH_MKL=true
WITH_MKL_SBLAS=true

COMMON_CXX_SRC = main.cpp SpMatrix.cpp LinearSolverC.cpp
COMMON_CXX_SRCDIR = ../standalone
COMMON_C_SRC = Matrix.c
COMMON_C_SRCDIR = ../../common

CXXFLAGS = -g -D_STANDALONE
#CXXFLAGS += -DSORT_CSR
CFLAGS = -g -D_STANDALONE
FFLAGS = -g -D_STANDALONE
LDFLAGS = -mlcmain=main -lstdc++

BUILD_ENV = fccpx

include ../standalone/Make.$(BUILD_ENV)

#######
SOLVERDIR = .

CMDS = ssl2
SOLVER_SRC = ssl2.c
SOLVER_OBJ = ssl2.o

all: $(CMDS)

solver.o: solver.F90
	$(FC) $(FFLAGS) -c $^

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

ssl2: $(COMMON_OBJS) ssl2.o solver.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
