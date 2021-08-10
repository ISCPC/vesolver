WITH_MKL=true
WITH_MKL_SBLAS=true
WITH_ONEAPI=true

CXXFLAGS = -g
CFLAGS = -g

include ../standalone/Make.icc
#include ../standalone/Make.gcc
#include ../standalone/Make.SXAT
#include ../standalone/Make.fcc

#######
SOLVERDIR = .

CMDS = cg pcg
SOLVER_SRC = cg.c pcg.c
SOLVER_OBJ = cg.o pcg.o

all: $(CMDS)

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

cg: $(COMMON_OBJS) cg.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

pcg: $(COMMON_OBJS) pcg.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
