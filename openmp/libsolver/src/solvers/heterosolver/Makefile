WITH_MKL=true
WITH_MKL_SBLAS=true

CXXFLAGS = -g
CFLAGS = -g
LDFLAGS = -lheterosolver_openmp -lblas_sequential -fopenmp

BUILD_ENV = SXAT

include ../standalone/Make.$(BUILD_ENV)

#######
SOLVERDIR = .

CMDS = hs
SOLVER_SRC = hs.c
SOLVER_OBJ = hs.o

all: $(CMDS)

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

hs: $(COMMON_OBJS) hs.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
