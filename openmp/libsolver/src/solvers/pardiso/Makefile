WITH_MKL=true
WITH_MKL_SBLAS=true
WITH_ONEAPI=true

CXXFLAGS = -g
#CXXFLAGS += -DSORT_CSR
CFLAGS = -g

BUILD_ENV = icc

include ../standalone/Make.$(BUILD_ENV)

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
