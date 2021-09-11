WITH_MKL=true
WITH_MKL_SBLAS=true
WITH_ONEAPI=true

LIS_ROOT = ../../../../../lib
#LIS_ROOT = $(HOME)/local

CXXFLAGS = -g -I$(LIS_ROOT)/include
#CXXFLAGS += -DSORT_CSR
CFLAGS = -g -I$(LIS_ROOT)/include
#CFLAGS += -DUSE_SXAURORA_SBLAS
#LDFLAGS = -L$(LIS_ROOT)/lib -llis -nofor-main -lstdc++
LDFLAGS = -L$(LIS_ROOT)/lib -llis

BUILD_ENV = gcc

include ../standalone/Make.$(BUILD_ENV)

#######
SOLVERDIR = .

CMDS = lis_csr lis_csc lis_ell lis_jad
SOLVER_SRC = lis_csr.c lis_csc.c lis_ell.c lis_jad.c
SOLVER_OBJ = lis_csr.o lis_csc.o lis_ell.o lis_jad.o

all: $(CMDS)

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

lis_csr: $(COMMON_OBJS) lis_csr.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

lis_csc: $(COMMON_OBJS) lis_csc.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

lis_ell: $(COMMON_OBJS) lis_ell.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

lis_jad: $(COMMON_OBJS) lis_jad.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
