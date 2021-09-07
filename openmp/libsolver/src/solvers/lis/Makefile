WITH_MKL=true
WITH_MKL_SBLAS=true
WITH_ONEAPI=true

LIS_ROOT = ../../../../../../dist

CXXFLAGS = -g -I$(LIS_ROOT)/include
#CXXFLAGS += -DSORT_CSR
CFLAGS = -g -I$(LIS_ROOT)/include
LDFLAGS = -L$(LIS_ROOT)/lib -llis -nofor-main -lstdc++

BUILD_ENV = gcc

include ../standalone/Make.$(BUILD_ENV)

#######
SOLVERDIR = .

CMDS = lis lis_csc lis_ell lis_jad
SOLVER_SRC = lis.c lis_csc.c lis_ell.c lis_jad.c
SOLVER_OBJ = lis.o lis_csc.o lis_ell.o lis_jad.o

all: $(CMDS)

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

lis: $(COMMON_OBJS) lis.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

lis_csc: $(COMMON_OBJS) lis_csc.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

lis_ell: $(COMMON_OBJS) lis_ell.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

lis_jad: $(COMMON_OBJS) lis_jad.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
