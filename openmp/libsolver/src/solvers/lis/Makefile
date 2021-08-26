WITH_MKL=true
WITH_MKL_SBLAS=true
WITH_ONEAPI=true

#LIS_ROOT = /home/u/u906/local
LIS_ROOT = ../../../../../../dist

CXXFLAGS = -g -I$(LIS_ROOT)/include
#CXXFLAGS += -DSORT_CSR
CFLAGS = -g -I$(LIS_ROOT)/include
LDFLAGS = -L$(LIS_ROOT)/lib -llis -nofor-main -lstdc++

BUILD_ENV = gcc

include ../standalone/Make.$(BUILD_ENV)

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
