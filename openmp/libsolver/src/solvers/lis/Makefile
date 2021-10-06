WITH_MKL=true
WITH_MKL_SBLAS=true
WITH_ONEAPI=true

LIS_ROOT = ../../../../../lib
#LIS_ROOT = $(HOME)/local
#LIS_ROOT = $(HOME)/local/ve

CXXFLAGS = -g -I$(LIS_ROOT)/include
#CXXFLAGS += -DSORT_CSR
CFLAGS = -g -I$(LIS_ROOT)/include
#LDFLAGS = -L$(LIS_ROOT)/lib -llis -nofor-main -lstdc++
LDFLAGS = -L$(LIS_ROOT)/lib -llis

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
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
