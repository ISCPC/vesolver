#
# Makefile for solverBench
#
WITH_AURORA=false

OPTLEVEL = -O3
DBGFLAGS = -g

#####################################################
# In General, You don't have to touch blow.
#####################################################
CXXFLAGS = $(DBGFLAGS) $(OPTLEVEL) -Wall -std=c++11 -fPIC -fopenmp -D_TIMELOG
CFLAGS = $(DBGFLAGS) $(OPTLEVEL) -Wall -std=c99 -fPIC -fopenmp -D_TIMELOG

ifeq ($(WITH_AURORA),true)
CXX=mpinc++
CC=mpincc
FC=mpinfort
CXXFLAGS += -vh
CFLAGS += -vh
FFLAGS += -vh
else
CXX=mpic++
CC=mpicc
FC=mpifort
endif

SRCDIR = ../src
CFLAGS += -I$(SRCDIR)
CXXFLAGS += -I$(SRCDIR)
CXXSRCS = vesolver_api.cpp veserver_api.cpp SpMatrix.cpp
CSRCS = timelog.c

###################
LIBS = libvesolver.so

all: $(LIBS)

#.cpp.o :
%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $<

#.c.o :
%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $<

OBJS = $(CXXSRCS:.cpp=.o)
OBJS += $(CSRCS:.c=.o)

libvesolver.so: $(OBJS)
	$(CXX) $(CXXFLAGS) -shared -o $@ $^

clean:
	rm *.o $(LIBS) 2> /dev/null || true
