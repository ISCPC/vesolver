#
# Makefile for solverBench
#
WITH_AURORA=true

OPTLEVEL = -O3
DBGFLAGS = -g

#####################################################
# In General, You don't have to touch blow.
#####################################################
CXXFLAGS = $(DBGFLAGS) $(OPTLEVEL) -Wall -std=c++11 -D_TIMELOG
CFLAGS = $(DBGFLAGS) $(OPTLEVEL) -Wall -std=c99 -D_TIMELOG
FFLAGS = $(DBGFLAGS) $(OPTLEVEL) -Wall -cpp
LDFLAGS =

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
INCLUDEDIR = ../include
LIBDIR = ../lib
CFLAGS += -I$(SRCDIR) -I$(INCLUDEDIR)
CXXFLAGS += -I$(SRCDIR) -I$(INCLUDEDIR)
FFLAGS += -J.
LDFLAGS += -L$(LIBDIR) -lvesolver

###################
CXXCMDS = test2 test3 test4
FCMDS = test1f test2f

all: $(CXXCMDS) $(FCMDS)

CXXOBJS = timelog.o SpMatrix.o
FOBJS = VESolverAPI.o

timelog.o: $(SRCDIR)/timelog.c
	$(CC) $(CFLAGS) -c $<

SpMatrix.o: $(SRCDIR)/SpMatrix.cpp
	$(CXX) $(CXXFLAGS) -c $<

VESolverAPI.o: $(SRCDIR)/VESolverAPI.F90
	$(FC) $(FFLAGS) -c $<

$(CXXCMDS): % : %.cpp $(CCXOBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(FCMDS): % : %.F90 $(FOBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o *.mod $(CXXCMDS) $(FCMDS) 2> /dev/null || true
