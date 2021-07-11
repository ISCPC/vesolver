DISTSDIR = ../../dists

INCLUDEDIR = $(DISTSDIR)/include
LIBDIR = $(DISTSDIR)/lib

CC = gcc
CFLAGS = -O3 -g -fopenmp -I. -I$(INCLUDEDIR) -D_TIMELOG

CXX = g++
CXXFLAGS = -O3 -g -fopenmp -I. -I$(INCLUDEDIR) -D_TIMELOG

LDFLAGS = -fopenmp
LDFLAGS += -L/opt/nec/ve/veos/lib64 -Wl,-rpath=/opt/nec/ve/veos/lib64 -lveo
LDFLAGS += -fopenmp -L$(LIBDIR) -lvesolver_api 

CSRCS = 
CXXSRCS = SpMatrix.cpp vestest.cpp

.c.o :
	$(CC) $(CFLAGS) -c $<
.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

COBJS = $(CSRCS:.c=.o)
CXXOBJS = $(CXXSRCS:.cpp=.o)
FOBJS = $(FSRCS:.c=.o)
OBJS = $(COBJS) $(CXXOBJS) $(FOBJS)

vesolvertest: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

clean:
	-rm vesolvertest

