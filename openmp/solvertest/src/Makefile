CC = ncc
CFLAGS = -O3 -fopenmp -I. -I../../../include -D_TIMELOG

CXX = nc++
CXXFLAGS = -O3 -fopenmp -I. -I../../../include -D_TIMELOG

LDFLAGS = -fopenmp -L../../../velib -lsolver -lcblas -lblas_openmp -lsblas_openmp
LDFLAGS += -lheterosolver_openmp -lblas_sequential

CSRCS = 
CXXSRCS = SpMatrix.cpp solvertest.cpp

.c.o :
	$(CC) $(CFLAGS) -c $<
.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

COBJS = $(CSRCS:.c=.o)
CXXOBJS = $(CXXSRCS:.cpp=.o)
FOBJS = $(FSRCS:.c=.o)
OBJS = $(COBJS) $(CXXOBJS) $(FOBJS)

solvertest: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

clean:
	-rm solvertest

