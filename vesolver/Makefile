#
# Makefile for solverBench
#
WITH_AURORA=false

#
# Enable/Disable PARDISO
#
WITH_MKL=false
WITH_MKL_STATIC=true

WITH_MPI=false

OPTLEVEL = -O3
DBGFLAGS = -g
#DBGFLAGS = -g -DDISABLE_INLINE

#####################################################
# In General, You don't have to touch blow.
#####################################################
CXXFLAGS = $(DBGFLAGS) -Wall $(OPTLEVEL) -fopenmp -D_TIMELOG
CFLAGS = $(DBGFLAGS) $(OPTLEVEL) -Wall -O3 -std=c99 -fopenmp -D_TIMELOG
FFLAGS = $(DBGFLAGS) $(OPTLEVEL) -fopenmp
LDFLAGS = $(DBGFLAGS) -fopenmp

ifeq ($(WITH_AURORA),true)
CXX=mpinc++
CC=mpincc
FC=mpinfort
CXXFLAGS += -ve
CFLAGS += -ve
FFLAGS += -ve
else
CXX=mpic++
CC=mpicc
FC=mpifort
endif

SRCDIR = ../src
DISTDIR = ../dist/ve
CFLAGS += -I$(SRCDIR)
#CSRCS = $(SRCDIR)/common.c $(SRCDIR)/solver.c $(SRCDIR)/vesolver.c
#CXXSRCS = $(SRCDIR)/veserver.cpp $(SRCDIR)/vesolver.cpp
CXXSRCS = veserver.cpp vesolver.cpp itersolver.cpp

###################
# Common Libraries
###################
#COMMON_SRCDIR=../../common
COMMON_SRCDIR=.
CFLAGS += -I $(COMMON_SRCDIR)

CSRCS += timelog.c
CXXSRCS += SpMatrix.cpp

###################
# for HeteroSolver
###################
ifeq ($(WITH_AURORA),true)
CFLAGS += -DHETEROSOLVER -DSXAT
CXXFLAGS += -DHETEROSOLVER -DSXAT
LDFLAGS += -lheterosolver_mpi_openmp -lblas_sequential
LDFLAGS += -lsblas_openmp
LDFLAGS += -lcblas
endif

###################
# for Pardiso
###################
ifeq ($(WITH_AURORA),false)
ifeq ($(WITH_MKL),true)
CFLAGS += -DWITH_PARDISO -DMKL_SBLAS
CXXFLAGS += -DWITH_PARDISO -DMKL_SBLAS

INTEL_LIB = /opt/intel/lib/intel64
MKL_LIB = /opt/intel/mkl/lib/intel64
MKL_INCLUDE = /opt/intel/mkl/include
ifeq ($(WITH_MKL_STATIC),true)
	CMPLR_PATH = /opt/intel/compilers_and_libraries/linux/lib/intel64
else
	LDFLAGS += -L$(INTEL_LIB) -L$(MKL_LIB)
	LDFLAGS += -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_intel_thread
	LDFLAGS += -lmkl_core -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64
	LDFLAGS += -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_sequential
	LDFLAGS += -lmkl_intel_lp64 -liomp5 -lmkl_rt
#LDFLAGS += -llapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
endif

#LDFLAGS += -L$(INTEL_LIB) -L$(MKL_LIB)
#LDFLAGS += -llapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

CFLAGS += -I $(MKL_INCLUDE)
CXXFLAGS += -I $(MKL_INCLUDE)
FFLAGS += -I $(MKL_INCLUDE)
else
LDFLAGS += -lopenblas
endif # WITH_MKL
endif # WITH_AURORA

###################
CMDS = vesolver

all: $(CMDS)

#.cpp.o :
%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $<

#.c.o :
%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $<

OBJS = $(CXXSRCS:.cpp=.o)
OBJS += $(CSRCS:.c=.o)

libelmerapi.so: $(SRCDIR)/elmerSolver.c $(SRCDIR)/elmerSolverAPI.f90
	$(CC) $(CFLAGS) -fPIC -c $(SRCDIR)/elmerSolver.c; \
	$(FC) $(FFLAGS) -fPIC -c $(SRCDIR)/elmerSolverAPI.f90; \
	$(FC) $(LDFLAGS) -shared -fPIC -o $@ elmerSolver.o elmerSolverAPI.o

vesolver: $(OBJS) $(LIBS)
ifeq ($(WITH_AURORA),true)
	mpinfort $(FCFLAGS) $(LDOPT) $(LDFLAGS) -cxxlib -o $@ $^
else
ifeq ($(WITH_MKL),true)
ifeq ($(WITH_MKL_STATIC),true)
	$(CXX) $(CXXFLAGS) -o $@ $^ \
	-Wl,--start-group \
	$(INTEL_LIB)/libiomp5.a \
	$(MKL_LIB)/libmkl_gf_lp64.a \
	$(MKL_LIB)/libmkl_gnu_thread.a \
	$(MKL_LIB)/libmkl_intel_thread.a \
	$(MKL_LIB)/libmkl_core.a \
	$(MKL_LIB)/libmkl_scalapack_lp64.a \
	$(MKL_LIB)/libmkl_blacs_openmpi_lp64.a \
	$(MKL_LIB)/libmkl_lapack95_lp64.a \
	$(MKL_LIB)/libmkl_blas95_lp64.a \
	$(MKL_LIB)/libmkl_sequential.a \
	$(MKL_LIB)/libmkl_intel_lp64.a \
	-Wl,--end-group \
	-L$(CMPLR_PATH) -lpthread -lm -ldl -lc $(LDFLAGS)
else
	$(CXX) $(CXXFLAGS) -o $@ $(LDFLAGS) $^ \
	-lpthread -lm -ldl -lc
endif
else
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
endif
endif

dist:
	cp $(CMDS) $(DISTDIR)/bin; \
	cp libelmerapi.so $(DISTDIR)/lib

clean:
	rm *.o $(CMDS) $(LIBS) 2> /dev/null || true
