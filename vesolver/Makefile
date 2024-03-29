#
# Makefile for solverBench
#
WITH_AURORA=true

#
# Enable/Disable ElmerFEM solver support
#
WITH_ELMERSOLVER=true

#
# Enable/Disable PARDISO
#
WITH_MKL=true

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
CXXSRCS = veserver.cpp vesolver.cpp

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
CFLAGS += -DHETEROSOLVER
CXXFLAGS += -DHETEROSOLVER
LDFLAGS += -lheterosolver_mpi_openmp -lblas_sequential
LDFLAGS += -lsblas_openmp
endif

###################
# for Pardiso
###################
ifeq ($(WITH_MKL),true)
ifeq ($(WITH_AURORA),false)
CFLAGS += -DWITH_PARDISO -DMKL_SBLAS
CXXFLAGS += -DWITH_PARDISO -DMKL_SBLAS

INTEL_LIB = /opt/intel/lib/intel64
MKL_LIB = /opt/intel/mkl/lib/intel64
MKL_INCLUDE = /opt/intel/mkl/include
CMPLR_PATH = /opt/intel/compilers_and_libraries/linux/lib/intel64

#LDFLAGS += -L$(INTEL_LIB) -L$(MKL_LIB)
#LDFLAGS += -llapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

CFLAGS += -I $(MKL_INCLUDE)
CXXFLAGS += -I $(MKL_INCLUDE)
FFLAGS += -I $(MKL_INCLUDE)
endif # WITH_AURORA
endif # WITH_MKL

###################
# for ElmerSolver
###################
ifeq ($(WITH_ELMERSOLVER),true)
CSRCS += elmerSolver.c

CFLAGS += -DELMERSOLVER
CXXFLAGS += -DELMERSOLVER

FSRCS += elmerSolverAPI.f90
LIBS += libelmerapi.so

ELMER_LIB_PATH = ../dist/ve/lib/elmersolver
ELMER_MOD_PATH = ../elmer/build/fmodules

LDFLAGS += -L$(ELMER_LIB_PATH) -lelmersolver
LDFLAGS += -lfhuti -lmatc -larpack

ifeq ($(WITH_AURORA),true)
FFLAGS += -module $(ELMER_MOD_PATH)
ifeq ($(WITH_MPI),true)
LDFLAGS += -lparpack -lnc++
else
#LDFLAGS += -lmpi_stubs
endif
LDFLAGS += -lcblas -llapack -lblas_openmp
#
else # WIH_AURORA
FFLAGS += -J $(ELMER_MOD_PATH)
endif # WITH_AURORA

endif # WITH_ELMERSOLVER

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
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
endif
endif

dist:
	cp $(CMDS) $(DISTDIR)/bin; \
	cp libelmerapi.so $(DISTDIR)/lib

clean:
	rm *.o $(CMDS) $(LIBS) 2> /dev/null || true
