#!/bin/bash

WITH_VESOLVER=true
WITH_MKL=false
WITH_MPI=true
WITH_OpenMP=true

#INSTALL_PATH=../install
INSTALL_PATH=../vesolver
VESOLVER_LIBRARY_DIR=~/local/lib

#
#
#
### Add Solver Options
OPTIONS="-DWITH_ELMERGUI:BOOL=FALSE"
OPTIONS="${OPTIONS} -DWITH_TIMELOG:BOOL=TRUE"

# If enable VESolver, MPI and OpenMPI must be enabled.
if [ "$WITH_VESOLVER" == "true" ]; then
    OPTIONS="${OPTIONS} -DWITH_OpenMP:BOOL=TRUE"
    OPTIONS="${OPTIONS} -DWITH_MPI:BOOL=TRUE"
    OPTIONS="${OPTIONS} -DCMAKE_TOOLCHAIN_FILE=../elmerfem/cmake/Toolchains/Elmer-linux-gcc-vesolver.cmake"
    OPTIONS="${OPTIONS} -DWITH_VESOLVER:BOOL=TRUE"
    OPTIONS="${OPTIONS} -DVESOLVER_LIBRARIES=${VESOLVER_LIBRARY_DIR}"
else
    if [ "$WITH_OpenMP" == "true" ]; then
        OPTIONS="${OPTIONS} -DWITH_OpenMP:BOOL=TRUE"
        OPTIONS="${OPTIONS} -DCMAKE_TOOLCHAIN_FILE=../elmerfem/cmake/Toolchains/Elmer-linux-gcc-openmp.cmake"
    else
        OPTIONS="${OPTIONS} -DWITH_OpenMP:BOOL=FALSE"
        OPTIONS="${OPTIONS} -DCMAKE_TOOLCHAIN_FILE=../elmerfem/cmake/Toolchains/Elmer-linux-gcc.cmake"
    fi

    if [ "$WITH_MPI" = "true" ]; then
        OPTIONS="${OPTIONS} -DWITH_MPI:BOOL=TRUE"
    else
        OPTIONS="${OPTIONS} -DWITH_MPI:BOOL=FALSE"
    fi
fi


# with MKL
if [ "$WITH_MKL" = "true" ]; then
    export MKLROOT=/opt/intel/mkl
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/lib/intel64:/opt/intel/mkl/lib/intel64:/opt/elmerfem/lib

    OPTIONS="${OPTIONS} -DWITH_MKL:BOOL=TRUE"
    OPTIONS="${OPTIONS} -DWITH_Mumps:BOOL=FALSE"
else
    OPTIONS="${OPTIONS} -DWITH_MKL:BOOL=FALSE"
    OPTIONS="${OPTIONS} -DWITH_Mumps:BOOL=TRUE"
fi

# with Hypre
#export HYPRE_INCLUDE_DIR=/usr/include/hypre
#export HYPRE_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu

#OPTIONS="${OPTIONS} -DWITH_Hypre:BOOL=TRUE -DHypre_INCLUDE_DIR=${HYPRE_INCLUDE_DIR} -DHypre_LIBRARIES=${HYPRE_LIBRARY_DIR}"
#OPTIONS="${OPTIONS} -DWITH_Hypre:BOOL=TRUE"

# with Trilinos
#OPTIONS="${OPTIONS} -DWITH_Trilinos:BOOL=TRUE"

# with feti4i
#export FETI4I_ROOT=/home/uno/elmer/elmerfem/feti4i
#OPTIONS="${OPTIONS} -DWITH_FETI4I:BOOL=TRUE -DFETI4I_ROOT=${FETI4I_ROOT}"


cmake ${OPTIONS} -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} ../elmerfem

if [ $? -eq 0 ]; then
    make -j install
fi
