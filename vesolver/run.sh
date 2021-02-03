#!/bin/bash

nprocs=1
RUN_ON_VE=1

# Solver
export VESOLVER=BICGSTAB2
#export VESOLVER=HS
#export VESOLVER=DUMMY

# Parallel Run mode
#export VES_MODE=GATHER_ON_VH
export VES_MODE=GATHER_ON_VE
#export VES_MODE=SYMMETRIC

if [ "$1" != "" ]
then
    nprocs=$1
fi

ELMERPATH=../dist/ve/lib/elmersolver

export LD_LIBRARY_PATH=/opt/intel/lib/intel64:/opt/intel/mkl/lib/intel64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${ELMERPATH}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=.:${LD_LIBRARY_PATH}

if [ $RUN_ON_VE -eq 1 ]; then
    #export VE_LD_LIBRARY_PATH=/opt/nec/ve/mpi/2.5.0/lib64/ve:/opt/nec/ve/nlc/2.0.0/lib
    export VE_LD_LIBRARY_PATH=/opt/nec/ve/mpi/2.13.0/lib64/ve:/opt/nec/ve/nlc/2.2.0/lib
    export VE_LD_LIBRARY_PATH=${ELMERPATH}:${VE_LD_LIBRARY_PATH}
    export VE_LD_LIBRARY_PATH=.:${VE_LD_LIBRARY_PATH}

    mpirun -np ${nprocs} ./vesolver
else
    if [ $nprocs -eq 1 ]; then
        ./vesolver
        #gdb ./vesolver
    else
        mpirun -np ${nprocs} ./vesolver
    fi
fi
