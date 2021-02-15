#!/bin/bash

nprocs=1
nthreads=8
RUN_ON_VE=0

# Usage: run.sh [solver] [run_mode] {data_path}
#
# 1st parameter: solver
#
case "$1" in
cg)
    export VESOLVER=BICGSTAB2
    ;;

hs)
    export VESOLVER=HS
    ;;
    
dum*)
    export VESOLVER=DUMMY
    ;;

*)
    export VESOLVER=BICGSTAB2
    ;;
esac

#
# get the number of parallels
#
if [ "$3" != "" ]
then
    export VESOLVER_DATA_PATH=$3
else
    export VESOLVER_DATA_PATH=.
fi
nmat=`ls ${VESOLVER_DATA_PATH}/a*.bin | wc -l`

if [ ${nmat} -le 0 ]; then
    echo "ERROR: Matrix data not found"
    exit
fi

#
# 2nd parameter: Rum mode
#
case "$2" in
ve)
    export VES_MODE=GATHER_ON_VE
    export VESOLVER_NPARA=${nmat}
    ;;
    
vh)
    export VES_MODE=GATHER_ON_VH
    ;;

sym*)
    export VES_MODE=SYMMETRIC
    nprocs=${nmat}
    nthreads=1
    ;;

*)
    export VES_MODE=GATHER_ON_VE
    export VESOLVER_NPARA=${nmat}
    ;;
esac


ELMERPATH=../dist/ve/lib/elmersolver

export LD_LIBRARY_PATH=/opt/intel/lib/intel64:/opt/intel/mkl/lib/intel64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${ELMERPATH}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=.:${LD_LIBRARY_PATH}

export OMP_NUM_THREADS=${nthreads}

if [ $RUN_ON_VE -eq 1 ]; then
    #export VE_LD_LIBRARY_PATH=/opt/nec/ve/mpi/2.5.0/lib64/ve:/opt/nec/ve/nlc/2.0.0/lib
    export VE_LD_LIBRARY_PATH=/opt/nec/ve/mpi/2.13.0/lib64/ve:/opt/nec/ve/nlc/2.2.0/lib
    export VE_LD_LIBRARY_PATH=${ELMERPATH}:${VE_LD_LIBRARY_PATH}
    export VE_LD_LIBRARY_PATH=.:${VE_LD_LIBRARY_PATH}

    mpirun -np ${nprocs} -ve 0-3 -env OMP_NUM_THREADS ${nthreads} ./vesolver
else
    if [ $nprocs -eq 1 ]; then
        ./vesolver
        #gdb ./vesolver
    else
        mpirun -np ${nprocs} ./vesolver
    fi
fi
