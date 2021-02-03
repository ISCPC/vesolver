#!/bin/bash

host_threads=1
host_procs=4
solver_threads=8
solver_procs=1
prog=./test3
RUN_ON_VE=1

if [ $# -ge 1 ]; then
    host_procs=$1    
    solver_procs=$1    
fi 

if [ $# -ge 2 ]; then
    solver_threads=$2    
fi 

if [ $# -ge 3 ]; then
    prog=$3
fi 

VESPATH=..
VESOLVER=${VESPATH}/vesolver/vesolver
export LD_LIBRARY_PATH=${VESPATH}/lib:${LD_LIBRARY_PATH}

ELMERPATH=../dist/ve

if [ $RUN_ON_VE -eq 1 ]; then
    HOST_OPTS="-np ${host_procs} -vh -env OMP_NUM_THREADS ${host_threads}"
    SOLVER_OPTS="-np ${solver_procs} -ve 3 -env OMP_NUM_THREADS ${solver_threads}"
    export VE_LD_LIBRARY_PATH=${VESPATH}/vesolver:${ELMERPATH}/lib/elmersolver:${VE_LD_LIBRARY_PATH}
else
    HOST_OPTS="-np ${host_procs} -x OMP_NUM_THREADS=${host_threads}"
    SOLVER_OPTS="-np ${solver_procs} -x OMP_NUM_THREADS=${solver_threads} -x MKL_NUM_THREADS=${solver_threads}"
    export LD_LIBRARY_PATH=${VESPATH}/vesolver:${ELMERPATH}/lib/elmersolver:${LD_LIBRARY_PATH}
fi

mpirun ${OPTS} ${HOST_OPTS} ${prog} : ${SOLVER_OPTS} ${VESOLVER}
