#!/bin/sh

DISTSDIR=${HOME}/Hybrid/vesolver_dev/openmp/dists
#DISTSDIR=${HOME}/local

DATADIR=${HOME}/work/solverBench

#export VE_NODE_NUMBER=3
export SOLVER_TYPE="ITER_CG"
#export SOLVER_TYPE="ITER_BICGSTAB2"

export VE_LD_LIBRARY_PATH=${DISTSDIR}/ve/lib:${VE_LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${DISTSDIR}/lib:${LD_LIBRARY_PATH}

#export CCX_VEO_LIBRARY_PATH=${DISTSDIR}/ve/lib/libvesolver.so
export VESOLVER_PATH=${DISTSDIR}/ve/lib/libvesolver.so

./vesolvertest ${DATADIR}/$1
