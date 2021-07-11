#!/bin/sh

DATADIR=${HOME}/work/solverBench

export VE_LD_LIBRARY_PATH=~/Hybrid/vesolver_dev/openmp/dists/ve/lib:${VE_LD_LIBRARY_PATH}

./solvertest ${DATADIR}/$1
