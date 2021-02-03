/*
MIT License

Copyright (c) 2021 Shunji Uno <shunji_uno@iscpc.jp>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "vesolver_api.h"

int main(int argc, char** argv) {
    MPI_Comm comm;
    int rank, err;
    double *A;
    int nz;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &comm);

    err = vesolver_init(comm);
    if (err) {
        printf("ERROR: vesolver_init failed.\n");
        MPI_Finalize();
        return -1;
    }

    switch(rank) {
        case 0:
            nz = 10;
            A = (double*)calloc(nz, sizeof(double));
            break;

        case 1:
            nz = 6;
            A = (double*)calloc(nz, sizeof(double));
            break;

        default:
            nz = 8;
            A = (double*)calloc(nz, sizeof(double));
            break;
    }
    for(int i=0; i<nz; i++) {
        A[i] = ((double)(rank+1))/10.0;
    }

    vesolver_activate(comm, 1);
    //vesolver_solve(A, nz);
    vesolver_psolve(A, nz, VES_MODE_GATHER_ON_VE);
    vesolver_deactivate();
    vesolver_fini();

    MPI_Finalize();
}
