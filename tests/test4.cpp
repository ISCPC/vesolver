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
#include "SpMatrix.hpp"
#include "vesolver_api.h"
#include "timelog.h"

int main(int argc, char** argv) {
    VESolverAPI vesolver;
    MPI_Comm comm;
    int rank, size, err;
    char infile[256];
    SpDistMatrix A;
    DistVector b, x;
    const char *prefix_a="a";
    const char *prefix_b="b";
    double res = 1.0e-13;
    int nloop=1;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &comm);
    MPI_Comm_size(comm, &size);

    err = vesolver.init(comm);
    if (err) {
        printf("ERROR: vesolver_init failed.\n");
        MPI_Finalize();
        return -1;
    }

    /*
     * read matrix/vector
     */
    /* Load Matrix A */
    snprintf(infile, 256, "%s%d.bin", prefix_a, rank);
    if (A.load_file(infile) < 0) {
        fprintf(stderr, "ERROR: cannot load matrix A from %s\n", infile);
        exit(1);
    }

    /* Load Vector b */
    snprintf(infile, 256, "%s%d.bin", prefix_b, rank);
    if (b.load_file(infile) < 0) {
        fprintf(stderr, "ERROR: cannot load vector b from %s\n", infile);
        exit(1);
    }

    x.alloc(A.neq);

    float tot1_ve=0.0, tot1_vh=0.0;
    float tot2_ve=0.0, tot2_vh=0.0;
    for(int i=0; i<nloop; i++) {
        TIMELOG(tl1);
        TIMELOG(tl2);
        float t1, t2;

        /* Test Gather On VE */
        MPI_Barrier(comm);
        TIMELOG_START(tl1);
        vesolver.activate(comm, 1);
        TIMELOG_START(tl2);
        vesolver.solve(VESOLVER_BICGSTAB2, A, b, x, res, VES_MODE_GATHER_ON_VE);
        //vesolver.solve(VESOLVER_HS, A, b, x, res, VES_MODE_GATHER_ON_VE);
        TIMELOG_GETTIME(t2, tl2);
        vesolver.deactivate();
        MPI_Barrier(comm);
        TIMELOG_GETTIME(t1, tl1);

        if (rank == 0) {
            printf("Loop %d : Gather_On_VE : %f [sec], %f [sec]\n", i, t2, t1);
            tot1_ve += t1;
            tot2_ve += t2;
            fflush(stdout);
        }

        /* Test Gather On VE */
        MPI_Barrier(comm);
        TIMELOG_START(tl1);
        vesolver.activate(comm, 1);
        TIMELOG_START(tl2);
        vesolver.solve(VESOLVER_BICGSTAB2, A, b, x, res, VES_MODE_GATHER_ON_VH);
        //vesolver.solve(VESOLVER_HS, A, b, x, res, VES_MODE_GATHER_ON_VH);
        TIMELOG_GETTIME(t2, tl2);
        vesolver.deactivate();
        MPI_Barrier(comm);
        TIMELOG_GETTIME(t1, tl1);

        if (rank == 0) {
            printf("Loop %d : Gather_On_VH : %f [sec], %f [sec]\n", i, t2, t1);
            tot1_vh += t1;
            tot2_vh += t2;
            fflush(stdout);
        }
    }

    vesolver.fini();
    if (rank == 0) {
        /* Print the solution vector */
        printf("%s\n", "******** Solution ********");
        for (int i=0; i<5; i++) {
            printf("x[%d] = %14.12f\n", i, x.value[i]);
        }
        printf("%s\n", "********** End ***********");

        printf("AVERAGE:Gather_On_VE : %f [sec], %f [sec]\n", tot2_ve/nloop, tot1_ve/nloop);
        printf("AVERAGE:Gather_On_VH : %f [sec], %f [sec]\n", tot2_vh/nloop, tot1_vh/nloop);
    }

    MPI_Finalize();
}
