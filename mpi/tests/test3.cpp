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

int main(int argc, char** argv) {
    VESolverAPI vesolver;
    MPI_Comm comm;
    int rank, size, err;
    char infile[256];
    SpMatrix A;
    Vector b, x;
    const char *prefix_a="a";
    const char *prefix_b="b";
    double res = 1.0e-13;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &comm);
    MPI_Comm_size(comm, &size);

    if (size != 1) {
        printf("ERROR: The size of test3 must be 1.\n");
        MPI_Finalize();
        return -1;
    }

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
    snprintf(infile, 256, "%s.bin", prefix_a);
    if (A.load_file(infile) < 0) {
        fprintf(stderr, "ERROR: cannot load matrix A from %s\n", infile);
        exit(1);
    }

    /* Load Vector b */
    snprintf(infile, 256, "%s.bin", prefix_b);
    if (b.load_file(infile) < 0) {
        fprintf(stderr, "ERROR: cannot load vector b from %s\n", infile);
        exit(1);
    }

    x.alloc(A.ncol);

    vesolver.activate(comm, 1);
    vesolver.solve(VESOLVER_BICGSTAB2, A, b, x, res);
    //vesolver.solve(VESOLVER_HS, A, b, x, res);
    vesolver.deactivate();
    vesolver.fini();

    /* Print the solution vector */
    printf("%s\n", "******** Solution ********");
    for (int i=0; i<5; i++) {
        printf("x[%d] = %14.12f\n", i, x.value[i]);
    }
    printf("%s\n", "********** End ***********");

    MPI_Finalize();
}
