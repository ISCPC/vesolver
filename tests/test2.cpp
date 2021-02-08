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
#include <string.h>
#include <mpi.h>
#include "SpMatrix.hpp"
#include "vesolver_api.h"

int main(int argc, char** argv) {
    VESolverAPI vesolver;
    MPI_Comm App, Step;
    int rank, size, err;
    SpDistMatrix A;
    double *b=NULL, *x=NULL;
    double res = 1.0e-13;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &App);
    MPI_Comm_size(App, &size);

    err = vesolver_init(App);
    if (err) {
        printf("ERROR: vesolver_init failed.\n");
        MPI_Finalize();
        return -1;
    }

    /*
     * create subgroup
     */
    MPI_Group MPI_GROUP_WORLD, mgroup;
    MPI_Comm_group( MPI_COMM_WORLD, &MPI_GROUP_WORLD );

    int mranks[2] = { 2, 3 };
    MPI_Group_incl( MPI_GROUP_WORLD, 2, mranks, &mgroup );
    MPI_Comm_create( App, mgroup, &Step );

    /* Input Preparation */
    A.nrow = A.ncol = A.neq = 7;
    if (rank == 2) {
        /* The row numbers 1, 2, 3, and 4 are assigned to the MPI process 0. */
        A.mingind=1;
        A.maxgind=4;
        A.nrow=A.maxgind-A.mingind+1;
        A.value  = (double*)   malloc(sizeof(double)*8);
        A.indice = (int*) malloc(sizeof(int)*8);
        A.pointers = (int*) malloc(sizeof(int)*(A.maxgind-A.mingind+2));
        A.rorder = (int*) malloc(sizeof(int)*7);
        b = (double*)   malloc(sizeof(double)*(A.maxgind-A.mingind+1));
        x = (double*)   malloc(sizeof(double)*7);
        double   aval_local[]  = { 1.1, 1.2, 2.2, 2.3, 3.3, 3.4, 4.1, 4.4 };
        int iaind_local[] = { 1, 2, 2, 3, 3, 4, 1, 4 };
        int iaptr_local[] = { 1, 3, 5, 7, 9 };
        int rorder_local[] = { 1, 2, 3, 4, 0, 0, 0 };
        double   b_local[] = { 3.5, 11.3, 33.9, 10.0 };
        memcpy (A.value,  aval_local,  sizeof(aval_local) );
        memcpy (A.indice, iaind_local, sizeof(iaind_local));
        memcpy (A.pointers, iaptr_local, sizeof(iaptr_local));
        memcpy (A.rorder, rorder_local, sizeof(rorder_local));
        memcpy (b,  b_local,     sizeof(b_local)    );
    } else if(rank == 3) {
        /* The row numbers 3, 4, 5, 6 and 7 are assigned to the MPI process 1. */
        A.mingind=3;
        A.maxgind=7;
        A.nrow=A.maxgind-A.mingind+1;
        A.value  = (double*)   malloc(sizeof(double)*11);
        A.indice = (int*) malloc(sizeof(int)*11);
        A.pointers = (int*) malloc(sizeof(int)*(A.maxgind-A.mingind+2));
        A.rorder = (int*) malloc(sizeof(int)*7);
        b = (double*)   malloc(sizeof(double)*(A.maxgind-A.mingind+1));
        x = (double*)   malloc(sizeof(double)*7);
        double   aval_local[]  = { 3.5, 3.7, 5.3, 5.5, 5.6, 6.5, 6.6, 6.7, 7.3, 7.6, 7.7 };
        int iaind_local[] = { 5, 7, 3, 5, 6, 5, 6, 7, 3, 6, 7 };
        int iaptr_local[] = { 1, 3, 3, 6, 9, 12 };
        int rorder_local[] = { 0, 0, 1, 2, 3, 4, 5 };
        double   b_local[] = { 33.0, 11.7, 77.0, 119.0, 121.4 };
        memcpy (A.value,  aval_local,  sizeof(aval_local) );
        memcpy (A.indice, iaind_local, sizeof(iaind_local));
        memcpy (A.pointers, iaptr_local, sizeof(iaptr_local));
        memcpy (A.rorder, rorder_local, sizeof(rorder_local));
        memcpy (b,  b_local,     sizeof(b_local)    );
    }


    if (Step != MPI_COMM_NULL) {
        for (int i=0; i<10; i++) {
            vesolver_activate(Step, 2);
            vesolver_psolve2(VES_MODE_SYMMETRIC, VESOLVER_HS, 
                A.neq, A.nrow, A.mingind, A.maxgind, A.pointers, A.indice, A.value, A.rorder,
                b, x, res);
            //vesolver_psolve_dcsr(VES_MODE_SYMMETRIC, VESOLVER_HS,
            //    A.neq, A.pointers, A.indice, A.value, A.mingind, A.maxgind, b, x, res);
            vesolver_deactivate();

            if (rank == 2) {
                /* Print the solution vector */
                printf("%s\n", "******** Solution ********");
                for (int i=0; i<7; i++) {
                    printf("x[%d] = %14.12f\n", i, x[i]);
                }
                printf("%s\n", "********** End ***********");
            }
        }

        MPI_Comm_free(&Step);
    }

    MPI_Barrier(App);
    vesolver_fini();

    MPI_Finalize();
}
