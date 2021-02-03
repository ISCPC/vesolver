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
#include "elmerSolver.h"

static int getPreCond() {
    char* str = getenv("SOLVERBENCH_PRECONDITIONING");

    if (!str) {
        return ELMER_PC_NONE;
    }

    if (strcmp(str, "Diagonal") == 0) {
        return ELMER_PC_DIAG;
    } else if (strcmp(str, "ILU0") == 0) {
        return ELMER_PC_ILU0;
    } else if (strcmp(str, "ILUT") == 0) {
        return ELMER_PC_ILUT;
    } else if (strcmp(str, "Multigrid") == 0) {
        return ELMER_PC_MG;
    }
    return ELMER_PC_NONE;
}

int elmersolver(INT_T neq, INT_T nnz, INT_T *pointers, INT_T *indice, double *value,
    double* b, double* x, INT_T solverId)
{
    int preCond = getPreCond();

    if(neq < 0) {
        return -1;
    }

    printf(" Solving the system of equations using the elmer solver\n");
    elmersolverapi_(&neq, &nnz, pointers, indice, value, b, x, &solverId, &preCond);

    return 0;
}

int elmersolver_distributed(INT_T neq, INT_T nnz, INT_T *pointers, INT_T *indice, double *value,
    double* b, double* x, INT_T solverId)
{
    printf("ERROR: elmersolver for MPI is not implemented yet.\n");
    return -1;
}

#if 0
int run_elmersolver(SpMatrix_t* A, double* b, double* x, INT_T solverId) {
    int nthread = 1;
    int preCond = getPreCond();

    if(A->nrow < 0) {
        return -1;
    }

    if (SPMATRIX_IS_SYMMETRIC(A)) {
        printf(" Solving the system of equations using the symmetric elmer solver\n");
    } else {
        printf(" Solving the system of equations using the unsymmetric elmer solver\n");
    }

    /* set MKL_NUM_THREADS */
    char *env = getenv("OMP_NUM_THREADS");
    if (env) {
        nthread = atoi(env);
        if (nthread < 1) {nthread = 1;}
    }
    printf(" number of threads =% d\n\n",nthread);
#if 1
    if (SPMATRIX_IS_SYMMETRIC(A)) {
        printf("INFO: extracting Symmetric Matrix...\n");
        SpMatrix_t* A0 = SpMatrix_extract_symmetric(A);
        elmersolverapi_(&(A0->nrow), &(A0->ndim), A0->pointers, A0->indice, A0->value, b, x, &solverId, &preCond);
        SpMatrix_free(A0);
    } else if (SPMATRIX_IS_CSC(A)) {
        printf("INFO: converting to CSR format...\n");
        SpMatrix_t* A0 = SpMatrix_transpose(A);
        elmersolverapi_(&(A0->nrow), &(A0->ndim), A0->pointers, A0->indice, A0->value, b, x, &solverId, &preCond);
        SpMatrix_free(A0);
    } else {
        elmersolverapi_(&(A->nrow), &(A->ndim), A->pointers, A->indice, A->value, b, x, &solverId, &preCond);
    }
#endif

    return 0;
}
#endif
