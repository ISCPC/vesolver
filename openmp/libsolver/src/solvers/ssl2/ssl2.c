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
#include <strings.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include "Matrix.h"
#include "PluginAPI.h"
#include "timelog.h"

typedef struct ssl2_info {
    int32_t nw;
    double* COEF;
    int32_t* ICOL;
    double* factor;
} ssl2_info_t;


void ssl2_vcge_(int32_t* k, int32_t* nw, int32_t* n, 
    double* a, int32_t* icol, const double* b, double* x, double* eps,
    int32_t* ierr);

void ssl2_vbcse_(int32_t* k, int32_t* nw, int32_t* n, 
    double* a, int32_t* icol, const double* b, double* x, double* eps,
    int32_t* ierr);

static int csr2ellpack_sym(Matrix_t* A) {
    ssl2_info_t* info = (ssl2_info_t*)malloc(sizeof(ssl2_info_t));
    int32_t neq = A->NROWS;
    int32_t* iaptr = A->pointers;
    int32_t* iaind = A->indice;
    double* aval = A->values;
    int ncol=0;

    printf("*** 1st step: Normalize ***\n");
    info->factor = (double*)calloc(sizeof(double), neq);
    if (info->factor == NULL) {
        return -1;
    }

    /*  extract diagonal vector from matrix A  */
    for (int i=0; i<neq; i++) {
        info->factor[i] = 1.0/sqrt(aval[iaptr[i]]);
    }

    /*  scale matrix A  */
    for (int i=0; i<neq; i++) {
        for (int j=iaptr[i]; j<iaptr[i+1]; j++) {
            aval[j] *= (info->factor[i]) * (info->factor[iaind[j]]);
        }
    }

    printf("*** 2nd step: Analyze ***\n");
    int* nrows = (int*)calloc(sizeof(int), neq);
    for(int i=0; i<neq; i++) {
        nrows[i] = 0;
    }
    for(int i=0; i<neq; i++) {
        int n = iaptr[i+1]-iaptr[i]-1;
        //printf("num[%d] = %d\n", i, n);
        if (n > ncol) ncol = n;
        for(int j=iaptr[i]+1; j<iaptr[i+1]; j++) {
            nrows[iaind[j]]++;
        }
    }
    int maxrow=0;
    for(int i=0; i<neq; i++) {
        if (nrows[i] > maxrow) maxrow = nrows[i];
    }
    printf("ncol = %d, maxrow = %d\n", ncol, maxrow);
    if (maxrow > ncol) ncol = maxrow;
    free(nrows);
    int32_t nw = info->nw = ncol*2;

    printf("*** 3rd step: Fill matrix *** \n");
    info->COEF = (double*)calloc(sizeof(double), neq*nw);
    info->ICOL = (int32_t*)calloc(sizeof(int32_t), neq*nw);

    // Fill matrix
    int* ku = (int*)calloc(sizeof(int), neq);
    int* kl = (int*)calloc(sizeof(int), neq);
    for(int i=0; i<neq; i++) {
        ku[i] = 0;
        kl[i] = 0;
        for(int j=0; j<nw; j++) {
            info->COEF[j*neq+i] = 0.0f;
            info->ICOL[j*neq+i] = i+1;
        }
    }

    for(int i=0; i<neq; i++) {
        for(int j=iaptr[i]+1; j<iaptr[i+1]; j++) {
            // Fill Upper Matrix
            // COEFF(i, ku[i]) = aval[j]
            // ICOL(i, ku[i]) = iaind[j]
            int k = ku[i]*neq + i;
            info->COEF[k] = aval[j];
            info->ICOL[k] = iaind[j]+1; // one-based index
            ku[i]++;

            // Fill Lower Matrix
            // COEFF(iaind[j], kl[iaind[j]]) = aval[j]
            // ICOL(iaind[j], kl[iaind[j]]) = i
            k = (kl[iaind[j]]+ncol)*neq + iaind[j];
            info->COEF[k] = aval[j];
            info->ICOL[k] = i+1; // one-based index
            kl[iaind[j]]++;
        }
    }
    free(ku);
    free(kl);

    A->info = (void*)info;
    return 0;
}

static int csr2ellpack_asym(Matrix_t* A) {
    ssl2_info_t* info = (ssl2_info_t*)malloc(sizeof(ssl2_info_t));
    int32_t neq = A->NROWS;
    int32_t* iaptr = A->pointers;
    int32_t* iaind = A->indice;
    double* aval = A->values;
    int32_t nw=0;
    
    /* Convert Asymmetric CSR Matrix to ELLPACK(FORTRAN) */
    printf("*** step 1: Analyze. ***\n");
    for(int i=0; i<neq; i++) {
        int n = iaptr[i+1]-iaptr[i];
        if (n > nw) nw = n;
    }
    info->nw = nw;
    printf("nw = %d\n", info->nw);

    printf("*** step 2: Fill matrix *** \n");
    info->COEF = (double*)calloc(sizeof(double), neq*nw);
    info->ICOL = (int32_t*)calloc(sizeof(int32_t), neq*nw);
    info->factor = NULL;

    // Fill matrix
    int* ku = (int*)calloc(sizeof(int), neq);
    for(int i=0; i<neq; i++) {
        ku[i] = 0;
        for(int j=0; j<nw; j++) {
            info->COEF[j*neq+i] = 0.0f;
            info->ICOL[j*neq+i] = i+1;
        }
    }

    for(int i=0; i<neq; i++) {
        for(int j=iaptr[i]; j<iaptr[i+1]; j++) {
            // Fill Matrix
            // COEFF(i, ku[i]) = aval[j]
            // ICOL(i, ku[i]) = iaind[j]
            int k = ku[i]*neq + i;
            info->COEF[k] = aval[j];
            info->ICOL[k] = iaind[j]+1; // one-based index
            ku[i]++;
        }
    }
    free(ku);

    A->info = (void*)info;
    return 0;
}

/*
 * API
 */
static int solve_pre(Matrix_t* A) {
    Matrix_convert_index(A, 0);

    if (MATRIX_IS_SYMMETRIC(A)) {
        csr2ellpack_sym(A);
    } else {
        csr2ellpack_asym(A);
    }
    return 0;
}

static int solve(Matrix_t *A, const double* b, double* x, const double tolerance) {
    int32_t neq = A->NROWS;
    ssl2_info_t* info = (ssl2_info_t*)(A->info);
    int32_t ierr;
    double eps = tolerance;

    printf("*** Call Solver ****\n");
    if (MATRIX_IS_SYMMETRIC(A)) {
        double *b1 = (double*)calloc(sizeof(double), neq);

        for (int i=0; i<neq; i++) {
            b1[i] = b[i] * info->factor[i];
        }
        ssl2_vcge_(&neq, &info->nw, &neq, info->COEF, info->ICOL, b1, x, &eps, &ierr);
        for (int i=0; i<neq; i++) {
            x[i] *= info->factor[i];
        }

        free(b1);
    } else {
        ssl2_vbcse_(&neq, &info->nw, &neq, info->COEF, info->ICOL, b, x, &eps, &ierr);
        //ssl2_vtfqe_(&neq, &info->nw, &neq, info->COEF, info->ICOL, b, x, &eps, &ierr);
    }

    return ierr;
}


static int solve_post(Matrix_t* A) {
    ssl2_info_t* info = (ssl2_info_t*)(A->info);

    if (info) {
        if(info->COEF) free(info->COEF);
        if(info->ICOL) free(info->ICOL);
        if(info->factor) free(info->factor);
        free(A->info);
        A->info = NULL;
    }

    //Matrix_free(A);
    return 0;
}


static void solve_free(SolverPlugin_t* solver) {
	return;
}

/*
 * Solver Plugin Interface
 */
#ifdef _STANDALONE
SolverPlugin_t* solver_init() {
#else
SolverPlugin_t* ssl2_init() {
#endif
	SolverPlugin_t* solver = (SolverPlugin_t*)malloc(sizeof(SolverPlugin_t));

	solver->set_option = NULL;
	solver->solve_pre = solve_pre;
	solver->solve = solve;
	solver->solve_post = solve_post;
	solver->free = solve_free;

	return solver;
}
