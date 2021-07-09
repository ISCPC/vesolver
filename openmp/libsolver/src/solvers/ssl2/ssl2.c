#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include "Matrix.h"
#include "timelog.h"

typedef struct ellpack_info {
    int32_t nw;
    double* COEF;
    int32_t* ICOL;
    double* factor;
} ellpack_info_t;


void ssl2_vcge_(int32_t* k, int32_t* nw, int32_t* n, 
    double* a, int32_t* icol, const double* b, double* x, double* eps,
    int32_t* ierr);

void ssl2_vbcse_(int32_t* k, int32_t* nw, int32_t* n, 
    double* a, int32_t* icol, const double* b, double* x, double* eps,
    int32_t* ierr);

static int csr2ellpack_sym(Matrix_t* A) {
    ellpack_info_t* info = (ellpack_info_t*)malloc(sizeof(ellpack_info_t));
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
    ellpack_info_t* info = (ellpack_info_t*)malloc(sizeof(ellpack_info_t));
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

Matrix_t* PreProcess(const Matrix_t* A0) {
    Matrix_t* A = Matrix_duplicate(A0);
    Matrix_convert_index(A, 0);

    if (MATRIX_IS_SYMMETRIC(A)) {
        csr2ellpack_sym(A);
    } else {
        csr2ellpack_asym(A);
    }
    return A;
}

int PostProcess(Matrix_t* A) {
    ellpack_info_t* info = (ellpack_info_t*)(A->info);

    if (info) {
        if(info->COEF) free(info->COEF);
        if(info->ICOL) free(info->ICOL);
        if(info->factor) free(info->factor);
        free(A->info);
        A->info = NULL;
    }

    Matrix_free(A);
    return 0;
}

int LinearSolve(const Matrix_t *A, const double *b, double *x, const double tolerance) {
    int32_t neq = A->NROWS;
    ellpack_info_t* info = (ellpack_info_t*)(A->info);
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
    }

    return ierr;
}
