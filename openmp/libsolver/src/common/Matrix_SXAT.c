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

#ifdef SXAT
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"

/****************************************************************
 * Library with SX-Aurora TSUBASA NLC
 ****************************************************************/
void Matrix_init(Matrix_t *A) {
    A->pointers = NULL;
    A->indice = NULL;
    A->values = NULL;
    A->info = NULL;
    A->optimized = 0;
}

void Matrix_free(Matrix_t *A) {
    /* Destruction of the handle */
    if (A->optimized == 1) {
        sblas_destroy_matrix_handle(A->hdl);
        A->optimized = 0;
    }

    Matrix_free_generic(A);
}

int Matrix_optimize(Matrix_t *A) {
    /* Creation of a handle from CSR format */
    int ierr = sblas_create_matrix_handle_from_csr_rd(
        A->NROWS, A->NROWS, A->pointers, A->indice, A->values,
        (MATRIX_INDEX_TYPE(A) == 0) ? SBLAS_INDEXING_0 : SBLAS_INDEXING_1,
        MATRIX_IS_SYMMETRIC(A) ? SBLAS_SYMMETRIC : SBLAS_GENERAL,
        &A->hdl);
    if (ierr != SBLAS_OK) {
        printf("ERROR: sblas_create_matrix_handle_from_csr_rd() failed with %d\n", ierr);
        exit(1);
    }

    /* Analysis of the sparse matrix A */
    ierr = sblas_analyze_mv_rd(SBLAS_NON_TRANSPOSE, A->hdl);
    if (ierr != SBLAS_OK) {
        printf("ERROR: sblas_analyze_mv_rd() failed with %d\n", ierr);
        exit(1);
    }

    A->optimized = 1;
    return 0;
}

int Matrix_MV(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z) {
    for(int i=0; i<A->NROWS; i++) {
        z[i] = y[i];
    }
    int ierr = sblas_execute_mv_rd(SBLAS_NON_TRANSPOSE, A->hdl, alpha, x, beta, z);
    if (ierr != SBLAS_OK) {
        printf("ERROR: sblas_execute_mv_rd() failed with %d\n", ierr);
        exit(1);
    }

    return 0;
}
#endif /* SXAT */
