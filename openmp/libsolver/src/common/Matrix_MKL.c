#ifdef MKL
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"

/****************************************************************
 * Library with MKL
 ****************************************************************/
// override
void Matrix_init(Matrix_t *A) {
    A->pointers = NULL;
    A->indice = NULL;
    A->values = NULL;
    A->info = NULL;
    A->optimized = 0;

    if (mkl_enable_instructions(MKL_ENABLE_AVX512) != 0) {
        printf("INFO: AVX-512 enabled.\n");
    }
}

// override
void Matrix_free(Matrix_t *A) {
    if (A->optimized == 1) {
        mkl_sparse_destroy(A->hdl);
        A->optimized = 0;
    }

    Matrix_free_generic(A);
}

// override
int Matrix_optimize(Matrix_t *A) {
    struct matrix_descr descr;

    if (A->optimized == 1) {
        mkl_sparse_destroy(A->hdl);
        A->optimized = 0;
    }

    //printf("Matrix_optimize:MKL: Optimizing coefficient matrix.\n");
    //fflush(stdout);
    descr.type = MATRIX_IS_SYMMETRIC(A) ? SPARSE_MATRIX_TYPE_SYMMETRIC : SPARSE_MATRIX_TYPE_GENERAL;
    descr.mode = MATRIX_IS_LOWER(A) ? SPARSE_FILL_MODE_LOWER : SPARSE_FILL_MODE_UPPER;
    descr.diag = MATRIX_IS_UNIT(A) ? SPARSE_DIAG_UNIT : SPARSE_DIAG_NON_UNIT;

    sparse_status_t ierr = mkl_sparse_d_create_csr(&(A->hdl),
            (MATRIX_INDEX_TYPE(A) == 0) ? SPARSE_INDEX_BASE_ZERO : SPARSE_INDEX_BASE_ONE,
            A->NROWS, A->NROWS, A->pointers, &(A->pointers[1]), A->indice, A->values);
    if (ierr != SPARSE_STATUS_SUCCESS) {
        printf("ERROR: mkl_sparse_create_csr() failed with %d\n", ierr);
        exit(1);
    }

    /* Analysis of the sparse matrix A */
    ierr = mkl_sparse_set_mv_hint(A->hdl, SPARSE_OPERATION_NON_TRANSPOSE, descr, 1000);
    if (ierr != SPARSE_STATUS_SUCCESS) {
        printf("ERROR: mkl_sparse_set_mv_hint() failed with %d\n", ierr);
        exit(1);
    }

    ierr = mkl_sparse_optimize(A->hdl);
    if (ierr != SPARSE_STATUS_SUCCESS) {
        printf("ERROR: mkl_sparse_optimize() failed with %d\n", ierr);
        exit(1);
    }

    A->optimized = 1;
    return 0;
}

// override
int Matrix_MV(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z) {
    struct matrix_descr descr;

    if (A->optimized == 0) {
        printf("ERROR: Matrix is not optimized.\n");
        exit(1);
    }

    descr.type = MATRIX_IS_SYMMETRIC(A) ? SPARSE_MATRIX_TYPE_SYMMETRIC : SPARSE_MATRIX_TYPE_GENERAL;
    descr.mode = MATRIX_IS_LOWER(A) ? SPARSE_FILL_MODE_LOWER : SPARSE_FILL_MODE_UPPER;
    descr.diag = MATRIX_IS_UNIT(A) ? SPARSE_DIAG_UNIT : SPARSE_DIAG_NON_UNIT;

    for(int i=0; i<A->NROWS; i++) {
        z[i] = y[i];
    }

    sparse_status_t ierr = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A->hdl, descr, x, beta, z);
    if (ierr != SPARSE_STATUS_SUCCESS) {
        printf("ERROR: mkl_sparse_mv() failed with %d\n", ierr);
        exit(1);
    }
    return 0;
}
#endif /* MKL */

