#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"

/****************************************************************
 * Common Functions
 ****************************************************************/
void Matrix_setMatrixCSR(Matrix_t *A, const int nrows, const int nnz, const int *aptr,
    const int *aind, const double *aval, const uint32_t flags) {
    A->NROWS = nrows;
    A->NNZ = nnz;
    A->flags = flags;

    A->pointers = (int*)calloc(sizeof(int), A->NROWS+1);
    A->indice = (int*)calloc(sizeof(int), A->NNZ);
    A->values = (double*)calloc(sizeof(double), A->NNZ);

    bcopy(aptr, A->pointers, sizeof(int)*(A->NROWS+1));
    bcopy(aind, A->indice, sizeof(int)*(A->NNZ));
    bcopy(aval, A->values, sizeof(double)*(A->NNZ));

    A->optimized = 0;
}

static inline void Matrix_free_common(Matrix_t *A) {
    if (A->pointers) {
        free(A->pointers);
        A->pointers = NULL;
    }

    if (A->indice) {
        free(A->indice);
        A->indice = NULL;
    }
    if (A->values) {
        free(A->values);
        A->values = NULL;
    }
    if (A->info) {
        free(A->info);
        A->info = NULL;
    }
}

Matrix_t* Matrix_duplicate(const Matrix_t* A) {
    Matrix_t* D = malloc(sizeof(Matrix_t));
    D->NROWS = A->NROWS;
    D->NNZ = A->NNZ;
    D->flags = A->flags;

    D->pointers = (int*)calloc(sizeof(int), A->NROWS+1);
    D->indice = (int*)calloc(sizeof(int), A->NNZ);
    D->values = (double*)calloc(sizeof(double), A->NNZ);

    bcopy(A->pointers, D->pointers, sizeof(int)*(A->NROWS+1));
    bcopy(A->indice, D->indice, sizeof(int)*(A->NNZ));
    bcopy(A->values, D->values, sizeof(double)*(A->NNZ));

    return D;
}

int Matrix_convert_index(Matrix_t* A, int base) {
    if ((MATRIX_INDEX_TYPE(A) == 1)&&(base == 0)) {
        /* Convert from one-based index to zero-based index */
        for (int i=0; i<=A->NROWS; i++) {
            A->pointers[i] = A->pointers[i]-1;
        }

        for (int i=0; i<A->NNZ; i++) {
            A->indice[i] = A->indice[i]-1;
        }

        A->flags = ((A->flags)&(~0xf)) | MATRIX_TYPE_INDEX0;
    } else if ((MATRIX_INDEX_TYPE(A) == 0)&&(base == 1)) {
        /* Convert from zero-based index to one-based index */
        for (int i=0; i<=A->NROWS; i++) {
            A->pointers[i] = A->pointers[i]+1;
        }

        for (int i=0; i<A->NNZ; i++) {
            A->indice[i] = A->indice[i]+1;
        }

        A->flags = ((A->flags)&(~0xf)) | MATRIX_TYPE_INDEX1;
    }

    return 0;
}

int Matrix_transpose(Matrix_t* A) {
    int* aptr = (int*)calloc(sizeof(int), A->NROWS+1);
    int* aind = (int*)calloc(sizeof(int), A->NNZ);
    double* aval = (double*)calloc(sizeof(double), A->NNZ);
#if 0
    /* phase0: */
    for(int i=0; i<=A->NROWS; i++) {
        aptr[i] = 0;
    }
#endif
    /* phase1 */
    for(int i=0; i<A->NNZ; i++) {
        int j = A->indice[i];
        aptr[j]++;
    }
    for(int i=1; i<=A->NROWS; i++) {
        aptr[i] += aptr[i-1];
    }

    /* phase2 */
    for (int i=A->NROWS-1; i>=0; i--) {
        for (int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
            int jj = --(aptr[A->indice[j]]);
            aind[jj] = i;
            aval[jj] = A->values[j];
        }
    }

    free(A->pointers);
    free(A->indice);
    free(A->values);
    A->pointers = aptr;
    A->indice = aind;
    A->values = aval;

    A->flags ^= MATRIX_TYPE_LOWER;

    return 0;
}

/****************************************************************
 * Base Library
 ****************************************************************/
#if (!defined(MKL) && !defined(SXAT))
void Matrix_init(Matrix_t *A) {
    A->pointers = NULL;
    A->indice = NULL;
    A->values = NULL;
    A->info = NULL;
}

void Matrix_free(Matrix_t *A) {
    A->optimized = 0;
    Matrix_free_common(A);
}

int Matrix_optimize(Matrix_t *A) {
    A->optimized = 1;
    return 0;
}

int Matrix_MV(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z) {
    if (MATRIX_INDEX_TYPE(A) == 1) {
        /*  matrix vector product  */
        if (MATRIX_IS_SYMMETRIC(A)) {
            for (int i=A->NROWS-1; i>=0; i--) {
                z[i] = (beta * y[i]) + (A->values[A->pointers[i]-1]* x[i]);
                for (int j=A->pointers[i]; j<A->pointers[i+1]-1; j++) {
                    int k = A->indice[j]-1;
                    z[i] += A->values[j] * x[k];
                    z[k] += A->values[j] * x[i];
                }
            }
        } else {
            #pragma omp parallel for
            for(int i=0; i<A->NROWS; i++) {
                z[i] = beta * y[i];
                for(int j=A->pointers[i]-1; j<A->pointers[i+1]-1; j++) {
                    z[i] += alpha * A->values[j] * x[A->indice[j]-1];
                }
            }
        }
    } else {
        /*  matrix vector product  */
        if (MATRIX_IS_SYMMETRIC(A)) {
            for (int i=A->NROWS-1; i>=0; i--) {
                z[i] = (beta * y[i]) + (A->values[A->pointers[i]]* x[i]);
                for (int j=A->pointers[i]+1; j<A->pointers[i+1]; j++) {
                    int k = A->indice[j];
                    z[i] += A->values[j] * x[k];
                    z[k] += A->values[j] * x[i];
                }
            }
        } else {
            #pragma omp parallel for
            for(int i=0; i<A->NROWS; i++) {
                z[i] = beta * y[i];
                for(int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
                    z[i] += alpha * A->values[j] * x[A->indice[j]];
                }
            }
        }
    }
    return 0;
}
#endif



/****************************************************************
 * Library with MKL
 ****************************************************************/
#ifdef MKL
void Matrix_init(Matrix_t *A) {
    A->pointers = NULL;
    A->indice = NULL;
    A->values = NULL;
    A->info = NULL;

    if (mkl_enable_instructions(MKL_ENABLE_AVX512) != 0) {
        printf("INFO: AVX-512 enabled.\n");
    }
}

void Matrix_free(Matrix_t *A) {
    mkl_sparse_destroy(A->hdl);
    A->optimized = 0;

    Matrix_free_common(A);
}

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


/****************************************************************
 * Library with SX-Aurora TSUBASA NLC
 ****************************************************************/
#ifdef SXAT
void Matrix_init(Matrix_t *A) {
    A->pointers = NULL;
    A->indice = NULL;
    A->values = NULL;
    A->info = NULL;
}

void Matrix_free(Matrix_t *A) {
    /* Destruction of the handle */
    sblas_destroy_matrix_handle(A->hdl);

    Matrix_free_common(A);
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
