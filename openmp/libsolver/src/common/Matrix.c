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

Matrix_t* Matrix_duplicate(const Matrix_t* A) {
    Matrix_t* D = malloc(sizeof(Matrix_t));
    D->NROWS = A->NROWS;
    D->NNZ = A->NNZ;
    D->flags = A->flags;
    D->info = NULL;

    D->pointers = (int*)calloc(sizeof(int), A->NROWS+1);
    D->indice = (int*)calloc(sizeof(int), A->NNZ);
    D->values = (double*)calloc(sizeof(double), A->NNZ);

    bcopy(A->pointers, D->pointers, sizeof(int)*(A->NROWS+1));
    bcopy(A->indice, D->indice, sizeof(int)*(A->NNZ));
    bcopy(A->values, D->values, sizeof(double)*(A->NNZ));

#ifdef MKL
    D->hdl = 0;
#endif
#ifdef SXAT
    D->hdl = 0;
#endif
#ifdef SSL2
    D->ax = D->w = NULL;
    D->iw = NULL;
#endif

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

    if (MATRIX_IS_SYMMETRIC(A)) {
        A->flags ^= MATRIX_TYPE_LOWER;
    } else {
        A->flags ^= MATRIX_TYPE_CSR;
    }

    return 0;
}

/****************************************************************
 * Base Library
 ****************************************************************/
#pragma weak Matrix_init = Matrix_init_generic
void Matrix_init_generic(Matrix_t *A) {
    A->pointers = NULL;
    A->indice = NULL;
    A->values = NULL;
    A->info = NULL;
}

#pragma weak Matrix_free = Matrix_free_generic
void Matrix_free_generic(Matrix_t *A) {
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

    A->optimized = 0;
}

#pragma weak Matrix_optimize = Matrix_optimize_generic
int Matrix_optimize_generic(Matrix_t *A) {
    A->optimized = 1;
    return 0;
}

#pragma weak Matrix_MV = Matrix_MV_generic
int Matrix_MV_generic(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z) {
    if (MATRIX_INDEX_TYPE(A) == 1) {
        /*  matrix vector product  */
        if (MATRIX_IS_SYMMETRIC(A)) {
            for (int i=A->NROWS-1; i>=0; i--) {
                z[i] = (beta * y[i]) + (A->values[A->pointers[i]-1]* x[i]);
                for (int j=A->pointers[i]; j<A->pointers[i+1]-1; j++) {
                    int k = A->indice[j]-1;
                    z[i] += alpha * A->values[j] * x[k];
                    z[k] += alpha * A->values[j] * x[i];
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
