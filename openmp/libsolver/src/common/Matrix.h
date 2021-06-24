#pragma once
#include <stdint.h>
#include <math.h>
#include <float.h>

#ifdef MKL
#include <mkl.h>
#include <mkl_spblas.h>
#else
#include <cblas.h>
#endif /* MKL */

#ifdef SXAT
#include <sblas.h>
#endif

#define DNRM2(N,X) cblas_dnrm2(N,X,1)
#define DDOT(N,X,Y) cblas_ddot(N,X,1,Y,1)

#define MATRIX_TYPE_INDEX0     (0)
#define MATRIX_TYPE_INDEX1     (1)
#define MATRIX_TYPE_ASYMMETRIC (0<<4)
#define MATRIX_TYPE_SYMMETRIC  (1<<4)
#define MATRIX_TYPE_LOWER      (1<<5)
#define MATRIX_TYPE_UPPER      (0)
#define MATRIX_TYPE_UNIT       (1<<6)
#define MATRIX_TYPE_NON_UNIT   (0)

#define MATRIX_INDEX_TYPE(A)    (((A)->flags)&0xf)
#define MATRIX_IS_SYMMETRIC(A)  (((A)->flags)&MATRIX_TYPE_SYMMETRIC)
#define MATRIX_IS_LOWER(A)      (((A)->flags)&MATRIX_TYPE_LOWER)
#define MATRIX_IS_UNIT(A)       (((A)->flags)&MATRIX_TYPE_UNIT)

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Matrix {
    int NROWS;
    int NNZ;
    uint32_t flags;
    int* pointers;
    int* indice;
    double* values;
    void* info;
#ifdef MKL
    sparse_matrix_t hdl;
#endif
#ifdef SXAT
    sblas_handle_t hdl;
#endif
    int optimized;
} Matrix_t;

void Matrix_init(Matrix_t *A);
void Matrix_setMatrixCSR(Matrix_t *A, const int nrows, const int nnz, const int *aptr,
    const int *aind, const double *aval, const uint32_t flags);
void Matrix_free(Matrix_t *A);
int Matrix_optimize(Matrix_t *A);
int Matrix_MV(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z);

Matrix_t* Matrix_duplicate(const Matrix_t* A);
int Matrix_convert_index(Matrix_t* A, int base);
int Matrix_transpose(Matrix_t* A);

#ifdef __cplusplus
}
#endif
