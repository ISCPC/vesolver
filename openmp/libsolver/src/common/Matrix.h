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

#pragma once
#include <stdint.h>
#include <math.h>
#include <float.h>

#ifdef MKL
#include <mkl.h>
#ifdef MKL_SBLAS
#include <mkl_spblas.h>
#endif /* MKL_SBLAS */
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
#define MATRIX_TYPE_CSC        (0<<8)
#define MATRIX_TYPE_CSR        (1<<8)
#define MATRIX_TYPE_ELLPACK    (2<<8)
#define MATRIX_TYPE_JAD        (3<<8)
#define MATRIX_TYPE_DIA        (4<<8)
#define MATRIX_TYPE_DCSC       (8<<8)
#define MATRIX_TYPE_DCSR       (9<<8)
#define MATRIX_TYPE_MASK       (0xf<<8)

#define MATRIX_INDEX_TYPE(A)    (((A)->flags)&0xf)
#define MATRIX_IS_SYMMETRIC(A)  (((A)->flags)&MATRIX_TYPE_SYMMETRIC)
#define MATRIX_IS_LOWER(A)      (((A)->flags)&MATRIX_TYPE_LOWER)
#define MATRIX_IS_UNIT(A)       (((A)->flags)&MATRIX_TYPE_UNIT)
#define MATRIX_IS_CSC(A)        ((((A)->flags)&MATRIX_TYPE_MASK) == MATRIX_TYPE_CSC)
#define MATRIX_IS_CSR(A)        ((((A)->flags)&MATRIX_TYPE_MASK) == MATRIX_TYPE_CSR)

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ellpack_info {
    int32_t nw;
    double* COEF;
    int32_t* ICOL;
} ellpack_info_t;

typedef struct jad_info {
    int32_t* perm;
    int32_t* ptr;
    int32_t* index;
    double* value;
    int32_t maxnzr;
} jad_info_t;

typedef struct Matrix {
    int NROWS;
    int NNZ;
    uint32_t flags;
    int* pointers;
    int* indice;
    double* values;
    ellpack_info_t _ellpack;
    jad_info_t _jad;
    void* info;
#ifdef MKL
    sparse_matrix_t hdl;
#endif
#ifdef SXAT
    sblas_handle_t hdl;
#endif
#ifdef SSL2
    double *w;
    int *iw;
#endif
    int optimized;
    int allocated;
} Matrix_t;

/*
 * Generic Matrix operations
 */
void Matrix_init_generic(Matrix_t *A);
void Matrix_free_generic(Matrix_t *A);
int Matrix_optimize_generic(Matrix_t *A);
int Matrix_MV_generic(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z);

/*
 * Architecture dependent Matrix operations (prototype for override)
 */
void Matrix_init(Matrix_t *A);
void Matrix_free(Matrix_t *A);
int Matrix_optimize(Matrix_t *A);
int Matrix_MV(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z);

/*
 * Architecture independent functions
 */
void Matrix_setMatrixCSR(Matrix_t *A, const int nrows, const int nnz, const int *aptr,
    const int *aind, const double *aval, const uint32_t flags);
Matrix_t* Matrix_duplicate(const Matrix_t* A);
int Matrix_convert_index(Matrix_t* A, int base);
int Matrix_transpose(Matrix_t* A);
int Matrix_extract_symmetric(Matrix_t* A);
int Matrix_create_ellpack(Matrix_t* A);
int Matrix_create_jad(Matrix_t* A);

#ifdef __cplusplus
}
#endif
