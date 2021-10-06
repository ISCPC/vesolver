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

#ifndef VESOLVER_H_
#define VESOLVER_H_
#include <stdint.h>

typedef int32_t INT_T;

/*
 * From Matrix.h
 */
typedef struct vesolver_option {
    int64_t solver;
    double res;
} vesolver_option_t;

typedef struct vesolver_desc {
    int32_t ves_handle;
    INT_T neq;
    INT_T nnz;
    uint32_t flags;
    uint64_t dptr;
    uint64_t pointers;
    uint64_t indice;
    uint64_t values;
    uint64_t bptr;
    uint64_t xptr;
    vesolver_option_t options;
} vesolver_desc_t;

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

#if 0
#define MATRIX_INDEX_TYPE(A)    (((A)->flags)&0xf)
#define MATRIX_IS_SYMMETRIC(A)  (((A)->flags)&MATRIX_TYPE_SYMMETRIC)
#define MATRIX_IS_LOWER(A)      (((A)->flags)&MATRIX_TYPE_LOWER)
#define MATRIX_IS_UNIT(A)       (((A)->flags)&MATRIX_TYPE_UNIT)
#endif

/*
 * From LinearSolver.h
 */
#define VESOLVER_OPTION_SOLVER       1

#define VESOLVER_ITER_CG             1
#define VESOLVER_ITER_BICGSTAB2      2
#define VESOLVER_DIRECT_HS           3
#define VESOLVER_DIRECT_PARDISO      4

/*
 * VE solver API
 */
#if 0
typedef int vesolver_handle_t;

#ifdef __cplusplus
extern "C" {
#endif

int vesolver_init();
void vesolver_finalize();
vesolver_handle_t vesolver_activate();
int vesolver_deactivate(vesolver_handle_t hdl);
int vesolver_set_option(vesolver_handle_t hdl, int id, int value);
#if 0
int vesolver_set_matrix_csr(vesolver_handle_t hdl, INT_T neq,
	INT_T* pointers, INT_T* indice, double* value, uint32_t flags);
int vesolver_set_matrix_csc(vesolver_handle_t hdl, INT_T neq,
	INT_T* pointers, INT_T* indice, double* value, uint32_t flags);
#else
matrix_desc_t* vesolver_alloc_matrix(vesolver_handle_t hdl, INT_T neq,
	INT_T nnz, uint32_t flags);
int vesolver_set_matrix(vesolver_handle_t hdl, matrix_desc_t *desc);
#endif
int vesolver_free_matrix(vesolver_handle_t hdl, matrix_desc_t* desc);
int vesolver_solve_sync(vesolver_handle_t hdl, double* b, double* x, double res);
//int vesolver_solve_async();
//int vesolver_solve_wait();

#ifdef __cplusplus
}
#endif
#endif

#endif /* VESOLVER_H_ */
