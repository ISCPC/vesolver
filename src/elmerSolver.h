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
#ifndef __ELMERSOLVER__
#define __ELMERSOLVER__
#ifdef __cplusplus
extern "C" {
#endif
#include <stdint.h>

typedef int32_t INT_T;
typedef uint32_t UINT_T;

#define ELMER_BICGSTAB      0
#define ELMER_BICGSTAB2     1
#define ELMER_BICGSTABL     2
#define ELMER_CG            3
#define ELMER_CGS           4
#define ELMER_TFQMR         5
#define ELMER_GMRES         6
#define ELMER_SGS           7
#define ELMER_GCR           8
#define ELMER_IDRS          9

#define ELMER_PC_NONE       0
#define ELMER_PC_DIAG       1
#define ELMER_PC_ILUT       2
#define ELMER_PC_MG         3
#define ELMER_PC_ILU0      10

int elmersolver(INT_T neq, INT_T nnz, INT_T *pointers, INT_T *indice, double *value,
    double* b, double* x, INT_T solverId);

int elmersolver_distributed(INT_T neq, INT_T nnz, INT_T *pointers, INT_T *indice, double *value,
    double* b, double* x, INT_T solverId);

void elmersolverapi_(INT_T *neq, INT_T *ndim, INT_T *pointers, INT_T *indice, double *value, double *b, double *x, INT_T *solverId, INT_T *);

#ifdef __cplusplus
}
#endif
#endif /* __ELMERSOLVER__ */
