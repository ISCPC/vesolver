#pragma once

typedef int32_t INT_T;

extern int itersolver(INT_T neq, INT_T nnz, INT_T *pointers, INT_T *indice, double *value,
    double* b, double* x, double res);

extern int itersolver_coo(INT_T neq, INT_T nnz, INT_T *pointers, INT_T *indice, double *value,
    double* b, double* x, double res);
