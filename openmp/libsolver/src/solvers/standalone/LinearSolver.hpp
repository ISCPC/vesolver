#pragma once
#include <stdint.h>

typedef int32_t INT_T;

#define LS_CONVERGENCE       0
#define LS_BICGSTAB_2_RHO    -1
#define LS_DIVERGENCE        -2
#define LS_MAXITER           -3

class LinearSolver {
public:
    LinearSolver() {}
    ~LinearSolver() {}

    virtual int init() {return 0;}
    virtual int finalize() {return 0;}

    virtual int setMatrixCSR(const INT_T neq, const INT_T nnz,
        const INT_T *pointers, const INT_T *indice, const double *value,
        const uint32_t flags) {return -1;}
    virtual int optimize() {return -1;}
    virtual int solve(const double* b, double* x, const double res) {return -1;}
    virtual double residual(const double* b, const double* x, const int mode) {return 0.0f;}
};

extern LinearSolver* getSolver();
