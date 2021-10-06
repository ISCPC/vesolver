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
