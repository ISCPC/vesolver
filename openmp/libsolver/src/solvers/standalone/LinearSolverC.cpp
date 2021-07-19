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
#include <cstdio>
#include <stdlib.h>
#include "LinearSolver.hpp"
#include "Matrix.h"
#include "PluginAPI.h"

extern "C" {
	SolverPlugin_t* solver_init();
}

/*
 * Definition IterSolver
 */
class LinearSolverC : public LinearSolver {
public:
    LinearSolverC() {
   		solver = solver_init();
    }

    ~LinearSolverC() {
    	solver->free(solver);
    	free(solver);
    }

    int init() override {
        A = (Matrix_t*)malloc(sizeof(Matrix_t));
        Matrix_init(A);
        return 0;
    };

    int finalize() override {
        if (D != NULL) {
            if (solver->solve_post(D) != 0) {
                printf("WARNING: PostProcess() failed.\n");
            }
            Matrix_free(D);
        }

        Matrix_free(A);
        return 0;
    };

    int setMatrixCSR(const INT_T neq, const INT_T nnz, const INT_T *pointers,
        const INT_T *indice, const double *value, const uint32_t flags) override {
        Matrix_setMatrixCSR(A, neq, nnz, pointers, indice, value, flags);
        return 0;
    }

    int optimize() override {
        D = Matrix_duplicate(A);
    	return solver->solve_pre(D);
    }

    int solve(const double* b, double* x, const double res)  override {
        return (D != NULL) ? solver->solve(D, b, x, res) : -1;
    }

    double residual(const double* b, const double* x, const int mode) override {
        double *z = (double*)calloc(sizeof(double), A->NROWS);
        double res = 0.0f;

        int cc = Matrix_optimize(A);
        Matrix_MV(A, 1.0, x, -1.0, b, z);
        res = DNRM2(A->NROWS, z);

        if (mode == 1) {
            res /= DNRM2(A->NROWS, b);
        }

        free(z);
        return res;
    }

private:
    Matrix_t* A;
    Matrix_t* D=NULL;
    SolverPlugin_t* solver;
};

LinearSolver* getSolver() {
    return (LinearSolver*)new LinearSolverC();
}
