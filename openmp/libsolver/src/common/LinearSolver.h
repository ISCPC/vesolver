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

#ifndef COMMON_LINEARSOLVER_H_
#define COMMON_LINEARSOLVER_H_
#include "Matrix.h"

typedef int32_t INT_T;
typedef int SolverHandle_t;

#ifdef __cplusplus
extern "C" {
#endif

int solver_init();
void solver_finalize();

SolverHandle_t solver_create_handle();
int solver_free_handle(SolverHandle_t hdl);

int solver_set_option(SolverHandle_t hdl, int id, int value);
int solver_set_matrix_csr(SolverHandle_t hdl,
		const INT_T neq, const INT_T nnz, const INT_T *pointers,
		const INT_T *indice, const double *value, const uint32_t flags);
int solver_solve(SolverHandle_t hdl, const double* b, double* x, const double res);
//int solver_solve_sync(SolverHandle_t hdl, const double* b, double* x, const double res);
//int solver_solve_wait(SolverHandle_t hdl);
double solver_calc_residual(SolverHandle_t hdl, const double* b, const double* x, const int mode);


//
// Option
//
#define SOLVER_OPTION_SOLVER        1

#define SOLVER_ITER_CG_SYM          1
#define SOLVER_ITER_CG_ASYM         2
#define SOLVER_ITER_BICGSTAB2       3

#define SOLVER_DIRECT_HS         1000
#define SOLVER_DIRECT_PARDISO    1001

#endif /* COMMON_LINEARSOLVER_H_ */

#ifdef __cplusplus
}
#endif
