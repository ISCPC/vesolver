/*
 * LinearSolver.h
 *
 *  Created on: Jun 7, 2021
 *      Author: uno
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
#define SOLVER_OPTION_SOLVER       1

#define SOLVER_ITER_CG             1
#define SOLVER_ITER_BICGSTAB2      2
#define SOLVER_DIRECT_HS           3
#define SOLVER_DIRECT_PARDISO      4

#endif /* COMMON_LINEARSOLVER_H_ */

#ifdef __cplusplus
}
#endif
