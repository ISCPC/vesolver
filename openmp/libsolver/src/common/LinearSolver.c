/*
 * api.c
 *
 *  Created on: Jun 6, 2021
 *      Author: uno
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "Matrix.h"
#include "PluginAPI.h"
#include "LinearSolver.h"
#include "timelog.h"

static Matrix_t* A;
static Matrix_t* D=NULL;

#define MAX_HANDLES 16

static SolverPlugin_t* solver;

/*
 * Internal functions
 */
SolverPlugin_t* cg_init();
SolverPlugin_t* elmer_cg_init();
SolverPlugin_t* bicgstab2_init();
#ifdef SXAT
SolverPlugin_t* hs_init();
#endif
#ifdef MKL
SolverPlugin_t* pardiso_init();
#endif

static SolverPlugin_t* get_solver(int solverId) {
	switch(solverId) {
	case SOLVER_ITER_CG_SYM:
		return cg_init();

	case SOLVER_ITER_CG_ASYM:
		return elmer_cg_init();

	case SOLVER_ITER_BICGSTAB2:
		return bicgstab2_init();

#ifdef SXAT
	case SOLVER_DIRECT_HS:
		return hs_init();
#endif

#ifdef MKL
	case SOLVER_DIRECT_PARDISO:
		return pardiso_init();
#endif
	}

	return NULL;
}

/*
 * API Definition for libLinearSolver
 */
int solver_init() {
    solver = NULL;
	return 0;
};

void solver_finalize() {
    if (solver != NULL) {
	    solver->free(solver);
        free(solver);
    }
	return;
};

SolverHandle_t solver_create_handle() {
	A = (Matrix_t*)malloc(sizeof(Matrix_t));
	Matrix_init(A);

	return 0;
}

int solver_free_handle(SolverHandle_t hdl) {
	if (solver->solve_post(D) != 0) {
		printf("WARNING: PostProcess() failed.\n");
	}
	Matrix_free(A);
	return 0;
}

int solver_set_option(SolverHandle_t hdl, int id, int value) {
	switch(id) {
	case SOLVER_OPTION_SOLVER:
		solver = get_solver(value);
		if (solver == NULL) {
			fprintf(stderr, "ERROR: Unsupported solver (id=%d)\n", value);
			return -1;
		}
		break;
	}
	return 0;
}

int solver_set_matrix_csr(SolverHandle_t hdl,
		const INT_T neq, const INT_T nnz, const INT_T *pointers,
		const INT_T *indice, const double *value, const uint32_t flags) {
    TIMELOG(tl);

	TIMELOG_START(tl);
    Matrix_setMatrixCSR(A, neq, nnz, pointers, indice, value, flags);
    TIMELOG_END(tl, "setMatrix");

#if 0 /* DEBUG */
    printf("INFO: Matrix_A: %s, %s, INDEX_%d, nrows=%ld, nnz=%ld, flags=0x%08x\n", 
        MATRIX_IS_CSR(A) ? "CSR" : "CSC", MATRIX_IS_SYMMETRIC(A) ? "SYMMETRIC" : "ASYMMETRIC",
        MATRIX_INDEX_TYPE(A), neq, nnz, flags);
#endif

	// Optimize Coefficient Matrix internally
	TIMELOG_START(tl);
    D = Matrix_duplicate(A);
	int cc = solver->solve_pre(D);
    TIMELOG_END(tl, "preProcess");

    return cc;
}

int solver_solve(SolverHandle_t hdl, const double* b, double* x, const double res) {
	//res=1.e-4;
	return (D != NULL) ? (solver->solve)(D, b, x, res) : -1;
}

//int solver_solve_sync(SolverHandle_t hdl, const double* b, double* x, const double res);
//int solver_solve_wait(SolverHandle_t hdl);

double solver_calc_residual(SolverHandle_t hdl, const double* b, const double* x, const int mode) {
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
