#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include <heterosolver.h>
#include "Matrix.h"
#include "PluginAPI.h"
#include "timelog.h"

#ifdef FTRACE
#include <ftrace.h>

#define FTRACE_REGION_BEGIN(s) (void)ftrace_region_begin(s)
#define FTRACE_REGION_END(s)   (void)ftrace_region_end(s)
#else
#define FTRACE_REGION_BEGIN(s)
#define FTRACE_REGION_END(s)
#endif

#define   NRHS   1   /* The number of right-hand side vectors */

typedef int32_t INT_T;

typedef struct hs_info {
    HS_handle_t hnd;
} hs_info_t;

/*
 * API
 */
static Matrix_t* solve_pre(const Matrix_t* A0) {
    int ierr;
    TIMELOG(tl);

    Matrix_t* A = Matrix_duplicate(A0);
    //Matrix_convert_index(A, 1);

    hs_info_t* info = (hs_info_t*)malloc(sizeof(hs_info_t));

    HS_int_t isym = MATRIX_IS_SYMMETRIC(A) ? HS_SYMMETRIC : HS_UNSYMMETRIC;
    HS_int_t iformat = MATRIX_IS_CSR(A) ? HS_CSR : HS_CSC;
    
    ierr = HS_init_handle(&(info->hnd), A->NROWS, A->NROWS, isym, iformat);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: HS_init_handle failed with %d.\n", ierr);
        return NULL;
    }

    ierr = HS_set_option(info->hnd, HS_ORDP, HS_ORDP_METIS);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: HS_set_option(HS_ORDP) failed with %d.\n", ierr);
        return NULL;
    }

    if (MATRIX_INDEX_TYPE(A) == 1) {
        ierr = HS_set_option(info->hnd, HS_INDEXING, HS_INDEXING_1);
        if (ierr != HS_RESULT_OK) {
            fprintf(stderr, "ERROR: HS_set_option(HS_INDEXING) failed with %d.\n", ierr);
            return NULL;
        }
    }

    FTRACE_REGION_BEGIN("HS_preprocess");
    TIMELOG_START(tl);
    ierr = HS_preprocess_rd(info->hnd, A->pointers, A->indice, A->values);
    TIMELOG_END(tl, "hs_preprocess");
    FTRACE_REGION_END("HS_preprocess");
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: HS_preprocess_rd failed with %d.\n", ierr);
        return NULL;
    }

    FTRACE_REGION_BEGIN("HS_factorize");
    TIMELOG_START(tl);
    ierr = HS_factorize_rd(info->hnd, A->pointers, A->indice, A->values);
    TIMELOG_END(tl, "hs_factorize");
    FTRACE_REGION_END("HS_factorize");
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: HS_factorize_rd failed with %d.\n", ierr);
        return NULL;
    }

    A->info = (void*)info;
    return A;
}

static int solve(const Matrix_t *A, const double* b, double* x, const double tolerance) {
    int ierr;
    TIMELOG(tl);

    hs_info_t* info = (hs_info_t*)(A->info);

    /* Solution Phase */
    //double res = 1.0e-13;  
    double res = tolerance;

    /* Solution Phase */
    FTRACE_REGION_BEGIN("HS_solve");
    TIMELOG_START(tl);
    ierr = HS_solve_rd(info->hnd, A->pointers, A->indice, A->values, NRHS, (double*)b, x, &res);
    TIMELOG_END(tl, "hs_solve");
    FTRACE_REGION_END("HS_solve");
    if (ierr != HS_RESULT_OK) {
        if (ierr == HS_ERROR_ACCURACY) {
            fprintf(stderr, "ERROR: HS_solve_rd failed with HS_ERROR_ACCURACY (res=%e).\n", res);
        } else {
            fprintf(stderr, "ERROR: HS_solve_rd failed with %d.\n", ierr);
        }
        return -1;
    }

    /* Handle Finalization */
    ierr = HS_finalize_handle(info->hnd);

    return 0;
}

static int solve_post(Matrix_t* A) {
    Matrix_free(A);

    return 0;
}

static void solve_free(SolverPlugin_t* solver) {
	return;
}

/*
 * Solver Plugin Interface
 */
#ifdef _STANDALONE
SolverPlugin_t* solver_init() {
#else
SolverPlugin_t* hs_init() {
#endif
	SolverPlugin_t* solver = (SolverPlugin_t*)malloc(sizeof(SolverPlugin_t));

	solver->set_option = NULL;
	solver->solve_pre = solve_pre;
	solver->solve = solve;
	solver->solve_post = solve_post;
	solver->free = solve_free;

	return solver;
}

