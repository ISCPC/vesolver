#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"
#include "timelog.h"
#include "PluginAPI.h"
#include "lis.h"

typedef int32_t INT_T;

typedef struct lis_info {
    LIS_MATRIX  Matrix;
	LIS_SOLVER  solver;
} lis_info_t;

#define STR_BUFSIZE 256

static int solve_pre(Matrix_t* A) {
    LIS_INT ierr;
    char string[STR_BUFSIZE];

    Matrix_convert_index(A, 0);

    lis_info_t* info = (lis_info_t*)malloc(sizeof(lis_info_t));
    int nn = A->NROWS;

    // Parameters
    int SOLVER_number = 4; /* BiCGStab */
    int PRECOND_number = 1; /* Jacobi */
    //int PRECOND_number = 7; /* SA-AMG */
    float LIS_AMG_THETA = 0.01;
    int itrmax = 1000;
    double err0 = 1.0e-6;

    // Matrix for Lis solver
	ierr = lis_matrix_create(0, &(info->Matrix));
	ierr = lis_matrix_set_size(info->Matrix, 0, nn);

    // Setup Lis Matrix from CSR (Compressed Row Storage)
    ierr = lis_matrix_set_csr(A->NNZ, A->pointers, A->indice, A->values, info->Matrix);
    ierr = lis_matrix_assemble(info->Matrix);

    // Lis Solver
    ierr = lis_solver_create(&(info->solver));

    // Lis options
    // solver number
    snprintf(string, STR_BUFSIZE, "-i %d", SOLVER_number);
    ierr = lis_solver_set_option(string,info->solver);
    // precondition number
    snprintf(string, STR_BUFSIZE, "-p %d", PRECOND_number);
    ierr = lis_solver_set_option(string,info->solver);
    // AMG options
    if (PRECOND_number == 7) {
        snprintf(string, STR_BUFSIZE, "-saamg_unsym true -saamg_theta %f", LIS_AMG_THETA);
        ierr = lis_solver_set_option(string,info->solver);
    }
    // Iterations, Tolerance
    snprintf(string, STR_BUFSIZE, "-maxiter %d -tol %f", itrmax, err0);
    ierr = lis_solver_set_option(string,info->solver);

    A->info = (void*)info;
    return 0;
}

static int solve(Matrix_t *A, const double* b, double* x, const double tolerance) {
	LIS_VECTOR lx,lb;
    LIS_INT itr;
    LIS_REAL rr;
    LIS_INT ierr;
    TIMELOG(tl1);

    lis_info_t* info = (lis_info_t*)(A->info);
    int nn = A->NROWS;
    int isw = 0;  // Temporary

    // Vector for Lis solver
    ierr = lis_vector_create(0, &lb);
    ierr = lis_vector_create(0, &lx);
    ierr = lis_vector_set_size(lb, 0, nn);
    ierr = lis_vector_set_size(lx, 0, nn);

    // Set Values to Lis Vector
    if(isw == 0) {
        for(int i=0; i<nn; i++) {
          // INITIALIZE VECTOR x(i)
          ierr = lis_vector_set_value(LIS_INS_VALUE, i, 0.0, lx);
          ierr = lis_vector_set_value(LIS_INS_VALUE, i, b[i], lb);
        }
    } else {
        for(int i=0; i<nn; i++) {
          // Use INPUT x(i) as INITIAL VECTOR
          ierr = lis_vector_set_value(LIS_INS_VALUE, i, x[i], lx);
          ierr = lis_vector_set_value(LIS_INS_VALUE, i, b[i], lb);
        }
    }

    // SOLVE
    ierr = lis_solve(info->Matrix, lb, lx, info->solver);
    if (ierr != 0) {
        printf("ERROR: lis_solve (error=%d)\n", ierr);
        return -1;
    }

    // get iterations, final residue
    ierr = lis_solver_get_iter(info->solver, &itr);
    ierr = lis_solver_get_residualnorm(info->solver, &rr);
    printf("INFO: lis_solve done. (itr=%d, residual=%le)\n", itr, rr);

    // FINAL RESULT
    for(int i=0; i<nn; i++) {
        ierr = lis_vector_get_value(lx, i, &x[i]);
    }

    ierr = lis_vector_destroy(lb);
    ierr = lis_vector_destroy(lx);

    return 0;
}


static int solve_post(Matrix_t* A) {
    LIS_INT ierr;

    lis_info_t* info = (lis_info_t*)(A->info);

    // Release memories
    ierr = lis_solver_destroy(info->solver);
    ierr = lis_matrix_destroy(info->Matrix);
    A->pointers = A->indice = NULL;
    A->values = NULL;

    return 0;
}

static void solve_free(SolverPlugin_t* solver) {
	lis_finalize();

	return;
}

/*
 * Solver Plugin Interface
 */
#ifdef _STANDALONE
SolverPlugin_t* solver_init() {
#else
SolverPlugin_t* lis_init() {
#endif
    int argc=0;
    char** argv=NULL;

	SolverPlugin_t* solver = (SolverPlugin_t*)malloc(sizeof(SolverPlugin_t));

	solver->set_option = NULL;
	solver->solve_pre = solve_pre;
	solver->solve = solve;
	solver->solve_post = solve_post;
	solver->free = solve_free;

	lis_initialize(&argc, &argv);

	return solver;
}

