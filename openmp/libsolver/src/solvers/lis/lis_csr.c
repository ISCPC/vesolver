#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"
#include "timelog.h"
#include "PluginAPI.h"

#undef _COMPLEX // workaround to avoid conflict with cblas.h
#include "lis_config.h"
#include "lis.h"

typedef int32_t INT_T;

typedef struct lis_info {
    LIS_MATRIX  Matrix;
	LIS_SOLVER  solver;
} lis_info_t;

#define STR_BUFSIZE 1024

static int solve_pre(Matrix_t* A) {
    LIS_INT ierr;
    char string[STR_BUFSIZE];

    Matrix_extract_symmetric(A);
    Matrix_convert_index(A, 0);

    lis_info_t* info = (lis_info_t*)malloc(sizeof(lis_info_t));
    int nn = A->NROWS;

    // Parameters
    //int SOLVER_number = 1; /* CG */
    //int SOLVER_number = 2; /* BiCG */
    //int SOLVER_number = 3; /* CGS */
    int SOLVER_number = 4; /* BiCGStab */
    //int SOLVER_number = 5; /* BiCGStab(l) */
    //int SOLVER_number = 6; /* GPBiCG */
    //int SOLVER_number = 7; /* TFQMR */
    //int SOLVER_number = 8; /* Orthomin(m) */
    //int SOLVER_number = 9; /* GMRES(m) */
    //int SOLVER_number = 10; /* Jacobi */
    //int SOLVER_number = 11; /* Gauss-Seidel */
    //int SOLVER_number = 12; /* SOR */
    //int SOLVER_number = 13; /* BiCGSafe */
    //int SOLVER_number = 14; /* CR */
    //int SOLVER_number = 15; /* BiCR */
    //int SOLVER_number = 16; /* CRS */
    //int SOLVER_number = 17; /* BiCRSTAB */
    //int SOLVER_number = 18; /* GPBiCR */
    //int SOLVER_number = 19; /* BiCRSafe */
    //int SOLVER_number = 20; /* FGMRES(m) */
    //int SOLVER_number = 21; /* IDR(s) */
    //int SOLVER_number = 22; /* IDR(1) */
    //int SOLVER_number = 23; /* MINRES */
    //int SOLVER_number = 24; /* COCG */
    //int SOLVER_number = 25; /* COCR */

    //int PRECOND_number = 0; /* None */
    int PRECOND_number = 1; /* Jacobi */
    //int PRECOND_number = 2; /* ILU */
    //int PRECOND_number = 3; /* SSOR */
    //int PRECOND_number = 4; /* Hybrid */
    //int PRECOND_number = 5; /* I+S */
    //int PRECOND_number = 6; /* SAINV */
    //int PRECOND_number = 7; /* SA-AMG */
    //int PRECOND_number = 8; /* Cront ILU */
    //int PRECOND_number = 9; /* ILUT */
    float LIS_AMG_THETA = 0.01;
    int itrmax = 5000;
    double err0 = 1.0e-6;

    // Matrix for Lis solver
	ierr = lis_matrix_create(0, &(info->Matrix));
	ierr = lis_matrix_set_size(info->Matrix, 0, nn);

    // Setup Lis Matrix from CSR (Compressed Row Storage)
    //ierr = lis_matrix_set_csr(A->NNZ, A->pointers, A->indice, (LIS_SCALAR*)A->values, info->Matrix);
    ierr = lis_matrix_set_csr(A->NNZ, A->pointers, A->indice, A->values, info->Matrix);
    if (ierr != 0) {
        fprintf(stderr,"ERROR: lis_matrixs_set_csr() fails with ierr=%d\n", ierr);
        return -1;
    }
    ierr = lis_matrix_assemble(info->Matrix);
    if (ierr != 0) {
        fprintf(stderr,"ERROR: lis_matrixs_set_assemble() fails with ierr=%d\n", ierr);
        return -1;
    }

#ifdef USE_SXAURORA_SBLAS
    // Setup SBLAS
    ierr = lis_matrix_set_sblas( nn, nn, A->pointers, A->indice, A->values);
    if (ierr != 0) {
        fprintf(stderr,"ERROR: lis_matrixs_set_sblas() fails with ierr=%d\n", ierr);
        return -1;
    }
    if (ierr != 0) fprintf(stderr,"lis.c#49: ierr=%d\n", ierr);
#endif

    // Lis Solver
    ierr = lis_solver_create(&(info->solver));
    if (ierr != 0) {
        fprintf(stderr,"ERROR: lis_solver_create() fails with ierr=%d\n", ierr);
        return -1;
    }

    // Lis options
    sprintf(string, "-i %d -p %d -maxiter %d -tol %f",
         SOLVER_number, PRECOND_number, itrmax, err0);
    // AMG options
    if (PRECOND_number == 7) {
        snprintf(string, STR_BUFSIZE, "%s -saamg_unsym true -saamg_theta %f", string, LIS_AMG_THETA);
    }

    printf("INFO: Lis solver with option \"%s\"\n", string);
    ierr = lis_solver_set_option(string,info->solver);
    if (ierr != 0) {
        fprintf(stderr,"WARNING: lis_solver_set_option() failed with ierr=%d\n", ierr);
        return -1;
    }

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

