/*
 ============================================================================
 Name        : solvertest.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <stdint.h>
#include <stdlib.h>
#include "SpMatrix.hpp"
#include "timelog.h"
#include "LinearSolver.h"


#define MAX_PATH_LEN 1024
#define DATA_SUFFIX ".bin"

uint32_t load_matrix(SpMatrix &A, char* file_a) {
    uint32_t flags=0;

    /* Load Matrix A */
    bool oldflag = getenv("SPMATRIX_FORMAT_OLD") ? true : false;
    if (A.load_file(file_a, oldflag) < 0) {
        fprintf(stderr, "WARNING: cannot load matrix A. Trying old format.\n");
        if (A.load_file(file_a, true) < 0) {
            fprintf(stderr, "ERROR: cannot load matrix A from %s\n", file_a);
            exit(1);
        }
    }

    if (A.type & SPMATRIX_TYPE_SYMMETRIC) {
        flags |= MATRIX_TYPE_SYMMETRIC;
        printf("INFO: Matrix_A: SYMMETRIC, ");
    } else {
        flags |= MATRIX_TYPE_ASYMMETRIC;
        printf("INFO: Matrix_A: ASYMMETRIC, ");
    }
    if (A.type & SPMATRIX_TYPE_INDEX1) {
        flags |= MATRIX_TYPE_INDEX1;
        printf("INDEX_1, ");
    } else {
        flags |= MATRIX_TYPE_INDEX0;
        printf("INDEX_0, ");
    }
    printf("nrows=%d, nnz=%d, neq=%d, type=0x%08x\n", A.nrow, A.ndim, A.neq, A.type);

    return flags;
}

void load_vector(Vector &b, char* file_b) {
    if (b.load_file(file_b) < 0) {
        fprintf(stderr, "ERROR: cannot load vector b from %s\n", file_b);
        exit(1);
    }
}

int main(int argc, char** argv) {
    SpMatrix A, A0;
    Vector b, b0, x;
    char file_a[MAX_PATH_LEN], file_b[MAX_PATH_LEN];
    //double res = 1.0e-10;
    double res = 1.0e-4;
    //double res = 1.0e-1;
    int cc;
    TIMELOG(tl);

    if (argc > 1) {
        snprintf(file_a, MAX_PATH_LEN, "%s/a" DATA_SUFFIX, argv[1]);
        snprintf(file_b, MAX_PATH_LEN, "%s/b" DATA_SUFFIX, argv[1]);
    } else {
        snprintf(file_a, MAX_PATH_LEN, "a" DATA_SUFFIX);
        snprintf(file_b, MAX_PATH_LEN, "b" DATA_SUFFIX);
    }

    /* Load Matrix A */
    uint32_t flags = load_matrix(A, file_a);

    /* Load Vector b */
    load_vector(b, file_b);

    printf("A: nrows=%d, ndim=%d neq=%d\n", A.nrow, A.ndim, A.neq);
    x.alloc(A.ncol);

    cc = solver_init();
    if (cc != 0) {
        printf("ERROR: init failed with %d\n", cc);
        exit(1);
    }

    SolverHandle_t handle = solver_create_handle();
    //cc = solver_set_option(handle, SOLVER_OPTION_SOLVER, SOLVER_ITER_BICGSTAB2);
    //cc = solver_set_option(handle, SOLVER_OPTION_SOLVER, SOLVER_ITER_CG);
    cc = solver_set_option(handle, SOLVER_OPTION_SOLVER, SOLVER_DIRECT_HS);
    //cc = solver_set_option(handle, SOLVER_OPTION_SOLVER, SOLVER_DIRECT_PARDISO);
    if (cc != 0) {
        printf("ERROR: set_option(SOLVER_OPTION_SOLVER) failed with %d\n", cc);
        exit(1);
    }

    TIMELOG_START(tl);
    cc = solver_set_matrix_csr(handle, A.nrow, A.ndim, A.pointers, A.indice, A.value, flags);
    if (cc != 0) {
        printf("ERROR: setMatrix failed with %d\n", cc);
        exit(1);
    }
    TIMELOG_END(tl, "setMatrix");

    TIMELOG_START(tl);
    cc = solver_solve(handle, b.value, x.value, res);
    if (cc != 0) {
        printf("ERROR: solve failed with %d\n", cc);
        exit(1);
    }
    TIMELOG_END(tl, "solve");

    /* Print the solution vector */
    printf("%s\n", "******** Solution ********");
    for (int i=0; i<5; i++) {
        printf("x[%d] = %14.12f\n", i, x.value[i]);
    }
    printf("%s\n", "********** End ***********");

    printf("TRUERESIDUAL: %e\n", solver_calc_residual(handle, b.value, x.value, 0));
    printf("TRUERESIDUAL(SCALED): %e\n", solver_calc_residual(handle, b.value, x.value, 1));

    cc = solver_free_handle(handle);

    solver_finalize();

    return 0;
}
