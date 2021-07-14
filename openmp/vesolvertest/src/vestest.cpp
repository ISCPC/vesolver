/*
 * vestest.cpp
 *
 *  Created on: Jun 17, 2021
 *      Author: uno
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <strings.h>
#include "SpMatrix.h"
#include "timelog.h"
#include "vesolver.h"

/*
 * Main function
 */
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

    if (A.type & SPMATRIX_TYPE_CSR) {
        flags |= MATRIX_TYPE_CSR;
        printf("INFO: Matrix_A: CSR, ");
    } else {
        flags |= MATRIX_TYPE_CSC;
        printf("INFO: Matrix_A: CSC, ");
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

static inline void set_matrix(matrix_desc_t* desc, INT_T* pointers, INT_T* indice, double* values) {
    bcopy(pointers, desc->pointers, sizeof(INT_T)*(desc->neq+1));
    bcopy(indice, desc->indice, sizeof(INT_T)*(desc->nnz));
    bcopy(values, desc->values, sizeof(double)*(desc->nnz));
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

    vesolver_init();
    vesolver_handle_t hdl = vesolver_activate();
    if (hdl < 0) {
        printf("ERROR: vesolver_activate() failed.\n");
        exit(1);
    }

    vesolver_set_option(hdl, VESOLVER_OPTION_SOLVER, VESOLVER_ITER_CG);
    //vesolver_set_option(hdl, VESOLVER_OPTION_SOLVER, VESOLVER_ITER_BICGSTAB2);

    TIMELOG_START(tl);
    matrix_desc_t* desc = vesolver_alloc_matrix(hdl, A.nrow, A.ndim, flags);
    set_matrix(desc, A.pointers, A.indice, A.value);
    cc = vesolver_set_matrix(hdl, desc);
    TIMELOG_END(tl, "setMatrix");
    if (cc != 0) {
        printf("ERROR: vesolver_set_matrix() failed with %d\n", cc);
        exit(1);
    }

    TIMELOG_START(tl);
    cc = vesolver_solve_sync(hdl, b.value, x.value, res);
    TIMELOG_END(tl, "solve");
    if (cc != 0) {
        printf("ERROR: vesolver_solve_sync() failed with %d\n", cc);
        exit(1);
    }

    /* Print the solution vector */
    printf("%s\n", "******** Solution ********");
    for (int i=0; i<5; i++) {
        printf("x[%d] = %14.12f\n", i, x.value[i]);
    }
    printf("%s\n", "********** End ***********");

    //printf("TRUERESIDUAL: %e\n", solver_calc_residual(handle, b.value, x.value, 0));
    //printf("TRUERESIDUAL(SCALED): %e\n", solver_calc_residual(handle, b.value, x.value, 1));

    vesolver_free_matrix(hdl);

    vesolver_deactivate(hdl);
    vesolver_finalize();

    return 0;
}


