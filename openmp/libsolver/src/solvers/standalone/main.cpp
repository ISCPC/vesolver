#include <cstdio>
#include <stdint.h>
#include <stdlib.h>
#include "SpMatrix.hpp"
#include "timelog.h"
#include "Matrix.h"
#include "LinearSolver.hpp"

#define MAX_PATH_LEN 1024
#define DATA_SUFFIX ".bin"

uint32_t load_matrix(SpMatrix &A, char* file_a) {
    uint32_t flags=0;
    TIMELOG(tl1);

    /* Load Matrix A */
    bool oldflag = getenv("SPMATRIX_FORMAT_OLD") ? true : false;
    if (A.load_file(file_a, oldflag) < 0) {
        fprintf(stderr, "WARNING: cannot load matrix A. Trying old format.\n");
        if (A.load_file(file_a, true) < 0) {
            fprintf(stderr, "ERROR: cannot load matrix A from %s\n", file_a);
            exit(1);
        }
    }

#ifdef SORT_CSR
    TIMELOG_START(tl1);
    A.sort();
    TIMELOG_END(tl1, "Sorting");
#endif

    if ((A.type & 0xf) == SPMATRIX_TYPE_CSR) {
        flags |= MATRIX_TYPE_CSR;
        printf("INFO: Matrix_A: CSR, ");
    } else {
        flags |= MATRIX_TYPE_CSC;
        printf("INFO: Matrix_A: CSC, ");
    }
    if (A.type & SPMATRIX_TYPE_SYMMETRIC) {
        flags |= MATRIX_TYPE_SYMMETRIC;
        printf("SYMMETRIC, ");
    } else {
        flags |= MATRIX_TYPE_ASYMMETRIC;
        printf("ASYMMETRIC, ");
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
    LinearSolver* solver = getSolver();
    SpMatrix A, A0;
    Vector b, b0, x;
    char file_a[MAX_PATH_LEN], file_b[MAX_PATH_LEN];
    //double res = 1.0e-10;
    //double res = 1.0e-8;
    double res = 1.0e-6;
    //double res = 1.0e-4;
    //double res = 1.0e-1;
    int cc;
    TIMELOG(tl1);
    TIMELOG(tl2);
    unsigned int loop_num = 1;

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

    cc = solver->init();
    if (cc != 0) {
        printf("ERROR: init failed with %d\n", cc);
        exit(1);
    }

    char* env = getenv("NLOOPS");
    if (env != NULL) {
        loop_num = atoi(env);
        if (loop_num > 10) {
            loop_num = 10;
        }
    }
    for (int i=0; i<loop_num; i++) {
        TIMELOG_START(tl2);
        TIMELOG_START(tl1);
        cc = solver->setMatrixCSR(A.nrow, A.ndim, A.pointers, A.indice, A.value, flags);
        if (cc != 0) {
            printf("ERROR: setMatrix failed with %d\n", cc);
            exit(1);
        }
        TIMELOG_END(tl1, "setMatrix");

        TIMELOG_START(tl1);
        cc = solver->optimize();
        if (cc != 0) {
            printf("ERROR: optimize failed with %d\n", cc);
            exit(1);
        }
        TIMELOG_END(tl1, "optimize");

        printf("INFO: solving... (res=%le)\n", res);
        TIMELOG_START(tl1);
        cc = solver->solve(b.value, x.value, res);
        if (cc != 0) {
            printf("ERROR: solve failed with %d\n", cc);
            exit(1);
        }
        TIMELOG_END(tl1, "solve");
        TIMELOG_END(tl2, "overall");
    }

    /* Print the solution vector */
    printf("%s\n", "******** Solution ********");
    for (int i=0; i<5; i++) {
        printf("x[%d] = %14.12f\n", i, x.value[i]);
    }
    printf("%s\n", "********** End ***********");

    printf("TRUERESIDUAL: %e\n", solver->residual(b.value, x.value, 0));
    printf("TRUERESIDUAL(SCALED): %e\n", solver->residual(b.value, x.value, 1));
    cc = solver->finalize();

    return 0;
}
