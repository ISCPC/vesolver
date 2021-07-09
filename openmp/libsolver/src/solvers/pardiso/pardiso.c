#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"
#include "timelog.h"
#include "PluginAPI.h"

typedef int32_t INT_T;

int pardiso_(long long *pt, INT_T *maxfct, INT_T *mnum, INT_T *mtype, INT_T *phase,
        INT_T *n, double *a, INT_T *ia, INT_T *ja, INT_T *perm, INT_T *nrhs,
        INT_T *iparm, INT_T *msglvl, const double *b, double *x, INT_T *error);

static Matrix_t* solve_pre(const Matrix_t* A0) {
    Matrix_t* A = Matrix_duplicate(A0);
    Matrix_convert_index(A, 1);
    return A;
}

static int solve(const Matrix_t *A, const double* b, double* x, const double tolerance) {
    INT_T neq = A->NROWS;
    INT_T maxfct=1, mnum=1, nrhs=1, *perm=NULL, msglvl=0, error=0;
    INT_T mtype = MATRIX_IS_SYMMETRIC(A) ? -2 : 11;
    long long pt[64];
    INT_T iparm[64];
    TIMELOG(tl1);

    /* 
     * Factorize
     */
    INT_T phase=12;
    iparm[0]=0;
    iparm[1]=3;
    //iparm[34]=1;

    for(int i=0; i<64; i++) { pt[i]=0; }

    TIMELOG_START(tl1);
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &neq,
        A->values, A->pointers, A->indice, perm, &nrhs,
        iparm, &msglvl, b, x, &error);
    TIMELOG_END(tl1, "paradiso_factorize");

    if (error != 0) {
        printf("ERROR: pardiso_factorize (error=%d)\n", error);
        return -1;
    }

    /* 
     * Solve
     */
    phase=33;
    iparm[1]=3;
    iparm[5]=0;
    double* buf = (double*)calloc(sizeof(double), neq);

    TIMELOG_START(tl1);
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &neq,
        A->values, A->pointers, A->indice, perm, &nrhs,
        iparm, &msglvl, b, buf, &error);
    TIMELOG_END(tl1, "paradiso_solve");

    bcopy(buf, x, sizeof(double)*neq);
    free(buf);

    if (error != 0) {
        printf("ERROR: pardiso_solve (error=%d)\n", error);
        return -1;
    }

    /* 
     * Clean up
     */
    phase=-1;
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &neq,
        A->values, A->pointers, A->indice, perm, &nrhs,
        iparm, &msglvl, b, x, &error);
    if (error != 0) {
        printf("ERROR: pardiso_cleanup (error=%d)\n", error);
        return -1;
    }

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
SolverPlugin_t* pardiso_init() {
#endif
	SolverPlugin_t* solver = (SolverPlugin_t*)malloc(sizeof(SolverPlugin_t));

	solver->set_option = NULL;
	solver->solve_pre = solve_pre;
	solver->solve = solve;
	solver->solve_post = solve_post;
	solver->free = solve_free;

	return solver;
}
