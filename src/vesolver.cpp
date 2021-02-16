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
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <cmath>
#include "vesolver_api.h"
#include "vesolver.hpp"
#include "timelog.h"
#ifdef ELMERSOLVER
#include "elmerSolver.h"
#endif
#ifdef HETEROSOLVER
#include <omp.h>
#include <heterosolver.h>
#endif
#ifdef WITH_PARDISO
#include <mkl.h>
#include <mkl_cluster_sparse_solver.h>
#endif

#include <sys/mman.h>
#include <errno.h>

//#define SANITY_CHECK 1

/*
 * functions for HeteroSolver
 */
//static int vesolver_hs_solve(VElib_Handle_t* hdl, SpMatrix_t *A, double* b, double* x, double* res)
int VESolver::hs_solve(SpMatrix& A, Vector& b, Vector& x, double res) {
#ifdef HETEROSOLVER
    /* Local variables */
    HS_int_t nrows = A.nrow;
    int ierr;

    /* Handle Initialization */
    HS_handle_t hnd; 

    const HS_int_t* pointers = A.pointers;
    const HS_int_t* indice = A.indice;
    const double* value = A.value;

#if 0
    HS_int_t isym = SPMATRIX_IS_SYMMETRIC(A) ? HS_SYMMETRIC : HS_UNSYMMETRIC;
    HS_int_t iformat = SPMATRIX_IS_CSR(A) ? HS_CSR : HS_CSC;
#else
    HS_int_t isym = HS_UNSYMMETRIC;
    HS_int_t iformat = HS_DCSR;
#endif
    
    printf("INFO:VESolver[%d]:Solving the system of equations using the Heterosolver on VE.\n\n", myrank);

    MPI_Bcast(&nrows, 1, MPI_INTEGER, 0, solver_comm);
    //printf("INFO: A: type=0x%x, isym=%d, iformat=%d, nrows=%d\n", isym, iformat, A.type, nrows);

    ierr = PHS_init_handle(&hnd, nrows, nrows, isym, iformat, HS_MASTER, 1, A.nrow, solver_comm);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: HS_init_handle failed with %d.\n", ierr);
        exit(1);
    }

    /* Specifying the number of OpenMP threads */
    //omp_set_num_threads(1);

    /* Specifying One-based indexing */
    //if (SPMATRIX_INDEX_TYPE(A) == 1) {
        ierr = PHS_set_option(hnd, HS_INDEXING, HS_INDEXING_1);
        if (ierr != HS_RESULT_OK) {
            fprintf(stderr, "ERROR: HS_set_option failed with %d.\n", ierr);
            exit(1);
        }
    //}

    /*
    PHS_set_option(A % HShandle, HS_INDEXING, HS_INDEXING_1, ierror)
    PHS_set_option(A % HShandle, HS_OUTPUT, HS_OUTPUT_B, ierror)
    PHS_set_option(A % HShandle, HS_DUPLICATE, HS_DUPLICATE_YES, ierror)
    PHS_set_option(A % HShandle, HS_ORDP, HS_ORDP_METIS, ierror)
    PHS_set_option(A % HShandle, HS_ORDQ, HS_ORDQ_STATIC, ierror)
    !PHS_set_option(A % HShandle, HS_PERTURBATION, 10, ierror)
    PHS_set_option(A % HShandle, HS_SCALING, HS_SCALING_ON, ierror)
    !PHS_set_option(A % HShandle, HS_NORM, HS_NORM_L2, ierror)
    */

    /* Preprocessing Phase */
    ierr = PHS_preprocess_rd(hnd, pointers, indice, value);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: PHS_preprocess_rd failed with %d.\n", ierr);
        exit(1);
    }

    /* Numeric Factorization Phase */
    ierr = PHS_factorize_rd(hnd, pointers, indice, value);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: PHS_factorize_rd failed with %d.\n", ierr);
        exit(1);
    }

    /* Solution Phase */
    ierr = PHS_solve_rd(hnd, pointers, indice, value, 1, b.value, x.value, &res);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: PHS_solve_rd failed with %d.\n", ierr);
        exit(1);
    }

    /* Handle Finalization */
    ierr = PHS_finalize_handle(hnd);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: PHS_finalize_handle failed with %d.\n", ierr);
        exit(1);
    }

    return 0;
#else
    printf("ERROR: HeterSolver is not enabled.\n");
    return -1;
#endif
}

int VESolver::pardiso_solve(SpMatrix& A, Vector& b, Vector& x, double res) {
#ifdef WITH_PARDISO
    long long pt[64];
    const MKL_INT maxfct=1;
    const MKL_INT mnum=1;
    const MKL_INT mtype=11;
    MKL_INT phase;
    MKL_INT perm=0;
    const MKL_INT nrhs=1;
    MKL_INT iparm[64];
    const MKL_INT msglvl=0;
    MKL_INT ierr;
    int comm = MPI_Comm_c2f(solver_comm);

    printf("INFO:VESolver[%d]:Solving the system of equations using the PARDISO.\n\n", myrank);
    //printf("INFO: A: type=0x%x, nrow=%d, neq=%d\n", A.type, A.nrow, A.neq);

    // Set up parameters explicitly
    for(int i=0;i<64;i++) { pt[i]=0; }
    for(int i=0; i<64; i++) { iparm[i] = 0; }

    /* Perform analysis */
#if 1
    iparm[0]=0;
    iparm[1]=3;
    iparm[5]=0;
    iparm[39]=0;

    phase = 11;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &comm, &ierr);

    /* Perform factorization */
    phase = 22;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &comm, &ierr);

    /* Perform solve */
    phase = 33;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &comm, &ierr);

    /* Finalization */
    phase = -1;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &comm, &ierr);
#else
    iparm[0]=1;
    iparm[1]=2;
    iparm[5]=0;
    iparm[10]=0;
    iparm[11]=0;
    iparm[17]=-1;
    iparm[26]=1;
    iparm[34]=0;

    phase = 11;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &ierr);

    phase = 22;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &ierr);

    phase = 33;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &ierr);

    phase = -1;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &ierr);
#endif

    return 0;
#else /* WITH_PARDISO */
    printf("ERROR: Pardiso is not enabled.\n");
    return -1;
#endif /* WITH_PARDISO */
}

/*
 * functions for Elmer Solver
 */
int VESolver::elmer_solve(SpMatrix& A, Vector& b, Vector& x, double res) {
#ifdef ELMERSOLVER
    printf("INFO:VESolver[%d]: Solving the system of equations using the elmer iterative solver on VE.\n\n", myrank);

    return elmersolver(A.ncol, A.ndim, A.pointers, A.indice, A.value, b.value, x.value, ELMER_BICGSTAB2);
#else
    printf("ERROR: elmerSolver is not enabled.\n");
    return -1;
#endif
}

/*
 * functions for Dummy Solver (for testing)
 */
int VESolver::dummy_solve(SpMatrix& A, Vector& b, Vector& x, double res) {
    printf("INFO:VESolver[%d]: Dummy solver called.\n", myrank);

#if 0
    for(int i=0; i<5; i++) {
    //for(int i=0; i<A.nrow; i++) {
        for(int j=A.pointers[i]-1; j<A.pointers[i+1]-1; j++) {
            printf("A[%10d, %10d] = %20.16le (ptr:%10d)\n", i+1, A.indice[j], A.value[j], A.pointers[i+1]);
        }
    }

    for(int i=0; i<5; i++) {
    //for(int i=0; i<A.neq; i++) {
        printf("b[%10d] = %20.16le\n", i+1, b.value[i]);
    }
#endif
    for(int i=0; i<A.ncol; i++) {
        x.value[i] = 0.0;
    }

    return 0;
}

int VESolver::solve(int solver, SpMatrix& A, Vector& b, Vector& x, double res) {
    int cc;
    TIMELOG(tl);

    TIMELOG_START(tl);
    switch(solver) {
        case VESOLVER_HS:
#ifdef HETEROSOLVER
            cc = hs_solve(A, b, x, res);
#else
            cc = pardiso_solve(A, b, x, res);
#endif
            break;

        case VESOLVER_BICGSTAB2:
            cc = elmer_solve(A, b, x, res);
            break;

        case VESOLVER_DUMMY:
            cc = dummy_solve(A, b, x, res);
            break;

        default:
            fprintf(stderr, "ERROR: Invalid solver type (solver=%d)\n", solver);
            cc = -1;
    }
    TIMELOG_END(tl, "vesolver_solve");

    return cc;
}

int VESolver::solve(int solver, SpDistMatrix **An, DistVector **bn, int n, Vector& x, double res) {
    int cc=-1;
    TIMELOG(tl);

    TIMELOG_START(tl);
    SpMatrix *A = SpDistMatrix::gather(An, n);
    TIMELOG_END(tl, "gather_A");

#ifdef SANITY_CHECK
    printf("INFO:VESolver[%d]:solve: A: nrow=%d, ncol=%d, neq=%d, ndim=%d, %d, %d\n", 
        myrank, A->nrow, A->ncol, A->neq, A->ndim, A->pointers[A->neq-1], A->pointers[A->neq]);
    for (int i=0; i<A->neq; i++) {
        if ((A->pointers[i] <= 0)||(A->pointers[i] > A->ndim+1)) {
            printf("INFO:VESolver[%d]:solve: Invalid pointer A[%d]=%d\n", myrank, i, A->pointers[i]);
        }

        int k = A->pointers[i]-1;
        int prev = A->indice[k];
        for (int j=A->pointers[i]; j<A->pointers[i+1]-1; j++) {
            if (A->indice[j] <= prev) {
                printf("INFO:VESolver[%d]:solve: Invalid indice A[%d, %d]=%d < %d\n", myrank, i, j, A->indice[j], prev);
            }
            prev = A->indice[j];
        }
    }
    for (int i=0; i<A->ndim; i++) {
        if ((A->indice[i] <= 0)||(A->indice[i] > A->ncol)) {
            printf("INFO:VESolver[%d]:solve: Invalid indice A[%d]=%d\n", myrank, i, A->indice[i]);
        }
    }
#endif

    TIMELOG_START(tl);
    Vector *b = DistVector::gather(bn, n, An);
    TIMELOG_END(tl, "gather_B");

    cc = solve(solver, *A, *b, x, res);

    /*
     * free Matrix/Vector
     */
    delete A;
    delete b;

    return cc;
}

int VESolver::hs_solve(SpDistMatrix& A, DistVector& b, Vector& x, double res) {
#ifdef HETEROSOLVER
    /* Local variables */
    int ierr;

    /* Handle Initialization */
    HS_handle_t hnd; 

    const HS_int_t* pointers = A.pointers;
    const HS_int_t* indice = A.indice;
    const double* value = A.value;

#if 0
    HS_int_t isym = SPMATRIX_IS_SYMMETRIC(A) ? HS_SYMMETRIC : HS_UNSYMMETRIC;
    HS_int_t iformat = SPMATRIX_IS_CSR(A) ? HS_CSR : HS_CSC;
#else
    HS_int_t isym = HS_UNSYMMETRIC;
    HS_int_t iformat = HS_DCSR;
#endif
    
    printf("INFO:VESovler[%d]: Solving the system of equations using the Heterosolver on VE.\n\n", myrank);
    //printf("INFO: A: type=0x%x, nrow=%d, neq=%d, mingind=%d, maxgind=%d\n", A.type, A.nrow, A.neq, A.mingind, A.maxgind);

    ierr = PHS_init_handle(&hnd, A.neq, A.neq, isym, iformat, HS_DIST_AB, A.mingind, A.maxgind, solver_comm);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: HS_init_handle failed with %d.\n", ierr);
        exit(1);
    }

    /* Specifying the number of OpenMP threads */
    //omp_set_num_threads(1);

    /* Specifying One-based indexing */
    //if (SPMATRIX_INDEX_TYPE(A) == 1) {
        ierr = PHS_set_option(hnd, HS_INDEXING, HS_INDEXING_1);
        if (ierr != HS_RESULT_OK) {
            fprintf(stderr, "ERROR: HS_set_option failed with %d.\n", ierr);
            exit(1);
        }
    //}

    ierr = PHS_set_option(hnd, HS_OUTPUT, HS_OUTPUT_X);
    ierr = PHS_set_option(hnd, HS_DUPLICATE, HS_DUPLICATE_YES);
    ierr = PHS_set_option(hnd, HS_ORDP, HS_ORDP_METIS);
    ierr = PHS_set_option(hnd, HS_ORDQ, HS_ORDQ_STATIC);
    //ierr = PHS_set_option(hnd, HS_PERTURBATION, 10);
    ierr = PHS_set_option(hnd, HS_SCALING, HS_SCALING_ON);
    //ierr = PHS_set_option(hnd, HS_NORM, HS_NORM_L2);

    /* Preprocessing Phase */
    ierr = PHS_preprocess_rd(hnd, pointers, indice, value);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: PHS_preprocess_rd failed with %d.\n", ierr);
        exit(1);
    }

    /* Numeric Factorization Phase */
    ierr = PHS_factorize_rd(hnd, pointers, indice, value);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: PHS_factorize failed with %d.\n", ierr);
        exit(1);
    }

    /* Solution Phase */
    ierr = PHS_solve_rd(hnd, pointers, indice, value, 1, b.value, x.value, &res);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: PHS_solve_rd failed with %d (res=%le).\n", ierr, res);
        exit(1);
    }

    /* Handle Finalization */
    ierr = PHS_finalize_handle(hnd);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: PHS_finalize_handle failed with %d.\n", ierr);
        exit(1);
    }

    return 0;
#else
    printf("ERROR: HeterSolver is not enabled.\n");
    return -1;
#endif
}

/*
 void cluster_sparse_solver (_MKL_DSS_HANDLE_t pt, const MKL_INT *maxfct, const MKL_INT *mnum, const MKL_INT *mtype, const MKL_INT *phase, const MKL_INT *n, const void *a, const MKL_INT *ia, const MKL_INT *ja, MKL_INT *perm, const MKL_INT *nrhs, MKL_INT *iparm, const MKL_INT *msglvl, void *b, void *x, const int *comm, MKL_INT *error);
*/
int VESolver::cpardiso_solve(SpDistMatrix& A, DistVector& b, Vector& x, double res) {
#ifdef WITH_PARDISO
    long long pt[64];
    const MKL_INT maxfct=1;
    const MKL_INT mnum=1;
    const MKL_INT mtype=11;
    MKL_INT phase;
    MKL_INT perm=0;
    const MKL_INT nrhs=1;
    MKL_INT iparm[64];
    const MKL_INT msglvl=0;
    MKL_INT ierr;
    int comm = MPI_Comm_c2f(solver_comm);

    printf("INFO:VESolver[%d]: Solving the system of equations using the Cluster PARDISO.\n\n", myrank);
    //printf("INFO: A: type=0x%x, nrow=%d, neq=%d, mingind=%d, maxgind=%d\n", A.type, A.nrow, A.neq, A.mingind, A.maxgind);

    // Set up parameters explicitly
    for(int i=0;i<64;i++) { pt[i]=0; }

    for(int i=0; i<64; i++) {
        iparm[i] = 0;
    }
#if 1
    iparm[0]=1;             // Do not use = 1 / use = 0 solver default parameters
    iparm[1]=2;             // Minimize fill-in with OpenMP nested dissection
    iparm[4]=0;             // No user input permutation
    iparm[5]=0;             // Write solution vector to x
    iparm[7]=0;             // Number of iterative refinement steps
    iparm[9]=13;            // Perturbation value 10^-iparm(10) in case of small pivots
    iparm[10]=1;            // Use scalings from symmetric weighted matching
    iparm[12]=1;            // Use permutations from nonsymmetric weighted matching

    //iparm[17]=-1;         // Output: Number of nonzeros in the factor LU
    //iparm[18]=-1;         // Output: Mflops for LU factorization
    iparm[20]=1;            // Do not use Bunch Kaufman pivoting
    iparm[26]=0;            // Do not check sparse matrix representation
    iparm[27]=0;            // Use double precision
    iparm[34]=0;            // Use Fortran indexing

    /* CPardiso matrix input format */
    iparm[39] = 3;          // Distributed solution phase, distributed solution vector
    iparm[40] = A.mingind;  // Beginning of solution domain
    iparm[41] = A.maxgind;  // End of solution domain
#else
    iparm[0]=0;
    iparm[1]=3;
    iparm[17]=-1;         // Output: Number of nonzeros in the factor LU
    iparm[39]=3;
    iparm[40]=A.mingind;
    iparm[41]=A.maxgind;
#endif

    /* Perform analysis */
    phase = 11;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &comm, &ierr);
    if (ierr != 0) {
        printf("WARNING: cluster_sparse_solver failed (phase=%d, ierr=%d)\n",
            phase, ierr);
    }

    /* Perform factorization */
    phase = 22;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &comm, &ierr);
    if (ierr != 0) {
        printf("WARNING: cluster_sparse_solver failed (phase=%d, ierr=%d)\n",
            phase, ierr);
    }

    /* Perform solve */
    phase = 33;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &comm, &ierr);
    if (ierr != 0) {
        printf("WARNING: cluster_sparse_solver failed (phase=%d, ierr=%d)\n",
            phase, ierr);
    }
#if 1
    /* Finalization */
    phase = -1;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &A.neq, A.value, A.pointers, A.indice,
        &perm, &nrhs, iparm, &msglvl,
        b.value, x.value, &comm, &ierr);
    if (ierr != 0) {
        printf("WARNING: cluster_sparse_solver failed (phase=%d, ierr=%d)\n",
            phase, ierr);
    }
#endif
    return 0;
#else /* WITH_PARDISO */
    printf("ERROR: CPardiso is not enabled.\n");
    return -1;
#endif /* WITH_PARDISO */
}

int VESolver::elmer_solve(SpDistMatrix& A, DistVector& b, Vector& x, double res) {
#ifdef ELMERSOLVER
    printf("INFO:VESolver[%d]: Solving the system of equations using the elmer iterative solver on VE.\n\n", myrank);

    return elmersolver(A.ncol, A.ndim, A.pointers, A.indice, A.value, b.value, x.value, ELMER_BICGSTAB2);
#else
    printf("ERROR: elmerSolver is not enabled.\n");
    return -1;
#endif
}

int VESolver::dummy_solve(SpDistMatrix& A, DistVector& b, Vector& x, double res) {
    printf("INFO:vesolver[%d]: Dummy solver called.\n", myrank);

#if 0
    for(int i=0; i<5; i++) {
    //for(int i=0; i<A.nrow; i++) {
        for(int j=A.pointers[i]-1; j<A.pointers[i+1]-1; j++) {
            printf("%d:A[%10d, %10d] = %20.16le (ptr:%10d)\n", rank, i+1, A.indice[j], A.value[j], A.pointers[i+1]);
        }
    }

    for(int i=0; i<5; i++) {
    //for(int i=0; i<A.neq; i++) {
        printf("%d:b[%10d] = %20.16le\n", rank, i+1, b.value[i]);
    }
#endif
    for(int i=0; i<A.neq; i++) {
        x.value[i] = 0.0;
    }

    return 0;
}

int VESolver::solve(int solver, SpDistMatrix& A, DistVector& b, Vector& x, double res) {
    int cc;
    TIMELOG(tl);

#ifdef SANITY_CHECK
    for (int i=0; i<A.ndim; i++) {
        if (std::isnan(A.value[i])) {
            printf("WARNING:VESolver[%d]::solve: Detect A[%d] == NaN (%lf)\n", myrank, i, A.value[i]);
            break;
        }
    }
    for (int i=0; i<b.size; i++) {
        if (std::isnan(b.value[i])) {
            printf("WARNING:VESolver[%d]::solve: Detect b[%d] == NaN (%lf)\n", myrank, i, b.value[i]);
            break;
        }
    }
#endif

    TIMELOG_START(tl);
    switch(solver) {
        case VESOLVER_HS:
            A.ConvertToDCSR();
#ifdef HETEROSOLVER
            cc = hs_solve(A, b, x, res);
#else
            cc = cpardiso_solve(A, b, x, res);
#endif
            break;

        case VESOLVER_BICGSTAB2:
            cc = elmer_solve(A, b, x, res);
            break;

        case VESOLVER_DUMMY:
            cc = dummy_solve(A, b, x, res);
            break;

        default:
            fprintf(stderr, "ERROR: Invalid solver type (solver=%d)\n", solver);
            cc = -1;
    }
    TIMELOG_END(tl, "vesolver_solve");

#ifdef SANITY_CHECK
    for (int i=0; i<x.size; i++) {
        if (std::isnan(x.value[i])) {
            printf("WARNING:VESolver[%d]::solve: Detect x[%d] == NaN (%lf)\n", myrank, i, x.value[i]);
            break;
        }
    }
#endif
    return cc;
}

/*
 *  procss
 */
inline int VESolver::receive_matrix_info(int src, SpMatrix* A) {
    vesolver_matrix_info_t info;
    MPI_Status status;
    int cc;

    cc = MPI_Recv(&info, sizeof(vesolver_matrix_info_t), MPI_BYTE, src, VES_MATINFO_TAG, ves_comm, &status);
    A->ndim = info.nval;
    A->type = info.type;
    A->nrow = info.nrow;
    A->ncol = info.ncol;
    A->neq  = info.neq;
    A->mingind = info.nl;
    A->maxgind = info.nt;
    printf("INFO:VESolver[%d]::receive_matrix_info: Receive matrix info from Rank#%d (0x%04x %d %d %d %d %d)\n",
        myrank, src, info.type, info.nrow, info.nval, info.neq, info.nl, info.nt);

    return cc;
}

inline int VESolver::receive_matrix_info(int src, SpDistMatrix* A) {
    return receive_matrix_info(src, (SpMatrix*)A);
}

inline int VESolver::receive_matrix_data(int src, SpMatrix* A, Vector* b) {
    int cc;
    MPI_Status status;

    printf("INFO:VESolver[%d]::vesolver_receive_matrix_data: Receive request from rank %d with %d %d %d %d\n",
        myrank, src, A->nrow, A->ncol, A->ndim, A->neq);

    if (A->alloc(A->nrow, A->ndim)) {
        printf("ERROR: Memory allocation error.(%s:%d)\n", __FILE__, __LINE__);
        return -1;
    }
    if (b->alloc(A->nrow)) {
        printf("ERROR: Memory allocation error.(%s:%d)\n", __FILE__, __LINE__);
        return -1;
    }

    // Receive Matrix A
    cc = MPI_Recv(A->value, A->ndim, MPI_DOUBLE, src, VES_MATDATA_TAG, ves_comm, &status);
    if (cc != 0) {
        printf("ERROR: MPI_Recv failed.(%s:%d)\n", __FILE__, __LINE__);
    }
    //printf("INFO:Receive:value status(count=%d, cancelled=%d, source=%d, tag=%d, error=%d\n",
    //    status.count, status.cancelled, status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);

    cc = MPI_Recv(A->indice, A->ndim, MPI_INTEGER, src, VES_MATDATA_TAG, ves_comm, &status);
    if (cc != 0) {
        printf("ERROR: MPI_Recv failed.(%s:%d)\n", __FILE__, __LINE__);
    }
    //printf("INFO:Receive:indice status(count=%d, cancelled=%d, source=%d, tag=%d, error=%d\n",
    //    status.count, status.cancelled, status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);

    cc = MPI_Recv(A->pointers, (A->nrow)+1, MPI_INTEGER, src, VES_MATDATA_TAG, ves_comm, &status);
    if (cc != 0) {
        printf("ERROR: MPI_Recv failed.(%s:%d)\n", __FILE__, __LINE__);
    }
    //printf("INFO:Receive:pointers status(count=%d, cancelled=%d, source=%d, tag=%d, error=%d\n",
    //    status.count, status.cancelled, status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);

    // Recieve Matrix b
    cc = MPI_Recv(b->value, A->nrow, MPI_DOUBLE, src, VES_MATDATA_TAG, ves_comm, &status);
    if (cc != 0) {
        printf("ERROR: MPI_Recv failed.(%s:%d)\n", __FILE__, __LINE__);
    }
    //printf("INFO:Receive:b status(count=%d, cancelled=%d, source=%d, tag=%d, error=%d\n",
    //    status.count, status.cancelled, status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);

#if 0
    for(int i=0; i<5; i++) {
        printf("INFO:Receive:Ai[%d] = %d %d %lf %lf\n", i, A->pointers[i],
            A->indice[i], A->value[i], b->value[i]);
    }
    for(int i=A->nrow-5; i<A->nrow; i++) {
        printf("INFO:Receive:Ai[%d] = %d %d %lf %lf\n", i, A->pointers[i],
            A->indice[i], A->value[i], b->value[i]);
    }
    for(int i=A->ndim-5; i<A->ndim; i++) {
        printf("INFO:Receive:Aj[%d] = (%d, %lf)\n", i, A->indice[i], A->value[i]);
    }
#endif

    return cc;
}

inline int VESolver::receive_matrix_data(int src, SpDistMatrix* A, DistVector* b) {
    int cc;
    MPI_Status status;

    if (A->alloc(A->nrow, A->ndim)) {
        printf("ERROR: Memory allocation error.(%s:%d)\n", __FILE__, __LINE__);
        return -1;
    }

    A->rorder = (INT_T*)calloc(A->neq, sizeof(INT_T));
    if (A->rorder == NULL) {
        printf("ERROR: Memory allocation error.(%s:%d)\n", __FILE__, __LINE__);
        return -1;
    }

    if (b->alloc(A->nrow)) {
        printf("ERROR: Memory allocation error.(%s:%d)\n", __FILE__, __LINE__);
        return -1;
    }

    // Receive Matrix A
    printf("INFO:VESolver[%d]::vesolver_receive_matrix_data: Receive request from rank %d with %d %d %d %d\n",
        myrank, src, A->nrow, A->ncol, A->ndim, A->neq);
    cc = MPI_Recv(A->value, A->ndim, MPI_DOUBLE, src, VES_MATDATA_TAG, ves_comm, &status);
    if (cc != 0) {
        printf("ERROR: MPI_Recv failed.(%s:%d)\n", __FILE__, __LINE__);
    }
    cc = MPI_Recv(A->indice, A->ndim, MPI_INTEGER, src, VES_MATDATA_TAG, ves_comm, &status);
    if (cc != 0) {
        printf("ERROR: MPI_Recv failed.(%s:%d)\n", __FILE__, __LINE__);
    }
    cc = MPI_Recv(A->pointers, A->nrow+1, MPI_INTEGER, src, VES_MATDATA_TAG, ves_comm, &status);
    if (cc != 0) {
        printf("ERROR: MPI_Recv failed.(%s:%d)\n", __FILE__, __LINE__);
    }
    cc = MPI_Recv(A->rorder, A->neq, MPI_INTEGER, src, VES_MATDATA_TAG, ves_comm, &status);
    if (cc != 0) {
        printf("ERROR: MPI_Recv failed.(%s:%d)\n", __FILE__, __LINE__);
    }

#if 0
    for(int i=0; i<5; i++) {
        printf("INFO:A[%d] = (%d, %d, %lf)\n", i, A->pointers[i], A->indice[i], A->value[i]);
    }
    fflush(stdout);
#endif
    // Recieve Matrix b
    cc = MPI_Recv(b->value, A->nrow, MPI_DOUBLE, src, VES_MATDATA_TAG, ves_comm, &status);

    return cc;
}

inline int VESolver::send_result(Vector& x) {
    return MPI_Send(x.value, x.size, MPI_DOUBLE, 0, VES_MATDATA_TAG, ves_comm);
}

inline int VESolver::send_result(int dest, Vector& x) {
    return MPI_Send(x.value, x.size, MPI_DOUBLE, dest, VES_MATDATA_TAG, ves_comm);
}

int VESolver::process_gather_on_vh(VES_request_t& req) {
    SpMatrix A;
    Vector b, x;
    int cc = 0;
    vesolver_request_t *payload = (vesolver_request_t*)VES_REQ_GETPAYLOAD(req);

    printf("INFO:VESolver[%d]::process: solver=%d, type=0x%08x, nprocs=%d, res=%lf.\n",
        myrank, payload->solver, payload->type, payload->nprocs, payload->res);

    receive_matrix_info(0, &A);
    receive_matrix_data(0, &A, &b);
    if (x.alloc(A.ncol) < 0) {
        printf("ERROR:VESolver::process_gather_on_vh: allocation error. [%s:%d]\n", __FILE__, __LINE__);
        return -1;
    } 

#ifdef SANITY_CHECK
    int nerr=0;
    printf("INFO:VESolver[%d]:process_gather_on_vh: A: nrow=%d, ncol=%d, neq=%d, ndim=%d, %d, %d\n", 
        myrank, A.nrow, A.ncol, A.neq, A.ndim, A.pointers[A.neq-1], A.pointers[A.neq]);
    for (int i=0; i<A.neq; i++) {
        if ((A.pointers[i] <= 0)||(A.pointers[i] > A.ndim+1)) {
            printf("INFO:VESolver[%d]:process_gather_on_vh: Invalid pointer A[%d]=%d\n", myrank, i, A.pointers[i]);
            nerr++;
        }

        int k = A.pointers[i]-1;
        int prev = A.indice[k];
        for (int j=A.pointers[i]; j<A.pointers[i+1]-1; j++) {
            if (A.indice[j] <= prev) {
                printf("INFO:VESolver[%d]:process_gather_on_vh: Invalid indice A[%d, %d]=%d < %d\n", myrank, i, j, A.indice[j], prev);
                nerr++;
            }
            prev = A.indice[j];
        }
        if (nerr > 100) {
            send_result(x);
            return -1;
        }
    }
    for (int i=0; i<A.ndim; i++) {
        if ((A.indice[i] <= 0)||(A.indice[i] > A.ncol)) {
            printf("INFO:VESolver[%d]:process_gather_on_vh: Invalid indice A[%d]=%d\n", myrank, i, A.indice[i]);
            nerr++;
        }
    }
#endif

    cc = solve(payload->solver, A, b, x, payload->res);
    send_result(x);
    return cc;
}

int VESolver::process_gather_on_ve(VES_request_t& req) {
    int cc=0;
    SpDistMatrix **An;
    DistVector **bn;
    Vector x;
    vesolver_request_t *payload = (vesolver_request_t*)VES_REQ_GETPAYLOAD(req);
    int nprocs = payload->nprocs;

    An = (SpDistMatrix**)calloc(sizeof(SpMatrix*), nprocs);
    bn = (DistVector**)calloc(sizeof(Vector*), nprocs);

    for(int i=0; i<nprocs; i++) {
        An[i] = new SpDistMatrix();
        receive_matrix_info(i, An[i]);

        bn[i] = new DistVector();
        receive_matrix_data(i, An[i], bn[i]);

        printf("INFO:%d: size: %d x %d, nonzero: %d, offset: %d, neq: %d, size_b: %ld\n",
            i, An[i]->ncol, An[i]->nrow, An[i]->ndim, An[i]->offset, An[i]->neq, bn[i]->size);
    }

    x.alloc(An[0]->neq);
    cc = solve(payload->solver, An, bn, payload->nprocs, x, payload->res);
    send_result(x);

    /*
     * free Matrix/Vector
     */
    for(int i=0; i<nprocs; i++) {
        delete An[i]; 
        delete bn[i]; 
    }
    free(An);
    free(bn);

    return cc;
}

int VESolver::process_symmetric(VES_request_t& req) {
    int cc=0;
    SpDistMatrix A;
    DistVector b;
    Vector x;
    vesolver_request_t *payload = (vesolver_request_t*)VES_REQ_GETPAYLOAD(req);

    printf("INFO:VESolver[%d]::process: solver=%d, type=0x%08x, nprocs=%d, res=%lf.\n",
        myrank, payload->solver, payload->type, payload->nprocs, payload->res);

    receive_matrix_info(myrank, &A);
    receive_matrix_data(myrank, &A, &b);
    if (x.alloc(A.neq) < 0) {
        printf("ERROR:VESolver[%d]::process: Memory allocation failed\n", myrank);
        return -1;
    }
    cc = solve(payload->solver, A, b, x, payload->res);
    send_result(myrank, x);

    return cc;
}

int VESolver::process(VES_request_t& req) {
    int cc=-1;
    vesolver_request_t *payload = (vesolver_request_t*)VES_REQ_GETPAYLOAD(req);

    /* Call solver */
    switch(payload->mode) {
        case VES_MODE_GATHER_ON_VH:
            cc = process_gather_on_vh(req);
            break;

        case VES_MODE_GATHER_ON_VE:
            cc = process_gather_on_ve(req);
            break;

        case VES_MODE_SYMMETRIC:
            cc = process_symmetric(req);
            break;
    }
    fflush(stdout);

    return cc;
}

/*
 * Functions for Stand-Alone test
 */
int test_gather_on_vh(VESolver& server, Vector& x, double res, int solver) {
    SpMatrix A;
    Vector b;
    char infile[256];
    const char* file_a = "a.bin";
    const char* file_b = "b.bin";
    int cc = -1;
    int rank;

    const char* mat_path = getenv("VESOLVER_DATA_PATH");
    if (mat_path == NULL) {
        mat_path = ".";
    }
    printf("INFO:test_gather_on_vh: the path to matrix and vector: %s\n", mat_path);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        snprintf(infile, 256, "%s/%s", mat_path, file_a);
        if (A.load_file(infile) < 0) {
            fprintf(stderr, "ERROR: cannot load matrix A from %s\n", file_a);
            exit(1);
        }
        snprintf(infile, 256, "%s/%s", mat_path, file_b);
        if (b.load_file(infile) < 0) {
            fprintf(stderr, "ERROR: cannot load vector b from %s\n", file_b);
            exit(1);
        }
        x.alloc(A.ncol);
    }
    cc = server.solve(solver, A, b, x, res);
    return cc;
}

int test_gather_on_ve(VESolver& server, Vector& x, double res, int solver) {
    int cc=0;
    SpDistMatrix **An;
    DistVector **bn;
    char infile[256];
    int nprocs=8;
    const char *prefix_a="a";
    const char *prefix_b="b";

    const char* mat_path = getenv("VESOLVER_DATA_PATH");
    if (mat_path == NULL) {
        mat_path = ".";
    }
    printf("INFO:test_gather_on_ve: the path to matrix and vector: %s\n", mat_path);

    char* env = getenv("VESOLVER_NPARA");
    if (env != NULL) {
        nprocs = atoi(env);
        printf("INFO:test_gather_on_ve: set the number of parallels = %d\n", nprocs);
    }

    An = (SpDistMatrix**)calloc(sizeof(SpMatrix*), nprocs);
    bn = (DistVector**)calloc(sizeof(Vector*), nprocs);
    for(int i=0; i<nprocs; i++) {
        An[i] = new SpDistMatrix();
        snprintf(infile, 256, "%s/%s%d.bin", mat_path, prefix_a, i);
        if (An[i]->load_file(infile) < 0) {
            fprintf(stderr, "ERROR: cannot load matrix A%d from %s\n", i, infile);
            exit(1);
        }
        An[i]->reordering();

        bn[i] = new DistVector();
        snprintf(infile, 256, "%s/%s%d.bin", mat_path, prefix_b, i);
        if (bn[i]->load_file(infile) < 0) {
            fprintf(stderr, "ERROR: cannot load vector b%d from %s\n", i, infile);
            exit(1);
        }

        printf("INFO:%d:%s: size: %d x %d, nonzero: %d, offset: %d, neq: %d, size_b: %ld\n",
            i, infile, An[i]->ncol, An[i]->nrow, An[i]->ndim, An[i]->offset, An[i]->neq, bn[i]->size);
    }

    /* Allocate Vector x */
    x.alloc(An[0]->neq);

    cc = server.solve(solver, An, bn, nprocs, x, res);

    /*
     * free Matrix/Vector
     */
    for(int i=0; i<nprocs; i++) {
        delete An[i]; 
        delete bn[i]; 
    }
    free(An);
    free(bn);

    return cc;
}

int test_symmetric(VESolver& server, Vector& x, double res, int solver) {
    int rank, cc=0;
    char infile[256];
    const char *prefix_a="a";
    const char *prefix_b="b";
    SpDistMatrix A;
    DistVector b;

    const char* mat_path = getenv("VESOLVER_DATA_PATH");
    if (mat_path == NULL) {
        mat_path = ".";
    }
    printf("INFO:test_symmetric: the path to matrix and vector: %s\n", mat_path);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Load Matrix A */
    snprintf(infile, 256, "%s/%s%d.bin", mat_path, prefix_a, rank);
    if (A.load_file(infile) < 0) {
        fprintf(stderr, "ERROR: cannot load matrix A from %s\n", infile);
        exit(1);
    }
    A.reordering();

    /* Load Vector b */
    snprintf(infile, 256, "%s/%s%d.bin", mat_path, prefix_b, rank);
    if (b.load_file(infile) < 0) {
        fprintf(stderr, "ERROR: cannot load vector b from %s\n", infile);
        exit(1);
    }

    /* Allocate Vector x */
    if (x.alloc(A.neq) < 0) {
        printf("ERROR:VESolver[%d]::process: Memory allocation failed\n", rank);
        exit(1);
    }

    printf("INFO:test_symmetric[%d]: A.nrow=%d A.neq=%d A.mingind=%d A.maxgind=%d\n", 
        rank, A.nrow, A.neq, A.mingind, A.maxgind);
#if 0
    for(int i=0; i<5; i++) {
        for(int j=A.pointers[i]-1; j<A.pointers[i+1]-1; j++) {
            printf("INFO:test_symmetric[%d]: A(%d, %d) = %lf\n", rank, i, A.indice[j-1], A.value[j-1]);
        }
    }
#endif
    server.solve(solver, A, b, x, res);

    return cc;
}

void test(VESolver &server) {
    int cc=-1;
    char* env;
    double res = 1.0e-13;
    int solver = VESOLVER_BICGSTAB2;
    int mode = VES_MODE_GATHER_ON_VH;
    Vector x;

    if ((env = getenv("VESOLVER")) != NULL) {
        if(strcmp(env, "BICGSTAB2") == 0) {
            solver = VESOLVER_BICGSTAB2;
            printf("Using BICGSTAB2 solver.\n");
        } else if(strcmp(env, "HS") == 0) {
            solver = VESOLVER_HS;
            printf("Using HeteroSolver.\n");
        } else {
            solver = VESOLVER_DUMMY;
            printf("Using DUMMY solver.\n");
        }
    }

    if ((env = getenv("VES_MODE")) != NULL) {
        if(strcmp(env, "SYMMETRIC") == 0) {
            mode = VES_MODE_SYMMETRIC;
            printf("Run as testing mode(VES_MODE_SYMMETRIC)...\n");
        } else if(strcmp(env, "GATHER_ON_VE") == 0) {
            mode = VES_MODE_GATHER_ON_VE;
            printf("Run as testing mode(VES_MODE_GETHER_ON_VE)...\n");
        } else {
            mode = VES_MODE_GATHER_ON_VH;
            printf("Run as testing mode(VES_MODE_GETHER_ON_VH)...\n");
        }
    }
        

    /* Call solver */
    switch(mode) {
        case VES_MODE_GATHER_ON_VH:
            cc = test_gather_on_vh(server, x, res, solver);
            break;

        case VES_MODE_GATHER_ON_VE:
            cc = test_gather_on_ve(server, x, res, solver);
            break;

        case VES_MODE_SYMMETRIC:
            cc = test_symmetric(server, x, res, solver);
            break;
    }

    /* Show results */
    if (cc == 0) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            /* Print the solution vector */
            printf("%s\n", "******** Solution ********");
            for (int i=0; i<5; i++) {
                printf("  %s%1d%s%14.12f\n", "x[", i , "] = ", x.value[i] );
            }
            printf("%s\n", "********** End ***********");
        }
    } else {
        printf("ERROR: vesolver_test failed.\n");
    }
}

/*
 * Main routine
 */
int main(int argc, char** argv) {
    VESolver server;
    MPI_Comm comm;
    int rank, err;
    VES_request_t req;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_split(MPI_COMM_WORLD, VES_COLOR_SERVER, rank, &comm);

    err = server.init(comm);
    if (err < 0) {
        test(server);
        MPI_Finalize();
        return 0;
    }
    printf("INFO: VESolver server initialized.\n");

    /*
     * Calculation loop
     */
    while (server.wait_activate() == 0) {
        while(server.wait_request(req) == 0) {
            server.process(req);
        }
        server.deactivate();
    }

    server.fini();
    MPI_Finalize();
    return 0;
}
