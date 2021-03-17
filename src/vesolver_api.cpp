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
#include <mpi.h>
#include <cstdio>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include "SpMatrix.hpp"
#include "veserver.hpp"
#include "vesolver.hpp"
#include "vesolver_api.h"
#include "timelog.h"

#include <sys/mman.h>
#include <errno.h>

//#define SANITY_CHECK 1

/*
 * VEsolver client methods
 */
/*
 * Send Matrix Info
 */
int VESolverAPI::send_matrix_info(SpMatrix& A) {
    vesolver_matrix_info_t info;

    info.type = A.type;
    info.nrow = A.nrow;
    info.ncol = A.ncol;
    info.nval = A.ndim;
    info.neq = A.neq;
    printf("INFO:VESolverAPI::send_matrix_info: Sending matrix info to Rank#%d (0x%04x %d %d)\n",
        0, info.type, info.nrow, info.nval);

    return MPI_Send(&info, sizeof(vesolver_matrix_info_t), MPI_BYTE, 0, VES_MATINFO_TAG, ves_comm);
}

int VESolverAPI::send_matrix_info(int rank, SpDistMatrix& A) {
    vesolver_matrix_info_t info;

    info.type = A.type;
    info.nrow = A.nrow;
    info.ncol = A.ncol;
    info.nval = A.ndim;
    info.neq = A.neq;
    info.nl = A.mingind;
    info.nt = A.maxgind;
    printf("INFO:VESolverAPI::send_matrix_info: Sending dist matrix info to Rank#%d (0x%04x %d %d %d %d %d)\n",
        rank, info.type, info.nrow, info.nval, info.neq, info.nl, info.nt);

    return MPI_Send(&info, sizeof(vesolver_matrix_info_t), MPI_BYTE, rank, VES_MATINFO_TAG, ves_comm);
}

/*
 * Send Matrix Data
 */
int VESolverAPI::send_matrix_data(SpMatrix& A, Vector& b) {
    int cc;

#if 0
    for(int i=0; i<5; i++) {
        printf("INFO:Send:Ai[%d] = %d %d %lf %lf\n", i, A.pointers[i],
            A.indice[i], A.value[i], b.value[i]);
    }
    for(int i=A.nrow-5; i<A.nrow; i++) {
        printf("INFO:Send:Ai[%d] = %d %d %lf %lf\n", i, A.pointers[i],
            A.indice[i], A.value[i], b.value[i]);
    }
    for(int i=A.ndim-5; i<A.ndim; i++) {
        printf("INFO:Send:Aj[%d] = (%d, %lf)\n", i, A.indice[i], A.value[i]);
    }
#endif

#if 0
    // Send Matrix A
    cc = MPI_Send(A.value, A.ndim, MPI_DOUBLE, 0, VES_MATDATA_TAG, ves_comm);
    cc = MPI_Send(A.indice, A.ndim, MPI_INTEGER, 0, VES_MATDATA_TAG, ves_comm);
    cc = MPI_Send(A.pointers, (A.nrow)+1, MPI_INTEGER, 0, VES_MATDATA_TAG, ves_comm);

    // Send Matrix b
    cc = MPI_Send(b.value, A.nrow, MPI_DOUBLE, 0, VES_MATDATA_TAG, ves_comm);
#else
    //if (mlock(pointers, datasz)) {
    //    perror("WARNING: mlock failed.");
    //}

    // Send value of Matrix A
    bcopy(A.value, buffer, sizeof(double)*(A.ndim));
    cc = MPI_Send(buffer, A.ndim, MPI_DOUBLE, 0, VES_MATDATA_TAG, ves_comm);

    // Send indice of Matrix A
    bcopy(A.indice, buffer, sizeof(INT_T)*(A.ndim));
    cc = MPI_Send(buffer, A.ndim, MPI_INTEGER, 0, VES_MATDATA_TAG, ves_comm);
    
    // Send pointers of Matrix A
    printf("INFO:send_matrix_data: data location of A.pointers=0x%08lx\n", A.pointers);
    bcopy(A.pointers, buffer, sizeof(INT_T)*(A.nrow+1));
    cc = MPI_Send(buffer, (A.nrow)+1, MPI_INTEGER, 0, VES_MATDATA_TAG, ves_comm);

    // Send Vector b
    bcopy(b.value, buffer, sizeof(double)*(A.nrow));
    cc = MPI_Send(buffer, A.nrow, MPI_DOUBLE, 0, VES_MATDATA_TAG, ves_comm);

    //munlock(pointers, datasz);
#endif

    return cc;
}

int VESolverAPI::send_matrix_data(int rank, SpDistMatrix& A, DistVector& b) {
    int cc;

#if 0
    // Send Matrix A
    cc = MPI_Send(A.value, A.ndim, MPI_DOUBLE, rank, VES_MATDATA_TAG, ves_comm);
    cc = MPI_Send(A.indice, A.ndim, MPI_INTEGER, rank, VES_MATDATA_TAG, ves_comm);
    cc = MPI_Send(A.pointers, A.nrow+1, MPI_INTEGER, rank, VES_MATDATA_TAG, ves_comm);
    cc = MPI_Send(A.rorder, A.neq, MPI_INTEGER, rank, VES_MATDATA_TAG, ves_comm);

    // Send Matrix b
    cc = MPI_Send(b.value, A.nrow, MPI_DOUBLE, rank, VES_MATDATA_TAG, ves_comm);
#else
    // Send value of Matrix A
    bcopy(A.value, buffer, sizeof(double)*(A.ndim));
    cc = MPI_Send(buffer, A.ndim, MPI_DOUBLE, rank, VES_MATDATA_TAG, ves_comm);

    // Send indice of Matrix A
    bcopy(A.indice, buffer, sizeof(INT_T)*(A.ndim));
    cc = MPI_Send(buffer, A.ndim, MPI_INTEGER, rank, VES_MATDATA_TAG, ves_comm);
    
    // Send pointers of Matrix A
    bcopy(A.pointers, buffer, sizeof(INT_T)*(A.nrow+1));
    cc = MPI_Send(buffer, (A.nrow)+1, MPI_INTEGER, rank, VES_MATDATA_TAG, ves_comm);

    // Send rorder of Matrix A
    //bcopy(A.rorder, buffer, sizeof(INT_T)*(A.neq));
    //cc = MPI_Send(buffer, A.neq, MPI_INTEGER, rank, VES_MATDATA_TAG, ves_comm);
    bcopy(A.order, buffer, sizeof(INT_T)*(A.nrow));
    cc = MPI_Send(buffer, A.nrow, MPI_INTEGER, rank, VES_MATDATA_TAG, ves_comm);

    // Send Vector b
    bcopy(b.value, buffer, sizeof(double)*(A.nrow));
    cc = MPI_Send(buffer, A.nrow, MPI_DOUBLE, rank, VES_MATDATA_TAG, ves_comm);
#endif

    return cc;
}

/*
 * Send Matrix Data
 */
inline int VESolverAPI::receive_result(Vector& x) {
    MPI_Status status;

    return MPI_Recv(x.value, x.size, MPI_DOUBLE, 0, VES_MATDATA_TAG, ves_comm, &status);
}

inline int VESolverAPI::receive_result(int src, Vector& x) {
    MPI_Status status;

    return MPI_Recv(x.value, x.size, MPI_DOUBLE, src, VES_MATDATA_TAG, ves_comm, &status);
}

inline void VESolverAPI::send_solve_request(int mode, int solver, int type, int nprocs, double res) {
    VES_request_t req;

    req.req_type = VES_REQ_LINEARSOLVER;

    vesolver_request_t *payload = (vesolver_request_t*)VES_REQ_GETPAYLOAD(req);
    payload->mode = mode;
    payload->solver = solver;
    payload->type = type;
    payload->nprocs = nprocs;
    payload->res = res;

    send_request(req);
}

/*
 * Solver Interface for non-MPI client
 */
int VESolverAPI::solve(int solver, SpMatrix& A, Vector& b, Vector& x, double res) {
    int cc;

    if (cl_comm == MPI_COMM_NULL) {
        return -1;
    }

#ifdef SANITY_CHECK
    int nerr=0;
    printf("INFO:VESolverAPI:solve: A: nrow=%d, ncol=%d, neq=%d, ndim=%d, %d, %d\n", 
        A.nrow, A.ncol, A.neq, A.ndim, A.pointers[A.neq-1], A.pointers[A.neq]);
    for (int i=0; i<A.neq; i++) {
        if ((A.pointers[i] <= 0)||(A.pointers[i] > A.ndim+1)) {
            printf("INFO:VESolverAPI:solve: Invalid pointer A[%d]=%d\n", i, A.pointers[i]);
            nerr++;
        }

        int k = A.pointers[i]-1;
        int prev = A.indice[k];
        for (int j=A.pointers[i]; j<A.pointers[i+1]-1; j++) {
            if (A.indice[j] <= prev) {
                printf("INFO:VESolverAPI:solve: Invalid indice A[%d, %d]=%d < %d\n", i, j, A.indice[j], prev);
                nerr++;
            }
            prev = A.indice[j];
        }
        if (nerr > 100) {
            return -1;
        }
    }
    for (int i=0; i<A.ndim; i++) {
        if ((A.indice[i] <= 0)||(A.indice[i] > A.ncol)) {
            printf("INFO:VESolverAPI:solve: Invalid indice A[%d]=%d\n", i, A.indice[i]);
            nerr++;
        }
    }
#endif

    /* send request */
    send_solve_request(VES_MODE_GATHER_ON_VH, solver, A.type, 1, res);
    send_matrix_info(A);
    send_matrix_data(A, b);
    cc = receive_result(x);

    return cc;
}

/*
 * Solver Interface for GATHER_ON_VH
 */
typedef struct SpMatrix_info {
    INT_T nrow;
    INT_T ndim;
    INT_T neq;
    INT_T pad;
} SpMatrix_info_t;

int VESolverAPI::solve(int solver, SpDistMatrix **An, DistVector **bn, int n, Vector& x, double res) {
    int cc=-1;
    TIMELOG(tl);

    /* gathter */
    TIMELOG_START(tl);
    SpMatrix *A = SpDistMatrix::gather(An, n);
    TIMELOG_END(tl, "gather_A");

    TIMELOG_START(tl);
    Vector *b = DistVector::gather(bn, n, An);
    TIMELOG_END(tl, "gather_B");

    /* solve on vesolver */
    solve(solver, *A, *b, x, res);

    /*
     * free Matrix/Vector
     */
    delete A;
    delete b;

    return cc;
}

int VESolverAPI::solve_gather_on_vh(int solver, SpDistMatrix& A, DistVector& b, Vector& x, double res) {
    int rank, nprocs;
    int cc=0;
    SpDistMatrix **An;
    DistVector **bn;
    MPI_Status status;

    if (cl_comm == MPI_COMM_NULL) {
        return -1;
    }

    MPI_Comm_rank(cl_comm, &rank);
    MPI_Comm_size(cl_comm, &nprocs);
    //printf("INFO:VESolverAPI[%d]::solve_gather_on_vh is called.\n", rank);

    /* Get matrix size */
    SpMatrix_info_t* info = (SpMatrix_info_t*)calloc(nprocs, sizeof(SpMatrix_info_t));
    info[rank].nrow = A.nrow;
    info[rank].ndim = A.ndim;
    info[rank].neq = A.neq;
    MPI_Allgather(&info[rank], sizeof(SpMatrix_info_t), MPI_BYTE, info, sizeof(SpMatrix_info_t), MPI_BYTE, cl_comm);
    
    /* Gathering data on master process */
    if (rank == 0) {
        for(int i=0; i<nprocs; i++) {
            printf("INFO:VESolverAPI::A[%d] nrow=%d, ndim=%d, neq=%d.\n",
                i, info[i].nrow, info[i].ndim, info[i].neq);
        }

        An = (SpDistMatrix**)calloc(sizeof(SpMatrix*), nprocs);
        bn = (DistVector**)calloc(sizeof(Vector*), nprocs);

        An[0] = &A;
        bn[0] = &b;
        // Receive Matrix/Vector
        for(int i=1; i<nprocs; i++) {
            // Memory Allocation
            An[i] = new SpDistMatrix();
            An[i]->alloc(info[i].nrow, info[i].ndim);
            An[i]->rorder = (INT_T*)calloc(info[i].neq, sizeof(INT_T));
            An[i]->nrow = An[i]->ncol = info[i].nrow;
            An[i]->ndim = info[i].ndim;
            An[i]->neq = info[i].neq;

            bn[i] = new DistVector();
            bn[i]->alloc(info[i].nrow);

            // Receive Matrix A
            MPI_Recv(An[i]->value, info[i].ndim, MPI_DOUBLE, i, VES_MATGATHER_TAG, cl_comm, &status);
            MPI_Recv(An[i]->indice, info[i].ndim, MPI_INTEGER, i, VES_MATGATHER_TAG, cl_comm, &status);
            MPI_Recv(An[i]->pointers, info[i].nrow+1, MPI_INTEGER, i, VES_MATGATHER_TAG, cl_comm, &status);
            MPI_Recv(An[i]->rorder, info[i].neq, MPI_INTEGER, i, VES_MATGATHER_TAG, cl_comm, &status);

            // Receive Vector b
            MPI_Recv(bn[i]->value, info[i].nrow, MPI_DOUBLE, i, VES_MATGATHER_TAG, cl_comm, &status);
        }

        /* Allocate Vector x */
        solve(solver, An, bn, nprocs, x, res);

        /*
         * free Matrix/Vector
         */
        for(int i=1; i<nprocs; i++) {
            delete An[i];
            delete bn[i];
        }
        free(An);
        free(bn);
    } else {
        // Send Matrix A
        MPI_Send(A.value, A.ndim, MPI_DOUBLE, 0, VES_MATGATHER_TAG, cl_comm);
        MPI_Send(A.indice, A.ndim, MPI_INTEGER, 0, VES_MATGATHER_TAG, cl_comm);
        MPI_Send(A.pointers, A.nrow+1, MPI_INTEGER, 0, VES_MATGATHER_TAG, cl_comm);
        MPI_Send(A.rorder, A.neq, MPI_INTEGER, 0, VES_MATGATHER_TAG, cl_comm);

        // Send Vector b
        MPI_Send(b.value, A.nrow, MPI_DOUBLE, 0, VES_MATGATHER_TAG, cl_comm);
    }
    free(info);

    /* deliver results */
    MPI_Bcast(x.value, A.neq, MPI_DOUBLE_PRECISION, 0, cl_comm);

    return cc;
}

/*
 * Solver Interface for GATHER_ON_VE
 */
int VESolverAPI::solve_gather_on_ve(int solver, SpDistMatrix& A, DistVector& b, Vector& x, double res) {
    int rank, nprocs;
    int cc = 0;

    if (cl_comm == MPI_COMM_NULL) {
        return -1;
    }

    MPI_Comm_rank(cl_comm, &rank);
    MPI_Comm_size(cl_comm, &nprocs);
    printf("INFO:VESolverAPI[%d]::solve_gather_on_ve is called.\n", rank);

    send_solve_request(VES_MODE_GATHER_ON_VE, solver, A.type, nprocs, res);
    send_matrix_info(0, A);
    send_matrix_data(0, A, b);

    if (rank == 0) {
        cc = receive_result(0, x);
    }

    /* deliver results */
    MPI_Bcast(x.value, A.neq, MPI_DOUBLE_PRECISION, 0, cl_comm);

    return cc;
}

/*
 * Solver Interface for SYMMETRIC
 */
int VESolverAPI::solve_symmetric(int solver, SpDistMatrix& A, DistVector& b, Vector& x, double res) {
    int rank, nprocs;
    int cc = 0;

    if (cl_comm == MPI_COMM_NULL) {
        return -1;
    }

    MPI_Comm_rank(cl_comm, &rank);
    MPI_Comm_size(cl_comm, &nprocs);
    printf("INFO:VESolverAPI[%d]::solve_symmetric is called.\n", rank);

    send_solve_request(VES_MODE_SYMMETRIC, solver, A.type, nprocs, res);
    send_matrix_info(rank, A);
    send_matrix_data(rank, A, b);
    cc = receive_result(rank, x);

    return cc;
}

/*
 * Solver Interface for MPI-Parallelized client
 */
int VESolverAPI::solve(int solver, SpDistMatrix& A, DistVector& b, Vector& x, double res, int mode) {
    int cc = -1; 

    if (cl_comm == MPI_COMM_NULL) {
        return -1;
    }

    A.reordering();

    /* Call solver */
    switch(mode) {
        case VES_MODE_GATHER_ON_VH:
            cc = solve_gather_on_vh(solver, A, b, x, res);
            break;

        case VES_MODE_GATHER_ON_VE:
            cc = solve_gather_on_ve(solver, A, b, x, res);
            break;

        case VES_MODE_SYMMETRIC:
            cc = solve_symmetric(solver, A, b, x, res);
            break;
    }

    return cc;
}

/*
 * C-API
 */
static VESolverAPI vesolver;

int vesolver_init(MPI_Comm comm) {
    return vesolver.init(comm);
}

int vesolver_init_f(int comm) {
    return vesolver.init(MPI_Comm_f2c(comm));
}

int vesolver_fini() {
    vesolver.fini();
    return 0;
}

int vesolver_activate(MPI_Comm comm, int nprocs) {
    return vesolver.activate(comm, nprocs);
}

int vesolver_activate_f(int comm, int nprocs) {
    return vesolver.activate(MPI_Comm_f2c(comm), nprocs);
}

int vesolver_deactivate() {
    vesolver.deactivate();
    return 0;
}

int vesolver_solve(int solver, int32_t mtype, int32_t neq, int32_t *pointers, int32_t *indice, double *value, double *b, double *x, double res) {
    SpMatrix A;
    Vector b0, x0;

    A.pointers = pointers;
    A.indice = indice;
    A.value = value;
    A.ndim = pointers[neq]-1;
    A.nrow = A.ncol = A.neq = neq;
    A.type = SPMATRIX_TYPE_CSR | SPMATRIX_TYPE_INDEX1 | SPMATRIX_TYPE_ASYMMETRIC;

    b0.size = neq;
    b0.value = b;
     
    x0.size = neq;
    x0.value = x;
     
    return vesolver.solve(solver, A, b0, x0, res);
}

int vesolver_psolve(int32_t mode, int32_t solver, int32_t neq, int32_t nrows, int32_t *pointers, int32_t *indice, double *value, int32_t *order, double *b, double *x, double res) {
    SpDistMatrix A;
    DistVector b0;
    Vector x0;

    A.pointers = pointers;
    A.indice = indice;
    A.value = value;
    A.order = order;
    A.ndim = pointers[nrows]-1;
    A.nrow = A.ncol = nrows;
    A.neq = neq;
    A.type = SPMATRIX_TYPE_CSR | SPMATRIX_TYPE_INDEX1 | SPMATRIX_TYPE_ASYMMETRIC | SPMATRIX_TYPE_DISTRIBUTE;

    b0.size = nrows;
    b0.value = b;
     
    x0.size = neq;
    x0.value = x;
     
    return vesolver.solve(solver, A, b0, x0, res, mode);
}
#if 0
int vesolver_psolve2(int32_t mode, int32_t solver, int32_t neq, int32_t nrows, int32_t nl, int32_t nt, int32_t *pointers, int32_t *indice, double *value, int32_t *rorder, double *b, double *x, double res) {
    SpDistMatrix A;
    DistVector b0;
    Vector x0;

    A.pointers = pointers;
    A.indice = indice;
    A.value = value;
    A.rorder = rorder;
    A.ndim = pointers[nrows]-1;
    A.nrow = A.ncol = nrows;
    A.neq = neq;
    A.mingind = nl;
    A.maxgind = nt;
    A.type = SPMATRIX_TYPE_CSR | SPMATRIX_TYPE_INDEX1 | SPMATRIX_TYPE_ASYMMETRIC | SPMATRIX_TYPE_DISTRIBUTE;

    b0.size = nrows;
    b0.value = b;
     
    x0.size = neq;
    x0.value = x;
     
    return vesolver.solve(solver, A, b0, x0, res, mode);
}

int vesolver_psolve_dcsr(int32_t mode, int32_t solver, int32_t neq, int32_t *pointers, int32_t *indice, double *value, int32_t nl, int32_t nt, double *b, double *x, double res) {
    SpDistMatrix A;
    DistVector b0;
    Vector x0;

    int nrows=nt-nl+1;

    A.pointers = pointers;
    A.indice = indice;
    A.value = value;
    A.ndim = pointers[nrows]-1;
    A.nrow = nrows;
    A.ncol = A.neq = neq;
    A.mingind = nl;
    A.maxgind = nt;
    A.type = SPMATRIX_TYPE_DCSR | SPMATRIX_TYPE_INDEX1 | SPMATRIX_TYPE_ASYMMETRIC | SPMATRIX_TYPE_DISTRIBUTE;

    b0.size = nrows;
    b0.value = b;
     
    x0.size = neq;
    x0.value = x;
     
    return vesolver.solve(solver, A, b0, x0, res, mode);
}
#endif
