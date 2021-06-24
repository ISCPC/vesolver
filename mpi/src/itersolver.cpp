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
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <cblas.h>
#include <stdint.h>
#include "itersolver.hpp"
#include <omp.h>
#ifdef SXAT
#include <sblas.h>
#endif
#include "timelog.h"

#define HUTI_MAXIT          1000
#define HUTI_DBUGLVL        10
//#define HUTI_TOLERANCE      1.0e-8
#define HUTI_MAXTOLERANCE   1.0e+5

#define HUTI_CONVERGENCE       0
#define HUTI_BICGSTAB_2_RHO    -1
#define HUTI_DIVERGENCE        -2
#define HUTI_MAXITER           -3

#define DNRM2(N,X) cblas_dnrm2(N,X,1)
#define DDOT(N,X,Y) cblas_ddot(N,X,1,Y,1)

class Matrix {
public:
    Matrix() {
    }

    ~Matrix() {
        if (pointers) { free(pointers); }
        if (indice) { free(indice); }
        if (values) { free(values); }

        if (iarow) { free(iarow); }

#ifdef SXAT
        /* Destruction of the handle */
        sblas_destroy_matrix_handle(handle);
#endif /* SXAT */
    }

    /* Matrix-vector multiplication z = alpha * A * x + beta * y*/
    virtual int MatVec(double alpha, double* x, double beta, double* y, double* z) {
#ifdef SXAT
        for(int i=0; i<NROWS; i++) {
            z[i] = y[i];
        }
        int ierr = sblas_execute_mv_rd(SBLAS_NON_TRANSPOSE, handle, alpha, x, beta, z);
        if (ierr != SBLAS_OK) {
            printf("ERROR: sblas_execute_mv_rd() failed with %d\n", ierr);
            exit(1);
        }
#else /* SXAT */
        if (pointers) {
            // CSR format
            #pragma omp parallel for
            for(int i=0; i<NROWS; i++) {
                z[i] = beta * y[i];
                for(int j=pointers[i]-1; j<pointers[i+1]-1; j++) {
                    z[i] += alpha * values[j] * x[indice[j]-1];
                }
            }
        } else {
            // COO format
            #pragma omp parallel for
            for(int i=0; i<NROWS; i++) {
                z[i] = beta * y[i];
            }
            for(int j=0; j<NNZ; j++) {
                int i = iarow[j]-1;
                z[i] += alpha * values[j] * x[indice[j]-1];
            }
        }
#endif /* SXAT */

        return 0;
    }

    void setMatrixCSR(int nrows, int nnz, int *aptr, int *aind, double *aval) {
        NROWS = nrows;
        NNZ = nnz;
        pointers = (int*)calloc(sizeof(int), NROWS+1);
        indice = (int*)calloc(sizeof(int), NNZ);
        values = (double*)calloc(sizeof(double), NNZ);

        bcopy(aptr, pointers, sizeof(int)*(NROWS+1));
        bcopy(aind, indice, sizeof(int)*NNZ);
        bcopy(aval, values, sizeof(double)*NNZ);

        /* create COO format */
        iarow = (int*)calloc(sizeof(int), NNZ);

        for(int i=0; i<nrows; i++) {
            for(int j=pointers[i]-1; j<pointers[i+1]-1; j++) {
                iarow[j] = i+1;
            }
        }

#ifdef SXAT
        /* Creation of a handle from CSR format */
        int ierr = sblas_create_matrix_handle_from_csr_rd(
            NROWS, NROWS, pointers, indice, values,
            SBLAS_INDEXING_1, SBLAS_GENERAL, &handle);
        if (ierr != SBLAS_OK) {
            printf("ERROR: sblas_create_matrix_handle_from_csr_rd() failed with %d\n", ierr);
            exit(1);
        }
       
        /* Analysis of the sparse matrix A */
        ierr = sblas_analyze_mv_rd(SBLAS_NON_TRANSPOSE, handle);
        if (ierr != SBLAS_OK) {
            printf("ERROR: sblas_analyze_mv_rd() failed with %d\n", ierr);
            exit(1);
        }
#endif
    }

    void setMatrixCOO(int nrows, int nnz, int *rows, int *cols, double *aval) {
        NROWS = nrows;
        NNZ = nnz;
        pointers = NULL;
        iarow = (int*)calloc(sizeof(int), NNZ);
        indice = (int*)calloc(sizeof(int), NNZ);
        values = (double*)calloc(sizeof(double), NNZ);

        bcopy(rows, iarow, sizeof(int)*NNZ);
        bcopy(cols, indice, sizeof(int)*NNZ);
        bcopy(aval, values, sizeof(double)*NNZ);

#ifdef SXAT
        /* Creation of a handle from COO format */
        int ierr = sblas_create_matrix_handle_from_coo_rd(
            NROWS, NROWS, NNZ, iarow, indice, values,
            SBLAS_INDEXING_1, SBLAS_GENERAL, &handle);
        if (ierr != SBLAS_OK) {
            printf("ERROR: sblas_create_matrix_handle_from_coo_rd() failed with %d\n", ierr);
            exit(1);
        }
       
        /* Analysis of the sparse matrix A */
        ierr = sblas_analyze_mv_rd(SBLAS_NON_TRANSPOSE, handle);
        if (ierr != SBLAS_OK) {
            printf("ERROR: sblas_analyze_mv_rd() failed with %d\n", ierr);
            exit(1);
        }
#endif
    }

    int getNRows() {
        return NROWS;
    }

private:
    int NROWS;
    int NNZ;
    int* pointers = NULL;
    int* indice = NULL;
    double* values = NULL;
    int* iarow = NULL;

#ifdef SXAT
    sblas_handle_t handle; 
#endif
};

static int solve(Matrix& A, double *b, double *x, double tolerance) {
    // Local variables
    double rho, oldrho, alpha, beta, omega1, omega2;
    double tau, delta, myy;
    int iter_count;
    double residual, rhsnorm;
    int ndim = A.getNRows();
    int debug=1;

    // work vectors
    double* work = (double*)calloc(sizeof(double), ndim*7);
    double* R = &work[0];
    double* S = &work[ndim];
    double* T = &work[ndim*2];
    double* U = &work[ndim*3];
    double* V = &work[ndim*4];
    double* W = &work[ndim*5];
    double* RTLD = &work[ndim*6];

    // Norms of right-hand side vector are used in convergence tests
    rhsnorm = DNRM2( ndim, b );

    // Generate vector X if needed
    for(int i=0; i<ndim; i++) {
        //x[i] = rand();
        x[i] = 1.0e-8;
        //x[i] = 1.0;
    }

    A.MatVec(-1.0, x, 1.0, b, R); //R = b - A * x;
    bcopy(R, RTLD, sizeof(double)*ndim);
    for(int i=0; i<ndim; i++) { U[i] = 0; }
    oldrho = 1; omega2 = 1; alpha = 0;

    /*
     * This is where the loop starts (that is we continue from here after
     * the first iteration)
     */
    for(iter_count=1; iter_count<HUTI_MAXIT; iter_count++) {
        oldrho = -omega2 * oldrho;

        /*
         * This is the even BiCG step
         */
        rho = DDOT( ndim, RTLD, R );
        if ( rho == 0 ) {
            free(work);
            return HUTI_BICGSTAB_2_RHO;
        }

        beta = ( rho * alpha ) / oldrho;
        oldrho = rho;
        //U = R - beta * U;
        for(int i=0; i<ndim; i++) {
            U[i] = R[i] - beta * U[i];
        }

        A.MatVec(1.0, U, 0.0, V, V); // V = A * U;

        alpha = oldrho / DDOT( ndim, RTLD, V );
        //R = R - alpha * V;
        for(int i=0; i<ndim; i++) {
            R[i] -= alpha * V[i];
        }

        A.MatVec(1.0, R, 0.0, S, S); // S = A * R
        //x = x + alpha * U
        for(int i=0; i<ndim; i++) {
            x[i] += alpha * U[i];
        }

        /*
         * This is the odd BiCG step
         */
        rho = DDOT( ndim, RTLD, S );
        if ( rho == 0 ) {
            free(work);
            return HUTI_BICGSTAB_2_RHO;
        }

        beta = ( rho * alpha ) / oldrho;
        oldrho = rho;
        //V = S - beta * V;
        for(int i=0; i<ndim; i++) {
            V[i] = S[i] - beta * V[i];
        }

        A.MatVec(1.0, V, 0.0, W, W); // W = A * V;

        alpha = oldrho / DDOT( ndim, RTLD, W );
        //U = R - beta * U;
        //R = R - alpha * V;
        //S = S - alpha * W;
        for(int i=0; i<ndim; i++) {
            U[i] = R[i] - beta * U[i];
            R[i] -= alpha * V[i];
            S[i] -= alpha * W[i];
        }

        A.MatVec(1.0, S, 0.0, T, T); // T = A * S;

        /*
         * This is the GCR(2) part
         */
        omega1 = DDOT( ndim, R, S );
        myy = DDOT( ndim, S, S );
        delta = DDOT( ndim, S, T );
        tau = DDOT( ndim, T, T );
        omega2 = DDOT( ndim, R, T );

        tau = tau - ( delta * delta ) / myy;
        omega2 = ( omega2 - ( delta * omega1 ) / myy ) / tau;
        omega1 = ( omega1 - delta * omega2 ) / myy;

        //x = x + omega1 * R + omega2 * S + alpha * U;
        //R = R - omega1 * S - omega2 * T;
        for(int i=0; i<ndim; i++) {
            x[i] += omega1 * R[i] + omega2 * S[i] + alpha * U[i];
            R[i] -= omega1 * S[i] + omega2 * T[i];
        }

        /*
         * Check the convergence against selected stopping criterion
         */
        A.MatVec(1.0, x, -1.0, b, S); // S = A * x - b;
        residual = DNRM2( ndim, S ) / rhsnorm;

        /*
         * Print debugging output if desired
         */
        if (debug) {
            if ( iter_count % HUTI_DBUGLVL == 0 ) {
                printf("%d %le\n", iter_count, residual);
            }
        }

        if ( residual < tolerance ) {
            free(work);
            return HUTI_CONVERGENCE;
        }

        if ( isnan(residual) || (residual > HUTI_MAXTOLERANCE )) {
            free(work);
            return HUTI_DIVERGENCE;
        }
    
        //U = U - omega1 * V - omega2 * W;
        for(int i=0; i<ndim; i++) {
            U[i] -= omega1 * V[i] + omega2 * W[i];
            // printf("U[%d] = %lf\n", i, U[i]);
        }
    }

    free(work);
    return HUTI_MAXITER;
}

int itersolver(INT_T neq, INT_T nnz, INT_T *pointers, INT_T *indice, double *value,
    double* b, double* x, double res)
{
    Matrix A;
    TIMELOG(tl);
    
    A.setMatrixCSR(neq, nnz, pointers, indice, value);
    TIMELOG_START(tl);
    int cc = solve(A, b, x, res);
    TIMELOG_END(tl, "solve");
    return cc;
}

int itersolver_coo(INT_T neq, INT_T nnz, INT_T *rows, INT_T *cols, double *value,
    double* b, double* x, double res)
{
    Matrix A;
    TIMELOG(tl);
    
    A.setMatrixCOO(neq, nnz, rows, cols, value);
    TIMELOG_START(tl);
    int cc = solve(A, b, x, res);
    TIMELOG_END(tl, "solve");
    return cc;
}
