#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"
#include "PluginAPI.h"

#define HUTI_MAXIT          50000
#define HUTI_DBUGLVL        100
//#define HUTI_TOLERANCE      1.0e-8
#define HUTI_MAXTOLERANCE   1.0e+5

#define HUTI_CONVERGENCE       0
#define HUTI_BICGSTAB_2_RHO    -1
#define HUTI_DIVERGENCE        -2
#define HUTI_MAXITER           -3

//#define SCALING 1
//#define DIAGONAL 1
//#define CALCULIX_NORM 1

typedef struct elmer_info {
    double* factor;
    double* diag;
} elmer_info_t;

static Matrix_t* solve_pre(const Matrix_t* A0) {
    Matrix_t* A = Matrix_duplicate(A0);

    elmer_info_t* info = (elmer_info_t*)malloc(sizeof(elmer_info_t));
#ifdef SCALING
    double* factor = (double*)calloc(sizeof(double), A->NROWS);
    if (factor == NULL) {
        return NULL;
    }

    Matrix_convert_index(A, 0);
    /*  extract diagonal vector from matrix A  */
    for (int i=0; i<A->NROWS; i++) {
        factor[i] = 1.0/sqrt(A->values[A->pointers[i]]);
    }

    /*  scale matrix A  */
    for (int i=0; i<A->NROWS; i++) {
        for (int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
            A->values[j] *= (factor[i]) * (factor[A->indice[j]]);
        }
    }

    info->factor = factor;
#endif
#ifdef DIAGONAL
    double* diag = (double*)calloc(sizeof(double), A->NROWS);
    if (diag == NULL) {
        return NULL;
    }

    Matrix_convert_index(A, 0);
    /*  extract diagonal vector from matrix A  */
    if (MATRIX_IS_SYMMETRIC(A)) {
        if(MATRIX_IS_LOWER(A)) {
            for (int i=0; i<A->NROWS; i++) {
                diag[i] = A->values[A->pointers[i+1]-1];
            }
        } else {
            printf("get diagonal elements of upper symmetric matrix.\n");
            for (int i=0; i<A->NROWS; i++) {
                diag[i] = A->values[A->pointers[i]];
            }
        }
    } else {
        printf("get diagonal elements of asymmetric matrix.\n");
        for (int i=0; i<A->NROWS; i++) {
            for (int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
                if (A->indice[j] == i) {
                    diag[i] = A->values[j];
                    break;
                }
            }
        }
    }

    info->diag = diag;
#endif /* DIAGONAL */
    A->info = (void*)info;
    return A;
}

static int solve_post(Matrix_t* A) {
    elmer_info_t* info = (elmer_info_t*)(A->info);
#ifdef SCALING
    if (info->factor) { free(info->factor); }
#endif
#ifdef DIAGONAL
    if (info->diag) { free(info->diag); }
#endif
    Matrix_free(A);
    return 0;
}

static inline void pcondl(const Matrix_t* A, double* u, double* v) {
    int ndim = A->NROWS;
#ifdef DIAGONAL
    elmer_info_t* info = (elmer_info_t*)(A->info);
    double *diag = info->diag;

    for (int i=0; i<ndim; i++) {
        if  ( fabs(diag[i]) > DBL_EPSILON ) {
            u[i] = v[i] / diag[i];
        } else {
            u[i] = v[i];
        } 
    }
#else
    if (u != v) {
        for (int i=0; i<ndim; i++) {
            u[i] = v[i];
        }
    }
#endif
}

static int solve(const Matrix_t *A, const double *b, double *x, const double tolerance) {
    // Local variables
    double rho, oldrho, alpha, beta, omega1, omega2;
    double tau, delta, myy;
    int iter_count;
    double residual, rhsnorm;
    int ndim = A->NROWS;
    int debug=1;

    // work vectors
    double* work = (double*)calloc(sizeof(double), ndim*8);
    double* R = &work[0];
    double* S = &work[ndim];
    double* T = &work[ndim*2];
    double* U = &work[ndim*3];
    double* V = &work[ndim*4];
    double* W = &work[ndim*5];
    double* RTLD = &work[ndim*6];
    double* T1 = &work[ndim*7];

    elmer_info_t* info = (elmer_info_t*)(A->info);
#ifdef SCALING
    double *factor = info->factor;
    double *b1 = (double*)calloc(sizeof(double), ndim);
    /*  scale right hand side (Ax=b -> Ax+b=0: negative sign)  */
    for (int i=0; i<ndim; i++) {
        //b[i] *= -factor[i];
        b1[i] = b[i] * factor[i];
    }
    b = (const double*)b1;
#endif

#ifdef CALCULIX_NORM
    tolerance = 0.0;
    for (int i=0; i<ndim; i++) {
        double err;
        err = fabs(b[i]);
        if (err > 1.0e-20) {
            tolerance += err;
        }
    }
    tolerance = (tolerance / ndim) * 0.005;
    printf("tolerance = %e\n", tolerance);
#endif

    // Norms of right-hand side vector are used in convergence tests
    rhsnorm = DNRM2( ndim, b );

    // Generate vector X if needed
    for(int i=0; i<ndim; i++) {
        //x[i] = rand();
        x[i] = 1.0e-8;
        //x[i] = 1.0;
    }

    Matrix_MV(A, -1.0, x, 1.0, b, U); // U = b - A * x;
    pcondl(A, R, U);
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

        Matrix_MV(A, 1.0, U, 0.0, T1, T1); // V = A * U;
        pcondl(A, V, T1);

        alpha = oldrho / DDOT( ndim, RTLD, V );
        //R = R - alpha * V;
        for(int i=0; i<ndim; i++) {
            R[i] -= alpha * V[i];
        }

        Matrix_MV(A, 1.0, R, 0.0, T1, T1); // S = A * R
        pcondl(A, S, T1);
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

        Matrix_MV(A, 1.0, V, 0.0, T1, T1); // W = A * V;
        pcondl(A, W, T1);

        alpha = oldrho / DDOT( ndim, RTLD, W );
        //U = R - beta * U;
        //R = R - alpha * V;
        //S = S - alpha * W;
        for(int i=0; i<ndim; i++) {
            U[i] = R[i] - beta * U[i];
            R[i] -= alpha * V[i];
            S[i] -= alpha * W[i];
        }

        Matrix_MV(A, 1.0, S, 0.0, T1, T1); // T = A * S;
        pcondl(A, T, T1);

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
#ifdef CALCULIX_NORM
        residual = 0.0;
        for (int i=0; i<ndim; i++) {
            double err = fabs(R[i]);
            if(err > residual) residual = err;
        }
#else
        //Matrix_MV(A, 1.0, x, -1.0, b, S); // S = A * x - b;
        //pcondl(A, S, S);
        //residual = DNRM2( ndim, S ) / rhsnorm;
        //residual = DNRM2( ndim, S );
        //residual = DNRM2( ndim, R );
        residual = DNRM2( ndim, R ) / rhsnorm;
#endif

        /*
         * Print debugging output if desired
         */
        if (debug) {
            if ( iter_count % HUTI_DBUGLVL == 0 ) {
                printf("%d %le\n", iter_count, residual);
                fflush(stdout);
            }
        }

        if ( residual < tolerance ) {
            printf("# of iterations = %d\n",iter_count);
            free(work);
#ifdef SCALING
            /*  Backscaling of the solution vector  */
            for (int i=0; i<ndim; i++) {
                x[i] *= factor[i];
            }
            free(b1);
#endif
            return HUTI_CONVERGENCE;
        }

        if ( isnan(residual) || (residual > HUTI_MAXTOLERANCE )) {
            free(work);
#ifdef SCALING
            free(b1);
#endif
            return HUTI_DIVERGENCE;
        }
    
        //U = U - omega1 * V - omega2 * W;
        for(int i=0; i<ndim; i++) {
            U[i] -= omega1 * V[i] + omega2 * W[i];
            // printf("U[%d] = %lf\n", i, U[i]);
        }
    }

    free(work);
#ifdef SCALING
    free(b1);
#endif
    return HUTI_MAXITER;
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
SolverPlugin_t* bicgstab2_init() {
#endif
	SolverPlugin_t* solver = (SolverPlugin_t*)malloc(sizeof(SolverPlugin_t));

	solver->set_option = NULL;
	solver->solve_pre = solve_pre;
	solver->solve = solve;
	solver->solve_post = solve_post;
	solver->free = solve_free;

	return solver;
}
