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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"
#include "PluginAPI.h"

#define HUTI_MAXIT          50000
#define HUTI_DBUGLVL        1000
//#define HUTI_TOLERANCE      1.0e-8
#define HUTI_MAXTOLERANCE   1.0e+5

#define HUTI_CONVERGENCE       0
#define HUTI_CG_RHO            -1
#define HUTI_DIVERGENCE        -2
#define HUTI_MAXITER           -3

#define SCALING 1
//#define DIAGONAL 1
//#define CALCULIX_NORM 1

typedef struct elmer_info {
    double* factor;
    double* diag;
} elmer_info_t;

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

/*
 * API
 */
static int solve_pre(Matrix_t* A) {
    elmer_info_t* info = (elmer_info_t*)malloc(sizeof(elmer_info_t));
#ifdef SCALING
    double* factor = (double*)calloc(sizeof(double), A->NROWS);
    if (factor == NULL) {
        return -1;
    }

    Matrix_convert_index(A, 0);
    if (MATRIX_IS_SYMMETRIC(A)) {
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
    } else {
        /*  extract diagonal vector from matrix A  */
        for (int i=0; i<A->NROWS; i++) {
            for (int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
                if (A->indice[j] == i) {
                    factor[i] = 1.0/sqrt(A->values[j]);
                    break;
                }
            }
        }

        /*  scale matrix A  */
        for (int i=0; i<A->NROWS; i++) {
            for (int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
                A->values[j] *= (factor[i]) * (factor[A->indice[j]]);
            }
        }
    }

    info->factor = factor;
#endif
#ifdef DIAGONAL
    double* diag = (double*)calloc(sizeof(double), A->NROWS);
    if (diag == NULL) {
        return -1;
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

    if (Matrix_optimize(A) < 0) {
        free(A);
        return -1;
    }

    A->info = (void*)info;
    return 0;
}

static int solve(Matrix_t *A, const double* b, double* x, const double tolerance) {
    // Local variables
    double rho, oldrho, alpha, beta;
    int iter_count;
    double residual, rhsnorm;
    int ndim = A->NROWS;
    int debug=0;

    // work vectors
    double* work = (double*)calloc(sizeof(double), ndim*4);
    double* P = &work[0];
    double* Q = &work[ndim];
    double* R = &work[ndim*2];
    double* Z = &work[ndim*3];

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

    Matrix_MV(A, -1.0, x, 1.0, b, R); //R = b - A * x;
    oldrho = 1; alpha = 0;

    /*
     * This is where the loop starts (that is we continue from here after
     * the first iteration)
     */
    for(iter_count=1; iter_count<HUTI_MAXIT; iter_count++) {
        pcondl(A, Z, R);
        rho = DDOT( ndim, R, Z );
        if ( rho == 0 ) {
            free(work);
            return HUTI_CG_RHO;
        }

        if (iter_count == 1) {
            // P = Z;
            for(int i=0; i<ndim; i++) {
                P[i] = Z[i];
            }
        } else {
            beta = rho / oldrho;
            // P = Z + beta * P
            for(int i=0; i<ndim; i++) {
                P[i] = Z[i] + beta * P[i];
            }
        }

        Matrix_MV(A, 1.0, P, 0.0, Q, Q); // Q = A * P;

        alpha = rho / DDOT( ndim, P, Q );

        for(int i=0; i<ndim; i++) {
            x[i] = x[i] + alpha * P[i];
            R[i] = R[i] - alpha * Q[i];
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
        //Matrix_MV(A, 1.0, x, -1.0, b, Q); // S = A * x - b;
        //residual = DNRM2( ndim, Q ) / rhsnorm;
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
    
        oldrho = rho;
    }

    free(work);
#ifdef SCALING
    free(b1);
#endif
    return HUTI_MAXITER;
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

static void solve_free(SolverPlugin_t* solver) {
	return;
}

/*
 * Solver Plugin Interface
 */
#ifdef _STANDALONE
SolverPlugin_t* solver_init() {
#else
SolverPlugin_t* elmer_cg_init() {
#endif
	SolverPlugin_t* solver = (SolverPlugin_t*)malloc(sizeof(SolverPlugin_t));

	solver->set_option = NULL;
	solver->solve_pre = solve_pre;
	solver->solve = solve;
	solver->solve_post = solve_post;
	solver->free = solve_free;

	return solver;
}

#if 0
  recursive subroutine  huti_dcgsolv  ( ndim, wrkdim, xvec, rhsvec, &
       ipar, dpar, work, matvecsubr, pcondlsubr, &
       pcondrsubr, dotprodfun, normfun, stopcfun )

    use huti_interfaces
    implicit none

    procedure( mv_iface_d ), pointer :: matvecsubr
    procedure( pc_iface_d ), pointer :: pcondlsubr
    procedure( pc_iface_d ), pointer :: pcondrsubr
    procedure( dotp_iface_d ), pointer :: dotprodfun
    procedure( norm_iface_d ), pointer :: normfun
    procedure( stopc_iface_d ), pointer :: stopcfun

    ! Parameters

    integer :: ndim, wrkdim
    double precision, dimension(ndim) :: xvec, rhsvec
    integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
    double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar
    double precision, dimension(ndim,wrkdim) :: work

    ! Local variables

    double precision :: alpha, beta, rho, oldrho
    integer iter_count, i
    double precision :: residual, rhsnorm, precrhsnorm

    !
    ! End of variable declarations
    !*********************************************************************

    !*********************************************************************
    ! The actual CG begins here (look the pseudo code in the
    ! "Templates..."-book, page 15)
    !
    ! First the initialization part
    !

    iter_count = 1

    ! Norms of right-hand side vector are used in convergence tests

    if ( HUTI_STOPC .eq. HUTI_TRESID_SCALED_BYB .or. &
         HUTI_STOPC .eq. HUTI_PRESID_SCALED_BYB ) then
       rhsnorm = normfun( HUTI_NDIM, B, 1 )
    end if
    if ( HUTI_STOPC .eq. HUTI_PRESID_SCALED_BYPRECB ) then
       call pcondlsubr( P, B, ipar )
       precrhsnorm = normfun( HUTI_NDIM, P, 1 )
    end if

    ! The following applies for all matrix operations in this solver

    HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED

    ! Generate vector X if needed

    if ( HUTI_INITIALX .eq. HUTI_RANDOMX ) then
       call  huti_drandvec   ( X, ipar )
    else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
       X = 1
    end if

    call matvecsubr( X, R, ipar )

    R = B - R

    !
    ! This is where the loop starts (that is we continue from here after
    ! the first iteration)
    !

300 continue

    call pcondlsubr( Q, R, ipar )
    call pcondrsubr( Z, Q, ipar )

    rho = dotprodfun( HUTI_NDIM, R, 1, Z, 1 )
    if ( rho .eq. 0 ) then
       HUTI_INFO = HUTI_CG_RHO
       go to 1000
    end if

    if ( iter_count .eq. 1 ) then
       P = Z
     else
       beta = rho / oldrho
       P = Z + beta * P
    end if

    call matvecsubr( P, Q, ipar )

    alpha = rho / dotprodfun( HUTI_NDIM, P, 1, Q, 1 )

    X = X + alpha * P
    R = R - alpha * Q

    !
    ! Check the convergence against selected stopping criterion
    !

    select case (HUTI_STOPC)
    case (HUTI_TRUERESIDUAL)
       call matvecsubr( X, Z, ipar )
       Z = Z - B
       residual = normfun( HUTI_NDIM, Z, 1 )
    case (HUTI_TRESID_SCALED_BYB)
       call matvecsubr( X, Z, ipar )
       Z = Z - B
       residual = normfun( HUTI_NDIM, Z, 1 ) / rhsnorm
    case (HUTI_PSEUDORESIDUAL)
       residual = normfun( HUTI_NDIM, R, 1 )
    case (HUTI_PRESID_SCALED_BYB)
       residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
    case (HUTI_PRESID_SCALED_BYPRECB)
       residual = normfun( HUTI_NDIM, R, 1 ) / precrhsnorm
    case (HUTI_XDIFF_NORM)
       Z = alpha * P
       residual = normfun( HUTI_NDIM, Z, 1 )
    case (HUTI_USUPPLIED_STOPC)
       residual = stopcfun( X, B, R, ipar, dpar )
    case default
       call matvecsubr( X, Z, ipar )
       Z = Z - B
       residual = normfun( HUTI_NDIM, Z, 1 )
    end select

    !
    ! Print debugging output if desired
    !

    if ( HUTI_DBUGLVL .ne. HUTI_NO_DEBUG ) then
       if ( mod(iter_count, HUTI_DBUGLVL) .eq. 0 ) then
          write (*, '(I8, E11.4)') iter_count, residual
       end if
    end if

    if ( residual .lt. HUTI_TOLERANCE ) then
       HUTI_INFO = HUTI_CONVERGENCE
       go to 1000
    end if

    IF( residual /= residual .OR. residual > HUTI_MAXTOLERANCE ) THEN
      HUTI_INFO = HUTI_DIVERGENCE
      GOTO 1000
    END IF


    oldrho = rho

    !
    ! Return back to the iteration loop (without initialization)
    !

    iter_count = iter_count + 1
    if ( iter_count .gt. HUTI_MAXIT ) then
       HUTI_INFO = HUTI_MAXITER
       go to 1000
    end if

    go to 300

    !
    ! This is where we exit last time (after enough iterations or breakdown)
    !

1000 continue
    if ( HUTI_DBUGLVL .ne. HUTI_NO_DEBUG ) then
       write (*, '(I8, E11.4)') iter_count, residual
    end if

    HUTI_ITERS = iter_count
    return

    ! End of execution
    !*********************************************************************

  end subroutine  huti_dcgsolv
#endif
