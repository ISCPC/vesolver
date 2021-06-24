/*

  pcgsolver.c	--	The preconditioned conjugate gradient solver(s)
  
  				Implemented solver/preconditioner:
  				    preconditioned cg-solver with
  				       * diagonal scaling only
  				       * incomplete Cholesky on full matrix
  				
  				Most functions are based on:
  				H.R. Schwarz: FORTRAN-Programme zur Methode 
                                             der finiten Elemente,
  				              Teubner, 1981 
 
                               The present version is based on the c-code by
                               Martin Ruecker and Ernst Rank of the Technical
                               University of Munich 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"
#include "PluginAPI.h"

/*
 * CG Solver
 */
#define ITG int32_t
#define NNEW(a,b,c) a=(b *)calloc((c),sizeof(b))
#define RENEW(a,b,c) a=(b *)realloc((b *)(a),(c)*sizeof(b))
#define SFREE(a) free(a)


/*--(parallel) conjugate gradient solver-------------------------------------------	*/
/*---------------------------------------------------------------------------------	*/
/*--parameter:                                                                   --	*/
/*--           A       compact row oriented storage of lower left of the scaled  --	*/
/*--                   matrix A                                                  --	*/
/*--           b       right hand side                                           --	*/
/*--           x       solution vector                                           --	*/
/*--           eps     required accuracy -> residual                             --	*/
/*--           niter   maximum number of iterations -> number of iterations      --	*/
/*---------------------------------------------------------------------------------	*/
static void CG (const Matrix_t *A, const double *b, double *x, double *eps, ITG *niter) {
    ITG	i=0, k=0, iam;
    double ekm1=0.0,c1=0.005,qam,ram=0.,err;
    double rr=0.0, pz=0.0, qk=0.0, rro=0.0;
    double *r, *p, *z;

    ITG neq = A->NROWS;

    NNEW(r,double,neq);
    NNEW(p,double,neq);
    NNEW(z,double,neq);

    /* solving the system of equations using the iterative solver */
    printf("Solving the system of equations using the iterative solver\n\n");

    /*..initialize result, search and residual vectors.................................	*/
    qam=0.;iam=0;
    for (i=0; i<neq; i++) {
        x[i] =  0.0;					/*..start vector x0 = 0....................	*/
        r[i] =  b[i];					/*..residual vector r0 = Ax+b = b..........	*/
        p[i] = -r[i];					/*..initial search vector..................	*/
        err=fabs(r[i]);
        if(err>1.e-20) {
            qam+=err;
            iam++;
        }
    }

    if(iam>0) {
        qam=qam/neq;
    } else {
        *niter=0;
        SFREE(r);SFREE(p);SFREE(z);
        return;
    }

    /*..main iteration loop............................................................	*/
    for (k=1; k<=(*niter); k++) {
        if ((ram<=c1*qam) && k!=1) break;
        rr = DDOT(neq, r, r);
#if 1
        if ((k % 100) == 0) {
            printf("iteration= %d, error= %e, limit=%e\n",k,ram,c1*qam);
        }
#endif

        /*......new search vector..........................................................	*/
        if (k!=1) {
            ekm1 = rr/rro;
            for (i=0; i<neq; i++)	p[i] = ekm1*p[i]-r[i];
        }
        /*......matrix vector product A p = z..............................................	*/
        Matrix_MV(A, 1.0, p, 0.0, z, z);

        /*......inner product pT z.........................................................	*/
        pz = DDOT(neq, p, z);

        /*......step size..................................................................	*/
        qk = rr/pz;

        /*......approximated solution vector, residual vector..............................	*/
        ram=0.;
        for (i=0; i<neq; i++) {
            x[i] = x[i] + qk*p[i];
            r[i] = r[i] + qk*z[i];
            err=fabs(r[i]);
            if(err>ram) ram=err;
        }

        /*......store previous residual....................................................	*/
        rro = rr;
    }
    if(k==*niter) {
        printf("*ERROR in CG: no convergence;");
        exit(1);
    } 
    *eps = rr;						/*..return residual............................	*/
    *niter = k;					/*..return number of iterations................	*/
    /*..That's it......................................................................	*/

    SFREE(r);SFREE(p);SFREE(z);
    return;
}


/*
 * API
 */
static Matrix_t* solve_pre(const Matrix_t* A0) {
    Matrix_t* A = Matrix_duplicate(A0);
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

    A->info = (void*)factor;
    return A;
}

static int solve(const Matrix_t *A, const double* b, double* x, const double res) {
    ITG neq = A->NROWS;
    double eps = res;

    printf("INFO: cgsolver neq=%d, len=%d, res=%lf\n", neq, A->NNZ, res);

    double *factor = (double*)(A->info);
    double *b1 = (double*)calloc(sizeof(double), neq);

    /*  scale right hand side (Ax=b -> Ax+b=0: negative sign)  */
    for (ITG i=0; i<neq; i++) {
        b1[i] = b[i] * -factor[i];
    }

    /*  SOLVER/PRECONDITIONING TYPE  */
    /*  Conjugate gradient solver without preconditioning  */
    ITG niter=5000000;
    CG(A, b1, x, &eps, &niter);

    /*  Backscaling of the solution vector  */
    for (ITG i=0; i<neq; i++) {
        x[i] *= factor[i];
    }

    printf("# of iterations = %d , eps=%e\n", niter, eps);

    free(b1);
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
SolverPlugin_t* cg_init() {
	SolverPlugin_t* solver = (SolverPlugin_t*)malloc(sizeof(SolverPlugin_t));

	solver->set_option = NULL;
	solver->solve_pre = solve_pre;
	solver->solve = solve;
	solver->solve_post = solve_post;
	solver->free = solve_free;

	return solver;
}
