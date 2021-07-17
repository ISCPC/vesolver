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


/*--Preconditioning of the equation system Ax+b=0----------------------------------	*/
/*---------------------------------------------------------------------------------	*/
/*--Partial Cholesky decomposition of scaled matrix A. The off-diagonal matrix   --	*/
/*--elements are divided by 1/(1+alpha) until a decomposition of the positive    --	*/
/*--definit matrix A exists. C is obtained by ignoring the fill-in during the    --	*/
/*--decomposition. In case of successfull decomposition the preconditioning      --	*/
/*--matrix C is returned. Otherwise the function is called again with new alpha. --	*/
/*--alpha has to be chosen as small as possible, because preconditioning effect  --	*/
/*--decreases with increasing alpha.-----------------------------------------------	*/
/*---------------------------------------------------------------------------------	*/
/*--parameter:                                                                   --	*/
/*--           A       compact row oriented storage of lower left of matrix A    --	*/
/*--           neq     order of A, number of equations                           --	*/
/*--           ia      column indices                                            --	*/
/*--           iz      row indices (diagonal elements indices)                   --	*/
/*--           alpha   see above                                                 --	*/
/*---------------------------------------------------------------------------------	*/
/*--The function corresponds to function PARTCH() in H.R. Schwarz: FORTRAN-Pro-  --	*/
/*--gramme zur Methode der finiten Elemente, p.117, Teubner, 1981                --	*/
/*---------------------------------------------------------------------------------	*/
static ITG PreConditioning (const Matrix_t *A, const double alpha, Matrix_t *C) {
    double          factor = 1.0/(1.0+alpha);
    ITG neq = A->NROWS;
    ITG *iz = A->pointers;
    ITG *ia = A->indice;
    double *aval = A->values;
    double *cval = C->values;

    C->NROWS = A->NROWS;
    C->NNZ = A->NNZ;
    C->pointers = A->pointers;
    C->indice = A->indice;
#if 0
    /*.division of the off-diagonal elements of A by factor (1.0+alpha)............... */
    cval[0] = aval[0];
    for (int i=1; i<neq; i++) {
        ITG jlo = iz[i];       /*..first non-zero element in current row...... */
        ITG jup = iz[i+1]-1;   /*..diagonal element........................... */
        cval[jup] = aval[jup];            /*..copy of diagonal element................... */
        for (int j=jlo; j<jup; j++)     /*..all non-zero off-diagonal elements......... */
            cval[j] = aval[j] * factor;
    }

    /*..partial Cholesky decomposition of scaled matrix A.............................. */
    for (int i=1; i<neq; i++) {   
        ITG jlo = iz[i];       /*..first non-zero element in current row...... */
        ITG jup = iz[i+1]-1;   /*..diagonal element........................... */
        for (int j=jlo; j<jup; j++) {   /*..all non-zero off-diagonal elements......... */
            cval[j] /= cval[iz[ia[j]+1]-1];
            ITG klo = j+1;              /*..next non-zero element in current row....... */
            ITG kup = jup;              /*..diagonal element in current row............ */
            for (int k=klo; k<=kup; k++) {
                ITG m = ia[k];
                ITG llo = iz[m];
                ITG lup = iz[m+1]-1;
                for (int l=llo; l<=lup; l++) {
                    if (ia[l]>ia[j]) break;
                    if (ia[l]<ia[j]) continue;
                    cval[k] -= cval[j]*cval[l];
                    break;
                }
            } 
        }
        ITG id = iz[i+1]-1;
        if (cval[id]<1.0e-6){
            return 0;
        }
        cval[id] = sqrt(cval[id]);
    }
#else
    /*.division of the off-diagonal elements of A by factor (1.0+alpha)............... */
    for (int i=0; i<neq; i++) {
        ITG jlo = iz[i];                /*..first non-zero element in current row...... */
        ITG jup = iz[i+1]-1;            /*..diagonal element........................... */ 
        cval[jup] = aval[jup];          /*..copy of diagonal element................... */
        for (int j=jlo; j<jup; j++) { /*..all non-zero off-diagonal elements......... */
            cval[j] = aval[j] * factor;
        }
    }
    //printf("INFO:PreConditioning: 1st phase done.\n");
    //fflush(stdout);

    /*..partial Cholesky decomposition of scaled matrix A.............................. */
    for (int i=1; i<neq; i++) {   
        ITG jlo = iz[i];       /*..first non-zero element in current row...... */
        ITG jup = iz[i+1]-1;   /*..diagonal element........................... */
        for (int j=jlo; j<jup; j++) {   /*..all non-zero off-diagonal elements......... */
            cval[j] /= cval[iz[ia[j]+1]-1];
            ITG klo = j+1;              /*..next non-zero element in current row....... */
            ITG kup = jup;              /*..diagonal element in current row............ */
            for (int k=klo; k<=kup; k++) {
                ITG m = ia[k];
                ITG llo = iz[m];
                ITG lup = iz[m+1]-1;
                for (int l=llo; l<=lup; l++) {
                    if (ia[l] >= ia[j]) {
                        if (ia[l] == ia[j]) {
                            cval[k] -= cval[j]*cval[l];
                        }
                        break;
                    }
                }
            } 
        }
        ITG id = iz[i+1]-1;
        if (cval[id]<1.0e-6){
            //printf("INFO:PreConditioning: 2nd phase done with 0. C[%d]=%lf\n", id, cval[id]);
            //fflush(stdout);
            return 0;
        }
        cval[id] = sqrt(cval[id]);
    }
    //printf("INFO:PreConditioning: 2nd phase done with 1.\n");
    //fflush(stdout);
#endif
    return 1;
}

/*--Solution of the equation system  M rho = r-------------------------------------	*/
/*---------------------------------------------------------------------------------	*/
/*--The solution of the equation system M rho = r, with M = C*CT, where C is the --	*/
/*--partial Cholesky decomposition of matrix A, represent a preconditioning step.--	*/
/*--The equation  system is solved by forward/backward substitution.---------------	*/
/*---------------------------------------------------------------------------------	*/
/*--parameter:                                                                   --	*/
/*--           C       compact row oriented storage of preconditioned matrix A   --	*/
/*--           neq     order of A, number of equations                           --	*/
/*--           ia      column indices                                            --	*/
/*--           iz      row indices (diagonal elements indices)                   --	*/
/*--           r       residual vector                                           --	*/
/*---------------------------------------------------------------------------------	*/
/*--The function corresponds to function CRHOER() in H.R. Schwarz: FORTRAN-Pro-  --	*/
/*--gramme zur Methode der finiten Elemente, p.118, Teubner, 1981                --	*/
/*---------------------------------------------------------------------------------	*/
//void Mrhor (const double *C, const ITG neq, const ITG *aind, const ITG *aptr, const double *r, double *rho) {
static void Mrhor (const Matrix_t *C, const double *r, double *rho) {
    ITG neq = C->NROWS;
    ITG *iz = C->pointers;
    ITG *ia = C->indice;
    double *cval = C->values;

    /*..solve equation sytem by forward/backward substitution.......................... */
    rho[0] = r[0];
    for (int i=1; i<neq; i++) {
        double s = 0.0;
        int jlo = iz[i];            /*..first non-zero element in current row...... */
        int jup = iz[i+1]-1;        /*..diagonal element in current row............ */
        for (int j=jlo; j<jup; j++) { /*..all non-zero off-diagonal element.......... */
            s += cval[j] * rho[ia[j]];
        }
        rho[i] = (r[i]-s) / cval[jup];
    }

    for (int i=neq-1; i>0; i--) {
        int jlo = iz[i];            /*..first non-zero element in current row...... */
        int jup = iz[i+1]-1;              /*..diagonal element in current row............ */
        rho[i] /= cval[jup];
        for (int j=jlo; j<jup; j++)    /*..all non-zero off-diagonal element.......... */
            rho[ia[j]] -= cval[j] * rho[i];
    }
    return;
}

//void PCG (double *A, double *x, double *b, ITG neq, ITG len, ITG *ia, 
//	  ITG *iz,double *eps, ITG *niter, ITG precFlg,
//	  double *rho, double *r, double *g, double *C, double *z) {
static void PCG (const Matrix_t *A, double *b, double *x, double *eps, ITG *niter) {
    ITG     i=0, k=0, iam;
    double  alpha=0.0, ekm1=0.0, rrho=0.0;
    double  rrho1=0.0, gz=0.0, qk=0.0;
    double  c1=0.005,qam,err,ram=0;

    double *rho, *r, *g, *z;
    Matrix_t C;

    ITG neq = A->NROWS;
    ITG len = A->NNZ;

    NNEW(r, double, neq);

    /*  initialize result and residual vectors  */
    qam=0.;iam=0;
    for (i=0; i<neq; i++) {
        x[i] = 0.0;	  
        r[i] = b[i];
        err=fabs(r[i]);
        if(err>1.e-20) { qam += err; iam++; }
    }

    if(iam>0) {
        qam=qam/iam;
    } else {
        *niter=0;
        SFREE(r);
        return;
    }

    NNEW(rho, double, neq);
    NNEW(g, double, neq);
    NNEW(z, double, neq);
    NNEW(C.values, double, len);

    /*  preconditioning  */
    printf("Cholesky preconditioning\n\n");

    printf("alpha=%f\n",alpha);
    int cnt=0;
    while (PreConditioning(A, alpha, &C) == 0) {
        alpha = (alpha == 0.0) ? 0.005 : (alpha * 2);
        printf("alpha=%f\n",alpha);
        if (cnt > 10) {
            printf("ERROR: Abnormal end\n");
            exit(1);
        }
        cnt++;
    }

    /* solving the system of equations using the iterative solver */
    printf("Solving the system of equations using the iterative solver\n\n");

    /*  main iteration loop  */
    for (k=1; k<=*niter; k++) {
        if ((ram<=c1*qam) && (k!=1)) {
            printf("DONE: iteration= %d, error= %e\n", k, ram);
            break;
        }

        /*  solve M rho = r, M=C CT  */
        Mrhor(&C, r, rho);

        /*  inner product (r,rho)  */
        //InnerProduct(r,rho,neq,&rrho);
        rrho = DDOT(neq, r, rho);

        /*  If working with Penalty-terms for boundary conditions you can get 
            numerical troubles because RRHO becomes a very large value. 
            With at least two iterations the result may be better !!! */

        /*  convergency check */
#if 1
        if ((k % 100) == 0) {
            //printf("iteration= %d, error= %e, limit=%e\n", k, ram, c1*qam);
            printf("iteration= %d, error= %e, limit=%e, rrho=%lf, gz=%lf\n", k, ram, c1*qam, rrho, gz);
        }
#endif

        if (k!=1) {
            ekm1 = rrho/rrho1;
            for (i=0; i<neq; i++) g[i] = ekm1*g[i]-rho[i];
        } else {
            /*  g1 = -rho0  */
            for (i=0; i<neq; i++) g[i] = -rho[i];
        }

        /*  solve equation system z = A g_k  */
        //MatVecProduct(A,g,neq,ia,iz,z);
        Matrix_MV(A, 1.0, g, 0.0, z, z);

        /*  inner product (g,z)  */
        //InnerProduct(g,z,neq,&gz);
        gz = DDOT(neq, g, z);

        qk = rrho/gz;
        ram=0.;
        for (i=0; i<neq; i++) {
            x[i] += qk*g[i];
            r[i] += qk*z[i];
            err=fabs(r[i]);
            if(err>ram) ram=err;
        }
        rrho1 = rrho;
    }

    if(k==*niter) {
        printf("*ERROR in PCG: no convergence;");
        //FORTRAN(stop,());
        return;
    } 
    *eps = rrho;
    *niter = k;

    SFREE(rho); SFREE(r); SFREE(g); SFREE(C.values); SFREE(z);
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

    if (Matrix_optimize(A) < 0) {
        free(A);
        return NULL;
    }

    Matrix_transpose(A);
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
    PCG(A, b1, x, &eps, &niter);

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
#ifdef _STANDALONE
SolverPlugin_t* solver_init() {
#else
SolverPlugin_t* pcg_init() {
#endif
	SolverPlugin_t* solver = (SolverPlugin_t*)malloc(sizeof(SolverPlugin_t));

	solver->set_option = NULL;
	solver->solve_pre = solve_pre;
	solver->solve = solve;
	solver->solve_post = solve_post;
	solver->free = solve_free;

	return solver;
}
