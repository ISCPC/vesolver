#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"

/****************************************************************
 * Common Functions
 ****************************************************************/
void Matrix_setMatrixCSR(Matrix_t *A, const int nrows, const int nnz, const int *aptr,
    const int *aind, const double *aval, const uint32_t flags) {
    A->NROWS = nrows;
    A->NNZ = nnz;
    A->flags = flags;

    A->pointers = (int*)calloc(sizeof(int), A->NROWS+1);
    A->indice = (int*)calloc(sizeof(int), A->NNZ);
    A->values = (double*)calloc(sizeof(double), A->NNZ);

    bcopy(aptr, A->pointers, sizeof(int)*(A->NROWS+1));
    bcopy(aind, A->indice, sizeof(int)*(A->NNZ));
    bcopy(aval, A->values, sizeof(double)*(A->NNZ));

    A->optimized = 0;
}

Matrix_t* Matrix_duplicate(const Matrix_t* A) {
    Matrix_t* D = malloc(sizeof(Matrix_t));
    D->NROWS = A->NROWS;
    D->NNZ = A->NNZ;
    D->flags = A->flags;
    D->info = NULL;

    D->pointers = (int*)calloc(sizeof(int), A->NROWS+1);
    D->indice = (int*)calloc(sizeof(int), A->NNZ);
    D->values = (double*)calloc(sizeof(double), A->NNZ);

    bcopy(A->pointers, D->pointers, sizeof(int)*(A->NROWS+1));
    bcopy(A->indice, D->indice, sizeof(int)*(A->NNZ));
    bcopy(A->values, D->values, sizeof(double)*(A->NNZ));

#ifdef MKL
    D->hdl = 0;
#endif
#ifdef SXAT
    D->hdl = 0;
#endif
#ifdef SSL2
    D->w = NULL;
    D->iw = NULL;
#endif

    return D;
}

int Matrix_convert_index(Matrix_t* A, int base) {
    if ((MATRIX_INDEX_TYPE(A) == 1)&&(base == 0)) {
        /* Convert from one-based index to zero-based index */
        for (int i=0; i<=A->NROWS; i++) {
            A->pointers[i] = A->pointers[i]-1;
        }

        for (int i=0; i<A->NNZ; i++) {
            A->indice[i] = A->indice[i]-1;
        }

        A->flags = ((A->flags)&(~0xf)) | MATRIX_TYPE_INDEX0;
    } else if ((MATRIX_INDEX_TYPE(A) == 0)&&(base == 1)) {
        /* Convert from zero-based index to one-based index */
        for (int i=0; i<=A->NROWS; i++) {
            A->pointers[i] = A->pointers[i]+1;
        }

        for (int i=0; i<A->NNZ; i++) {
            A->indice[i] = A->indice[i]+1;
        }

        A->flags = ((A->flags)&(~0xf)) | MATRIX_TYPE_INDEX1;
    }

    return 0;
}

int Matrix_transpose(Matrix_t* A) {
    int* aptr = (int*)calloc(sizeof(int), A->NROWS+1);
    int* aind = (int*)calloc(sizeof(int), A->NNZ);
    double* aval = (double*)calloc(sizeof(double), A->NNZ);
#if 0
    /* phase0: */
    for(int i=0; i<=A->NROWS; i++) {
        aptr[i] = 0;
    }
#endif
    /* phase1 */
    for(int i=0; i<A->NNZ; i++) {
        int j = A->indice[i];
        aptr[j]++;
    }
    for(int i=1; i<=A->NROWS; i++) {
        aptr[i] += aptr[i-1];
    }

    /* phase2 */
    for (int i=A->NROWS-1; i>=0; i--) {
        for (int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
            int jj = --(aptr[A->indice[j]]);
            aind[jj] = i;
            aval[jj] = A->values[j];
        }
    }

    free(A->pointers);
    free(A->indice);
    free(A->values);
    A->pointers = aptr;
    A->indice = aind;
    A->values = aval;

    if (MATRIX_IS_SYMMETRIC(A)) {
        A->flags ^= MATRIX_TYPE_LOWER;
    } else {
        A->flags ^= MATRIX_TYPE_CSR;
    }

    return 0;
}

static int csr2ellpack_sym(Matrix_t* A) {
    ellpack_info_t* info = &(A->_ellpack);
    int32_t neq = A->NROWS;
    int32_t* iaptr = A->pointers;
    int32_t* iaind = A->indice;
    double* aval = A->values;
    int ncol=0;

    //printf("*** 1st step: Analyze ***\n");
    int* nrows = (int*)calloc(sizeof(int), neq);
    for(int i=0; i<neq; i++) {
        nrows[i] = 0;
    }
    for(int i=0; i<neq; i++) {
        int n = iaptr[i+1]-iaptr[i]-1;
        //printf("num[%d] = %d\n", i, n);
        if (n > ncol) ncol = n;
        for(int j=iaptr[i]+1; j<iaptr[i+1]; j++) {
            nrows[iaind[j]]++;
        }
    }
    int maxrow=0;
    for(int i=0; i<neq; i++) {
        if (nrows[i] > maxrow) maxrow = nrows[i];
    }
    //printf("ncol = %d, maxrow = %d\n", ncol, maxrow);
    if (maxrow > ncol) ncol = maxrow;
    free(nrows);
    int32_t nw = info->nw = ncol*2;

    //printf("*** 2nd step: Fill matrix *** \n");
    info->COEF = (double*)calloc(sizeof(double), neq*nw);
    info->ICOL = (int32_t*)calloc(sizeof(int32_t), neq*nw);

    // Fill matrix
    int* ku = (int*)calloc(sizeof(int), neq);
    int* kl = (int*)calloc(sizeof(int), neq);
    for(int i=0; i<neq; i++) {
        ku[i] = 0;
        kl[i] = 0;
        for(int j=0; j<nw; j++) {
            info->COEF[j*neq+i] = 0.0f;
            info->ICOL[j*neq+i] = i+1;
        }
    }

    for(int i=0; i<neq; i++) {
        for(int j=iaptr[i]+1; j<iaptr[i+1]; j++) {
            // Fill Upper Matrix
            // COEFF(i, ku[i]) = aval[j]
            // ICOL(i, ku[i]) = iaind[j]
            int k = ku[i]*neq + i;
            info->COEF[k] = aval[j];
            info->ICOL[k] = iaind[j]+1; // one-based index
            ku[i]++;

            // Fill Lower Matrix
            // COEFF(iaind[j], kl[iaind[j]]) = aval[j]
            // ICOL(iaind[j], kl[iaind[j]]) = i
            k = (kl[iaind[j]]+ncol)*neq + iaind[j];
            info->COEF[k] = aval[j];
            info->ICOL[k] = i+1; // one-based index
            kl[iaind[j]]++;
        }
    }
    free(ku);
    free(kl);

    return 0;
}

static int csr2ellpack_asym(Matrix_t* A) {
    ellpack_info_t* info = &(A->_ellpack);
    int32_t neq = A->NROWS;
    int32_t* iaptr = A->pointers;
    int32_t* iaind = A->indice;
    double* aval = A->values;
    int32_t nw=0;
    
    /* Convert Asymmetric CSR Matrix to ELLPACK(FORTRAN) */
    //printf("*** step 1: Analyze. ***\n");
    for(int i=0; i<neq; i++) {
        int n = iaptr[i+1]-iaptr[i];
        if (n > nw) nw = n;
    }
    info->nw = nw;
    //printf("INFO: nw = %d\n", info->nw);

    //printf("*** step 2: Fill matrix *** \n");
    info->COEF = (double*)calloc(sizeof(double), neq*nw);
    info->ICOL = (int32_t*)calloc(sizeof(int32_t), neq*nw);

    // Fill matrix
    int* ku = (int*)calloc(sizeof(int), neq);
    for(int i=0; i<neq; i++) {
        ku[i] = 0;
        for(int j=0; j<nw; j++) {
            info->COEF[j*neq+i] = 0.0f;
            info->ICOL[j*neq+i] = i+1;
        }
    }

    for(int i=0; i<neq; i++) {
        for(int j=iaptr[i]; j<iaptr[i+1]; j++) {
            // Fill Matrix
            // COEFF(i, ku[i]) = aval[j]
            // ICOL(i, ku[i]) = iaind[j]
            int k = ku[i]*neq + i;
            info->COEF[k] = aval[j];
            info->ICOL[k] = iaind[j]+1; // one-based index
            ku[i]++;
        }
    }
    free(ku);

    return 0;
}

int Matrix_create_ellpack(Matrix_t* A) {
    Matrix_convert_index(A, 0);

    if (MATRIX_IS_SYMMETRIC(A)) {
        csr2ellpack_sym(A);
    } else {
        csr2ellpack_asym(A);
    }
    return 0;
}

/*
 * Convert CSR to JAD
 */
struct sort_elem {
    struct sort_elem* next;
    int row;
};

struct sort_list {
    struct sort_elem* top;
    struct sort_elem* tail;
    int nelems;
};

int Matrix_create_jad(Matrix_t* A) {
    jad_info_t* jad = &(A->_jad);
    Matrix_convert_index(A, 0);

    /* Allocate temporary array */
    struct sort_list* list = (struct sort_list*)calloc(sizeof(struct sort_list), A->NROWS+1);
    struct sort_elem* elem = (struct sort_elem*)calloc(sizeof(struct sort_elem), A->NROWS);

    if (MATRIX_IS_SYMMETRIC(A)) {
        printf("ERROR: Sorry, Symmetric matrix is not support yet.");
        return -1;
    } else {
        /* Allocate JAD */
        jad->perm = (int32_t*)calloc(sizeof(int32_t), A->NROWS);
        jad->index = (int32_t*)calloc(sizeof(int32_t), A->NNZ);
        jad->value = (double*)calloc(sizeof(double), A->NNZ);

        int maxnzr = 0;
        for (int i=0; i<A->NROWS; i++) {
            elem[i].row = i;
            int ii = A->pointers[i+1]-A->pointers[i];

            (list[ii].nelems)++;
            if (list[ii].top == NULL) {
                list[ii].top = list[ii].tail = &(elem[i]);
            } else {
                list[ii].tail = list[ii].tail->next = &(elem[i]);
            }
            if (ii > maxnzr) {
                maxnzr = ii;
            }
        }

        int i=0;
        //printf("INFO: maxnzr = %d\n", maxnzr);
        for (int j=maxnzr; j>=0; j--) {
            //printf("INFO: j=%d, i=%d\n", j, i);
            for (struct sort_elem* e=list[j].top; e!=NULL; e=e->next) {
                jad->perm[i] = e->row;
                i++;
            }
        }

        int n = A->NROWS;
        int k = 0;
        jad->ptr = (int32_t*)calloc(sizeof(int32_t), maxnzr+1);
        for (int i=0; i<=maxnzr; i++) {
            n -= list[i].nelems;
            jad->ptr[i] = k;
            //printf("INFO: i, k = %d, %d -------\n", i, k);
            for (int j=0; j<n; j++) {
                int ii = A->pointers[jad->perm[j]]+i;
                jad->index[k] = A->indice[ii];
                jad->value[k] = A->values[ii];
                //printf("INFO: (%d, %d) = %lf\n", jad->perm[j], jad->index[k], jad->value[k]);
                k++;
            }
        }
        jad->maxnzr = maxnzr;
    }
    //fflush(stdout);

    free(elem);
    free(list);

    return 0;
}

/****************************************************************
 * Base Library
 ****************************************************************/
#pragma weak Matrix_init = Matrix_init_generic
void Matrix_init_generic(Matrix_t *A) {
    A->pointers = NULL;
    A->indice = NULL;
    A->values = NULL;
    A->info = NULL;
}

#pragma weak Matrix_free = Matrix_free_generic
void Matrix_free_generic(Matrix_t *A) {
    if (A->pointers) {
        free(A->pointers);
        A->pointers = NULL;
    }

    if (A->indice) {
        free(A->indice);
        A->indice = NULL;
    }
    if (A->values) {
        free(A->values);
        A->values = NULL;
    }
    if (A->info) {
        free(A->info);
        A->info = NULL;
    }

    A->optimized = 0;
}

#pragma weak Matrix_optimize = Matrix_optimize_generic
int Matrix_optimize_generic(Matrix_t *A) {
    A->optimized = 1;
    return 0;
}

#pragma weak Matrix_MV = Matrix_MV_generic
int Matrix_MV_generic(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z) {
    if (MATRIX_INDEX_TYPE(A) == 1) {
        /*  matrix vector product  */
        if (MATRIX_IS_SYMMETRIC(A)) {
            for (int i=A->NROWS-1; i>=0; i--) {
                z[i] = (beta * y[i]) + (A->values[A->pointers[i]-1]* x[i]);
                for (int j=A->pointers[i]; j<A->pointers[i+1]-1; j++) {
                    int k = A->indice[j]-1;
                    z[i] += alpha * A->values[j] * x[k];
                    z[k] += alpha * A->values[j] * x[i];
                }
            }
        } else {
            #pragma omp parallel for
            for(int i=0; i<A->NROWS; i++) {
                z[i] = beta * y[i];
                for(int j=A->pointers[i]-1; j<A->pointers[i+1]-1; j++) {
                    z[i] += alpha * A->values[j] * x[A->indice[j]-1];
                }
            }
        }
    } else {
        /*  matrix vector product  */
        if (MATRIX_IS_SYMMETRIC(A)) {
            for (int i=A->NROWS-1; i>=0; i--) {
                z[i] = (beta * y[i]) + (A->values[A->pointers[i]]* x[i]);
                for (int j=A->pointers[i]+1; j<A->pointers[i+1]; j++) {
                    int k = A->indice[j];
                    z[i] += A->values[j] * x[k];
                    z[k] += A->values[j] * x[i];
                }
            }
        } else {
#if 1
            #pragma omp parallel for
            for(int i=0; i<A->NROWS; i++) {
                z[i] = beta * y[i];
                for(int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
                    z[i] += alpha * A->values[j] * x[A->indice[j]];
                }
            }
#else
            if (beta != 0.0f) {
                if (alpha == -1.0f) {
                    #pragma omp parallel for
                    for(int i=0; i<A->NROWS; i++) {
                        z[i] = beta * y[i];
                        for(int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
                            z[i] -= A->values[j] * x[A->indice[j]];
                        }
                    }
                } else {
                    #pragma omp parallel for
                    for(int i=0; i<A->NROWS; i++) {
                        z[i] = beta * y[i];
                        for(int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
                            z[i] += alpha * A->values[j] * x[A->indice[j]];
                        }
                    }
                }
            } else {
                if (alpha == 1.0f) {
                    int j0;
                    double t0;
                    int* jj0 = A->indice;
                    double* vv0 = A->values;
                    #pragma omp parallel for private(j0,t0)
                    for(int i=0; i<A->NROWS; i++) {
                        t0 = 0.0f;
                        int is = A->pointers[i];
                        int ie = A->pointers[i+1];

                        for(int j=is; j<ie; j++) {
                            j0 = jj0[j];
                            t0 += vv0[j] * x[j0];
                            //t0 += vv0[j] * x[jj0[j]];
                        }
                        z[i] = t0;
                    }
                } else {
                    #pragma omp parallel for
                    for(int i=0; i<A->NROWS; i++) {
                        z[i] = 0.0f;
                        for(int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
                            z[i] += alpha * A->values[j] * x[A->indice[j]];
                        }
                    }
                }
            }
#endif
        }
    }
    return 0;
}
