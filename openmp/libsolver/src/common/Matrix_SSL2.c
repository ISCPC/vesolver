#ifdef SSL2
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"
#include "cssl.h"

#define ELLPACK 1

/****************************************************************
 * Library with Fujitsu SSL2
 ****************************************************************/
void Matrix_init(Matrix_t *A) {
    A->pointers = NULL;
    A->indice = NULL;
    A->values = NULL;
    A->info = NULL;
    A->optimized = 0;

    A->w = NULL;
    A->iw = NULL;
}

void Matrix_free(Matrix_t *A) {
    /* Destruction of the handle */
    if (A->optimized == 1) {
        A->optimized = 0;
    }

    Matrix_free_generic(A);

    if (A->w != NULL) {
        free(A->w);
        A->w = NULL;
    }
    if (A->iw != NULL) {
        free(A->iw);
        A->iw = NULL;
    }
}

#ifdef ELLPACK
int Matrix_optimize(Matrix_t *A) {
    Matrix_create_ellpack(A);

    A->optimized = 1;
    return 0;
}

int Matrix_MV(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z) {
    const ellpack_info_t* info = &(A->_ellpack);
    int icon, ierr=1;

    //int ierr = c_dm_vmvse(a, k, nw, n, icol, x, y, &icon);
    ierr = c_dm_vmvse(info->COEF, A->NROWS, info->nw, A->NROWS, info->ICOL, x, z, &icon);
    if (alpha != 1.0) {
        for(int i=0; i<A->NROWS; i++) {
            z[i] *= alpha;
        }
    }

    if ((y != NULL) && (beta != 0)) {
        for(int i=0; i<A->NROWS; i++) {
            z[i] += beta * y[i];
        }
    }

    return ierr ? -1 : 0;
}
#else
int Matrix_optimize(Matrix_t *A) {
    if (MATRIX_IS_CSR(A)) {
        Matrix_transpose(A);
    }

    Matrix_convert_index(A, 1);

    if (A->w == NULL) {
        A->w = (double*)calloc(A->NNZ, sizeof(double));
    } else {
        A->w = (double*)realloc((void*)A->w, (A->NNZ)*sizeof(double));
    }
    if (A->iw == NULL) {
        A->iw = (int*)calloc((A->NNZ)*2, sizeof(int));
    } else {
        A->iw = (int*)realloc((void*)A->iw, (A->NNZ)*2*sizeof(int));
    }

    if ((A->w == NULL)||(A->iw == NULL)) {
        printf("ERROR: Memory Allocation error.\n");
        return -1;
    }
    
    A->optimized = 1;
    return 0;
}

int Matrix_MV(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z) {
    int icon, ierr=1;

    //int ierr = c_dm_vmvscc(a, nz, nrow, nfcnz, n, x, y, w, iw, &icon);
    ierr = c_dm_vmvscc(A->values, A->NNZ, A->indice, A->pointers, A->NROWS, x, z, A->w, A->iw, &icon);
    if (alpha != 1.0) {
        for(int i=0; i<A->NROWS; i++) {
            z[i] *= alpha;
        }
    }

    if ((y != NULL) && (beta != 0)) {
        for(int i=0; i<A->NROWS; i++) {
            z[i] += beta * y[i];
        }
    }

    return ierr ? -1 : 0;
}
#endif /* ELLPACK */
#endif /* SSL2 */
