#ifdef SSL2
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdint.h>
#include "Matrix.h"

/****************************************************************
 * Library with Fujitsu SSL2
 ****************************************************************/
void Matrix_init(Matrix_t *A) {
    A->pointers = NULL;
    A->indice = NULL;
    A->values = NULL;
    A->info = NULL;
    A->optimized = 0;
}

void Matrix_free(Matrix_t *A) {
    /* Destruction of the handle */
    if (A->optimized == 1) {
        A->optimized = 0;
    }

    Matrix_free_generic(A);
}

int Matrix_optimize(Matrix_t *A) {
    A->optimized = 1;
    return 0;
}

int Matrix_MV(const Matrix_t *A, const double alpha, const double* x, const double beta, const double* y, double* z) {
    return 0;
}
#endif /* SSL2 */
