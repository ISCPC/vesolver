#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "SpMatrix.h"

#define WORK_BUFSZ ((A)->nrow)

void Matrix_sort(SpMatrix_t *A) {
    int max_ncol = 0;
    int* work_ind = calloc(sizeof(int), WORK_BUFSZ);
    int* work_val = calloc(sizeof(double), WORK_BUFSZ);

    for(int i=0; i<A->nrow; i++) {
        int ncol = A->pointers[i+1] - A->pointers[i];
        if (ncol == 1) continue;

        //printf("Sorting Row %d...", i);
        int min_col, max_col;
        min_col = max_col = A->indice[A->pointers[i]];
        for(int j=A->pointers[i]+1; j<A->pointers[i+1]; j++) {
            if (A->indice[j] > max_col) {
                max_col = A->indice[j];
            }
            if (A->indice[j] < min_col) {
                min_col = A->indice[j];
            }
        }

        /* claer buf */
        for (int k=0; k<(max_col-min_col+1); k++) {
            work_ind[k] = -1;
        }


        for(int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
            work_ind[A->indice[j]-min_col] = 1;
            work_val[A->indice[j]-min_col] = A->value[j];
        }

        int jj = A->pointers[i];
        for (int k=0; k<(max_col-min_col+1); k++) {
            if (work_ind[k] > 0) {
                A->indice[jj] = min_col+k;
                A->value[jj] = work_val[k];
                jj++;
            }
        }
        //printf("\n");
    }

    free(work_ind);
    free(work_val);
    //printf("Max Colmun = %d\n", max_ncol);
    printf("Sorting done\n");
    fflush(stdout);
}

int main(int argc, char** argv) {
    SpMatrix_t *A;
    char *infile;
    int* order;
    int base=0;
    
    if (argc > 1) {
        infile = argv[1];
    } else {
        printf("Usage: %s [infile]\n", argv[0]);
        exit(1);
    }

    A = SpMatrix_load_csr_matrix(infile);

    if (A == NULL) {
        return(1);
    }

    if(SPMATRIX_IS_SYMMETRIC(A)) {
        printf("Type: Symmetric, ");
    } else {
        printf("Type: Asymmetric, ");
    }

    if(SPMATRIX_INDEX_TYPE(A)) {
        printf("ONE-based index, ");
        base = 1;
    } else {
        printf("ZERO-based index, ");
        base = 0;
    }

    if(SPMATRIX_IS_CSR(A)) {
        printf("CSR format\n");
    } else {
        printf("CSC format\n");
    }

#if 0
    if(SPMATRIX_IS_DISTRIBUTE(A)) {
        order = (int*)calloc(sizeof(int), A->nrow);
        for (int i=0; i<A->neq; i++) {
            order[A->order[i]-1] = i;
        }
    }
#endif

    int npara = 48;
    int d = (A->nrow-1) / npara + 1;
    printf("Matrix size: %d x %d (=> %d para, rows per proc %d)\n",
        A->ncol, A->nrow, npara, d);
    printf("The number of non-zero elements: %10d\n", A->ndim);

    Matrix_sort(A);

    if (base == 0) {
        for(int n=0; n<npara; n++) {
            int nnz=0;
            int maxncols = 0;
            int maxdist = 0;
            int maxcol = 0;
            int mincol = A->nrow;
            for(int i=d*n; (i<d*(n+1))&&(i<A->nrow); i++) {
                int ncols = A->pointers[i+1]-A->pointers[i];
                int min = A->indice[A->pointers[i]];
                int max = A->indice[A->pointers[i+1]-1];
                nnz+=ncols;

                //printf("INFO: Row#%d  ncols=%d min=%d max=%d (%d)\n", 
                //    i, ncols, mincol, maxcol, maxcol-mincol);
                if (ncols > maxncols) {
                    maxncols = ncols;
                }
                if ((max-min) > maxdist) {
                    maxdist = max-min;
                }
                if (max>maxcol) {
                    maxcol = max;
                }
                if (min<mincol) {
                    mincol = min;
                }
#if 0
                for(int j=A->pointers[i]; j<A->pointers[i+1]; j++) {
                    printf("%10d (-> %10d): A[%10d, %10d] = %+20.16le\n",
                         j, A->pointers[i+1], i, A->indice[j], A->value[j]);
                }
#endif
            }
            printf("INFO: n=%d, nnz=%d, maxncols=%d, maxdist=%d (%d-%d)\n", n, nnz, maxncols, maxdist, mincol, maxcol);
        }
    } else {
        for(int i=0; i<A->nrow; i++) {
#if 0
            for(int j=A->pointers[i]-1; j<A->pointers[i+1]-1; j++) {
                printf("%10d (-> %10d): A[%10d, %10d] = %+20.16le\n",
                     j+1, A->pointers[i+1], i+1, A->indice[j], A->value[j]);
            }
#endif
        }
    }

    return 0;
}
