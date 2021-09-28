#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "SpMatrix.h"

int convert_matrix(char* infile, char* outfile_a, char* outfile_b, char* outfile_x) {
    FILE *fp;
    SpMatrix_t *A;
    Vector_t *b, *x;
    char dmy[16];
    int idx, idx0, idx1, nzs, nrows;
    int ptr, ind;
    int dmy0, dmy1;
    double val, val2;
    int i;
    int neq;

    if ((fp = fopen(infile, "r")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open file %s\n", infile);
        exit(1);
    }

    /*
     * Reading Matrix A
     */
    if (fscanf(fp, "%s%d%d%d", dmy, &idx0, &neq, &nzs) < 4) {
        fprintf(stderr, "ERROR: Invalid header for A.\n");
        exit(1);
    }

    A = (SpMatrix_t*)malloc(sizeof(SpMatrix_t));
    A->type = SPMATRIX_TYPE_CSR | SPMATRIX_TYPE_INDEX0 | SPMATRIX_TYPE_ASYMMETRIC;
    A->ncol = A->nrow = neq;
    A->ndim = nzs;
    A->offset = 0;
    A->neq = neq;
    printf("0x%08x %d %d %d %d\n", A->type, A->nrow, A->ndim, A->offset, neq);

    A->pointers = (int*)calloc(A->nrow+1, sizeof(int));
    A->indice = (int*)calloc(nzs, sizeof(int));
    A->value = (double*)calloc(nzs, sizeof(double));

    if ((fscanf(fp, "%s%d%d", dmy, &dmy0, &ptr) < 3) || (strncmp(dmy, "A_rows", 16))) {
        fprintf(stderr, "ERROR: Data corruption detected at reading A_rows\n");
        exit(1);
    }
    A->pointers[0] = ptr;
    for(i=1; i<=A->nrow; i++) {
        if (fscanf(fp, "%s%d%d", dmy, &dmy0, &ptr) < 3) {
            fprintf(stderr, "ERROR: Invalid data format for A_rows.\n");
            exit(1);
        }
        A->pointers[i] = ptr;
    }

    if ((fscanf(fp, "%s%d%d%lf", dmy, &dmy0, &ind, &val) < 4) || (strncmp(dmy, "A_cols", 16))) {
        fprintf(stderr, "ERROR: Data corruption detected at reading A_cols\n");
        exit(1);
    }
    A->indice[0] = ind;
    A->value[0] = val;
    for(i=1; i<A->ndim; i++) {
        if (fscanf(fp, "%s%d%d%lf", dmy, &dmy0, &ind, &val) < 4) {
            fprintf(stderr, "ERROR: Invalid data format for A_cols.\n");
            exit(1);
        }
        A->indice[i] = ind;
        A->value[i] = val;
    }

    b = (Vector_t*)malloc(sizeof(Vector_t));
    b->value = (double*)calloc(A->nrow, sizeof(double));
    b->size = A->nrow;
    if ((fscanf(fp, "%s%d%lf", dmy, &dmy0, &val) < 3) || (strncmp(dmy, "b", 16))) {
        fprintf(stderr, "ERROR: Data corruption detected at reading b.\n");
        exit(1);
    }
    b->value[0] = val;
    for(i=1; i<A->nrow; i++) {
        if (fscanf(fp, "%s%d%lf", dmy, &dmy0, &val) < 3) {
            fprintf(stderr, "ERROR: Invalid data format for b.\n");
            exit(1);
        }
        b->value[i] = val;
    }

    x = (Vector_t*)malloc(sizeof(Vector_t));
    x->value = (double*)calloc(A->nrow, sizeof(double));
    x->size = A->nrow;
    if ((fscanf(fp, "%s%d%lf", dmy, &dmy0, &val) < 3) || (strncmp(dmy, "x", 16))) {
        fprintf(stderr, "ERROR: Data corruption detected at reading x.\n");
        exit(1);
    }
    x->value[0] = val;
    for(i=1; i<A->nrow; i++) {
        if (fscanf(fp, "%s%d%lf", dmy, &dmy0, &val) < 3) {
            fprintf(stderr, "ERROR: Invalid data format for x.\n");
            exit(1);
        }
        x->value[i] = val;
    }

    SpMatrix_save_csr_matrix(outfile_a, A);
    SpMatrix_save_vector(outfile_b, b);
    SpMatrix_save_vector(outfile_x, x);

    fclose(fp);
    return 0;
}

int main(int argc, char** argv) {
    char infile[256], outfile_a[256], outfile_b[256];
    char* prefix = "mat.txt";
    int size=4;

    convert_matrix("mat.txt", "a.bin", "b.bin", "x.bin");

    return 0;
}
