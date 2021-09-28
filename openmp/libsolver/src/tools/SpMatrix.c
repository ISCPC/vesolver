/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */
/*              Copyright (C) 2020 Shunji Uno */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include "SpMatrix.h"
#ifdef WITH_NETCDF
#include <netcdf.h>
#endif /* WITH_NETCDF */


#ifdef WITH_NETCDF
/*
 * NetCDF CSR matrix IO library
 */
int nc_save_csr_matrix(char *filename, SpMatrix_t *A) {
    int ncid;
    int info_dimid, neq_dimid, nzs_dimid;
    int varid_info, varid_rows, varid_cols, varid_values;
    INT_T info[4];
    int retval;

    if ((retval = nc_create(filename, NC_CLOBBER, &ncid))) {
        return -1;
    }

    /* 
     * Define data
     */
    if ((retval = nc_def_dim(ncid, "info", 4, &info_dimid))) {
        return -1;
    }
    if ((retval = nc_def_var(ncid, "info", NC_INT, 1, &info_dimid, &varid_info))) {
        return -1;
    }

    /* A_rows */
    if ((retval = nc_def_dim(ncid, "neq", A->nrow+1, &neq_dimid))) {
        return -1;
    }
    if ((retval = nc_def_var(ncid, "rows", NC_INT, 1, &neq_dimid, &varid_rows))) {
        return -1;
    }

    /* A_cols & A_values*/
    if ((retval = nc_def_dim(ncid, "nzs", A->ndim, &nzs_dimid))) {
        return -1;
    }
    if ((retval = nc_def_var(ncid, "cols", NC_INT, 1, &nzs_dimid, &varid_cols))) {
        return -1;
    }
    if ((retval = nc_def_var(ncid, "values", NC_DOUBLE, 1, &nzs_dimid, &varid_values))) {
        return -1;
    }

    if ((retval = nc_enddef(ncid))) {
        return -1;
    }

    /* 
     * Put data
     */
    info[0] = A->nrow;
    info[1] = A->ndim;
    info[2] = A->offset;
    info[3] = A->type; // SPMATRIX_TYPE_CSR | SPMATRIX_TYPE_INDEX1 | SPMATRIX_TYPE_ASYMMETRIC;
    if ((retval = nc_put_var_int(ncid, varid_info, info))) {
        return -1;
    }

    if ((retval = nc_put_var_int(ncid, varid_rows, A->pointers))) {
        return -1;
    }

    if ((retval = nc_put_var_int(ncid, varid_cols, A->indice))) {
        return -1;
    }

    if ((retval = nc_put_var_double(ncid, varid_values, A->value))) {
        return -1;
    }

    if ((retval = nc_close(ncid))) {
        return -1;
    }

    return 0;
}

SpMatrix_t* nc_load_csr_matrix(char *filename) {
    int ncid;
    int varid_info, varid_rows, varid_cols, varid_values;
    int retval;
    INT_T info[4];
    SpMatrix_t *A=malloc(sizeof(SpMatrix_t));


    if ((retval = nc_open(filename, NC_NOWRITE, &ncid))) {
        return NULL;
    }

    /* Get Info */
    if ((retval = nc_inq_varid(ncid, "info", &varid_info))) {
        return NULL;
    }

    if ((retval = nc_get_var_int(ncid, varid_info, info))) {
        return NULL;
    }

    // printf("info: %d %d %d\n", info[0], info[1], info[2]);
    A->ncol = A->nrow = info[0];
    A->ndim = info[1];
    A->offset = info[2];
    A->type = info[3];

    A->pointers = calloc(A->nrow+1, sizeof(INT_T));
    A->indice = calloc(A->ndim, sizeof(INT_T));
    A->value = calloc(A->ndim, sizeof(double));

    /* A_rows */
    if ((retval = nc_inq_varid(ncid, "rows", &varid_rows))) {
        return NULL;
    }

    if ((retval = nc_get_var_int(ncid, varid_rows, A->pointers))) {
        return NULL;
    }

    /* A_cols & A_values*/
    if ((retval = nc_inq_varid(ncid, "cols", &varid_cols))) {
        return NULL;
    }

    if ((retval = nc_get_var_int(ncid, varid_cols, A->indice))) {
        return NULL;
    }

    if ((retval = nc_inq_varid(ncid, "values", &varid_values))) {
        return NULL;
    }

    if ((retval = nc_get_var_double(ncid, varid_values, A->value))) {
        return NULL;
    }

    return A;
}

int nc_save_vector(char *filename, Vector_t *v) {
    int ncid;
    int dimid, varid;
    int retval;

    if ((retval = nc_create(filename, NC_CLOBBER, &ncid))) {
        return -1;
    }

    /* 
     * Define data
     */
    if ((retval = nc_def_dim(ncid, "size", v->size, &dimid))) {
        return -1;
    }
    if ((retval = nc_def_var(ncid, "values", NC_DOUBLE, 1, &dimid, &varid))) {
        return -1;
    }

    if ((retval = nc_enddef(ncid))) {
        return -1;
    }

    /* 
     * Put data
     */
    if ((retval = nc_put_var_double(ncid, varid, v->value))) {
        return -1;
    }

    if ((retval = nc_close(ncid))) {
        return -1;
    }

    return 0;
}

Vector_t* nc_load_vector(char *filename) {
    int ncid, varid, ndims, dimid;
    int retval;
    size_t dimlen;
    Vector_t *v=malloc(sizeof(Vector_t));

    if ((retval = nc_open(filename, NC_NOWRITE, &ncid))) {
        return NULL;
    }

    if ((retval = nc_inq_varid(ncid, "values", &varid))) {
        return NULL;
    }

    if ((retval = nc_inq_varndims(ncid, varid, &ndims))) {
        return NULL;
    }

    if (ndims > 1) {
        return NULL;
    }

    if ((retval = nc_inq_vardimid(ncid, varid, &dimid))) {
        return NULL;
    }

    if ((retval = nc_inq_dimlen(ncid, dimid, &dimlen))) {
        return NULL;
    }

    v->size = dimlen;
    v->value = calloc(dimlen, sizeof(double));

    if ((retval = nc_get_var_double(ncid, varid, v->value))) {
        return NULL;
    }

    return v;
}
#endif /* WITH_NETCDF */

/*
 * Binary CSR matrix IO library
 */
int SpMatrix_save_csr_matrix(char *filename, SpMatrix_t *A) {
    int fd;
    INT_T info[8];
    size_t datasz;

    if ((fd = open(filename, O_CREAT|O_RDWR|O_TRUNC, 0644)) < 0) {
        return -1;
    }

    /* 
     * Put data
     */
    info[0] = A->nrow;
    info[1] = A->ndim;
    info[2] = A->offset;
    info[3] = A->type;
    info[4] = A->neq;

    datasz = sizeof(INT_T)*8;
    if (write(fd, info, datasz) != datasz) {
        return -1;
    }
    
    datasz = sizeof(INT_T)*(A->nrow+1);
    if (write(fd, A->pointers, datasz) != datasz) {
        printf("ERROR: cannot save pointers of A from %s\n", filename);
        return -1;
    }

    datasz = sizeof(INT_T)*(A->ndim);
    if (write(fd, A->indice, datasz) != datasz) {
        printf("ERROR: cannot save indice of A from %s\n", filename);
        return -1;
    }

    datasz = sizeof(double)*(A->ndim);
    if (write(fd, A->value, datasz) != datasz) {
        printf("ERROR: cannot save values of A from %s\n", filename);
        return -1;
    }

    if (SPMATRIX_IS_DISTRIBUTE(A)) {
        datasz = sizeof(int)*(A->nrow);
        if (write(fd, A->order, datasz) != datasz) {
            printf("ERROR: cannot save order of A from %s\n", filename);
            return -1;
        }
    }

    close(fd);

    return 0;
}

SpMatrix_t* SpMatrix_load_csr_matrix(char *filename) {
    int fd;
    INT_T info[8];
    SpMatrix_t *A=malloc(sizeof(SpMatrix_t));
    size_t datasz;

    if ((fd = open(filename, O_RDONLY)) < 0) {
        return NULL;
    }

    datasz = sizeof(INT_T)*8;
    if (read(fd, info, datasz) != datasz) {
        return NULL;
    }
    
    A->ncol = A->nrow = info[0];
    A->ndim = info[1];
    A->offset = info[2];
    A->type = info[3];

    A->pointers = calloc(A->nrow+1, sizeof(INT_T));
    A->indice = calloc(A->ndim, sizeof(INT_T));
    A->value = calloc(A->ndim, sizeof(double));

    datasz = sizeof(INT_T)*(A->nrow+1);
    if (read(fd, A->pointers, datasz) != datasz) {
        printf("ERROR: cannot load pointers of A from %s\n", filename);
        return NULL;
    }

    datasz = sizeof(INT_T)*(A->ndim);
    if (read(fd, A->indice, datasz) != datasz) {
        printf("ERROR: cannot load indice of A from %s\n", filename);
        return NULL;
    }

    datasz = sizeof(double)*(A->ndim);
    if (read(fd, A->value, datasz) != datasz) {
        printf("ERROR: cannot load values of A from %s\n", filename);
        return NULL;
    }

    if (SPMATRIX_IS_DISTRIBUTE(A)) {
        A->neq = info[4];
        A->order = (int*)calloc(sizeof(int), A->nrow);
        datasz = sizeof(int)*(A->nrow);
        if (read(fd, A->order, datasz) != datasz) {
            printf("ERROR: cannot load order of A from %s\n", filename);
            return NULL;
        }
    }
    close(fd);

    return A;
}

int SpMatrix_save_vector(char *filename, Vector_t *v) {
    int fd;
    INT_T num = v->size;
    size_t datasz;

    if ((fd = open(filename, O_CREAT|O_RDWR|O_TRUNC, 0644)) < 0) {
        return -1;
    }

    datasz = sizeof(INT_T);
    if (write(fd, &num, datasz) != datasz) {
        return -1;
    }
    
    datasz = sizeof(double)*(num);
    if (write(fd, v->value, datasz) != datasz) {
        printf("ERROR: cannot save vector from %s\n", filename);
        return -1;
    }

    close(fd);

    return 0;
}

Vector_t* SpMatrix_load_vector(char *filename) {
    int fd;
    Vector_t *v=malloc(sizeof(Vector_t));
    INT_T num;
    size_t datasz;

    if ((fd = open(filename, O_RDONLY)) < 0) {
        return NULL;
    }

    datasz = sizeof(INT_T);
    if (read(fd, &num, datasz) != datasz) {
        return NULL;
    }
    
    v->size = num;
    v->value = calloc(num, sizeof(double));

    datasz = sizeof(double)*(num);
    if (read(fd, v->value, datasz) != datasz) {
        printf("ERROR: cannot load vector from %s\n", filename);
        return NULL;
    }

    close(fd);

    return v;
}

/*
 * Basic Functions
 */
void Vector_free(Vector_t* v) {
    free(v->value);
    free(v);
}

void SpMatrix_free(SpMatrix_t* A) {
    free(A->pointers);
    free(A->indice);
    free(A->value);
    free(A);
}

SpMatrix_t* SpMatrix_transpose(const SpMatrix_t* const A0) {
    SpMatrix_t* A1;

    A1 = (SpMatrix_t*)malloc(sizeof(SpMatrix_t));
    A1->type = A0->type;
    A1->ncol = A0->ncol;
    A1->nrow = A0->nrow;
    A1->ndim = A0->ndim;
    A1->pointers = (INT_T*)calloc(sizeof(INT_T), A1->ncol+1);
    A1->indice = (INT_T*)calloc(sizeof(INT_T), A1->ndim);
    A1->value = (double*)calloc(sizeof(double), A1->ndim);

    /* phase0: */
    for(int i=0; i<=A1->ncol; i++) {
        A1->pointers[i] = 0;
    }

    /* phase1 */
    for(int i=0; i<A0->ndim; i++) {
        int j = A0->indice[i]-1;
        A1->pointers[j]++;
    }
    A1->pointers[0]++;
    for(int i=1; i<=A1->ncol; i++) {
        A1->pointers[i] += A1->pointers[i-1];
    }

    /* phase2 */
    for (int i=A0->nrow-1; i>=0; i--) {
        for (int j=A0->pointers[i]-1; j<A0->pointers[i+1]-1; j++) {
            int jj = --(A1->pointers[A0->indice[j]-1]);
            A1->indice[jj-1] = i+1;
            A1->value[jj-1] = A0->value[j];
        }
    }

    return A1;
}

SpMatrix_t* SpMatrix_extract_symmetric(const SpMatrix_t* const A0) {
    SpMatrix_t* A1;

    A1 = (SpMatrix_t*)malloc(sizeof(SpMatrix_t));
    A1->type = A0->type;
    A1->ncol = A0->ncol;
    A1->nrow = A0->nrow;
    A1->ndim = A0->ndim*2-A0->ncol;
    A1->pointers = (INT_T*)calloc(sizeof(INT_T), A1->ncol+1);
    A1->indice = (INT_T*)calloc(sizeof(INT_T), A1->ndim);
    A1->value = (double*)calloc(sizeof(double), A1->ndim);

    /* phase0: */
    for(int i=0; i<A1->ncol; i++) {
        A1->pointers[i] = A0->pointers[i+1]-A0->pointers[i]-1;
        //printf("Phase0: A1->pointers[%d]=%d\n", i, A1->pointers[i]);
    }
    A1->pointers[A1->ncol] = 0;

    /* phase1 */
    for(int i=0; i<A0->ndim; i++) {
        int j = A0->indice[i]-1;
        A1->pointers[j]++;
    }
    A1->pointers[0]++;
    for(int i=1; i<=A1->ncol; i++) {
        A1->pointers[i] += A1->pointers[i-1];
    }
    //for(int i=0; i<=A1->ncol; i++) {
    //    printf("Phase1: A1->pointers[%d]=%d\n", i, A1->pointers[i]);
    //}

    /* phase2 */
    for (int i=A0->nrow-1; i>=0; i--) {
        for (int j=A0->pointers[i+1]-2; j>A0->pointers[i]-2; j--) {
            int jj = --(A1->pointers[i]);
            A1->indice[jj-1] = A0->indice[j];
            A1->value[jj-1] = A0->value[j];
            //printf("Phase2-1: A1(%d, %d) = %lf (jj=%d)\n", i+1, A0->indice[j], A0->value[j], jj);
        }
        for (int j=A0->pointers[i]; j<A0->pointers[i+1]-1; j++) {
            int jj = --(A1->pointers[A0->indice[j]-1]);
            A1->indice[jj-1] = i+1;
            A1->value[jj-1] = A0->value[j];
            //printf("Phase2-2: A1(%d, %d) = %lf (jj=%d)\n", A0->indice[j], A1->indice[jj-1], A0->value[j], jj);
        }
        //printf("--\n");
    }

    return A1;
}

#ifdef WITH_MPI
#include <mpi.h>

SpMatrix_t* SpMatrix_distribute(const SpMatrix_t* const A0) {
    int rank, size;
    SpMatrix_t* A;
    int nl=1, nt=1;

    if ( (!SPMATRIX_IS_CSR(A0)) || (SPMATRIX_INDEX_TYPE(A0) != 1)) {
        return NULL;
    }

    A = (SpMatrix_t*)malloc(sizeof(SpMatrix_t));

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

#if 0
    int n = A0->nrow / size;
    nl = n*rank+1;
    nt = (rank == size-1) ? A0->nrow : n*(rank+1);
#else
    int n = A0->ndim / size;
    for(int i=0; i<=A0->nrow; i++) {
        if (n*rank < A0->pointers[i]) {
            nl = i+1;
            break;
        }
    }
    for(int i=nl-1; i<=A0->nrow; i++) {
        if (n*(rank+1) < A0->pointers[i]) {
            nt = i;
            break;
        }
    }
#endif

    A->offset = nl;
    A->nrow = nt-nl+1;
    A->ncol = A0->ncol;
    A->ndim = A0->pointers[nt] - A0->pointers[nl-1];
    A->pointers = (INT_T*)calloc(sizeof(INT_T), nt-nl+2);
    A->indice = (INT_T*)calloc(sizeof(INT_T), A->ndim);
    A->value = (double*)calloc(sizeof(double), A->ndim);
    printf("INFO:%d: row: %d - %d, nzs: %d\n", rank, nl, nt, A->ndim);

    int k = 0;
    for(int i=nl-1; i<nt; i++) {
        A->pointers[i-nl+1] = k+1;
        //printf("INFO:%d: copy Row %d => %d\n", rank, i, i-nl+1);
        for(int j=A0->pointers[i]-1; j<A0->pointers[i+1]-1; j++) {
            //printf("INFO:%d: copy A0[%d, %d] to A1[%d, %d] : %lf\n", rank, i, j, i-nl+1, k, A0->value[j]);
            A->indice[k] = A0->indice[j];
            A->value[k] = A0->value[j];
            k++;
        }
    }
    A->pointers[nt-nl+1] = k+1;

    MPI_Barrier(MPI_COMM_WORLD);

    return A;
}
#endif

Vector_t* create_vector(size_t size, double default_value) {
    Vector_t* vb;

    vb = (Vector_t*)malloc(sizeof(Vector_t));
    vb->size = size;
    vb->value = (double*)calloc(sizeof(double), vb->size);
    for(int i=0; i<vb->size; i++) {
        vb->value[i] = default_value;
    }
    return vb;
}

