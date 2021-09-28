#ifndef _SPMATRIX_H_
#define _SPMATRIX_H_
#include <stdint.h>

typedef int32_t INT_T;
typedef uint32_t UINT_T;

#define SPMATRIX_TYPE_CSC        (0)
#define SPMATRIX_TYPE_CSR        (1)
#define SPMATRIX_TYPE_INDEX0     (0<<4)
#define SPMATRIX_TYPE_INDEX1     (1<<4)
#define SPMATRIX_TYPE_ASYMMETRIC (0<<8)
#define SPMATRIX_TYPE_SYMMETRIC  (1<<8)
#define SPMATRIX_TYPE_DISTRIBUTE (1<<12)

#define SPMATRIX_IS_CSC(A)       ((((A)->type)&0xf) == SPMATRIX_TYPE_CSC)
#define SPMATRIX_IS_CSR(A)       ((((A)->type)&0xf) == SPMATRIX_TYPE_CSR)
#define SPMATRIX_INDEX_TYPE(A)   ((((A)->type)>>4)&0xf)
#define SPMATRIX_IS_SYMMETRIC(A) (((A)->type)&SPMATRIX_TYPE_SYMMETRIC)
#define SPMATRIX_IS_DISTRIBUTE(A) (((A)->type)&SPMATRIX_TYPE_DISTRIBUTE)

typedef struct {
    UINT_T    type;
    INT_T    *pointers;
    INT_T    *indice;
    double   *value;
    INT_T     ndim;
    INT_T     ncol;
    INT_T     nrow;
    INT_T     offset;
    INT_T     neq;
    INT_T    *order;
    INT_T    *rorder;
} SpMatrix_t;

typedef struct {
    double   *value;
    INT_T     size;
} Vector_t;

void SpMatrix_free(SpMatrix_t* A);
void Matrix_free(SpMatrix_t* b);

#ifdef WITH_NETCDF
/*
 * NetCDF CSR matrix IO library
 */
int nc_save_csr_matrix(char *filename, SpMatrix_t *A);
SpMatrix_t* nc_load_csr_matrix(char *filename);

int nc_save_vector(char *filename, Vector_t *v);
Vector_t* nc_load_vector(char *filename);
#endif /* WITH_NETCDF */

/*
 * Binary CSR matrix IO library
 */
int SpMatrix_save_csr_matrix(char *filename, SpMatrix_t *A);
SpMatrix_t* SpMatrix_load_csr_matrix(char *filename);

int SpMatrix_save_vector(char *filename, Vector_t *v);
Vector_t* SpMatrix_load_vector(char *filename);

/*
 * Basic Functions
 */
void Vector_free(Vector_t* v);
void SpMatrix_free(SpMatrix_t* A);
SpMatrix_t* SpMatrix_transpose(const SpMatrix_t* const A0);
SpMatrix_t* SpMatrix_extract_symmetric(const SpMatrix_t* const A0);
#ifdef WITH_MPI
SpMatrix_t* SpMatrix_distribute(const SpMatrix_t* const A0);
#endif

Vector_t* create_vector(size_t size, double default_value);
#endif /* _SPMATRIX_H_ */
