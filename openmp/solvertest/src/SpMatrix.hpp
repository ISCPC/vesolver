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

#pragma once
#include <stdint.h>

typedef int32_t INT_T;
typedef uint32_t UINT_T;

#define SPMATRIX_TYPE_CSC        (0)
#define SPMATRIX_TYPE_CSR        (1)
#define SPMATRIX_TYPE_DCSC       (2)
#define SPMATRIX_TYPE_DCSR       (3)
#define SPMATRIX_TYPE_INDEX0     (0<<4)
#define SPMATRIX_TYPE_INDEX1     (1<<4)
#define SPMATRIX_TYPE_ASYMMETRIC (0<<8)
#define SPMATRIX_TYPE_SYMMETRIC  (1<<8)
#define SPMATRIX_TYPE_DISTRIBUTE (1<<12)
#define SPMATRIX_TYPE_REORDERED  (1<<16)

#if 0
#define SPMATRIX_IS_CSC(A)        ((((A)->type)&0xf) == SPMATRIX_TYPE_CSC)
#define SPMATRIX_IS_CSR(A)        ((((A)->type)&0xf) == SPMATRIX_TYPE_CSR)
#define SPMATRIX_IS_DCSC(A)       ((((A)->type)&0xf) == SPMATRIX_TYPE_CSC)
#define SPMATRIX_IS_DCSR(A)       ((((A)->type)&0xf) == SPMATRIX_TYPE_CSR)
#define SPMATRIX_INDEX_TYPE(A)    ((((A)->type)>>4)&0xf)
#define SPMATRIX_IS_SYMMETRIC(A)  (((A)->type)&SPMATRIX_TYPE_SYMMETRIC)
#define SPMATRIX_IS_DISTRIBUTE(A) (((A)->type)&SPMATRIX_TYPE_DISTRIBUTE)
#else
#define SPMATRIX_IS_CSC()        (((type)&0xf) == SPMATRIX_TYPE_CSC)
#define SPMATRIX_IS_CSR()        (((type)&0xf) == SPMATRIX_TYPE_CSR)
#define SPMATRIX_IS_DCSC()       (((type)&0xf) == SPMATRIX_TYPE_DCSC)
#define SPMATRIX_IS_DCSR()       (((type)&0xf) == SPMATRIX_TYPE_DCSR)
#define SPMATRIX_INDEX_TYPE()    (((type)>>4)&0xf)
#define SPMATRIX_IS_SYMMETRIC()  (((type)&SPMATRIX_TYPE_SYMMETRIC)
#define SPMATRIX_IS_DISTRIBUTE() ((type)&SPMATRIX_TYPE_DISTRIBUTE)
#define SPMATRIX_IS_REORDERED()  ((type)&SPMATRIX_TYPE_REORDERED)
#endif

class SpMatrix {
public:
    SpMatrix();
    ~SpMatrix();

    int load_file(const char *filename, bool oldflag=false);
    int alloc(int64_t neq, int64_t nnz);
    int transpose();

    UINT_T    type;
    INT_T    *pointers;
    INT_T    *indice;
    double   *value;
    INT_T     ndim;
    INT_T     ncol;
    INT_T     nrow;
    INT_T     offset;
    INT_T     neq;
    INT_T     mingind;
    INT_T     maxgind;
    INT_T    *order;
    INT_T    *rorder;

protected:
    bool    selfAllocated;
};

class SpDistMatrix : public SpMatrix {
public:
    SpDistMatrix();
    ~SpDistMatrix();

    int reordering(); 
    static SpMatrix* gather(SpDistMatrix** An, int nprocs);
    int ConvertToDCSR();

private:
};

class Vector {
public:
    Vector();
    ~Vector();

    int load_file(const char *filename);
    int alloc(int64_t n);

    double   *value;
    int64_t   size;

private:
    bool    selfAllocated;
};

class DistVector : public Vector {
public:
    DistVector();
    ~DistVector();

    static Vector* gather(DistVector** bn, int nprocs, SpDistMatrix** An);

private:
};

#if 0
void SpMatrix_free(SpMatrix_t* A);
void Matrix_free(SpMatrix_t* b);

/*
 * Binary CSR matrix IO library
 */
SpMatrix_t* SpMatrix_load_csr_matrix(char *filename);
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
#endif
