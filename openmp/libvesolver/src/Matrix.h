#pragma once
#define MATRIX_TYPE_INDEX0     (0)
#define MATRIX_TYPE_INDEX1     (1)
#define MATRIX_TYPE_ASYMMETRIC (0<<4)
#define MATRIX_TYPE_SYMMETRIC  (1<<4)
#define MATRIX_TYPE_LOWER      (1<<5)
#define MATRIX_TYPE_UPPER      (0)
#define MATRIX_TYPE_UNIT       (1<<6)
#define MATRIX_TYPE_NON_UNIT   (0)

#define MATRIX_INDEX_TYPE(A)    (((A)->flags)&0xf)
#define MATRIX_IS_SYMMETRIC(A)  (((A)->flags)&MATRIX_TYPE_SYMMETRIC)
#define MATRIX_IS_LOWER(A)      (((A)->flags)&MATRIX_TYPE_LOWER)
#define MATRIX_IS_UNIT(A)       (((A)->flags)&MATRIX_TYPE_UNIT)
