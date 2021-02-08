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
#include <cstdio>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <strings.h>
#include <fcntl.h>
#include "SpMatrix.hpp"
#include "timelog.h"

/*
 * Functions for SpMatrix
 */
SpMatrix::SpMatrix() {
    pointers = NULL;
    indice = NULL;
    value = NULL;
    order = NULL;
    rorder = NULL;
    ncol = nrow = 0;
    ndim = 0;
    neq = 0;
    selfAllocated = false;
}

SpMatrix::~SpMatrix() {
    if (selfAllocated) {
        if (pointers) { free(pointers); }
        if (value) { free(value); }
        if (order) { free(order); }
    }
    if (selfAllocated||SPMATRIX_IS_REORDERED()) {
        if (indice) { free(indice); }
    }
    if (rorder) { free(rorder); }
    pointers = NULL;
    indice = NULL;
    value = NULL;
    order = NULL;
    rorder = NULL;
    ncol = nrow = 0;
    ndim = 0;
    neq = 0;
}

int SpMatrix::alloc(int64_t nrows, int64_t nnz) {
    selfAllocated = true;
    pointers = (INT_T*)calloc(sizeof(INT_T), nrows+1);
    indice = (INT_T*)calloc(sizeof(INT_T), nnz);
    value = (double*)calloc(sizeof(double), nnz);
    return (value) ? 0 : -1;
}

int SpMatrix::load_file(const char *filename) {
    int fd;
    int32_t info[8];
    ssize_t datasz;

    if ((fd = open(filename, O_RDONLY)) < 0) {
        return -1;
    }

    datasz = sizeof(INT_T)*8;
    if (read(fd, info, datasz) != datasz) {
        return -1;
    }
    
    ncol = nrow = info[0];
    ndim = info[1];
    offset = info[2];
    type = info[3];

    alloc(nrow, ndim);

    datasz = sizeof(INT_T)*(nrow+1);
    if (read(fd, pointers, datasz) != datasz) {
        printf("ERROR: cannot load pointers of A from %s\n", filename);
        return -1;
    }

    datasz = sizeof(INT_T)*(ndim);
    if (read(fd, indice, datasz) != datasz) {
        printf("ERROR: cannot load indice of A from %s\n", filename);
        return -1;
    }

    datasz = sizeof(double)*(ndim);
    if (read(fd, value, datasz) != datasz) {
        printf("ERROR: cannot load values of A from %s\n", filename);
        return -1;
    }

    if (SPMATRIX_IS_DISTRIBUTE()) {
        neq = info[4];
        order = (INT_T*)calloc(sizeof(INT_T), nrow);
        datasz = sizeof(int)*(nrow);
        if (read(fd, order, datasz) != datasz) {
            printf("ERROR: cannot load order of A from %s\n", filename);
            return -1;
        }
    } else {
        neq = nrow;
    }

    close(fd);

    return 0;
}

/*
 * Functions for SpDistMatrix
 */
SpDistMatrix::SpDistMatrix() {}
SpDistMatrix::~SpDistMatrix() {}

int SpDistMatrix::reordering() {
    INT_T* newind = (INT_T*)calloc(sizeof(INT_T), ndim);

    /* reordering indice */
    if (!SPMATRIX_IS_REORDERED()) {
        for(int i=0; i<ndim; i++) {
            if ((indice[i] > nrow)||(indice[i]<=0)) {
                printf("ERROR: invalid data in indice[%d] = %d.\n", i, indice[i]);
            } else {
                newind[i] = order[indice[i]-1];
            }
        }
        ncol = neq;
        if (selfAllocated) {
            free(indice);
        }
        indice = newind;
        type |= SPMATRIX_TYPE_REORDERED;
    }

    /* create reverse order table */
    rorder = (INT_T*)calloc(neq, sizeof(INT_T));
    for(int i=0; i<nrow; i++) {
        if ((order[i]-1 >= neq)||(order[i]<=0)) {
            printf("ERROR: invalid data in order[%d] = %d.\n", i, order[i]);
        } else {
            rorder[order[i]-1] = i+1;
        }
    }
    /* find max/min index */
    maxgind = order[0];
    mingind = order[0];
    for(int i=1; i<nrow; i++) {
        if (order[i] > maxgind) maxgind = order[i];
        if (order[i] < mingind) mingind = order[i];
    }

    return 0;
}

/* item list operation */
struct list_item {
    struct list_item* next;
    int   indice;
    double value;
};

inline struct list_item* item_new(int indice, double value) {
    struct list_item *item;

    item = (struct list_item*)malloc(sizeof(struct list_item));
    if (item) {
        item->indice = indice;
        item->value = value;
        item->next = NULL;
    } else {
        printf("ERROR:item_new: alloc failed\n");
        exit(1);
    }

    return item;
}

struct sorted_list {
    struct list_item *top;
};

inline struct sorted_list* list_new() {
    struct sorted_list *list;

    list = (struct sorted_list*)malloc(sizeof(struct sorted_list));
    if (list) {
        list->top = NULL;
    } else {
        printf("ERROR:list_new: alloc failed\n");
        exit(1);
    }

    return list;
}

inline void list_free(struct sorted_list *list) {
#if 0
    struct list_item* item;

    for (item = list->top; item != NULL; item = item->next) {
        free(item);
    }
#endif
    free(list);
}

inline struct list_item* list_top(struct sorted_list* list) {
    return list->top;
}

inline void list_insert(struct sorted_list *list, int indice, double value) {
    struct list_item* newitem;
    struct list_item* item = list->top;

    if (item == NULL) {
        newitem = item_new(indice, value);
        list->top = newitem;
        return;
    }

    if (indice < item->indice) {
        newitem = item_new(indice, value);
        newitem->next = list->top;
        list->top = newitem;
        return;
    } else if (indice == item->indice) {
        item->value += value;
        return;
    }

    for (; item->next != NULL; item = item->next) {
        if (indice < item->next->indice) {
            newitem = item_new(indice, value);
            newitem->next = item->next;
            item->next = newitem;
            return;
        } else if (indice == item->next->indice) {
            item->next->value += value;
            return;
        }
    }

    newitem = item_new(indice, value);
    item->next = newitem;
}

SpMatrix* SpDistMatrix::gather(SpDistMatrix** An, int nprocs) {
    TIMELOG(tl);

    SpMatrix* A = new SpMatrix();

    /*
     * preprocess
    */
    int ndim=0;
    int maxnelem=0;
    int neq = An[0]->neq;

    TIMELOG_START(tl);
    int* duplex = (int*)calloc(sizeof(int), neq);
    int* rank = (int*)calloc(sizeof(int), neq);
    for(int i=0; i<neq; i++) {
        int lastcol=0;
        int nelem = 0;
        for(int r=0; r<nprocs; r++) {
            int ii = An[r]->rorder[i];
            if (ii != 0) {
                nelem += An[r]->pointers[ii] - An[r]->pointers[ii-1];
                duplex[i]++;
                rank[i] = r;
                if (lastcol >= An[r]->indice[An[r]->pointers[ii-1]]) {
                    duplex[i] = -1000000;
                }
                lastcol = An[r]->indice[An[r]->pointers[ii]-1];
            }
        }
        if (nelem > maxnelem) {
            maxnelem = nelem;
        }
    }
    for(int r=0; r<nprocs; r++) {
        ndim += An[r]->ndim;
    }
    TIMELOG_END(tl, "gather_preprocess");
    printf("INFO: maxnelem = %d\n", maxnelem);

    /*
     * gathering Matrix A
    */
    TIMELOG_START(tl);
    A->alloc(neq, ndim);
    int k=0;
    for(int i=0; i<neq; i++) {
        A->pointers[i] = k+1;
        //printf("INFO: Overwrapped row (row=%d)\n", i);
        struct sorted_list *tmplist = list_new();
        for(int r=0; r<nprocs; r++) {
            int ii = An[r]->rorder[i];
            if (ii != 0) {
                for(int j=An[r]->pointers[ii-1]-1; j<An[r]->pointers[ii]-1; j++) {
                    list_insert(tmplist, An[r]->indice[j], An[r]->value[j]);
                    //printf("A%d(%d, %d) = %lf\n", r, i, An[r]->indice[j], An[r]->value[j]);
                }
            }
        }

        //for(struct list_item* item = list_top(tmplist); item != NULL; item = item->next) {
        for(struct list_item* item = list_top(tmplist); item != NULL;) {
            struct list_item* next;
            //printf("INFO: Insert A(%d, %d) = %lf\n", i, item->indice, item->value);
            A->indice[k] = item->indice;
            A->value[k] = item->value;
            next = item->next;
            free(item);
            item = next;
            k++;
        }
        list_free(tmplist);
    }
    A->pointers[neq] = k+1;
    printf("Gather: ndim=%d, neq=%d\n", k, neq);
    TIMELOG_END(tl, "gather_main");

    A->type = SPMATRIX_TYPE_CSR | SPMATRIX_TYPE_INDEX1 | SPMATRIX_TYPE_ASYMMETRIC;
    A->ncol = A->nrow = neq;
    A->ndim = k;
    A->offset = 1;
    A->neq = neq;

    free(duplex);
    free(rank);

    return A;
}

int SpDistMatrix::ConvertToDCSR() {
    TIMELOG(tl);

    TIMELOG_START(tl);
    /* sorting indice */
    nrow = maxgind-mingind+1;
    INT_T* newptr = (INT_T*)calloc(sizeof(INT_T), nrow+1);
    INT_T* newind = (INT_T*)calloc(sizeof(INT_T), ndim);
    double* newval = (double*)calloc(sizeof(double), ndim);
    int jj=0;
    for (int i=0; i<nrow; i++) {
        newptr[i] = jj+1;
        struct sorted_list *tmplist = list_new();
        if (rorder[mingind+i-1] != 0) {
            int k = rorder[mingind+i-1]-1;
            for (int j=pointers[k]-1; j<pointers[k+1]-1; j++) {
                list_insert(tmplist, indice[j], value[j]);
            }
            for(struct list_item* item = list_top(tmplist); item != NULL;) {
                struct list_item* next;
                //printf("INFO: Insert A(%d, %d) = %lf\n", i, item->indice, item->value);
                newind[jj] = item->indice;
                newval[jj] = item->value;
                next = item->next;
                free(item);
                item = next;
                jj++;
            }
            list_free(tmplist);
        }
    }
    newptr[nrow] = jj+1;
    if (selfAllocated) {
        free(pointers);
        free(indice);
        free(value);
    }
    pointers = newptr;
    indice = newind;
    value = newval;

    type = SPMATRIX_TYPE_DCSR | SPMATRIX_TYPE_INDEX1 | SPMATRIX_TYPE_ASYMMETRIC;
    TIMELOG_END(tl, "ConverToDCSR");

    return 0;
}

/*
 * Functions for Vector
 */
Vector::Vector() {
    size = 0;
    value = NULL;
    selfAllocated = false;
}

Vector::~Vector() {
    if (selfAllocated) {
        if (value) { free(value); }
    }
    size = 0;
    value = NULL;
}

int Vector::load_file(const char *filename) {
    int fd;
    int32_t num;
    ssize_t datasz;

    if ((fd = open(filename, O_RDONLY)) < 0) {
        return -1;
    }

    datasz = sizeof(int32_t);
    if (read(fd, &num, datasz) != datasz) {
        return -1;
    }
    
    selfAllocated = true;
    size = (int64_t)num;
    value = (double*)calloc(num, sizeof(double));

    datasz = sizeof(double) * num;
    if (read(fd, value, datasz) != datasz) {
        printf("ERROR: cannot load vector from %s\n", filename);
        return -1;
    }
    close(fd);

    return 0;
}

int Vector::alloc(int64_t n) {
    selfAllocated = true;
    size = n;
    value = (double*)calloc(n, sizeof(double));
    return (value) ? 0 : -1;
}

/*
 * Functions for DistVector
 */
DistVector::DistVector() {}
DistVector::~DistVector() {}

Vector* DistVector::gather(DistVector** bn, int nprocs, SpDistMatrix** An) {
    Vector* b = new Vector();
    int neq = An[0]->neq;

    b->alloc((int64_t)neq);

    /*
     * gathering vector b
     */
    for(int i=0; i<neq; i++) {
        for(int r=0; r<nprocs; r++) {
            int ii = An[r]->rorder[i];
            if (ii != 0) {
                b->value[i] += bn[r]->value[ii-1];
            }
        }
    }

    return b;
}
