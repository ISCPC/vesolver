/*
 * vesolver.c
 *
 *  Created on: Jun 11, 2021
 *      Author: uno
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "vesolver.h"

#ifdef SXAT
#include <ve_offload.h>
typedef struct veo_proc_handle veo_proc_handle;
typedef struct veo_thr_ctxt veo_thr_ctxt;
typedef struct veo_args veo_args;

#else
// Dummy VEO prototype definition to avoid build error
typedef void veo_proc_handle;
typedef void veo_thr_ctxt;
typedef void veo_args;

int veo_api_version(); 
veo_proc_handle* veo_proc_create(int venode) {return NULL;}
int veo_proc_destroy(veo_proc_handle* proc) {return -1;}
uint64_t veo_load_library(veo_proc_handle* proc, const char* libname) {return 0;}
veo_thr_ctxt* veo_context_open(veo_proc_handle* proc) {return NULL;}
int veo_context_close(veo_thr_ctxt* ctx) {return -1;}
int veo_alloc_mem(veo_proc_handle* h, uint64_t* addr, const size_t size) {return -1;}
int veo_free_mem(veo_proc_handle* h, uint64_t addr) {return -1;}
int veo_read_mem(veo_proc_handle* h, void* dst, uint64_t src, size_t size) {return -1;}
int veo_write_mem(veo_proc_handle* h, uint64_t dst, const void* src, size_t size) {return -1;}

int veo_call_sync(veo_proc_handle* h, uint64_t addr, veo_args* ca, uint64_t* result) {return -1;}
int veo_call_wait_result(veo_thr_ctxt* ctx, uint64_t reqid, uint64_t* retp) {return -1;}
uint64_t veo_call_async(veo_thr_ctxt* ctx, uint64_t addr, veo_args* args) {return 0;}

uint64_t veo_get_sym(veo_proc_handle* proc, uint64_t libhdl, const char* symname) {return 0;}
veo_args* veo_args_alloc(void) {return NULL;}
void veo_args_free(veo_args* ca) {return;}
int veo_args_set_u64(veo_args* ca, int argnum, int64_t val) {return -1;}
int veo_args_set_i32(veo_args* ca, int argnum, int32_t val) {return -1;}

#define VEO_REQUEST_ID_INVALID (~0UL)

enum veo_command_state {
    VEO_COMMAND_OK = 0,
    VEO_COMMAND_EXCEPTION,
    VEO_COMMAND_ERROR,
    VEO_COMMAND_UNFINISHED,
};
#endif /* SXAT */

/*
 * Main content
 */
static struct veo_func_table {
    char* symname;
    uint64_t vemva;
} veo_func_table[] = {
        /*
         * uint64_t vesolver_init();
         * uint64_t vesolver_finalize();
         * uint64_t vesolver_create_handle_by_desc(uint64_t vedescp);
         * uint64_t vesolver_free_handle(uint64_t vedescp);
         * uint64_t vesolver_solve(uint64_t vedescp);
         */
        {"vesolver_init", 0},
        {"vesolver_finalize", 0},
        {"vesolver_create_handle_by_desc", 0},
        {"vesolver_free_handle", 0},
        {"vesolver_solve", 0}
};

#define CCX_VEO_DEFAULT_LIBRARY_PATH "/opt/local/ve/lib/libLinearSolver.so"

#define VESOLVER_INIT           0
#define VESOLVER_FINALIZE       1
#define VESOLVER_CREATE_HANDLE  2
#define VESOLVER_FREE_HANDLE    3
#define VESOLVER_SOLVE          4
#define VESOLVER_NUM_FUNCS      5

typedef struct vesolver_option {
    int64_t solver;
    double res;
} vesolver_option_t;

typedef struct vesolver_desc {
    int32_t ves_handle;
    INT_T neq;
    INT_T nnz;
    uint32_t flags;
    uint64_t dptr;
    uint64_t pointers;
    uint64_t indice;
    uint64_t values;
    uint64_t bptr;
    uint64_t xptr;
    vesolver_option_t options;
} vesolver_desc_t;

typedef struct vesolver_instance {
    veo_proc_handle* proc;
    uint64_t veo_handle;
} vesolver_instance_t;

typedef struct vesolver_context {
    vesolver_instance_t* instance;
    veo_thr_ctxt* ctx;
    matrix_desc_t*   matdesc;
    vesolver_desc_t* vedesc;
    vesolver_option_t option;
} vesolver_context_t;

#define VESOLVER_NUM_INSTANCES 1
static vesolver_instance_t instances[VESOLVER_NUM_INSTANCES];
#define VESOLVER_NUM_CONTEXTS 4
static vesolver_context_t  contexts[VESOLVER_NUM_CONTEXTS];

/*
 * VEO-call wrapper functions
 */
#define VEO_ALLOC_MEM(ctx, addr, size) veo_alloc_mem_wrapper((ctx), (addr), (size), __LINE__)
#define VEO_FREE_MEM(ctx, addr) veo_free_mem_wrapper((ctx), (addr), __LINE__)
#define VEO_READ_MEM(ctx, dst, src, size) veo_read_mem_wrapper((ctx), (dst), (src), (size), __LINE__)
#define VEO_WRITE_MEM(ctx, dst, src, size) veo_write_mem_wrapper((ctx), (dst), (src), (size), __LINE__)
#define VEO_CALL_ASYNC(ctx, symid, argp) veo_call_async_wrapper((ctx), (symid), (argp), __LINE__)
#define VEO_CALL_WAIT(ctx, reqid) veo_call_wait_wrapper((ctx), (reqid), __LINE__)
#define VEO_CALL_SYNC(ctx, symid, argp) veo_call_sync_wrapper((ctx), (symid), (argp), __LINE__)

static inline void veo_alloc_mem_wrapper(const vesolver_context_t* context, uint64_t* addr, const size_t size, int line) {
    int rc = veo_alloc_mem(context->instance->proc, addr, size);
    if (rc != 0) {
        printf("ERROR: Cannot allocate memory on VE (line:%d)\n", line);
        exit(1);
    }
}

static inline void veo_free_mem_wrapper(const vesolver_context_t* context, uint64_t addr, int line) {
    int rc = veo_free_mem(context->instance->proc, addr);
    if (rc != 0) {
        printf("WARNING: Cannot free memory on VE (line:%d)\n", line);
    }
}

static inline void veo_read_mem_wrapper(const vesolver_context_t* context, void* dst, uint64_t src, size_t size, int line) {
    int rc = veo_read_mem(context->instance->proc, dst, src, size);
    if (rc != 0) {
        printf("ERROR: Cannot write memory to VE (line:%d)\n", line);
        exit(1);
    }
}

static inline void veo_write_mem_wrapper(const vesolver_context_t* context, uint64_t dst, const void* src, size_t size, int line) {
    int rc = veo_write_mem(context->instance->proc, dst, src, size);
    if (rc != 0) {
        printf("ERROR: Cannot write memory to VE (line:%d)\n", line);
        exit(1);
    }
}

static inline uint64_t veo_call_async_wrapper(const vesolver_context_t* context, int symid, veo_args* argp, int line) {
    uint64_t id = veo_call_async(context->ctx, veo_func_table[symid].vemva, argp);
    if (id == VEO_REQUEST_ID_INVALID) {
        printf("ERROR: veo_call_async(\"%s\") failed (line:%d)\n", veo_func_table[symid].symname, line);
        exit(1);
    }
    return id;
}

static inline uint64_t veo_call_wait_wrapper(const vesolver_context_t* context, uint64_t reqid, int line) {
    uint64_t retval;

    int rc = veo_call_wait_result(context->ctx, reqid, &retval);
    if (rc != VEO_COMMAND_OK) {
        printf("ERROR: veo_call_wait_result failed with rc=%d (line:%d)\n", rc, line);
        exit(1);
    }
    return retval;
}

static inline uint64_t veo_call_sync_wrapper(const vesolver_context_t* context, int symid, veo_args* argp, int line) {
    uint64_t retval;

    int rc = veo_call_sync(context->instance->proc, veo_func_table[symid].vemva, argp, &retval);
    if (rc != VEO_COMMAND_OK) {
        printf("ERROR: veo_call_sync(\"%s\") failed with rc=%d (line:%d)\n", veo_func_table[symid].symname, rc, line);
        exit(1);
    }
    return retval;
}

/*
 * API functions
 */
int vesolver_init() {
    int venum = 0;
    int version = veo_api_version();

    printf("INFO: SX-Aurora TSUBASA VE Offloading API version: %d\n", version);

    if (version != 9) {
        printf("WARNING: unsupported VE Offloading API.\n");
    }

    char* venum_str = getenv("VE_NODE_NUMBER");
    if (venum_str != NULL) {
        venum = atoi(venum_str);
    }
    printf("INFO: Using VE%d\n", venum);

    vesolver_instance_t* instance = &instances[0];

    instance->proc = veo_proc_create(venum);
    if (instance->proc == NULL) {
        printf("ERROR: Creating VE offload process failed.\n");
        fflush(stdout);
        return -1;
    }

    char* libpath = getenv("CCX_VEO_LIBRARY_PATH");
    if (libpath == NULL) {
        libpath = CCX_VEO_DEFAULT_LIBRARY_PATH;
    }

    instance->veo_handle = veo_load_library(instance->proc, libpath);
    if (instance->veo_handle == 0) {
        printf("ERROR: Loading %s failed.\n", libpath);
        fflush(stdout);
        return -1;
    }

    for (int i=0; i<VESOLVER_NUM_FUNCS; i++) {
        veo_func_table[i].vemva = veo_get_sym(instance->proc, instance->veo_handle, veo_func_table[i].symname);
        if (veo_func_table[i].vemva == 0) {
            printf("ERROR: Cannot find symbol '%s' in libLinearSolver.so.\n", veo_func_table[i].symname);
            fflush(stdout);
            return -1;
        }
    }

    return 0;
}

void vesolver_finalize() {
    vesolver_instance_t* instance = &instances[0];

    veo_proc_destroy(instance->proc);
    instance->proc = NULL;
    instance->veo_handle = 0;

    for (int i=0; i<VESOLVER_NUM_FUNCS; i++) {
        veo_func_table[i].vemva = 0;
    }
}

vesolver_handle_t vesolver_activate() {
    vesolver_handle_t hdl=0;
    vesolver_context_t* context = &contexts[hdl];
    context->instance = &instances[0];

    context->ctx = veo_context_open(context->instance->proc);
    if (context->ctx == NULL) {
        printf("ERROR: opening VE offloading context failed.\n");
        fflush(stdout);
        return -1;
    }

    veo_args *argp = veo_args_alloc();
    VEO_CALL_SYNC(context, VESOLVER_INIT, argp);
    veo_args_free(argp);

    return hdl;
}

int vesolver_deactivate(vesolver_handle_t hdl) {
    vesolver_context_t* context = &contexts[hdl];

    veo_args *argp = veo_args_alloc();
    VEO_CALL_SYNC(context, VESOLVER_FINALIZE, argp);
    veo_args_free(argp);

    if (context->ctx == NULL) {
        return -1;
    }

    return veo_context_close(context->ctx);
}


int vesolver_set_option(vesolver_handle_t hdl, int id, int value) {
    vesolver_context_t* context = &contexts[hdl];

    switch(id) {
    case VESOLVER_OPTION_SOLVER:
        context->option.solver = value;
        break;
    }
    return 0;
}

#if 0
int vesolver_set_matrix_csr(vesolver_handle_t hdl, INT_T neq,
        INT_T* pointers, INT_T* indice, double* value, uint32_t flags) {
    return -1;
}

int vesolver_set_matrix_csc(vesolver_handle_t hdl, INT_T neq,
        INT_T* pointers, INT_T* indice, double* value, uint32_t flags) {
    veo_args *argp = veo_args_alloc();
    uint64_t dptr, aptr, aind, val;
    matrix_desc_t desc;
    uint64_t ndim = pointers[neq];

    desc.flags = flags;
    desc.neq = neq;

    // send matrix description
    VEO_ALLOC_MEM(&dptr, sizeof(matrix_desc_t));
    VEO_WRITE_MEM(dptr, &desc, (neq+1)*sizeof(int));
    veo_args_set_u64(argp, 0, dptr);

    // Send matrix
    VEO_ALLOC_MEM(&aptr, (neq+1)*sizeof(int));
    VEO_WRITE_MEM(aptr, pointers, (neq+1)*sizeof(int));
    veo_args_set_u64(argp, 1, aptr);

    VEO_ALLOC_MEM(&aind, ndim*sizeof(int));
    VEO_WRITE_MEM(aind, indice, ndim*sizeof(int));
    veo_args_set_u64(argp, 2, aind);

    VEO_ALLOC_MEM(&val, ndim*sizeof(double));
    VEO_WRITE_MEM(val, value, ndim*sizeof(double));
    veo_args_set_u64(argp, 3, val);

    uint64_t id = VEO_CALL_ASYNC(VESOLVER_CREATE_HANDLE, argp);
    uint64_t retval = VEO_CALL_WAIT(id);
    veo_args_free(argp);
    return retval;
}
#endif

matrix_desc_t* vesolver_alloc_matrix(vesolver_handle_t hdl, INT_T neq,
        INT_T nnz, uint32_t flags) {
    vesolver_context_t* context = &contexts[hdl];

    matrix_desc_t* desc = (matrix_desc_t*)malloc(sizeof(matrix_desc_t));
    desc->neq = neq;
    desc->nnz = nnz;
    desc->flags = flags;
    desc->pointers = (INT_T*)calloc(neq+1, sizeof(INT_T));
    desc->indice = (INT_T*)calloc(nnz, sizeof(INT_T));
    desc->values = (double*)calloc(nnz, sizeof(double));
    context->matdesc = desc;

    // Create VE Solver descriptor & Allocate memory on VE
    vesolver_desc_t* vedesc = (vesolver_desc_t*)malloc(sizeof(vesolver_desc_t));

    vedesc->ves_handle = hdl;
    vedesc->neq = neq;
    vedesc->nnz = nnz;
    vedesc->flags = flags;
    vedesc->options.solver = context->option.solver;
    VEO_ALLOC_MEM(context, &(vedesc->dptr), sizeof(vesolver_desc_t));
    VEO_ALLOC_MEM(context, &(vedesc->pointers), (neq+1)*sizeof(int));
    VEO_ALLOC_MEM(context, &(vedesc->indice), nnz*sizeof(int));
    VEO_ALLOC_MEM(context, &(vedesc->values), nnz*sizeof(double));

    VEO_ALLOC_MEM(context, &(vedesc->bptr), neq*sizeof(double));
    VEO_ALLOC_MEM(context, &(vedesc->xptr), neq*sizeof(double));
    context->vedesc = vedesc;

    return desc;
}

int vesolver_set_matrix(vesolver_handle_t hdl, matrix_desc_t *desc) {
    vesolver_context_t* context = &contexts[hdl];
    vesolver_desc_t* vedesc = context->vedesc;
    veo_args *argp = veo_args_alloc();
    INT_T neq = desc->neq;
    INT_T ndim = desc->nnz;

    if (hdl != vedesc->ves_handle) {
        fprintf(stderr, "ERROR:vesolver_set_matrix: Invalid handle.\n");
        return -1;
    }
    // send matrix description & data
    VEO_WRITE_MEM(context, vedesc->dptr, context->vedesc, sizeof(vesolver_desc_t));
    VEO_WRITE_MEM(context, vedesc->pointers, desc->pointers, (neq+1)*sizeof(int));
    VEO_WRITE_MEM(context, vedesc->indice, desc->indice, ndim*sizeof(int));
    VEO_WRITE_MEM(context, vedesc->values, desc->values, ndim*sizeof(double));
    veo_args_set_u64(argp, 0, vedesc->dptr);

    uint64_t id = VEO_CALL_ASYNC(context, VESOLVER_CREATE_HANDLE, argp);
    uint64_t retval = VEO_CALL_WAIT(context, id);
    veo_args_free(argp);

    return (retval == vedesc->dptr) ? 0 : 1;
}

int vesolver_free_matrix(vesolver_handle_t hdl) {
    vesolver_context_t* context = &contexts[hdl];
    vesolver_desc_t* vedesc = context->vedesc;
    matrix_desc_t* desc = context->matdesc;

    veo_args *argp = veo_args_alloc();
    veo_args_set_u64(argp, 0, vedesc->dptr);
    VEO_CALL_SYNC(context, VESOLVER_FREE_HANDLE, argp);
    veo_args_free(argp);

    VEO_FREE_MEM(context, vedesc->dptr);
    VEO_FREE_MEM(context, vedesc->pointers);
    VEO_FREE_MEM(context, vedesc->indice);
    VEO_FREE_MEM(context, vedesc->values);
    VEO_FREE_MEM(context, vedesc->bptr);
    VEO_FREE_MEM(context, vedesc->xptr);
    free(vedesc);
    context->vedesc = NULL;

    free(desc->pointers);
    free(desc->indice);
    free(desc->values);
    free(desc);
    context->matdesc = NULL;
    return 0;
}

int vesolver_solve_sync(vesolver_handle_t hdl, double* b, double* x, double res) {
    vesolver_context_t* context = &contexts[hdl];
    vesolver_desc_t* vedesc = context->vedesc;

    veo_args *argp = veo_args_alloc();
    INT_T neq = vedesc->neq;

    vedesc->options.res = res;
    VEO_WRITE_MEM(context, vedesc->dptr, context->vedesc, sizeof(vesolver_desc_t));
    VEO_WRITE_MEM(context, vedesc->bptr, b, neq*sizeof(double));
    veo_args_set_u64(argp, 0, vedesc->dptr);

    uint64_t retval = VEO_CALL_SYNC(context, VESOLVER_SOLVE, argp);
    veo_args_free(argp);

    VEO_READ_MEM(context, x, vedesc->xptr, neq*sizeof(double));

    return retval;
}

#if 0
int vesolver_solve_async();
int vesolver_solve_wait();
#endif

#if 0
void aurora_hs_main(double *ad, double *au, double *adb, double *aub, double *sigma,
        double *b, ITG *icol, ITG *irow, ITG *neq, ITG *nzs){
    ITG *pointers=NULL;
    ITG *indice=NULL;
    double *value=NULL;
    long long ndim = (*neq)+(*nzs);

    if (veo_handle <= 0) {
        return;
    }

    mrow = (*neq);
    ncol = (*neq);

    // Call VE function named "factorize_solve"
    struct veo_args *argp1 = veo_args_alloc();
    uint64_t id, retval;
    uint64_t bptr, xptr;

    // Call VE function named "factorize"
    veo_args_set_i32(argp1, 0, mrow);
    veo_args_set_i32(argp1, 1, ncol);

    VEO_ALLOC_MEM(&aptr, (*neq+1)*sizeof(int));
    VEO_WRITE_MEM(aptr, pointers, (*neq+1)*sizeof(int));
    veo_args_set_u64(argp1, 2, aptr);

    VEO_ALLOC_MEM(&aind, ndim*sizeof(int));
    VEO_WRITE_MEM(aind, indice, ndim*sizeof(int));
    veo_args_set_u64(argp1, 3, aind);

    VEO_ALLOC_MEM(&val, ndim*sizeof(double));
    VEO_WRITE_MEM(val, value, ndim*sizeof(double));
    veo_args_set_u64(argp1, 4, val);

    id = VEO_CALL_ASYNC(LIBCCX_FACTORIZE, argp1);

    // Call VE function named "solve"
    struct veo_args *argp2 = veo_args_alloc();

    VEO_ALLOC_MEM(&bptr, mrow*sizeof(double));
    VEO_WRITE_MEM(bptr, b, mrow*sizeof(double));
    veo_args_set_u64(argp2, 0, bptr);

    VEO_ALLOC_MEM(&xptr, ncol*sizeof(double));
    veo_args_set_u64(argp2, 1, xptr);

    retval = VEO_CALL_WAIT(id);
    if (retval !=0) {
        printf("ERROR: %s:factorize failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    }
    id = VEO_CALL_ASYNC(LIBCCX_SOLVE, argp2);
    SFREE(pointers);
    SFREE(indice);
    SFREE(value);
    veo_args_clear(argp1);

    retval = VEO_CALL_WAIT(id);
    if (retval !=0) {
        printf("ERROR: %s:solve failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    }
    veo_read_mem(proc, b, xptr, ncol*sizeof(double));

    // Call VE function named "finalize"
    id = VEO_CALL_ASYNC(LIBCCX_FINALIZE, argp1);
    veo_args_free(argp2);
    veo_free_mem(proc, bptr);
    veo_free_mem(proc, xptr);
    retval = VEO_CALL_WAIT(id);
    if (retval !=0) {
        printf("ERROR: %s:finalize failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    }
    veo_args_free(argp1);

    veo_free_mem(proc, aptr);
    veo_free_mem(proc, aind);
    veo_free_mem(proc, val);

    return;
}

void aurora_cg_main(double *ad, double *au, double *adb, double *aub, double *sigma,
        double *b, ITG *icol, ITG *irow, ITG *neqp, ITG *nzsp){
    ITG *pointers=NULL;
    ITG *indice=NULL;
    double *value=NULL;
    ITG ndim = (*neqp)+(*nzsp);

    if (veo_handle <= 0) {
        return;
    }

    printf("Solving the system of equations using the iterative solver on VE\n\n");
    NNEW(pointers,ITG,*neqp+1);
    NNEW(indice,ITG,ndim);
    NNEW(value,double,ndim);
    set_matrixes(ad, au, adb, aub, sigma, icol, irow, neqp, nzsp, pointers, indice, value, ndim);
    ITG neq = *neqp;

    // Call VE function named "cgve_solve"
    // uint64_t cgve_solve(ITG neq, ITG len, uint64_t aptr, uint64_t aind, uint64_t val, uint64_t bptr, uint64_t xptr);
    struct veo_args *argp = veo_args_alloc();
    uint64_t retval;
    uint64_t bptr, xptr;

    veo_args_set_i32(argp, 0, neq);
    veo_args_set_i32(argp, 1, ndim);

    VEO_ALLOC_MEM(&aptr, (neq+1)*sizeof(int));
    VEO_WRITE_MEM(aptr, pointers, (neq+1)*sizeof(int));
    veo_args_set_u64(argp, 2, aptr);

    VEO_ALLOC_MEM(&aind, ndim*sizeof(int));
    VEO_WRITE_MEM(aind, indice, ndim*sizeof(int));
    veo_args_set_u64(argp, 3, aind);

    VEO_ALLOC_MEM(&val, ndim*sizeof(double));
    VEO_WRITE_MEM(val, value, ndim*sizeof(double));
    veo_args_set_u64(argp, 4, val);

    VEO_ALLOC_MEM(&bptr, neq*sizeof(double));
    VEO_WRITE_MEM(bptr, b, neq*sizeof(double));
    veo_args_set_u64(argp, 5, bptr);

    VEO_ALLOC_MEM(&xptr, neq*sizeof(double));
    veo_args_set_u64(argp, 6, xptr);
    retval = VEO_CALL_SYNC(LIBCCX_CGVE_SOLVE, argp);
    if (retval !=0) {
        printf("ERROR: %s failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    }
    SFREE(pointers);
    SFREE(indice);
    SFREE(value);
    veo_args_free(argp);
    veo_read_mem(proc, b, xptr, neq*sizeof(double));
    veo_free_mem(proc, bptr);
    veo_free_mem(proc, xptr);
    veo_free_mem(proc, aptr);
    veo_free_mem(proc, aind);
    veo_free_mem(proc, val);

    return;
}
#endif
