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

#include <stdio.h>
#include <stdint.h>
#include "timelog.h"
#include "vesolver.h"
#include "LinearSolver.h"

/*
uint64_t vesolver_init();
uint64_t vesolver_finalize();
uint64_t vesolver_create_handle_by_desc(uint64_t descp);
uint64_t vesolver_free_handle(uint64_t descp);
uint64_t vesolver_solve(uint64_t descp);
*/

typedef struct vesolver_instance {
    SolverHandle_t hdl;
} vesolver_instance_t;

#define VESOLVER_NUM_INSTANCES 1
static vesolver_instance_t instances[VESOLVER_NUM_INSTANCES];

#define VESOLVER_SUCCESS 0
#define VESOLVER_ERROR   1

uint64_t vesolver_init() {
    int cc = solver_init();
    return (cc<0) ? VESOLVER_ERROR : VESOLVER_SUCCESS;
}

uint64_t vesolver_finalize() {
    solver_finalize();
    return VESOLVER_SUCCESS;
}

uint64_t vesolver_create_handle_by_desc(uint64_t vedescp) {
    vesolver_desc_t* vedesc = (vesolver_desc_t*)vedescp;
    vesolver_instance_t* instance = &instances[vedesc->ves_handle];
    int cc;

    SolverHandle_t hdl = solver_create_handle();
    cc = solver_set_option(hdl, SOLVER_OPTION_SOLVER, (int)(vedesc->options.solver));
    if (cc<0) {
        fprintf(stderr, "ERROR:vesolver_create_handle_by_desc: Invalid solver with %d.\n", cc);
        return 0;
    }

    cc = solver_set_matrix_csr(hdl, vedesc->neq, vedesc->nnz,
            (INT_T*)(vedesc->pointers), (INT_T*)(vedesc->indice),
            (double*)(vedesc->values), vedesc->flags);
    if (cc<0) {
        fprintf(stderr, "ERROR:vesolver_create_handle_by_desc: setMatrix failed with %d.\n", cc);
        return 0;
    }

    instance->hdl = hdl;
    return vedescp;
}

uint64_t vesolver_free_handle(uint64_t vedescp) {
    vesolver_desc_t* vedesc = (vesolver_desc_t*)vedescp;
    vesolver_instance_t* instance = &instances[vedesc->ves_handle];

    int cc = solver_free_handle(instance->hdl);
    return (cc<0) ? VESOLVER_ERROR : VESOLVER_SUCCESS;
}

uint64_t vesolver_solve(uint64_t vedescp) {
    vesolver_desc_t* vedesc = (vesolver_desc_t*)vedescp;
    vesolver_instance_t* instance = &instances[vedesc->ves_handle];

    int cc = solver_solve(instance->hdl, (double*)(vedesc->bptr), (double*)(vedesc->xptr),
            vedesc->options.res);

    return (cc<0) ? VESOLVER_ERROR : VESOLVER_SUCCESS;
}
