/*
 * vesolver.c
 *
 *  Created on: Jun 17, 2021
 *      Author: uno
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
    vesolver_instance_t instance = &instances[vedesc->ves_handle];
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
    vesolver_instance_t instance = &instances[vedesc->ves_handle];

    int cc = solver_free_handle(instance->hdl);
    return (cc<0) ? VESOLVER_ERROR : VESOLVER_SUCCESS;
}

uint64_t vesolver_solve(uint64_t vedescp) {
    vesolver_desc_t* vedesc = (vesolver_desc_t*)vedescp;
    vesolver_instance_t instance = &instances[vedesc->ves_handle];

    int cc = solver_solve(instance->hdl, (double*)(vedesc->bptr), (double*)(vedesc->xptr),
            vedesc->options.res);

    return (cc<0) ? VESOLVER_ERROR : VESOLVER_SUCCESS;
}
