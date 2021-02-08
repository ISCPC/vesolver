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
#include <mpi.h>

#define VES_MODE_GATHER_ON_VH  0
#define VES_MODE_GATHER_ON_VE  1
#define VES_MODE_SYMMETRIC     2

#define VESOLVER_HS            0
#define VESOLVER_BICGSTAB2   100
#define VESOLVER_DUMMY      9999

#ifdef __cplusplus
extern "C" {
#endif

    int vesolver_init(MPI_Comm comm);
    int vesolver_init_f(int comm);
    int vesolver_fini();
    int vesolver_activate(MPI_Comm comm, int nprocs);
    int vesolver_activate_f(int comm, int nprocs);
    int vesolver_deactivate();
    int vesolver_solve(int solver, int32_t mtype, int32_t neq, int32_t *pointers, int32_t *indice, double *value, double *b, double *x, double res);
    int vesolver_psolve(int32_t mode, int32_t solver, int32_t neq, int32_t nrows, int32_t *pointers, int32_t *indice, double *value, int32_t *rorder, double *b, double *x, double res);
#if 0
    int vesolver_psolve2(int32_t mode, int32_t solver, int32_t neq, int32_t nrows, int32_t nl, int32_t nt, int32_t *pointers, int32_t *indice, double *value, int32_t *rorder, double *b, double *x, double res);
    int vesolver_psolve_dcsr(int32_t mode, int32_t solver, int32_t neq, int32_t *pointers, int32_t *indice, double *value, int32_t nl, int32_t nt, double *b, double *x, double res);
#endif

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#include "veserver.hpp"
#include "SpMatrix.hpp"

class VESolverAPI : public VEServerAPI {
public:
    VESolverAPI() {}
    ~VESolverAPI() {}

    int solve(int solver, SpMatrix& A, Vector& b, Vector& x, double res);
    int solve(int solver, SpDistMatrix **An, DistVector **bn, int n, Vector& x, double res);
    int solve(int solver, SpDistMatrix& A, DistVector& b, Vector& x, double res, int mode);
    int solve_gather_on_vh(int solver, int32_t neq, int32_t nrows, int32_t *pointers, int32_t *indice, double *value, int32_t *rorder, double *b, double *x, double res);

private:
    int send_matrix_info(SpMatrix& A);
    int send_matrix_info(int rank, SpDistMatrix& A);
    int send_matrix_data(SpMatrix& A, Vector& b);
    int send_matrix_data(int rank, SpDistMatrix& A, DistVector& b);
    inline int receive_result(Vector& x);
    inline int receive_result(int src, Vector& x);

    int solve_gather_on_vh(int solver, SpDistMatrix& A, DistVector& b, Vector& x, double res);
    int solve_gather_on_ve(int solver, SpDistMatrix& A, DistVector& b, Vector& x, double res);
    int solve_symmetric(int solver, SpDistMatrix& A, DistVector& b, Vector& x, double res);

    inline void send_solve_request(int mode, int solver, int type, int nprocs, double res);
};
#endif
