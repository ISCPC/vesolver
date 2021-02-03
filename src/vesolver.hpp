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
#include "SpMatrix.hpp"
#include "veserver.hpp"

#define VES_MATINFO_TAG     1000
#define VES_MATDATA_TAG     1001
#define VES_MATGATHER_TAG   1002

class VESolver : public VEServer {
public:
    VESolver() {};
    ~VESolver() {};

    int process(VES_request_t& req) override;
    int solve(int solver, SpMatrix& A, Vector& b, Vector& x, double res);
    int solve(int solver, SpDistMatrix **An, DistVector **bn, int n, Vector& x, double res);
    int solve(int solver, SpDistMatrix& A, DistVector& b, Vector& x, double res);

private:
    inline int receive_matrix_info(int src, SpMatrix* A);
    inline int receive_matrix_info(int src, SpDistMatrix* A);
    inline int receive_matrix_data(int src, SpMatrix* A, Vector* b);
    inline int receive_matrix_data(int src, SpDistMatrix* A, DistVector* b);
    inline int send_result(Vector& x); 
    inline int send_result(int dest, Vector& x); 

    int process_gather_on_vh(VES_request_t& req);
    int process_gather_on_ve(VES_request_t& req);
    int process_symmetric(VES_request_t& req);

    int hs_solve(SpMatrix& A, Vector& b, Vector& x, double res);
    int pardiso_solve(SpMatrix& A, Vector& b, Vector& x, double res); 
    int elmer_solve(SpMatrix& A, Vector& b, Vector& x, double res);
    int dummy_solve(SpMatrix& A, Vector& b, Vector& x, double res);
    int hs_solve(SpDistMatrix& A, DistVector& b, Vector& x, double res);
    int cpardiso_solve(SpDistMatrix& A, DistVector& b, Vector& x, double res);
    int elmer_solve(SpDistMatrix& A, DistVector& b, Vector& x, double res);
    int dummy_solve(SpDistMatrix& A, DistVector& b, Vector& x, double res);
};

typedef struct vesolver_request {
    int32_t  mode;
    int32_t  solver;
    uint32_t type;
    int32_t  nprocs;
    double res;
} vesolver_request_t;

typedef struct vesolver_matrix_info {
    uint32_t type;
    int32_t  nrow;
    int32_t  ncol;
    int32_t  nval;
    int32_t  neq;
    int32_t  nl;
    int32_t  nt;
} vesolver_matrix_info_t;
