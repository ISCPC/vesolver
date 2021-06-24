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
#include <mpi.h>

/*
 * Common definition
 */
#define VES_COLOR_HOST    0
#define VES_COLOR_SERVER  100

/*
 * Tags for communication
 */
#define VES_ACTIVATE_TAG 1
#define VES_REQUEST_TAG  2

/*
 * VEserver status
 */
typedef struct ves_info {
    int32_t status;       // server process status
#define VES_ST_INIT      0
#define VES_ST_INACTIVE  1
#define VES_ST_ACTIVE    2
    int32_t client_leader; // the world rank of client leader
    int32_t server_leader; // the world rank of server leader
    int32_t nprocs;    // the number of server processes
    int32_t nactives;  // the number of active server processes
} VES_info_t;

/*
 * VEserver request packet
 */
typedef struct {
    int32_t req_type;
#define VES_REQ_QUIT             0
#define VES_REQ_LINEARSOLVER     1
#define VES_REQ_NUM              2
    int32_t dummy;
#define VES_REQ_PAYLOAD_SIZE 4080
    uint8_t payload[VES_REQ_PAYLOAD_SIZE];
} VES_request_t;

#define VES_REQ_GETPAYLOAD(req) (void*)(&((req).payload))


/*
 * VEServer definition
 */
class VEServer {
public:
    VEServer();
    ~VEServer();

    int init(MPI_Comm comm);
    void fini();

    int wait_activate();
    void deactivate();

    int wait_request(VES_request_t& req);
    virtual int process(VES_request_t& req);

protected:
    MPI_Comm ves_comm;     /* intercommunicator */
    MPI_Comm solver_comm;  /* intracommunicator for active solver */
    VES_info_t ves_info;
    int myrank;

private:
    int *rankmap;
    MPI_Comm sv_comm;      /* intracommunicator for server world */
    MPI_Group sv_group, solver_group;
};

/*
 * API definition to commuicat VEServer
 */
class VEServerAPI {
public:
    VEServerAPI();
    ~VEServerAPI();

    int init(MPI_Comm comm);
    void fini();

    int activate(MPI_Comm comm, int nprocs);
    void deactivate();

    int send_request(VES_request_t& req);
    inline MPI_Comm get_ves_comm() { return ves_comm; }
    inline MPI_Comm get_cl_comm() { return cl_comm; }

protected:
    MPI_Comm ves_comm; /* intercommunicator */
    MPI_Comm cl_comm;  /* intracommunicator for client world */
    int myrank;

private:
    VES_info_t ves_info;
};
