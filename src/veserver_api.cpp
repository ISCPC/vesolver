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
#include <mpi.h>
#include "veserver.hpp"

VEServerAPI::VEServerAPI() {
    ves_info.status = VES_ST_INIT;
}

VEServerAPI::~VEServerAPI() {
    ves_info.status = VES_ST_INIT;
}

int VEServerAPI::init(MPI_Comm comm) {
    int size, ves_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_size(comm, &ves_rank);

    if (ves_rank >= size) {
        fprintf(stderr, "ERROR: vesolver doesn't work correctly.");
        return -1;
    }

    /* Get VEsolver info */
    MPI_Bcast(&ves_info, sizeof(VES_info_t), MPI_BYTE, ves_rank, MPI_COMM_WORLD);
    return 0;
}

void VEServerAPI::fini() {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf("INFO:VEServerAPI::fini: terminating... ");
        fflush(stdout);

        /* Send exit request */
        VES_info_t newinfo;
        newinfo.client_leader = -1;
        MPI_Send(&newinfo, sizeof(VES_info_t), MPI_BYTE, ves_info.server_leader, VES_ACTIVATE_TAG, MPI_COMM_WORLD);
        printf("Done.\n");
    }
    ves_info.status = VES_ST_INIT;
}

int VEServerAPI::activate(MPI_Comm comm, int nprocs) {
    int rank, size;

    cl_comm = comm;
    if (comm != MPI_COMM_NULL) {
        MPI_Comm_rank(cl_comm, &rank);
        MPI_Comm_size(cl_comm, &size);
        if (rank == 0) {
            /* Send Activate request */
            VES_info_t newinfo;
            MPI_Comm_rank(MPI_COMM_WORLD, &(newinfo.client_leader));
            newinfo.nactives = nprocs;
            MPI_Send(&newinfo, sizeof(VES_info_t), MPI_BYTE, ves_info.server_leader, VES_ACTIVATE_TAG, MPI_COMM_WORLD);
        }
        MPI_Intercomm_create(cl_comm, 0, MPI_COMM_WORLD, ves_info.server_leader, VES_ACTIVATE_TAG, &ves_comm);
        MPI_Bcast(&ves_info, sizeof(VES_info_t), MPI_BYTE, 0, ves_comm);
        //printf("INFO:VEServerAPI::activate: VEServer activated. (%d/%d processes)\n", ves_info.nactives, ves_info.nprocs);
    } else {
        ves_info.status = VES_ST_INACTIVE;
    }
    return 0;
}

void VEServerAPI::deactivate() {
    VES_request_t req;

    req.req_type = VES_REQ_QUIT;
    send_request(req);

    MPI_Comm_free(&ves_comm);
    ves_info.status = VES_ST_INIT;
};

int VEServerAPI::send_request(VES_request_t& req) {
    int rank;

    MPI_Comm_rank(cl_comm, &rank);
    if (rank == 0) {
        MPI_Bcast(&req, sizeof(VES_request_t), MPI_BYTE, MPI_ROOT, ves_comm);
    } else {
        MPI_Bcast(&req, sizeof(VES_request_t), MPI_BYTE, MPI_PROC_NULL, ves_comm);
    }

    return 0;
}
