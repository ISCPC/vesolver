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

/*
 *
 */
VEServer::VEServer() {
    ves_info.status = VES_ST_INIT;
}

VEServer::~VEServer() {
    ves_info.status = VES_ST_INIT;
}

int VEServer::init(MPI_Comm comm) {
    sv_comm = comm;

    /* Step 1: Send VESinfo to vesolver members to let them know world rank of 
        vesolver master process */
    MPI_Comm_size(comm, &(ves_info.nprocs));
    MPI_Comm_rank(MPI_COMM_WORLD, &ves_info.server_leader);
    MPI_Comm_rank(sv_comm, &myrank);

    MPI_Comm_group( sv_comm, &sv_group );
    rankmap = (int*)calloc(sizeof(int), ves_info.nprocs);
    for (int i=0; i<ves_info.nprocs; i++) {
        rankmap[i] = i;
    }

    ves_info.status = VES_ST_INACTIVE;
    MPI_Bcast(&ves_info, sizeof(VES_info_t), MPI_BYTE, 0, sv_comm);

    if (ves_info.server_leader == 0) {
        sv_comm = MPI_COMM_WORLD;
        solver_comm = MPI_COMM_WORLD;
        return -1; /* Standalone mode */
    }

    /* Step 2: Send VESinfo to all members */
    MPI_Bcast(&ves_info, sizeof(VES_info_t), MPI_BYTE, ves_info.server_leader, MPI_COMM_WORLD);
    return 0;
}

void VEServer::fini() {
    if (rankmap) {
        free(rankmap);
    }
    ves_info.status = VES_ST_INIT;
};

int VEServer::wait_activate() {
    int root;
    MPI_Status status;

    if (myrank == 0) {
        VES_info_t newinfo;
        MPI_Recv(&newinfo, sizeof(VES_info_t), MPI_BYTE, MPI_ANY_SOURCE, VES_ACTIVATE_TAG, MPI_COMM_WORLD, &status);
        ves_info.client_leader = newinfo.client_leader;
        ves_info.nactives = newinfo.nactives;
        root = MPI_ROOT;
    } else {
        root = MPI_PROC_NULL;
    }

    MPI_Bcast(&ves_info, sizeof(VES_info_t), MPI_BYTE, 0, sv_comm);
    if (ves_info.client_leader < 0) {
        printf("INFO:VEServer[%d]: got exit request.\n", myrank);
        return -2;
    }

    MPI_Group_incl( sv_group, ves_info.nactives, rankmap, &solver_group );
    MPI_Comm_create( sv_comm, solver_group, &solver_comm );

    if (myrank < ves_info.nactives) {
        MPI_Intercomm_create(solver_comm, 0, MPI_COMM_WORLD, ves_info.client_leader, VES_ACTIVATE_TAG, &ves_comm);
        ves_info.status = VES_ST_ACTIVE;
        MPI_Bcast(&ves_info, sizeof(VES_info_t), MPI_BYTE, root, ves_comm);
        printf("INFO:VEServer[%d]: activated. (leader=%d)\n", myrank, ves_info.client_leader);
    } else {
        printf("INFO:VEServer[%d]: inactivated.\n", myrank);
        ves_info.status = VES_ST_INACTIVE;
    }

    return 0;
}

void VEServer::deactivate() {
    if (ves_info.status == VES_ST_ACTIVE) {
        MPI_Comm_free(&ves_comm);
        MPI_Comm_free(&solver_comm);
    }
    MPI_Group_free(&solver_group);
    ves_info.status = VES_ST_INACTIVE;
};

int VEServer::wait_request(VES_request_t& req) {
    if (ves_info.status != VES_ST_ACTIVE) {
        return -1;
    }

    printf("INFO:VEServer[%d]: waiting request...\n", myrank);
    MPI_Bcast(&req, sizeof(VES_request_t), MPI_BYTE, 0, ves_comm);

    if (req.req_type == VES_REQ_QUIT) {
        return -2;
    }

    return 0;
}

int VEServer::process(VES_request_t& req) {
    return 0;
}
