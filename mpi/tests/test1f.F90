! MIT License
!
! Copyright (c) 2021 Shunji Uno <shunji_uno@iscpc.jp>
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

! This is a sample program that solves a system of linear equations Ax = b for x.
! Here, the matrix A has an unsymmetric pattern that is stored in
! CSR format with One-based indexing. 
! NOTICE: 
!         Here, error checking is omitted.
!         Practical programs must surely check errors.
#define VES_MODE_GATHER_ON_VH  0
#define VES_MODE_GATHER_ON_VE  1
#define VES_MODE_SYMMETRIC     2

#define VESOLVER_HS            0
#define VESOLVER_BICGSTAB2   100
#define VESOLVER_DUMMY      9999

program example_openmp_unsym_csr_indexing_1
    USE VESolverAPI
    implicit none
    include "mpif.h"
    !             Matrix A                      X        B
    ! [ 1.1, 1.2,   0,   0,   0,   0,   0 ]   [1.0]   [  3.5]
    ! [   0, 2.2, 2.3,   0,   0,   0,   0 ]   [2.0]   [ 11.3]
    ! [   0,   0, 3.3, 3.4, 3.5,   0, 3.7 ]   [3.0]   [ 66.9]
    ! [ 4.1,   0,   0, 4.4,   0,   0,   0 ] * [4.0] = [ 21.7]
    ! [   0,   0, 5.3,   0, 5.5, 5.6,   0 ]   [5.0]   [ 77.0]
    ! [   0,   0,   0,   0, 6.5, 6.6, 6.7 ]   [6.0]   [119.0]
    ! [   0,   0, 7.3,   0,   0, 7.6, 7.7 ]   [7.0]   [121.4]
    ! The number of rows of A
    integer, parameter ::  mrow=7
    ! The number of columns of A
    integer, parameter ::  ncol=7

    ! Values of the nonzero elements of A in row-major order 
    real(8), dimension(19) :: aval = (/ &
        1.1d0, 1.2d0, 2.2d0, 2.3d0, 3.3d0, &
        3.4d0, 3.5d0, 3.7d0, 4.1d0, 4.4d0, &
        5.3d0, 5.5d0, 5.6d0, 6.5d0, 6.6d0, &
        6.7d0, 7.3d0, 7.6d0, 7.7d0 &
        /)
    ! Column indices of nonzero elements
    integer, dimension(19) :: iaind = (/ &
        1, 2, 2, 3, 3, &
        4, 5, 7, 1, 4, &
        3, 5, 6, 5, 6, &
        7, 3, 6, 7 &
        /)
    ! Starting points of the rows of the arrays aval and iaind
    integer, dimension(mrow+1) :: iaptr = (/ 1, 3, 5, 9, 11, 14, 17, 20 /)
    ! The number of right-hand side vectors
    integer :: nrhs = 1
    ! The right-hand side vector b
    real(8), dimension(mrow) :: b = (/ 3.5d0, 11.3d0, 66.9d0, 21.7d0, 77.0d0, 119.0d0, 121.4d0 /)

    ! The solution vector x
    real(8), dimension(ncol) :: x
    ! The target accuracy for the residual norm
    real(8) :: res = 1.0e-13 
    ! Local variables
    integer :: ierr, comm, myrank, nprocs, i, j, lx

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_split(MPI_COMM_WORLD, 0, myrank, comm, ierr)
    call MPI_Comm_size(comm, nprocs, ierr)

    IF (nprocs.ne.1) THEN
        write(*,*) "ERROR: The size of test3 must be 1."
        call MPI_Finalize(ierr)
        STOP
    END IF

    call VESolver_Init(comm, ierr)
    IF (ierr.ne.0) THEN
        write(*,*) "ERROR: vesolver_init failed."
        call MPI_Finalize(ierr)
        STOP
    END IF

    call VESolver_Activate(comm, 1, ierr)
    call VESolver_Solve(VESOLVER_BICGSTAB2, mrow, aval, iaptr, iaind, b, x, res, ierr)
    call VESolver_Deactivate()
    call VESolver_Fini()

    ! Print the solution vector
    write (*,"(A26)") "******** Solution ********"
    do j = 1, nrhs
        lx = ncol
        do i = 1, ncol
            write (*,"(2X,A2,I1,A4,f14.12)") "x(",i,") = ",x( i + (j-1)*lx)
        end do
    end do
    write (*,"(A26)") "********** End ***********"

    call MPI_Finalize(ierr)
end program example_openmp_unsym_csr_indexing_1
