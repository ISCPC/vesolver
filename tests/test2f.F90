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
! DCSR format with One-based indexing. 
! NOTICE: 
!         Here, error checking is omitted.
!         Practical programs must surely check errors.
#define VES_MODE_GATHER_ON_VH  0
#define VES_MODE_GATHER_ON_VE  1
#define VES_MODE_SYMMETRIC     2
    
#define VESOLVER_HS            0
#define VESOLVER_BICGSTAB2   100
#define VESOLVER_DUMMY      9999

program example_hybrid_unsym_dcsr_distab_indexing_1
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
   !
   ! Here, both A and B are distributed into MPI process 0 and 1.
   !  + MPI process 0:
   !              Matrix A                             B
   !       [ 1.1, 1.2,   0,   0,   0,   0,   0 ]     [  3.5]
   !       [   0, 2.2, 2.3,   0,   0,   0,   0 ]     [ 11.3]
   !       [   0,   0, 3.3, 3.4,   0,   0,   0 ]     [ 66.9]
   !       [ 4.1,   0,   0, 4.4,   0,   0,   0 ]     [ 21.7]
   !       [   0,   0,   0,   0,   0,   0,   0 ]     
   !       [   0,   0,   0,   0,   0,   0,   0 ]     
   !       [   0,   0,   0,   0,   0,   0,   0 ]
   !
   !  + MPI process 1:
   !              Matrix A                             B
   !       [   0,   0,   0,   0,   0,   0,   0 ]
   !       [   0,   0,   0,   0,   0,   0,   0 ]
   !       [   0,   0,   0,   0, 3.5,   0, 3.7 ]     [ 66.9]
   !       [   0,   0,   0,   0,   0,   0,   0 ]     [ 21.7]
   !       [   0,   0, 5.3,   0, 5.5, 5.6,   0 ]     [ 77.0]
   !       [   0,   0,   0,   0, 6.5, 6.6, 6.7 ]     [119.0]
   !       [   0,   0, 7.3,   0,   0, 7.6, 7.7 ]     [121.4]
   ! The number of rows of A
   integer, parameter ::  neq=7
   integer :: nrow=0
   
   ! The number of columns of A
   integer, parameter ::  ncol=7
   
   ! Values of the nonzero elements of A in row-major order
   real(8), allocatable :: aval(:)
   ! Column indices of nonzero elements
   integer, allocatable :: iaind(:)
   ! Starting points of the rows of the arrays aval and iaind
   integer, allocatable :: iaptr(:)
   ! Starting points of the rows of the arrays aval and iaind
   integer, allocatable :: order(:)
   ! The minimum row index of non-zero elements
   !  that belongs to each MPI process (1 <= mingind <= mrow)
   integer :: mingind = 0
   ! The maximum row index of non-zero elements
   !  that belongs to each MPI process (1 <= maxgind <= mrow)
   integer :: maxgind = 0
   ! The number of right-hand side vectors
   integer :: nrhs = 1
   ! The right-hand side vector b
   real(8), allocatable :: b(:)
   
   ! The solution vector x
   real(8), allocatable :: x(:)
   ! The target accuracy for the residual norm
   real(8) :: res = 1.0e-13 
   
   ! Local variables
   integer :: ierr, comm, myrank, nprocs, i, j, lx
   ! macro definition

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
   call MPI_Comm_split(MPI_COMM_WORLD, 0, myrank, comm, ierr)
   call MPI_Comm_size(comm, nprocs, ierr)
   
   ! Input Preparation
   if(myrank == 0) then
      ! The row numbers 1, 2, 3, and 4 are assigned to the MPI process 0.
      mingind=1 
      maxgind=4
      nrow=maxgind-mingind+1
      allocate(aval(8))
      allocate(iaind(8))
      allocate(iaptr(nrow+1))
      allocate(order(nrow))
      allocate(b(nrow), x(ncol))
      aval  = (/ 1.1d0, 1.2d0, 2.2d0, 2.3d0, 3.3d0, 3.4d0, 4.1d0, 4.4d0 /)
      iaind = (/ 1, 2, 2, 3, 3, 4, 1, 4 /)
      iaptr = (/ 1, 3, 5, 7, 9 /)
      order = (/ 1, 2, 3, 4 /)
      b     = (/ 3.5d0, 11.3d0, 33.9d0, 10.0d0 /)
   else if(myrank == 1) then
      ! The row numbers 3, 4, 5, 6 and 7 are assigned to the MPI process 1.
      mingind=3
      maxgind=7
      nrow=maxgind-mingind+1
      allocate(aval(11))
      allocate(iaind(11))
      allocate(iaptr(nrow+1))
      allocate(order(nrow))
      allocate(b(nrow), x(ncol))
      aval  = (/ 3.5d0, 3.7d0, 5.3d0, 5.5d0, 5.6d0, 6.5d0, 6.6d0, 6.7d0, 7.3d0, 7.6d0, 7.7d0 /)
      iaind = (/ 3, 5, 1, 3, 4, 3, 4, 5, 1, 4, 5 /)
      iaptr = (/ 1, 3, 3, 6, 9, 12 /)
      order = (/ 3, 4, 5, 6, 7 /)
      b = (/ 33.0d0, 11.7d0, 77.0d0, 119.0d0, 121.4d0 /)
   end if

    call VESolver_Init(comm, ierr)
    IF (ierr.ne.0) THEN
        write(*,*) "ERROR: vesolver_init failed."
        call MPI_Finalize(ierr)
        STOP
    END IF

   call VESolver_Activate(comm, 1, ierr)
   !call VESolver_PSolve(VES_MODE_GATHER_ON_VE, VESOLVER_BICGSTAB2, &
   call VESolver_PSolve(VES_MODE_GATHER_ON_VH, VESOLVER_BICGSTAB2, &
        neq, nrow, aval, iaptr, iaind, order, b, x, res, ierr)
   call VESolver_Deactivate()
   call VESolver_Fini()
   
   ! Print the solution vector
   if (myrank == 0) then
      write (*,"(A26)") "******** Solution ********"
      do j = 1, nrhs
          lx = ncol
          do i = 1, ncol
              write (*,"(2X,A2,I1,A4,f14.12)") "x(",i,") = ",x( i + (j-1)*lx)
          end do
      end do
      write (*,"(A26)") "********** End ***********"
      deallocate(aval, iaind, iaptr, b, x)
   else
      deallocate(aval, iaind, iaptr, b)
   end if
   call mpi_finalize(ierr)
   
end program example_hybrid_unsym_dcsr_distab_indexing_1
