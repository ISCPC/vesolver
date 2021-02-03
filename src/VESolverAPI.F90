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

MODULE VESolverAPI

    IMPLICIT NONE
    public VESolver_Init, VESolver_Fini

    INTERFACE
        SUBROUTINE VESolver_Fini() bind (c,name="vesolver_fini")
            USE iso_c_binding
        END SUBROUTINE
    END INTERFACE

    INTERFACE
        SUBROUTINE VESolver_Deactivate() bind (c,name="vesolver_deactivate")
            USE iso_c_binding
        END SUBROUTINE
    END INTERFACE


CONTAINS

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_Init(comm, err)
!--------------------------------------------------------------------------------
        INTEGER, intent(in) :: comm
        INTEGER, intent(out) :: err

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Init_C(comm) bind (c,name="vesolver_init_f")
                USE iso_c_binding
                INTEGER(c_int), VALUE :: comm
            END FUNCTION
        END INTERFACE

        err = VESolver_Init_C(comm)
!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_Init
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_Activate(comm, nprocs, err)
!--------------------------------------------------------------------------------
        INTEGER, intent(in) :: comm, nprocs
        INTEGER, intent(out) :: err

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Activate_C(comm, nprocs) bind (c,name="vesolver_activate_f")
                USE iso_c_binding
                INTEGER(c_int), VALUE :: comm, nprocs
            END FUNCTION
        END INTERFACE

        err = VESolver_Activate_C(comm, nprocs)
!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_Activate
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_Solve(solver, neq, Values, Rows, Cols, b, x, res, err)
!--------------------------------------------------------------------------------
        INTEGER, intent(in) :: solver, neq
        REAL(8), DIMENSION(:) :: Values, b, x
        INTEGER, DIMENSION(:) :: Rows, Cols
        !REAL(8), DIMENSION(:), TARGET CONTIG :: Values, b, x
        !INTEGER, DIMENSION(:), TARGET CONTIG :: Rows, Cols
        REAL(8) :: res
        INTEGER :: err

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Solve_C(solver, mtype, neq, &
                Rows, Cols, Values, b, x, res) bind (c,name="vesolver_solve")
                USE iso_c_binding
                INTEGER(c_int), VALUE, intent(in) :: solver, mtype, neq
                REAL(c_double), DIMENSION(*), intent(in) :: Values, b, x
                INTEGER(c_int), DIMENSION(*), intent(in) :: Rows, Cols
                REAL(c_double), VALUE, intent(in) :: res
            END FUNCTION
        END INTERFACE

        !WRITE(*,*) 'INFO: Send Request to VE (solver, res) = (', solver, res, ').'
        err = VESolver_Solve_C(solver, 17, neq, Rows, Cols, Values, b, x, res)

!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_Solve
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_PSolve(mode, solver, neq, nRows, Values, Rows, Cols, Rorder, b, x, res, err)
!--------------------------------------------------------------------------------
        INTEGER, intent(in) :: mode, solver, neq, nRows
        REAL(8), DIMENSION(:), intent(in) :: Values, b, x
        INTEGER, DIMENSION(:), intent(in) :: Rows, Cols, Rorder
        !REAL(8), DIMENSION(:), TARGET CONTIG :: Values, b, x
        !INTEGER, DIMENSION(:), TARGET CONTIG :: Rows, Cols
        REAL(8) :: res
        INTEGER :: err

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_PSolve_C(mode, solver, neq, &
                nRows, Rows, Cols, Values, Rorder, b, x, res) bind (c,name="vesolver_psolve")
                USE iso_c_binding
                INTEGER(c_int), VALUE, intent(in) :: mode, solver, neq, nRows
                REAL(c_double), DIMENSION(*), intent(in) :: Values, b, x
                INTEGER(c_int), DIMENSION(*), intent(in) :: Rows, Cols, Rorder
                REAL(c_double), VALUE, intent(in) :: res
            END FUNCTION
        END INTERFACE

        !WRITE(*,*) 'INFO: Send Request to VE (solver, res) = (', solver, res, ').'
        err = VESolver_PSolve_C(mode, solver, neq, nRows, &
            Rows, Cols, Values, Rorder, b, x, res)

!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_PSolve
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_PSolve_dcsr(mode, solver, neq, Values, Rows, Cols, nl, nt, b, x, res, err)
!--------------------------------------------------------------------------------
        INTEGER, intent(in) :: mode, solver, neq, nl, nt
        REAL(8), DIMENSION(:), intent(in) :: Values, b, x
        INTEGER, DIMENSION(:), intent(in) :: Rows, Cols
        REAL(8) :: res
        INTEGER :: err

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_PSolve_dcsr_C(mode, solver, neq, &
                Rows, Cols, Values, nl, nt, b, x, res) bind (c,name="vesolver_psolve_dcsr")
                USE iso_c_binding
                INTEGER(c_int), VALUE, intent(in) :: mode, solver, neq, nl, nt
                REAL(c_double), DIMENSION(*), intent(in) :: Values, b, x
                INTEGER(c_int), DIMENSION(*), intent(in) :: Rows, Cols
                REAL(c_double), VALUE, intent(in) :: res
            END FUNCTION
        END INTERFACE

        !WRITE(*,*) 'INFO: Send Request to VE (solver, res) = (', 9999, res, ').'
        err = VESolver_PSolve_dcsr_C(mode, solver, neq, Rows, Cols, Values, nl, nt, b, x, res)

!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_PSolve_dcsr
!--------------------------------------------------------------------------------

END MODULE VESolverAPI
