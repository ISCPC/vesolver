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

SUBROUTINE ElmerSolverAPI(neq, ndim, pointers, indice, val, b, x, solverId, preCond)
    USE Types
    USE SolverUtils
    USE Messages

    INTEGER :: neq, ndim
    INTEGER, target :: pointers(neq+1),indice(ndim)
    double precision, target :: val(ndim)
    double precision :: x(neq),b(neq)
    INTEGER :: solverId
    INTEGER :: preCond
    INTEGER :: i,j
!    INTEGER, pointer :: pointers(:),indice(:)
!    double precision, pointer :: val(:)
!    double precision :: x(:),b(:)
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t) :: Solver
!    double precision :: Norm
!    INTEGER :: DOFs
!
    WRITE(*,'(A,I10,I10)') 'INFO: ', neq, ndim
!    DO i=1,10
!        WRITE(*,'(A,I10,I10,F20.6)') 'INFO:FORTRAN: ', pointers(i), indice(i), val(i)
!    END DO
!    call flush(6)
!
    allocate(A)
    allocate(A % Diag(neq))
    A % FORMAT = MATRIX_CRS
    A % NumberOfRows = neq
    A % Rows => pointers
    A % Cols => indice
    A % Values => val
    A % COMPLEX = .FALSE.

    DO i=1,neq
        DO j=A % Rows(i),A % Rows(i+1)-1
            IF (A % Cols(j).eq.i) THEN
                A % Diag(i) = j
                EXIT
            END IF
        END DO
    END DO

    !DO i=1,neq
    !  write(*,*) 'INFO:Matrix_A:Diag: ', i, A % Diag(i)
    !END DO
    
    !Norm = 0._dp
    !DOFs = 0

    !DO i=1,A % NumberOfRow
    !DO i=1,5
    !    DO j=A % Rows(i),A % Rows(i+1)-1
    !        WRITE(*,'(A,I10,I10,F20.6)') 'INFO:FORTRAN: ', A % Rows(i), A % Cols(j), A % Values(j)
    !    END DO
    !END DO

!    Linear System Max Iterations = 500
!    Linear System Convergence Tolerance = 1.0e-10
!    BiCGstabl polynomial degree = 2
!    Linear System Preconditioning = ILU0
!    Linear System ILUT Tolerance = 1.0e-3
!    Linear System Abort Not Converged = False
!    Linear System Residual Output = 10
!    Linear System Precondition Recompute = 1

    CALL ListAddNewInteger( Solver % Values, 'Linear System Max Iterations', 5000 )
    CALL ListAddNewConstReal( Solver % Values, 'Linear System Convergence Tolerance', 1.0D-10 )
    CALL ListAddNewInteger( Solver % Values, 'BiCGstabl polynomial degree', 2)
    CALL ListAddNewConstReal( Solver % Values, 'Linear System ILUT Tolerance', 1.0D-3 )
    CALL ListAddNewLogical( Solver % Values, 'Linear System Abort Not Converged', .FALSE. )
    CALL ListAddNewInteger( Solver % Values, 'Linear System Residual Output', 10 )
    CALL ListAddNewInteger( Solver % Values, 'Linear System Precondition Recompute', 100 )
    CALL ListAddNewLogical( Solver % Values, 'Linear System Complex', .FALSE. )
    !CALL ListAddNewLogical( Solver % Values, 'Linear System Symmetric', .TRUE. )
    MaxOutputLevel = 32

    SELECT CASE(solverId)
    CASE(0)
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'BiCGStab' )
    CASE(1)
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'BiCGStab2' )
    CASE(2)
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'BiCGStabl' )
    CASE(3)
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'CG' )
    CASE(4)
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'CGS' )
    CASE(5)
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'TFQMR' )
    CASE(6)
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'GMRES' )
    CASE(7)
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'SGS' )
    CASE(8)
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'GCR' )
    CASE(9)
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'IDRS' )
    CASE DEFAULT
        CALL ListAddNewString( Solver % Values, 'Linear System Iterative Method', 'BiCGStab' )
    END SELECT

    SELECT CASE(preCond)
    CASE(1)
        CALL ListAddNewString( Solver % Values, 'Linear System Preconditioning', 'Diagonal' )
    CASE(10)
        CALL ListAddNewString( Solver % Values, 'Linear System Preconditioning', 'ILU0' )
        !CALL ListAddNewLogical( Solver % Values, 'Edge Basis', .TRUE. )
    CASE(2)
        CALL ListAddNewString( Solver % Values, 'Linear System Preconditioning', 'ILUT' )
    CASE(3)
        CALL ListAddNewString( Solver % Values, 'Linear System Preconditioning', 'Multigrid' )
    END SELECT

    !CALL ListAddNewString( Solver % Values, 'Linear System Preconditioning', 'ILU1' )
    !CALL ListAddNewConstReal( Solver % Values, 'Linear System ILUT Tolerance', 1.0D-3 )
    !CALL ListAddNewInteger( Solver % Values, 'Linear System Precondition Recompute', 1 )

!    CALL SolveConstraintModesSystem(A, x, b, Solver)
!    CALL SolveLinearSystem(A, x, b, Norm, DOFs, Solver)
    CALL IterSolver(A, x, b, Solver)
END SUBROUTINE ElmerSolverAPI

SUBROUTINE pdseupd
  RETURN
END SUBROUTINE pdseupd

SUBROUTINE pdneupd
  RETURN
END SUBROUTINE pdneupd

SUBROUTINE pzneupd
  RETURN
END SUBROUTINE pzneupd

SUBROUTINE pdsaupd
  RETURN
END SUBROUTINE pdsaupd

SUBROUTINE pdnaupd
  RETURN
END SUBROUTINE pdnaupd

SUBROUTINE pznaupd
  RETURN
END SUBROUTINE pznaupd

