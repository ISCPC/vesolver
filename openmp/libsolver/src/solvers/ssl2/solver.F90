SUBROUTINE SSL2_VCGE(K, NW, N, A, ICOL, B, X, EPS, ICON)
    INTEGER :: K, NW, N
    REAL*8 :: A(K,NW)
    INTEGER :: ICOL(K,NW)
    REAL*8 :: B(N), X(N)
    REAL*8 :: EPS
!
    INTEGER :: IPC=1, ITMAX=10000, ISW=1, IGUSS=0
    INTEGER,PARAMETER :: MAXT=48
    REAL*8 :: OMEGA=0.1
    INTEGER :: ITER, ICON
    REAL*8 :: RZ
!    REAL*8 :: VW(K*NW+N*4), IVW(K*NW+N*4)
    REAL*8 :: VW(N+MAXT,7), IVW(MAXT,2)
!
    ICON=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Symmetric CG Solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef WITH_SSL2
    write(*,*) "INFO: call DVCGE() with IPC=",IPC
!    CALL DVCGE(A, K, NW, N, ICOL, B, IPC, ITMAX, ISW, OMEGA, &
!        EPS, IGUSS, X, ITER, RZ, VW, IVW, ICON)
    CALL DM_VCGE(A, K, NW, N, ICOL, B, IPC, ITMAX, ISW, OMEGA, &
        EPS, IGUSS, X, ITER, RZ, VW, IVW, ICON)
    IF (ICON.ne.0) THEN
        write(*,*) "ERROR: DVCGE() failed with CODE=",ICON
    ELSE
        write(*,*) "INFO: DVCGE() : (ITER, RZ)",ITER,RZ
    END IF
#else
    DO i=1,K
        write(*,*) (A(i,j),j=1,NW)
    END DO
    DO i=1,K
         write(*,*) (ICOL(i,j),j=1,NW)
    END DO
#endif
!
END SUBROUTINE SSL2_VCGE

SUBROUTINE SSL2_VBCSE(K, NW, N, A, ICOL, B, X, EPS, ICON)
    INTEGER :: K, NW, N
    REAL*8 :: A(K,NW)
    INTEGER :: ICOL(K,NW)
    REAL*8 :: B(N), X(N)
    REAL*8 :: EPS
!
    INTEGER :: ITMAX=10000, IGUSS=0
    INTEGER :: ITER, ICON, L=2
    REAL*8 :: VW(K*8)
!
    ICON=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ASymmetric BiCGStab Solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef WITH_SSL2
    write(*,*) "INFO: call DVBCSE()"
!    CALL DVBCSE(A, K, NW, N, ICOL, B, ITMAX, EPS, IGUSS, L, &
!        X, ITER, VW, ICON)
    CALL DM_VBCSE(A, K, NW, N, ICOL, B, ITMAX, EPS, IGUSS, L, &
        X, ITER, ICON)
    IF (ICON.ne.0) THEN
        write(*,*) "ERROR: DVBCSE() failed with CODE=",ICON
    ELSE
        write(*,*) "INFO: DVCGE() : ITER=",ITER
    END IF
#else
    DO i=1,K
        write(*,*) (A(i,j),j=1,NW)
    END DO
    DO i=1,K
         write(*,*) (ICOL(i,j),j=1,NW)
    END DO
#endif
!
END SUBROUTINE SSL2_VBCSE
