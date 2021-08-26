!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DM_VCGE : Symmetric CG solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
#ifdef WITH_SSL2
    write(*,*) "INFO: call DM_DVCGE() with IPC=",IPC
!    CALL DVCGE(A, K, NW, N, ICOL, B, IPC, ITMAX, ISW, OMEGA, &
!        EPS, IGUSS, X, ITER, RZ, VW, IVW, ICON)
    CALL DM_VCGE(A, K, NW, N, ICOL, B, IPC, ITMAX, ISW, OMEGA, &
        EPS, IGUSS, X, ITER, RZ, VW, IVW, ICON)
    IF (ICON.ne.0) THEN
        write(*,*) "ERROR: DM_DVCGE() failed with CODE=",ICON
    ELSE
        write(*,*) "INFO: DM_DVCGE() : (ITER, RZ)",ITER,RZ
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DM_VBCSE: ASymmetric BiCGStab Solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
#ifdef WITH_SSL2
    write(*,*) "INFO: call DM_DVBCSE()"
!    CALL DVBCSE(A, K, NW, N, ICOL, B, ITMAX, EPS, IGUSS, L, &
!        X, ITER, VW, ICON)
    CALL DM_VBCSE(A, K, NW, N, ICOL, B, ITMAX, EPS, IGUSS, L, &
        X, ITER, ICON)
    IF (ICON.ne.0) THEN
        write(*,*) "ERROR: DM_DVBCSE() failed with CODE=",ICON
    ELSE
        write(*,*) "INFO: DM_DVBCSE() : ITER=",ITER
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DM_VTFQE: ASymmetric TFQMR Solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SSL2_VTFQE(K, NW, N, A, ICOL, B, X, EPS, ICON)
    INTEGER :: K, NW, N
    REAL*8 :: A(K,NW)
    INTEGER :: ICOL(K,NW)
    REAL*8 :: B(N), X(N)
    REAL*8 :: EPS
!
    INTEGER :: ITMAX=10000, IGUSS=0
    INTEGER :: ITER, ICON
    REAL*8 :: VW(K*8)
!
    ICON=0
!
#ifdef WITH_SSL2
    write(*,*) "INFO: call DM_VTFQE()"
    CALL DM_VTFQE(A, K, NW, N, ICOL, B, ITMAX, EPS, IGUSS, &
        X, ITER, ICON)
    IF (ICON.ne.0) THEN
        write(*,*) "ERROR: DM_VTFQE() failed with CODE=",ICON
    ELSE
        write(*,*) "INFO: DM_VTFQE() : ITER=",ITER
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
END SUBROUTINE SSL2_VTFQE
