!==============================================================================
! DRIVER FOR ROSENBROCK CODE RODAS ON VAN DER POL EQUATION
!
! Compilation commands:
!   $ gfortran -O2 dr_rodas.f90 rodas.f90 decsol.f90 dc_decsol.f90 -o decrod
!   $ gfortran -O2 dr_rodas.f90 rodas.f90 lapack.f90 lapackc.f90 dc_lapack.f90 -o laprod
!
! NOTE: Even with second level compiler optimization, the 'decsol' SUBROUTINES
!       are faster than using 'lapack' SUBROUTINES.
!==============================================================================
PROGRAM MAIN
    IMPLICIT REAL(KIND=8) (A-H,O-Z)
    
    ! --- PARAMETERS FOR RODAS (FULL JACOBIAN)
    PARAMETER (ND = 2, LWORK = 2*ND*ND + 14*ND + 20, LIWORK = ND + 20)
    
    ! --- DECLARATIONS
    DIMENSION Y(ND), WORK(LWORK), IWORK(LIWORK)
    EXTERNAL FVPOL, JVPOL, SOLOUT
    
    ! --- PROBLEM PARAMETERS
    RPAR = 1.D-6       ! Van der Pol stiffness parameter
    N    = 2              ! Dimension of the system
    
    ! --- SOLVER CONFIGURATION
    IFCN = 0           ! Problem is autonomous
    IJAC = 1           ! Compute Jacobian analytically
    MLJAC = N          ! Jacobian is a full matrix
    IMAS = 0           ! Differential equation is in explicit form
    IOUT = 1           ! Use output routine during integration
    
    ! --- INITIAL CONDITIONS
    X = 0.0D0
    Y(1) = 2.0D0
    Y(2) = -0.66D0
    
    ! --- INTEGRATION PARAMETERS
    XEND = 2.0D0        ! Endpoint of integration
    RTOL = 1.0D-10      ! Relative tolerance
    ATOL = 1.0D-10      ! Absolute tolerance
    ITOL = 0
    H = 1.0D-6          ! Initial step size
    
    ! --- INITIALIZE WORK ARRAYS
    DO I = 1, 20
        IWORK(I) = 0
        WORK(I) = 0.D0
    END DO
    
    ! --- CALL THE SOLVER
    CALL RODAS(N, FVPOL, IFCN, X, Y, XEND, H, &
               RTOL, ATOL, ITOL, &
               JVPOL, IJAC, MLJAC, MUJAC, FVPOL, IDFX, &
               FVPOL, IMAS, MLMAS, MUMAS, &
               SOLOUT, IOUT, &
               WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)
    
    ! --- OUTPUT RESULTS
!    WRITE (6, 99) X, Y(1), Y(2)
!99  FORMAT(1X, 'X =', F6.3, '    Y =', 2E18.10)
    
    ! --- PRINT STATISTICS
    WRITE (6, 90) RTOL
90  FORMAT('       rtol=', D8.2)
    
    WRITE (6, 91) (IWORK(J), J=14,20)
91  FORMAT(' fcn=', I7, ' jac=', I6, ' step=', I6, &
              ' accpt=', I6, ' rejct=', I4, ' dec=', I6, &
              ' sol=', I7)
    
END PROGRAM MAIN

!==============================================================================
! OUTPUT ROUTINE - PRINTS SOLUTION AT REGULAR INTERVALS
! Fixed version that handles all intermediate points properly
!==============================================================================
SUBROUTINE SOLOUT(NR, XOLD, X, Y, CONT, LRC, N, RPAR, IPAR, IRTRN)
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION Y(N), CONT(LRC)
    COMMON /INTERN/ XOUT
    
    !DT = 1.D-03  ! Step size
    DT = 0.2D0
    IF (NR == 1) THEN
        ! First call - print initial condition
        WRITE (6, 99) X, Y(1), Y(2), NR-1
        XOUT = DT
    ELSE
        ! Print all intermediate points in the current step
        DO WHILE (XOUT <= X)
            Y1 = CONTRO(1, XOUT, CONT, LRC)
            Y2 = CONTRO(2, XOUT, CONT, LRC)
            WRITE (6, 99) XOUT, Y1, Y2, NR-1
            XOUT = XOUT + DT
        END DO
    END IF
    
99  FORMAT(1X, 'X =', F6.3, '    Y =', 2E18.10, '    NSTEP =', I6)

END SUBROUTINE SOLOUT


!==============================================================================
! RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
!==============================================================================
SUBROUTINE FVPOL(N, X, Y, F, RPAR, IPAR)
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION Y(N), F(N)
    
    F(1) = Y(2)
    F(2) = ((1 - Y(1)**2) * Y(2) - Y(1)) / RPAR
    
    RETURN
END SUBROUTINE FVPOL


!==============================================================================
! JACOBIAN OF VAN DER POL'S EQUATION
!==============================================================================
SUBROUTINE JVPOL(N, X, Y, DFY, LDFY, RPAR, IPAR)
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION Y(N), DFY(LDFY, N)
    
    DFY(1,1) = 0.0D0
    DFY(1,2) = 1.0D0
    DFY(2,1) = (-2.0D0 * Y(1) * Y(2) - 1.0D0) / RPAR
    DFY(2,2) = (1.0D0 - Y(1)**2) / RPAR
    
    RETURN
END SUBROUTINE JVPOL