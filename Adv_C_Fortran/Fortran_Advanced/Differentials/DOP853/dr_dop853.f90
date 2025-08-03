! * * * * * * * * * * * * * * * * * * * * * * * * *
! --- DRIVER FOR DOPRI5 ON VAN DER POL'S EQUATION
! * * * * * * * * * * * * * * * * * * * * * * * * *
! eh dr_dop853 dop853
program main
    ! Parameters
    integer, parameter :: NDGL = 2, NRD = 2
    integer, parameter :: LWORK = 11*NDGL + 8*NRD + 21
    integer, parameter :: LIWORK = NRD + 21
    double precision, parameter :: TOL = 1.0d-9
    
    ! Variables
    integer :: n, iout, itol, ipar, idid, i, j
    double precision :: x, xend, rpar, rtol, atol
    double precision, dimension(NDGL) :: y
    double precision, dimension(LWORK) :: work
    integer, dimension(LIWORK) :: iwork
    
    ! External procedures
    external :: fvpol, solout
    
    ! Problem setup
    n = 2                    ! Dimension of the system
    rpar = 1.0d-3           ! Parameter for Van der Pol equation
    iout = 3                ! Dense output option
    
    ! Initial conditions
    x = 0.0d0
    y(1) = 2.0d0
    y(2) = 0.0d0
    
    ! Integration endpoint
    xend = 2.0d0
    
    ! Error tolerances
    itol = 0                ! Both tolerances are scalars
    rtol = TOL
    atol = TOL
    
    ! Initialize work arrays
    do i = 1, 10
        iwork(i) = 0
        work(i) = 0.0d0
    end do
    
    iwork(5) = n            ! Number of components for dense output
    iwork(4) = 1            ! Stiffness test frequency

    ! --- CALL OF THE SUBROUTINE DOPRI8
    CALL DOP853(N,FVPOL,X,Y,XEND, &
        RTOL,ATOL,ITOL, &
        SOLOUT,IOUT, &
        WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
    
    ! --- PRINT FINAL SOLUTION
    WRITE (6,99) X,Y(1),Y(2)
99  FORMAT(1X,'X =',F5.2,'    Y =',2E18.10)
    
    ! --- PRINT STATISTICS
    WRITE (6,90) TOL
90  FORMAT('       tol=',D8.2)
    WRITE (6,91) (IWORK(J),J=17,20)
91  FORMAT(' fcn=',I5,' step=',I4,' accpt=',I4,' rejct=',I3)

    ! Print completion status
    select case(idid)
        case(1)
            write(*,*) 'Integration completed successfully'
        case(2)
            write(*,*) 'Integration interrupted by solout'
        case(-1)
            write(*,*) 'Input is not consistent'
        case(-2)
            write(*,*) 'Larger NMAX is needed'
        case(-3)
            write(*,*) 'Step size becomes too small'
        case(-4)
            write(*,*) 'Problem is probably stiff'
        case default
            write(*,*) 'Unknown return code:', idid
    end select

END PROGRAM main

SUBROUTINE SOLOUT (NR, XOLD, X, Y, N, CON, ICOMP, ND, &
    RPAR, IPAR, IRTRN, XOUT)

    ! --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
    ! --- BY USING "CONTD8", THE CONTINUOUS COLLOCATION SOLUTION
    implicit none
    
    ! Arguments
    integer, intent(in) :: nr, n, nd
    integer, intent(in) :: icomp(nd), ipar
    double precision, intent(in) :: xold, x, rpar
    double precision, intent(in) :: y(n), con(8*nd)
    integer, intent(inout) :: irtrn
    double precision, intent(inout) :: xout
    
    ! External function
    double precision, external :: contd8

    IF (NR == 1) THEN
        WRITE (6,99) X,Y(1),Y(2),NR-1
        XOUT=0.1D0
    ELSE
        10 CONTINUE
        IF (X >= XOUT) THEN
            WRITE (6,99) XOUT,CONTD8(1,XOUT,CON,ICOMP,ND), &
            CONTD8(2,XOUT,CON,ICOMP,ND),NR-1
            XOUT=XOUT+0.1D0
            GOTO 10
        END IF
    END IF
    99 FORMAT(1X,'X =',F5.2,'    Y =',2E18.10,'    NSTEP =',I4)
    
END SUBROUTINE SOLOUT

!======================================================================
! VAN DER POL EQUATION RIGHT-HAND SIDE
! System: y1' = y2
!         y2' = ((1-y1^2)*y2 - y1)/eps
!======================================================================
subroutine fvpol(n, x, y, f, rpar, ipar)
    implicit none
    
    ! Arguments
    integer, intent(in) :: n, ipar
    double precision, intent(in) :: x, rpar
    double precision, intent(in) :: y(n)
    double precision, intent(out) :: f(n)
    
    ! Local variables
    double precision :: eps
    
    eps = rpar
    f(1) = y(2)
    f(2) = ((1.0d0 - y(1)**2) * y(2) - y(1)) / eps
    
end subroutine fvpol