!   gfortran arenstrof_orbit_dopri5.f90 dopri5.f90 -o ao
program main
    implicit none
    integer, parameter :: ndgl=4,nrdens=2
    integer :: n
    double precision :: t, tstop
    double precision, dimension(2) :: rpar
    double precision, dimension(ndgl) :: y
    double precision, parameter :: tol = 1.0d-7

    double precision :: atol, rtol
    common /tolerances/ atol, rtol

    ! --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
    RTOL=1.0D-7
    ATOL=RTOL

    ! --- dimension of the system
    n=ndgl

    ! --- initial values and endpoint of integration
    rpar(1)=0.012277471d0
    rpar(2)=1.d0-rpar(1)
    t=0.0d0
    tstop=17.0652165601579625588917206249d0
    
    y(1)=0.994d0
    y(2)=0.0d0
    y(3)=0.0d0
    y(4)=-2.00158510637908252240537862224d0
    

    ! time series saved into 't'; variables saved into 'y' vectors
    call integrate_dopri5(n, nrdens, y, t, tstop, rpar)

end program main

! subroutine to calculate dopri5 array lengths
subroutine calculate_dopri5_arrays(neq, nrdens, lwork, liwork)
    implicit none
    integer, intent(in) :: neq, nrdens
    integer, intent(out) :: lwork, liwork
    
    ! for dense output
    lwork = 8 * neq + 5 * nrdens + 21
    liwork = nrdens + 21
    
    write(6, '(a,i0,a,i0)') 'calculated dopri5 array sizes: lwork = ', lwork, ', liwork = ', liwork
    
end subroutine calculate_dopri5_arrays

subroutine integrate_dopri5(neq, nrdens, y, x, tend, rpar)
    implicit none
    integer, intent(in) :: neq, nrdens
    double precision, intent(inout) :: x, tend
    double precision, dimension(neq), intent(inout) :: y
    double precision, dimension(2), intent(in) :: rpar

    
    integer :: i, j, lwork, liwork, itol, iout, ipar, idid
    integer, allocatable, dimension(:) :: iwork
    double precision, allocatable, dimension(:) :: work
    double precision :: atol, rtol
    common /tolerances/ atol, rtol

    external faren, solout

    ! error tolerances
    itol = 0    ! both tolerances are scalars
    
    ! --- output routine (and dense output) is used during integration
    iout=2
    
    ! user parameters
    ipar = 0

    ! error message parameter
    idid = 0
    
    call calculate_dopri5_arrays(neq, nrdens, lwork, liwork)
    allocate(work(lwork))
    allocate(iwork(liwork))

    ! --- DEFAULT VALUES FOR PARAMETERS
    DO I=1,20
        IWORK(I)=0
        WORK(I)=0.D0
    END DO

    ! in integrate_dop853 subroutine, replace initialization with:
    ! --- dense output is used for the two position coordinates 1 and 2
    iwork(5)=nrdens
    iwork(21)=1
    iwork(22)=2

    ! --- call of the subroutine dopri5
    call dopri5(neq, faren, x, y, tend, &
                rtol, atol, itol, solout,&
                iout, work,lwork,iwork,&
                liwork,rpar,ipar,idid)

    ! check integration status
    select case (idid)
    case (1)
        write(6, '(/a)') 'integration completed successfully.'
    case (2)
        write(6, '(/a)') 'integration completed - interrupted by solout.'
    case (-1)
        write(6, '(/a)') 'error: input is not consistent.'
    case (-2)
        write(6, '(/a)') 'error: larger nmax is needed.'
    case (-3)
        write(6, '(/a)') 'error: step size becomes too small.'
    case (-4)
        write(6, '(/a)') 'error: problem is probably stiff (use stiff solver).'
    case default
        write(6, '(/a,i0)') 'integration finished with status: ', idid
    end select

    ! --- print final solution
    write (6,99) y(1),y(2)
    99 format(1x,'x = xend     y =',2e18.10)
! --- print statistics
    write (6,91) rtol,(iwork(j),j=17,20)
    91 format('     tol=',d8.2,'   fcn=',i5,' step=',i4, &
    ' accpt=',i4,' rejct=',i3)

    deallocate(work, iwork)

end subroutine integrate_dopri5

!======================================================================
!
!======================================================================
subroutine faren(n,x,y,f,rpar,ipar)
    ! --- arenstorf orbit
    implicit none
    integer, intent(in) :: n, ipar
    double precision, intent(in) :: x
    double precision, dimension(2), intent(in) :: rpar
    double precision, dimension(n), intent(in) :: y
    double precision, dimension(n), intent(out) :: f
    
    double precision :: amu, amup, r1, r2

    ! silence unused parameter warnings
    if (.false.) then
        write(*,*) x, ipar
    end if

    amu=rpar(1)
    amup=rpar(2)
    
    f(1)=y(3)
    f(2)=y(4)
    
    r1=(y(1)+amu)**2+y(2)**2
    r1=r1*sqrt(r1)
    r2=(y(1)-amup)**2+y(2)**2
    r2=r2*sqrt(r2)
    
    f(3)=y(1)+2*y(4)-amup*(y(1)+amu)/r1-amu*(y(1)-amup)/r2
    f(4)=y(2)-2*y(3)-amup*y(2)/r1-amu*y(2)/r2

end subroutine faren

! output routine for simple solution monitoring
subroutine solout(nr, xold, x, y, n, con, icomp, nd, rpar, ipar, irtrn)
    implicit none
    integer, intent(in) :: nr, n, nd, ipar
    double precision, intent(in) :: xold, x
    double precision, dimension(n), intent(in) :: y
    double precision, dimension(5*nd), intent(in) :: con
    integer, dimension(nd), intent(in) :: icomp
    double precision, dimension(2), intent(in) :: rpar
    integer, intent(inout) :: irtrn
    double precision, parameter :: dt = 2.0D0

    ! Persistent variable to track next output point
    double precision :: xout
    save xout
    data xout /0.0d0/

    ! External interpolation function
    double precision, external :: contd5

    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) xold, rpar, ipar, irtrn
    endif

    if (nr == 1) then
        write (6,99) x, y(1), y(2), nr-1
        xout = x + dt
    else
        do while (x >= xout)
            write (6,99) xout, contd5(1, xout, con, icomp, nd), &
                          contd5(2, xout, con, icomp, nd), nr-1
            xout = xout + dt
        end do
    end if

99  format(1x, 'x =', f6.2, '    y =', 2e18.10, '    nstep =', i4)

end subroutine solout
