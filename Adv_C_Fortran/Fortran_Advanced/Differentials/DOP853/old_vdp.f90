!   gfortran old_vdp.f90 dop853.f90 -o old
!   gfortran old_vdp.f90 dop853.f -o old -std=legacy

! No errors solution with a bit of cheating
!   gfortran -O2 -finit-local-zero old_vdp.f90 dop853.f90 -o old 
program main
    implicit none
    integer, parameter :: neq = 2    
    integer :: n
    double precision :: t, tstop, rpar
    double precision, dimension(neq) :: y
    double precision, parameter :: TOL = 1.0d-9

    double precision :: atol, rtol
    common /tolerances/ atol, rtol

    ! Problem setup
    n = 2         ! Dimension of the system
    rpar = 1.0d-3 ! Parameter for Van der Pol equation
    
    ! Initial conditions
    y(1) = 2.0d0    
    y(2) = 0.0d0  
    t = 0.0d0    ! Start time
    tstop = 2.0d0   ! Integration endpoint
    
    ! Shared between integrate_DOP853() via common operators
    rtol = TOL  ! relative tolerance
    atol = TOL  ! absolute tolerance
    ! Before calling DOP853, add:
    if (rtol <= 0.0D0) rtol = 1.0D-9  ! Default relative tol
    if (atol <= 0.0D0) atol = 1.0D-9  ! Default absolute tol

    ! Time series saved into 't'; variables saved into 'y' vectors
    call integrate_DOP853(neq, y, t, tstop, rpar)

end program main

! Subroutine to calculate DOP853 array lengths
subroutine calculate_dop853_arrays(neq, lwork, liwork)
    implicit none
    integer, intent(in) :: neq
    integer, intent(out) :: lwork, liwork
    
    ! For dense output with all components
    lwork = 11 * neq + 8 * neq + 21
    liwork = neq + 21
    
    write(6, '(A,I0,A,I0)') 'Calculated DOP853 array sizes: LWORK = ', &
        lwork, ', LIWORK = ', liwork
    
end subroutine calculate_dop853_arrays

subroutine integrate_DOP853(neq, y, t, tend, rpar)
    implicit none
    integer, intent(in) :: neq
    double precision, intent(inout) :: t, tend
    double precision, dimension(neq), intent(inout) :: y
    double precision, intent(inout) :: rpar
    
    integer :: i, lwork, liwork, itol, iout, ipar, idid
    integer, allocatable, dimension(:) :: iwork
    double precision, allocatable, dimension(:) :: work
    double precision :: atol, rtol
    common /tolerances/ atol, rtol

    external fvpol, solout

    ! Error tolerances
    itol = 0    ! Both tolerances are scalars
    iout = 2    ! \
    
    ! User parameters
    ipar = 0

    ! Error message parameter
    idid = 0
    
    call calculate_dop853_arrays(neq, lwork, liwork)
    allocate(work(lwork))
    allocate(iwork(liwork))

    ! In integrate_DOP853 subroutine, replace initialization with:
    ! Initialize work arrays
    do i = 1, 10
        iwork(i) = 0
        work(i) = 0.0d0
    end do

    call DOP853(neq, fvpol, t, y, tend, rtol, atol, itol, solout, iout, &
                work, lwork, iwork, liwork, rpar, ipar, idid)

    ! Check integration status
    select case (idid)
    case (1)
        write(6, '(/A)') 'Integration completed successfully.'
    case (2)
        write(6, '(/A)') 'Integration completed - interrupted by solout.'
    case (-1)
        write(6, '(/A)') 'Error: Input is not consistent.'
    case (-2)
        write(6, '(/A)') 'Error: Larger NMAX is needed.'
    case (-3)
        write(6, '(/A)') 'Error: Step size becomes too small.'
    case (-4)
        write(6, '(/A)') 'Error: Problem is probably stiff (use stiff solver).'
    case default
        write(6, '(/A,I0)') 'Integration finished with status: ', idid
    end select

    deallocate(work, iwork)

end subroutine integrate_DOP853

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

    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) x, ipar
    end if
    
    eps = rpar
    f(1) = y(2)
    f(2) = ((1.0d0 - y(1)**2) * y(2) - y(1)) / eps
    
end subroutine fvpol

! Output routine for simple solution monitoring
subroutine solout(nr, told, t, y, n, con, icomp, nd, &
                 rpar, ipar, irtrn, tout)
    implicit none
    integer, intent(in) :: nr, n, nd, ipar
    integer, dimension(nd), intent(in) :: icomp
    double precision, intent(in) :: told, rpar
    double precision, intent(in) :: t
    double precision, dimension(n), intent(in) :: y
    double precision, dimension(8*nd), intent(in) :: con
    integer, intent(inout) :: irtrn
    double precision, intent(inout) :: tout
    double precision, save :: next_output = 0.5D0  ! First output time
    double precision, parameter :: dt = 0.5D0      ! Output interval
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) told, rpar, ipar, irtrn, icomp(1), con(1), nd, tout
    endif
    
    if (nr == 1) then
        ! Print initial condition
        write(6, '(1X,F10.6,2F14.6,I8)') t, y(1), y(2), nr-1
    else
        ! Print when we pass output times
        do while (told < next_output .and. next_output <= t)
            write(6, '(1X,F10.6,2F14.6,I8)') next_output, y(1), y(2), nr-1
            next_output = next_output + dt
        end do
    endif
    
end subroutine solout