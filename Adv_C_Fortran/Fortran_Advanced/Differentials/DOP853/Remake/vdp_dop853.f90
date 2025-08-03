! gfortran vdp_dop853.f90 dop853_REMADE.f90 -o vdp
program main
    use dop853_mod
    implicit none
    
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: neq = 2    
    integer :: n
    real(dp) :: t, tstop
    real(dp), dimension(neq) :: y
    real(dp), dimension(1) :: rpar  ! Changed to array
    real(dp), parameter :: TOL = 1.0e-9_dp

    !real(dp), dimension(neq) :: atol, rtol  ! Already fixed from previous error
    !common /tolerances/ atol, rtol
    ! Problem setup
    n = 2         ! Dimension of the system
    rpar(1) = 1.0e-3_dp ! Parameter for Van der Pol equation (now array element)
    
    ! Initial conditions
    y(1) = 2.0_dp    
    y(2) = 0.0_dp  
    t = 0.0_dp    ! Start time
    tstop = 2.0_dp   ! Integration endpoint
    
    ! Set tolerances
    !rtol = TOL
    !atol = TOL

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
    
    write(6, '(A,I0,A,I0)') 'Calculated DOP853 array sizes: LWORK = ', lwork, ', LIWORK = ', liwork
    
end subroutine calculate_dop853_arrays

subroutine integrate_DOP853(neq, y, t, tend, rpar)
    use dop853_mod
    implicit none
    
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: neq
    real(dp), intent(inout) :: t, tend
    real(dp), dimension(neq), intent(inout) :: y
    real(dp), dimension(*), intent(inout) :: rpar  ! Changed to assumed-size array
    integer, dimension(1) :: ipar  ! Changed to array
    
    integer :: lwork, liwork, itol, iout, idid
    integer, allocatable, dimension(:) :: iwork
    real(dp), allocatable, dimension(:) :: work
    real(dp), parameter :: TOL = 1.0e-9_dp  ! Change desired accuray here
    real(dp), dimension(neq) :: atol, rtol

    external fvpol, solout

    ! Error tolerances
    itol = 1    ! Both tolerances are arrays
    iout = 2    ! Dense output
    
    ! User parameters
    ipar = 0

    ! Error message parameter
    idid = 0
    
    call calculate_dop853_arrays(neq, lwork, liwork)
    allocate(work(lwork))
    allocate(iwork(liwork))

    ! Initialize arrays
    work = 0.0_dp
    iwork = 0

    ! Set tolerances
    rtol = TOL
    atol = TOL

    call dop853(neq, fvpol, t, y, tend, rtol, atol, itol, solout, iout, &
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
    
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    ! Arguments
    integer, intent(in) :: n, ipar
    real(dp), intent(in) :: x, rpar
    real(dp), intent(in) :: y(n)
    real(dp), intent(out) :: f(n)
    
    ! Local variables
    real(dp) :: eps

    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) x, ipar
    end if
    
    eps = rpar
    f(1) = y(2)
    f(2) = ((1.0_dp - y(1)**2) * y(2) - y(1)) / eps
    
end subroutine fvpol

! Output routine for simple solution monitoring
subroutine solout(nr, told, t, y, n, con, icomp, nd, &
                 rpar, ipar, irtrn, tout)
    implicit none
    
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer, intent(in) :: nr, n, nd, ipar
    integer, dimension(nd), intent(in) :: icomp
    real(dp), intent(in) :: told, rpar
    real(dp), intent(in) :: t
    real(dp), dimension(n), intent(in) :: y
    real(dp), dimension(8*nd), intent(in) :: con
    integer, intent(inout) :: irtrn
    real(dp), intent(inout) :: tout
    real(dp), save :: next_output = 0.5_dp  ! First output time
    real(dp), parameter :: dt = 0.5_dp      ! Output interval
    
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