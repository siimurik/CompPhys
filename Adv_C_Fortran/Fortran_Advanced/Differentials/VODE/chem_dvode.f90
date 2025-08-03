! Compile with:
! $ gfortran chem_dvode.f90 vode.f dgesl.f dgbsl.f dgbfa.f dgefa.f -o chem -lopenblas -std=legacy
! $ gfortran chem_dvode.f90 VODEPACK.f90 -o chem -lopenblas

program dvode_ex
    implicit none
    external fex, jex
    integer, parameter :: dp = kind(0.0d0)
    real(kind=dp), dimension(3) :: atol, y
    integer :: neq, mf, liw, lrw
    integer, allocatable, dimension(:) :: iwork
    real(kind=dp), allocatable, dimension(:) :: rwork
    real(kind=dp) :: rtol, t, tout
    real(kind=dp) :: rpar
    integer :: ipar

    ! Initialize parameters
    rpar = 0.0_dp
    ipar = 0

    ! Problem setup
    neq = 3
    mf = 21
    
    ! Calculate required array sizes automatically
    call calculate_vode_arrays(neq, mf, lrw, liw)
    
    ! Allocate work arrays
    allocate(rwork(lrw))
    allocate(iwork(liw))
    
    ! Initialize arrays
    iwork = 0
    rwork = 0.0_dp

    ! Initial conditions
    y(1) = 1.0_dp
    y(2) = 0.0_dp
    y(3) = 0.0_dp
    t = 0.0_dp
    tout = 0.4_dp

    ! Tolerance setup
    rtol = 1.0d-4
    atol(1) = 1.0d-8
    atol(2) = 1.0d-14
    atol(3) = 1.0d-6

    ! Run the main integration
    call run_dvode_integration(fex, jex, neq, y, t, tout, rtol, atol, &
                              rwork, lrw, iwork, liw, mf, rpar, ipar)

    ! Print final statistics
    call print_dvode_statistics(iwork, liw)

    ! Clean up
    deallocate(rwork)
    deallocate(iwork)

end program dvode_ex

! Subroutine to automatically calculate VODE array lengths
subroutine calculate_vode_arrays(neq, mf, lrw, liw)
    implicit none
    integer, intent(in) :: neq, mf
    integer, intent(out) :: lrw, liw
    integer, parameter :: ml = 5, mu = 5  ! Default bandwidths for banded case
    
    ! Calculate LRW based on method flag
    select case (mf)
    case (10)
        ! Nonstiff (Adams) method
        lrw = 20 + 16 * neq
    case (21, 22)
        ! Stiff (BDF) method with full Jacobian
        lrw = 22 + 9 * neq + 2 * neq**2
    case (24, 25)
        ! Stiff (BDF) method with banded Jacobian
        ! For this example, assuming ML=MU=5 (adjust as needed)
        ! In practice, you would pass ML and MU as parameters
        lrw = 22 + 11 * neq + (3 * ml + 2 * mu) * neq
    case default
        ! Default to full Jacobian case
        lrw = 22 + 9 * neq + 2 * neq**2
    end select
    
    ! Calculate LIW based on method flag
    select case (mf)
    case (10)
        ! Nonstiff method
        liw = 30
    case (21, 22, 24, 25)
        ! Stiff methods
        liw = 30 + neq
    case default
        ! Default to stiff case
        liw = 30 + neq
    end select
    
    write(6, '(A,I0,A,I0)') 'Calculated array sizes: LRW = ', lrw, ', LIW = ', liw
    
end subroutine calculate_vode_arrays

! Main integration subroutine
subroutine run_dvode_integration(f, jac, neq, y, t, tout, rtol, atol, &
                                rwork, lrw, iwork, liw, mf, rpar, ipar)
    implicit none
    external f, jac
    integer, parameter :: dp = kind(0.0d0)
    integer, intent(in) :: neq, lrw, liw, mf, ipar
    real(kind=dp), intent(in) :: rtol, rpar
    real(kind=dp), intent(inout) :: t, tout
    real(kind=dp), dimension(neq), intent(inout) :: y
    real(kind=dp), dimension(neq), intent(in) :: atol
    real(kind=dp), dimension(lrw), intent(inout) :: rwork
    integer, dimension(liw), intent(inout) :: iwork
    
    ! Local variables
    integer :: itol, itask, istate, iopt, iout

    ! Solver options
    itol = 2      ! atol is an array
    itask = 1     ! normal computation
    istate = 1    ! first call
    iopt = 0      ! no optional input

    write(6, '(A)') 'Starting DVODE integration...'
    write(6, '(A)') '      Time          Y(1)          Y(2)          Y(3)'
    write(6, '(A)') '----------------------------------------------------'

    ! Main integration loop
    do iout = 1, 12
        call dvode(f, neq, y, t, tout, itol, rtol, atol, itask, istate, &
                   iopt, rwork, lrw, iwork, liw, jac, mf, rpar, ipar)
        
        write(6, 20) t, y(1), y(2), y(3)
20      format(' ', d12.4, 3d14.6)
        
        if (istate < 0) then
            write(6, 90) istate
90          format(///' Error halt: ISTATE =', i3)
            select case (istate)
            case (-1)
                write(6, '(A)') 'Excess work done on this call. (Perhaps wrong MF.)'
            case (-2)
                write(6, '(A)') 'Excess accuracy requested. (Tolerances too small.)'
            case (-3)
                write(6, '(A)') 'Illegal input detected. (See printed message.)'
            case (-4)
                write(6, '(A)') 'Repeated error test failures. (Check all input.)'
            case (-5)
                write(6, '(A)') 'Repeated convergence failures. (Perhaps bad Jacobian or wrong MF/tolerances.)'
            case (-6)
                write(6, '(A)') 'Error weight became zero during problem. (Solution component vanished.)'
            end select
            stop 1
        else
            tout = tout * 10.0_dp
        endif
    enddo
    
    write(6, '(A)') 'Integration completed successfully.'

end subroutine run_dvode_integration

! Subroutine to print final statistics
subroutine print_dvode_statistics(iwork, liw)
    implicit none
    integer, intent(in) :: liw
    integer, dimension(liw), intent(in) :: iwork
    
    write(6, '(/A)') 'Final Integration Statistics:'
    write(6, '(A)') '============================'
    write(6, 60) iwork(11), iwork(12), iwork(13), iwork(19), &
                 iwork(20), iwork(21), iwork(22)
60  format(' No. steps =', i4, / &
           ' No. f-s =', i4, / &
           ' No. J-s =', i4, / &
           ' No. LU-s =', i4, / &
           ' No. nonlinear iterations =', i4, / &
           ' No. nonlinear convergence failures =', i4, / &
           ' No. error test failures =', i4, /)
    
    ! Additional useful statistics
    if (liw >= 22) then
        write(6, '(A,I0)') ' Method order used: ', iwork(14)
        write(6, '(A,I0)') ' Length of real work array used: ', iwork(17)
        write(6, '(A,I0)') ' Length of integer work array used: ', iwork(18)
    endif

end subroutine print_dvode_statistics

! Right-hand side function
subroutine fex(neq, t, y, ydot, rpar, ipar)
    implicit none
    integer, parameter :: dp = kind(0.0d0)
    integer, intent(in) :: neq, ipar
    real(dp), intent(in) :: t, rpar
    real(dp), intent(in), dimension(neq) :: y
    real(dp), intent(out), dimension(neq) :: ydot
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) t, rpar, ipar  ! Never executed
    endif
    
    ydot(1) = -0.04_dp * y(1) + 1.0d4 * y(2) * y(3)
    ydot(3) = 3.0d7 * y(2) * y(2)
    ydot(2) = -ydot(1) - ydot(3)
    
end subroutine fex

! Jacobian function
subroutine jex(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
    implicit none
    integer, parameter :: dp = kind(0.0d0)
    integer, intent(in) :: neq, ml, mu, nrpd, ipar
    real(dp), intent(in) :: t, rpar
    real(dp), intent(in), dimension(neq) :: y
    real(dp), intent(out), dimension(nrpd, neq) :: pd
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) t, ml, mu, rpar, ipar  ! Never executed
    endif
    
    pd(1,1) = -0.04_dp
    pd(1,2) = 1.0d4 * y(3)
    pd(1,3) = 1.0d4 * y(2)
    pd(2,1) = 0.04_dp
    pd(2,3) = -pd(1,3)
    pd(3,2) = 6.0d7 * y(2)
    pd(2,2) = -pd(1,2) - pd(3,2)
    
end subroutine jex