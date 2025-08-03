! Compile with:
! gfortran predator_prey_dlsoda.f90 odepack.f odepack_interface.f90 odepack_common.f90 odepack_sub1.f odepack_sub2.f odepack_mod.f90 -o predator_prey -lopenblas -std=legacy

program dlsoda_predator_prey
    use iso_c_binding
    use odepack_interface
    use odepack_common
    implicit none
    
    external f2, jdum
    integer, parameter :: dp = kind(0.0d0)
    integer, parameter :: neq = 2
    real(kind=dp), dimension(neq) :: atol, y
    integer :: iopt, iout, istate, itask
    integer :: itol, jt, liw, lrw
    integer, dimension(30) :: iwork
    real(kind=dp) :: rtol, t, tout
    real(kind=dp), dimension(100) :: rwork
    real(kind=dp) :: dt
    integer :: i_step, n_step
    type(odepack_common_data), target :: common_data
    
    ! Initialize arrays
    iwork = 0
    rwork = 0.0_dp
    
    ! Initialize common_data structure
    common_data%ierr = 0
    
    ! Initial conditions for predator-prey system
    ! y(1) = prey population
    ! y(2) = predator population
    y(1) = 1.0D0    ! Initial prey population
    y(2) = 1.0D0    ! Initial predator population
    
    t = 0.0_dp
    tout = 0.0_dp
    n_step = 1000
    dt = 0.01D0     ! Time step
    
    ! Tolerance setup
    itol = 1                ! Scalar tolerances
    rtol = 1.0d-8
    atol = 1.0d-10
    
    ! Solver options
    itask = 1
    istate = 1
    iopt = 0                ! No optional inputs
    lrw = size(rwork)
    liw = size(iwork)
    jt = 2                  ! Internal generated Jacobian
    
    write(*,'(a)') "PREDATOR-PREY SYSTEM USING DLSODA"
    write(*,'(a)') "dy1/dt = y1*(alpha - delta*y2)"
    write(*,'(a)') "dy2/dt = y2*(y1*gamma - beta)"
    write(*,'(a)') "alpha=2.0, beta=1.0, delta=1.0, gamma=3.0"
    write(*,'(a)') ""
    write(*,'(a)') "   TIME         PREY      PREDATOR"
    write(*,'(a)') "----------------------------------------"
    
    ! Open output file
    open(unit=11, file='predator_prey.dat')
    write(11, '(3g16.8)') t, y(1), y(2)
    write(*, '(f8.4, 2f12.6)') t, y(1), y(2)
    
    ! Integration loop
    do i_step = 1, n_step
        tout = t + dt
        
        call dlsoda(f2, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, &
                   rwork, lrw, iwork, liw, jdum, jt, common_data)
        
        if (istate < 0) then
            write(*,*) 'Integration failed with istate =', istate
            exit
        endif
        
        ! Output every 10th step for readability
        if (mod(i_step, 10) == 0) then
            write(*, '(f8.4, 2f12.6)') t, y(1), y(2)
        endif
        
        write(11, '(3g16.8)') t, y(1), y(2)
    end do
    
    close(11)
    
    ! Final statistics
    write(*,'(a)') ""
    write(*,'(a,i6)') 'Total integration steps: ', iwork(11)
    write(*,'(a,i6)') 'Function evaluations: ', iwork(12)
    write(*,'(a,i6)') 'Jacobian evaluations: ', iwork(13)
    
end program dlsoda_predator_prey

subroutine f2(neq, t, y, ydot, common_data)
    use odepack_common
    implicit none
    integer, parameter :: dp = kind(0.0d0)
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), intent(in), dimension(neq) :: y
    real(dp), intent(out), dimension(neq) :: ydot
    type(odepack_common_data), intent(inout) :: common_data
    
    real(dp) :: alpha, beta, delta, gamma

    ! Predator-prey parameters
    alpha = 2.0D0   ! prey population growth parameter
    beta = 1.0D0    ! predator population extinction parameter
    delta = 1.0D0   ! species interaction parameter 1
    gamma = 3.0D0   ! species interaction parameter 2
    
    ! Predator-prey equations (Lotka-Volterra)
    ydot(1) = y(1) * (alpha - delta * y(2))  ! prey equation
    ydot(2) = y(2) * (y(1) * gamma - beta)   ! predator equation
    
    common_data%ierr = 0
    
    ! Dummy use of t to avoid compiler warnings
    if (t /= t) then
        write(*,*) 'Time variable is NaN'
    endif
    
end subroutine f2

subroutine jdum(neq, t, y, ml, mu, pd, nrowpd, common_data)
    use odepack_common
    implicit none
    integer, parameter :: dp = kind(0.0d0)
    integer, intent(in) :: neq, ml, mu, nrowpd
    real(dp), intent(in) :: t
    real(dp), intent(in) :: y(neq)
    real(dp), intent(out) :: pd(nrowpd, neq)
    type(odepack_common_data), intent(inout) :: common_data
    
    common_data%ierr = 0
    pd = 0.0_dp
    
end subroutine jdum