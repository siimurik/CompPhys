!    gfortran -g -fcheck=all dlsoda_ex.f90 odepack.f odepack_interface.f90 odepack_common.f90 odepack_sub1.f odepack_sub2.f odepack_mod.f90 -o dex -lopenblas -std=legacy && ./dex 

program dlsoda_ex
    use iso_c_binding
    use odepack_interface
    use odepack_common
    implicit none
    
    external fex, jdum
    integer, parameter :: dp = kind(0.0d0)
    real(kind=dp), dimension(3) :: atol, y
    integer :: iopt, iout, istate, itask
    integer :: itol, jt, liw, lrw, neq
    integer, dimension(23) :: iwork
    real(kind=dp) :: rtol, t, tout
    real(kind=dp), dimension(70) :: rwork
    type(odepack_common_data), target :: common_data
    
    ! Initialize ALL arrays completely
    iwork = 0
    rwork = 0.0_dp
    
    ! Initialize common_data structure
    common_data%ierr = 0
    ! Add any other common_data initialization here if needed
    
    ! Problem setup
    neq = 3
    y(1) = 1.0_dp
    y(2) = 0.0_dp
    y(3) = 0.0_dp
    t = 0.0_dp
    tout = 0.4_dp
    
    ! Tolerance setup
    itol = 2
    rtol = 1.0d-4
    atol(1) = 1.0d-6
    atol(2) = 1.0d-10
    atol(3) = 1.0d-6
    
    ! Solver options
    itask = 1
    istate = 1
    iopt = 0
    lrw = 70
    liw = 23
    jt = 2
    
    ! Main integration loop
    do iout = 1, 12
        call dlsoda(fex, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, &
                   rwork, lrw, iwork, liw, jdum, jt, common_data)
        
        write(6, 99010) t, y(1), y(2), y(3)
        99010 format(' At t =', d12.4, '   Y =', 3d14.6)
        
        if (istate < 0) then
            write(6, 99020) istate
            99020 format(///' Error halt.. ISTATE =', i3)
            stop 1
        else
            tout = tout * 10.0_dp
        endif
    enddo
    
    ! Final statistics
    write(6, 99030) iwork(11), iwork(12), iwork(13), iwork(19), rwork(15)
    99030 format(/' No. steps =', i4, '  No. f-s =', i4, '  No. J-s =', i4/ &
                 ' Method last used =', i2, '   Last switch was at t =', d12.4)
    
end program dlsoda_ex

! Fixed FEX subroutine with common_data parameter
subroutine fex(neq, t, y, ydot, common_data)
    use odepack_common
    implicit none
    integer, parameter :: dp = kind(0.0d0)
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), intent(in), dimension(3) :: y
    real(dp), intent(out), dimension(3) :: ydot
    type(odepack_common_data), intent(inout) :: common_data
    
    ! Initialize error flag
    common_data%ierr = 0
    
    ! Chemical kinetics equations
    ydot(1) = -0.04_dp * y(1) + 1.0d4 * y(2) * y(3)
    ydot(3) = 3.0d7 * y(2) * y(2)
    ydot(2) = -ydot(1) - ydot(3)
    
end subroutine fex

! Fixed JDUM subroutine with common_data parameter
subroutine jdum(neq, t, y, ml, mu, pd, nrowpd, common_data)
    use odepack_common
    implicit none
    integer, parameter :: dp = kind(0.0d0)
    integer, intent(in) :: neq, ml, mu, nrowpd
    real(dp), intent(in) :: t
    real(dp), intent(in) :: y(neq)
    real(dp), intent(out) :: pd(nrowpd, neq)
    type(odepack_common_data), intent(inout) :: common_data
    
    ! Initialize error flag
    common_data%ierr = 0
    
    ! Dummy Jacobian - not used when jt=2
    ! Initialize pd array to avoid uninitialized memory
    pd = 0.0_dp
    
end subroutine jdum