!   gfortran ex.f90 odepack.f odepack_interface.f90 odepack_common.f90 ...
!   odepack_sub1.f odepack_sub2.f odepack_mod.f90 -o ex -llapack       ...
!   -lblas -std=legacy
program dlsoda_chemical_kinetics
    use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
    use odepack_interface, only: DINTDY, DROOTS, DSTODA
    use odepack_common
    implicit none
    
    ! Parameters based on the documentation example
    integer, parameter :: neq_size = 3
    integer, parameter :: lrw = 70    ! As specified in documentation
    integer, parameter :: liw = 23    ! As specified in documentation
    
    ! Variables for DLSODA
    external :: fex, jdum
    integer :: neq, itol, itask, istate, iopt, lrw_actual, liw_actual, jt
    double precision :: y(neq_size), t, tout, rtol, atol(neq_size)
    double precision :: rwork(lrw)
    integer :: iwork(liw)
    type(odepack_common_data), target :: common_data
    
    ! Problem-specific variables
    integer :: iout
    
    ! Initialize the chemical kinetics problem from documentation
    ! dy1/dt = -.04*y1 + 1.e4*y2*y3
    ! dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
    ! dy3/dt = 3.e7*y2**2
    
    neq = neq_size
    y(1) = 1.0d0         ! Initial condition: y1(0) = 1
    y(2) = 0.0d0         ! Initial condition: y2(0) = 0
    y(3) = 0.0d0         ! Initial condition: y3(0) = 0
    t = 0.0d0            ! Initial time
    tout = 0.4d0         ! First output time
    
    ! Set tolerances as in documentation example
    itol = 2                ! Vector tolerances
    rtol = 1.0d-4          ! Relative tolerance
    atol(1) = 1.0d-6       ! Absolute tolerance for y1
    atol(2) = 1.0d-10      ! Absolute tolerance for y2 (smaller because y2 has smaller values)
    atol(3) = 1.0d-6       ! Absolute tolerance for y3
    
    ! Set task and options
    itask = 1              ! Normal computation of output values
    istate = 1             ! First call
    iopt = 0               ! No optional inputs
    jt = 2                 ! Jacobian type: internally generated full Jacobian
    
    ! Set work array sizes
    lrw_actual = lrw
    liw_actual = liw
    
    write(*,'(a)') '# Chemical Kinetics Problem Solution using DLSODA'
    write(*,'(a)') '# Based on the example in DLSODA documentation'
    write(*,'(a)') '# t           y1          y2          y3'
    
    ! Main integration loop - 12 output points as in documentation
    do iout = 1, 12
        call dlsoda(fex, neq, y, t, tout, itol, rtol, atol, itask, &
                   istate, iopt, rwork, lrw_actual, iwork, liw_actual, &
                   jdum, jt, common_data)
        
        write(*,'(a,d12.4,a,3d14.6)') ' At t =', t, '   Y =', y(1), y(2), y(3)
        
        if (istate < 0) then
            write(*,'(a)') ''
            write(*,'(a)') ''
            write(*,'(a,i3)') ' Error halt.. ISTATE =', istate
            select case(istate)
            case(-1)
                write(*,'(a)') ' Excess work done on this call (perhaps wrong JT).'
            case(-2)
                write(*,'(a)') ' Excess accuracy requested (tolerances too small).'
            case(-3)
                write(*,'(a)') ' Illegal input detected.'
            case(-4)
                write(*,'(a)') ' Repeated error test failures (check all inputs).'
            case(-5)
                write(*,'(a)') ' Repeated convergence failures (perhaps bad Jacobian).'
            case(-6)
                write(*,'(a)') ' Error weight became zero during problem.'
            case(-7)
                write(*,'(a)') ' Work space insufficient to finish.'
                write(*,'(a,i0)') ' Required RWORK length: ', iwork(17)
                write(*,'(a,i0)') ' Required IWORK length: ', iwork(18)
            end select
            stop
        end if
        
        ! Increase tout by factor of 10 for next output
        tout = tout * 10.0d0
    end do
    
    ! Final statistics as in documentation example
    write(*,'(a,i4,a,i4,a,i4)') &
        ' No. steps =', iwork(11), '  No. f-s =', iwork(12), '  No. J-s =', iwork(13)
    write(*,'(a,i2,a,d12.4)') &
        ' Method last used =', iwork(19), '   Last switch was at t =', rwork(15)
    
end program dlsoda_chemical_kinetics

! Subroutine defining the chemical kinetics equations (FEX from documentation)
subroutine fex(neq, t, y, ydot, common_data)
    use odepack_common
    implicit none
    integer, intent(in) :: neq
    double precision, intent(in) :: t
    double precision, intent(in) :: y(neq)
    double precision, intent(out) :: ydot(neq)
    type(odepack_common_data), intent(inout) :: common_data
    
    ! Chemical kinetics equations from documentation:
    ! dy1/dt = -.04*y1 + 1.e4*y2*y3
    ! dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
    ! dy3/dt = 3.e7*y2**2
    
    ydot(1) = -0.04d0 * y(1) + 1.0d4 * y(2) * y(3)
    ydot(3) = 3.0d7 * y(2) * y(2)
    ydot(2) = -ydot(1) - ydot(3)
    
end subroutine fex

! Dummy Jacobian subroutine (JDUM from documentation)
subroutine jdum(neq, t, y, ml, mu, pd, nrowpd, common_data)
    use odepack_common
    implicit none
    integer, intent(in) :: neq, ml, mu, nrowpd
    double precision, intent(in) :: t
    double precision, intent(in) :: y(neq)
    double precision, intent(out) :: pd(nrowpd, neq)
    type(odepack_common_data), intent(inout) :: common_data
    
    ! This subroutine is not called when jt=2 (internal Jacobian)
    ! but must be present in the call - this is JDUM from documentation
    
end subroutine jdum