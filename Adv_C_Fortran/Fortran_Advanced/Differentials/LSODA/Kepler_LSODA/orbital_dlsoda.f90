! Compile with:
! gfortran orbital_dlsoda.f90 odepack.f odepack_interface.f90 odepack_common.f90 odepack_sub1.f odepack_sub2.f odepack_mod.f90 -o orbital -lopenblas -std=legacy

program dlsoda_orbital
    use iso_c_binding
    use odepack_interface
    use odepack_common
    implicit none
    
    external orbital_equations, jac_dummy_silent
    integer, parameter :: dp = kind(0.0d0)
    integer, parameter :: neq = 4
    double precision, dimension(neq) :: atol, y
    integer :: iopt, istate, itask
    integer :: itol, jt, liw, lrw
    integer, dimension(30) :: iwork
    double precision :: rtol, t, tout
    double precision, dimension(100) :: rwork
    double precision :: v_0, theta, dt
    integer :: i_step, n_step
    type(odepack_common_data), target :: common_data
    double precision :: earth_radius, altitude, initial_radius
    
    ! Initialize arrays
    iwork = 0
    rwork = 0.0D0
    
    ! Initialize common_data structure
    common_data%ierr = 0
    
    ! Orbital parameters
    earth_radius = 6371.0D0      ! Earth radius in km
    altitude = 200.0D0           ! Altitude above Earth surface in km
    initial_radius = earth_radius + altitude  ! Total distance from Earth center
    v_0 = 7.8D0                  ! Initial velocity in km/s
    theta = 90.7D0                ! Launch angle in degrees
    
    ! Initial conditions - object starts at (r, 0) with velocity at angle theta
    y(1) = initial_radius        ! x position (km)
    y(2) = 0.0D0                 ! y position (km)
    y(3) = v_0 * cosd(theta)     ! x velocity (km/s)
    y(4) = v_0 * sind(theta)     ! y velocity (km/s)
    
    t = 0.0D0
    tout = 0.0D0
    n_step = 10000               ! More steps for full orbit
    dt = 1.0D0                   ! 1 second time step
    
    ! Tolerance setup
    itol = 1                     ! Scalar tolerances
    rtol = 1.0d-8
    atol = 1.0d-10
    
    ! Solver options
    itask = 1
    istate = 1
    iopt = 0                     ! No optional inputs
    lrw = size(rwork)
    liw = size(iwork)
    jt = 2                       ! Internal generated Jacobian
    
    write(*,'(a)') "ORBITAL MECHANICS USING DLSODA"
    write(*,'(a)') "Object orbiting Earth in elliptical orbit"
    write(*,'(a,f6.1,a)') "Initial altitude: ", altitude, " km"
    write(*,'(a,f4.1,a)') "Initial velocity: ", v_0, " km/s"
    write(*,'(a,f4.1,a)') "Launch angle: ", theta, " degrees"
    write(*,'(a)') ""
    write(*,'(a)') "   TIME(s)      X(km)       Y(km)       VX(km/s)    VY(km/s)     R(km)"
    write(*,'(a)') "------------------------------------------------------------------------"
    
    ! Open output file
    open(unit=11, file='orbit.dat')
    write(11, '(6g16.8)') t, y(1), y(2), y(3), y(4), sqrt(y(1)**2 + y(2)**2)
    write(*, '(f8.1, 4f12.3, f10.1)') t, y(1), y(2), y(3), y(4), sqrt(y(1)**2 + y(2)**2)
    
    ! Integration loop
    do i_step = 1, n_step
        tout = t + dt
        
        call dlsoda(orbital_equations, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, &
                   rwork, lrw, iwork, liw, jac_dummy_silent, jt, common_data)
        
        if (istate < 0) then
            write(*,*) 'Integration failed with istate =', istate
            write(*,*) 'Last successful time:', rwork(15)
            exit
        endif
        
        ! Output every 200th step for readability (every 200 seconds)
        if (mod(i_step, 200) == 0) then
            write(*, '(f8.1, 4f12.3, f10.1)') t, y(1), y(2), y(3), y(4), sqrt(y(1)**2 + y(2)**2)
        endif
        
        write(11, '(6g16.8)') t, y(1), y(2), y(3), y(4), sqrt(y(1)**2 + y(2)**2)
        
        ! Stop if object crashes into Earth
        if (sqrt(y(1)**2 + y(2)**2) < earth_radius) then
            write(*,*) 'Object crashed into Earth at t =', t
            exit
        endif
    end do
    
    close(11)
    
    ! Final statistics
    write(*,'(a)') ""
    write(*,'(a,i6)') 'Total integration steps: ', iwork(11)
    write(*,'(a,i6)') 'Function evaluations: ', iwork(12)
    write(*,'(a,i6)') 'Jacobian evaluations: ', iwork(13)
    write(*,'(a,f8.1,a)') 'Final time: ', t, ' seconds'
    
end program dlsoda_orbital


subroutine orbital_equations(neq, t, y, ydot, common_data)
    use odepack_common
    implicit none
    integer, parameter :: dp = kind(0.0d0)
    integer, intent(in) :: neq
    double precision, intent(in) :: t
    double precision, intent(in), dimension(neq) :: y
    double precision, intent(out), dimension(neq) :: ydot
    type(odepack_common_data), intent(inout) :: common_data
    
    double precision :: GM, r, r3
    
    ! Earth's gravitational parameter GM in km³/s²
    ! GM = G * M_earth where G = 6.67430e-11 m³/kg/s² and M_earth = 5.972e24 kg
    ! Converting to km³/s²: GM = 3.986004418e14 m³/s² = 3.986004418e5 km³/s²
    GM = 3.986004418D5       ! km³/s²
    !GM = 6.67430D-11 * 5.972D24
    
    ! Calculate distance from Earth center
    r = sqrt(y(1)**2 + y(2)**2)
    !r = 6379.D0
    r3 = r**3
    
    ! Orbital mechanics equations
    ydot(1) = y(3)                    ! dx/dt = vx
    ydot(2) = y(4)                    ! dy/dt = vy
    ydot(3) = -GM * y(1) / r3         ! dvx/dt = -GM*x/r³
    ydot(4) = -GM * y(2) / r3         ! dvy/dt = -GM*y/r³
    
    common_data%ierr = 0
    
    ! Dummy use of t to avoid compiler warnings
    if (t /= t) then
        write(*,*) 'Time variable is NaN'
    endif
    
end subroutine orbital_equations

!subroutine jdum(neq, t, y, ml, mu, pd, nrowpd, common_data)
!    use odepack_common
!    implicit none
!    integer, parameter :: dp = kind(0.0d0)
!    integer, intent(in) :: neq, ml, mu, nrowpd
!    double precision, intent(in) :: t
!    double precision, intent(in) :: y(neq)
!    double precision, intent(out) :: pd(nrowpd, neq)
!    type(odepack_common_data), intent(inout) :: common_data
!    
!    common_data%ierr = 0
!    pd = 0.0D0
!    
!end subroutine jdum

!-----------------------------------------------------------------------
! For real/double precision scalars
!-----------------------------------------------------------------------
subroutine r8_fake_use(x)
    implicit none
    double precision, intent(in) :: x
    
    ! Do absolutely nothing with x, but compiler thinks it's used
    if (x .lt. -huge(x)) then
        continue  ! This will never execute but compiler sees x is "used"
    end if
end subroutine r8_fake_use

!-----------------------------------------------------------------------
! For real/double precision arrays (1D)
!-----------------------------------------------------------------------
subroutine r8_array_fake_use(x)
    implicit none
    double precision, intent(in) :: x(:)
    
    ! Use array in a way that will never execute
    if (size(x) .lt. 0) then
        continue  ! Never executes since size() is always >= 0
    end if
end subroutine r8_array_fake_use

!-----------------------------------------------------------------------
! For integers
!-----------------------------------------------------------------------
subroutine int_fake_use(i)
    implicit none
    integer, intent(in) :: i
    
    ! Do nothing with i but compiler thinks it's used
    if (i .lt. -huge(i)) then
        continue  ! Never executes
    end if
end subroutine int_fake_use

!-----------------------------------------------------------------------
! 'Silent' Version of the Jacobian subroutine:
!-----------------------------------------------------------------------
subroutine jac_dummy_silent(neq, t, y, ml, mu, pd, nrowpd, common_data)
    use odepack_common
    implicit none
    integer, intent(in) :: neq, ml, mu, nrowpd
    double precision, intent(in) :: t, y(neq)
    double precision, intent(out) :: pd(nrowpd, neq)
    type(odepack_common_data), intent(inout) :: common_data
    
    ! Silence unused variable warnings
    call int_fake_use(neq)
    call r8_fake_use(t)
    call r8_array_fake_use(y)
    call int_fake_use(ml)
    call int_fake_use(mu)
    call int_fake_use(nrowpd)
    
    ! Initialize outputs
    pd = 0.0d0
    common_data%ierr = 0

    contains

    subroutine r8_array_fake_use(x)
        implicit none
        double precision, intent(in) :: x(:)
        
        ! Use array in a way that will never execute
        if (size(x) .lt. 0) then
            continue  ! Never executes since size() is always >= 0
        end if
    end subroutine r8_array_fake_use
    
    ! Dummy routine - not called when JT=2
end subroutine jac_dummy_silent

