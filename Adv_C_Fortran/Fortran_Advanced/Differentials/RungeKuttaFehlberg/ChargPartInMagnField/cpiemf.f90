! ===============================================================================
! CHARGED PARTICLE MOTION IN EARTH'S MAGNETIC FIELD
! ===============================================================================
! 
! This program simulates the motion of a charged particle (proton) in Earth's
! dipole magnetic field using the Runge-Kutta-Fehlberg (RKF45) integration method.
!
! PHYSICAL SYSTEM:
! The motion is governed by the Lorentz force equation:
!   m * dv/dt = q * (v × B)
! 
! Where:
!   - m = particle mass (kg)
!   - q = particle charge (C)
!   - v = particle velocity vector (m/s)
!   - B = magnetic field vector (T)
!
! DIFFERENTIAL EQUATION SYSTEM:
! The 6 first-order ODEs represent position and velocity components:
!   dx/dt = vx    (position derivatives)
!   dy/dt = vy
!   dz/dt = vz
!   dvx/dt = (q/m) * (vy*Bz - vz*By)    (Lorentz force components)
!   dvy/dt = (q/m) * (vz*Bx - vx*Bz)
!   dvz/dt = (q/m) * (vx*By - vy*Bx)
!
! COORDINATE SYSTEM:
! - Origin at Earth's center
! - Units: kilometers for distance, km/s for velocity
! - Earth radius R_earth = 6378.137 km
!
! MAGNETIC FIELD MODEL:
! Earth's dipole field approximation:
!   B = (μ₀/4π) * (M/r³) * [3(M̂·r̂)r̂ - M̂]
! Simplified to: B₀ * magnetic_dipole_field(x,y,z)
! ===============================================================================

PROGRAM charged_particle_motion
    IMPLICIT NONE
    
    ! Variable declarations
    REAL    :: y(6), work(39), m, q, qm, B0, R_maa
    REAL    :: t, dt, tout, relerr, abserr
    INTEGER :: iwork(5), neqn, nt, i, iflag
    COMMON qm, R_maa
    EXTERNAL func
    
    ! Open output file for trajectory data
    OPEN(10, FILE="proton.dat")
    
    ! Physical constants
    q = 1.602e-19        ! Elementary charge (C)
    m = 1.673e-27         ! Proton mass (kg)
    qm = q/m             ! Charge-to-mass ratio (C/kg)
    B0 = 3.07e-5         ! Earth's magnetic dipole moment coefficient (T)
    R_maa = 6378.137     ! Earth's radius (km)
    
    ! Numerical integration parameters
    neqn = 6             ! Number of differential equations
    relerr = 1.e-8       ! Relative error tolerance
    abserr = 1.e-8       ! Absolute error tolerance
    t = 0.               ! Initial time (s)
    dt = 0.01            ! Time step (s)
    nt = 20000           ! Number of time steps
    iflag = 1            ! RKF45 control flag
    
    ! Initial conditions
    ! Position (km): x=30000, y=0, z=0 (above Earth's equator)
    y(1) = 30000.        ! x-position
    y(2) = 0.            ! y-position  
    y(3) = 0.            ! z-position
    
    ! Velocity (km/s): complex 3D motion
    y(4) = 0.            ! x-velocity
    y(5) = 10000.        ! y-velocity
    y(6) = 20000.        ! z-velocity
    
    ! Time integration loop
    DO i = 1, nt
        tout = t + dt
        CALL rkf45(func,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)

        if (iflag .ne. 2) then
            iflag = 2
        endif

        WRITE(10,*) tout, y
    enddo
    
    CLOSE(10)

END PROGRAM

! ===============================================================================
! SUBROUTINE: MAGNETIC FIELD CALCULATION
! ===============================================================================
! Calculates Earth's dipole magnetic field at position (x,y,z)
! 
! INPUT:  x, y, z - position coordinates (normalized by Earth radius)
! OUTPUT: B(3)    - magnetic field components [Bx, By, Bz] (Tesla)
!
! The dipole field formula assumes the magnetic dipole is aligned with z-axis:
!   Bx = -B₀ * 3xz/r⁵
!   By = -B₀ * 3yz/r⁵  
!   Bz = -B₀ * (2z² - x² - y²)/r⁵
! where r = √(x² + y² + z²)
! ===============================================================================
SUBROUTINE calc_magn_field(x, y, z, B)
    IMPLICIT NONE
    REAL x, y, z, B(3), r, B0
    
    ! Calculate distance from Earth's center
    r = SQRT(x**2 + y**2 + z**2)
    
    ! Earth's magnetic field strength at surface
    B0 = 3.07e-5
    
    ! Calculate dipole field components
    B(1) = -B0 * 3. * x * z / r**5    ! Bx component
    B(2) = -B0 * 3. * y * z / r**5    ! By component  
    B(3) = -B0 * (2.*z*z - x*x - y*y) / r**5  ! Bz component
      
END SUBROUTINE

! ===============================================================================
! SUBROUTINE: DIFFERENTIAL EQUATION SYSTEM
! ===============================================================================
! Defines the system of 6 first-order ODEs for charged particle motion
!
! INPUT:  t    - current time
!         y(6) - state vector [x, y, z, vx, vy, vz]
! OUTPUT: dy(6) - derivatives [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]
!
! The equations represent:
! - Position derivatives: dr/dt = v
! - Velocity derivatives: dv/dt = (q/m) * (v × B)
! ===============================================================================
SUBROUTINE func(t, y, dy)
    IMPLICIT NONE
    REAL :: t, y(6), dy(6), B(3), qm, R_maa
    COMMON qm, R_maa

    IF (.FALSE.) THEN
        WRITE(*,*) t
    END IF
    
    ! Position derivatives (velocity components)
    dy(1) = y(4)  ! dx/dt = vx
    dy(2) = y(5)  ! dy/dt = vy
    dy(3) = y(6)  ! dz/dt = vz
    
    ! Calculate magnetic field at current position (normalized coordinates)
    CALL calc_magn_field(y(1)/R_maa, y(2)/R_maa, y(3)/R_maa, B)
    
    ! Velocity derivatives (Lorentz force: F = q(v × B))
    dy(4) = qm * (y(5)*B(3) - y(6)*B(2))  ! dvx/dt
    dy(5) = qm * (y(6)*B(1) - y(4)*B(3))  ! dvy/dt  
    dy(6) = qm * (y(4)*B(2) - y(5)*B(1))  ! dvz/dt
      
END SUBROUTINE