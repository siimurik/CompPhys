program main
!==================================================================
! Integration of a function using Gauss 16 and 32 points 
! with doubling number of intervals 
! till  error = |I_32 - I_16| < eps
!==================================================================
    implicit none
    double precision    :: a, b, integral, eps
    double precision    :: t1, t2, time
    integer :: nint
    double precision, parameter :: pi = 4.0*atan(1.0)
    double precision, external  :: f
    
! Integration limits and accuraccy.
    a =  5.D0
    b = 15.D0
    eps = 1.0e-15

    write(*,99)

    call cpu_time(t1)
    call gauss2(f,a,b,eps,integral,nint)
    call cpu_time(t2)
    time = t2 - t1
    
    write(*,103) nint, integral
    write(*,104) time
     99 format ('Integration of a function using Gauss 16 and 32 points with doubling', /, &
                'number of intervals till  error = |I_32 - I_16| < eps.',/)
    !103 format ('   intervals = ', i8, /,'   integral  = ', f12.8)
    103 format (/'   intervals = ', i8, /,'   integral  = ', 1pe15.8)
    104 format (/'Calculations took ', 1pe10.4, ' seconds.')
    stop
end program main

function f(x)
!----------------------------------------
! Function for integration
!----------------------------------------
    implicit none
    double precision :: f, x

    !f = sin(x)
    !f = x*cos(10.0*x**2)/(x**2 + 1.0)
!    f = (x-1251.D0/65.D0)**2.D0*(18.D0*x-3.D0*x**2.D0)
    f = 0.1D0*exp(-0.1D0*x)
    return

end function f

subroutine gauss2(f,a,b,eps,integral,nint)
!============================================================
! Integration of f(x) on [a,b]
! Method: Gauss quadratures with doubling number of intervals  
!         till  error = |I_16 - I_8| < eps
! written by: Alex Godunov (October 2009)
!------------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a	  - Lower limit of integration
! b	  - Upper limit of integration
! eps - tolerance
! OUT:
! integral - Result of integration
! nint     - number of intervals to achieve accuracy
!============================================================

    implicit none 
    double precision    :: a, b, eps, integral, gauss16, gauss32
    double precision    :: s1, s2, h, ax, bx
    integer :: nint, n, i
    integer, parameter :: nmax = 16384   ! max number of intervals
    double precision, external     :: f

! loop over number of intervals (starting from 1 interval)
    n = 1
    write(*,100)
    100 format ('   i           ax                  bx              s1              s2')
    do while (n <= nmax)
        s1 = 0.0
        s2 = 0.0
        h  = (b-a)/dfloat(n)
        do i = 1, n
            ax =  a + h*dfloat(i-1)
            bx = ax + h
            s1 = s1 + gauss16( f,ax,bx)
            s2 = s2 + gauss32(f,ax,bx)
        end do
        write(*,101) i, ax, bx, s1, s2
        101 format (i5, 4f20.10)
        if (abs(s2-s1) <= eps .and. abs(s2-s1)/abs(s2+s1) <= eps) then
            integral = s2
            nint = n
            exit
        end if
        integral = s2
        nint = n
        n = n*2
    end do
    return
end subroutine gauss2

function gauss16(f,a,b)
    !==========================================================
    ! Integration of f(x) on [a,b]
    ! Method: Gauss 16 points  
    ! written by: Alex Godunov (October 2009)
    !----------------------------------------------------------
    ! IN:
    ! f   - Function to integrate (supplied by a user)
    ! a	  - Lower limit of integration
    ! b	  - Upper limit of integration
    ! OUT:
    ! gauss16 - Result of integration
    !==========================================================
    implicit none
    integer, parameter :: n = 8
    double precision    :: gauss16, f, a, b
    double precision    :: ti(n), ci(n)
    data ti/0.0950125098, 0.2816035507, 0.4580167776, 0.6178762444, &  
            0.7554044083, 0.8656312023, 0.9445750230, 0.9894009349/ 
    data ci/0.1894506104, 0.1826034150, 0.1691565193, 0.1495959888, &
            0.1246289712, 0.0951585116, 0.0622535239, 0.0271524594/ 
    double precision    :: r, m, c
    integer :: i
    
    r = 0.0;
    m = (b-a)/2.0;
    c = (b+a)/2.0;
    
    do i = 1, n 
        r = r + ci(i)*(f(m*(-1.0)*ti(i) + c) + f(m*ti(i) + c))
    end do
    gauss16 = r*m
    return
end function gauss16

function gauss32(f,a,b)
    !==========================================================
    ! Integration of f(x) on [a,b]
    ! Method: Gauss 32 points  
    !----------------------------------------------------------
    ! IN:
    ! f   - Function to integrate (supplied by a user)
    ! a	  - Lower limit of integration
    ! b	  - Upper limit of integration
    ! OUT:
    ! gauss32 - Result of integration
    !==========================================================
    implicit none
    integer, parameter :: n = 32
    double precision    :: gauss32, f, a, b
    double precision    :: x32(n), w32(n)
    ! Define the abscissas and weights for 32-point Gauss quadrature
    data x32 /  -0.0483076656877383D0,	0.0483076656877383D0,   -0.144471961582797D0,	0.144471961582797D0,    &	
                -0.239287362252137D0,	0.239287362252137D0	,   -0.331868602282128D0,	0.331868602282128D0,    &	
                -0.421351276130635D0,	0.421351276130635D0	,   -0.506899908932229D0,	0.506899908932229D0,    &	
                -0.587715757240762D0,	0.587715757240762D0	,   -0.663044266930215D0,	0.663044266930215D0,    &	
                -0.73218211874029D0	 ,  0.73218211874029D0	,   -0.794483795967942D0,	0.794483795967942D0,    &	
                -0.84936761373257D0	 ,  0.84936761373257D0	,   -0.896321155766052D0,	0.896321155766052D0,    &	
                -0.93490607593774D0	 ,  0.93490607593774D0	,   -0.964762255587506D0,	0.964762255587506D0,    &	
                -0.985611511545268D0,	0.985611511545268D0	,   -0.997263861849482D0,	0.997263861849482D0 /

    data w32 /  0.0965400885147278D0,	0.0965400885147278D0,	0.0956387200792749D0,	0.0956387200792749D0,	&
                0.0938443990808046D0,	0.0938443990808046D0,	0.0911738786957639D0,	0.0911738786957639D0,	&
                0.0876520930044038D0,	0.0876520930044038D0,	0.0833119242269467D0,	0.0833119242269467D0,	&
                0.0781938957870703D0,	0.0781938957870703D0,	0.0723457941088485D0,	0.0723457941088485D0,	&
                0.0658222227763618D0,	0.0658222227763618D0,	0.0586840934785355D0,	0.0586840934785355D0,	&
                0.0509980592623762D0,	0.0509980592623762D0,	0.0428358980222267D0,	0.0428358980222267D0,	&
                0.0342738629130214D0,	0.0342738629130214D0,	0.0253920653092621D0,	0.0253920653092621D0,	&
                0.0162743947309057D0,	0.0162743947309057D0,	0.0070186100094701D0,	0.0070186100094701D0/
    double precision    :: r, m, c
    integer :: i
    
    r = 0.0;
    m = (b-a)/2.0;
    c = (b+a)/2.0;
    
    do i = 1, n
        r =  r + w32(i) * f(m*x32(i) + c)
    end do
    r = r*m
    return
end function gauss32
