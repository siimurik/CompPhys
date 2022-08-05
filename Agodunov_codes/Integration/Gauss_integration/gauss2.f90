program main
!==================================================================
! Integration of a function using Gauss 8 and 16 points 
! with doubling number of intervals 
! till  error = |I_16 - I_8| < eps
!==================================================================
    implicit none
    real(8)    :: a, b, integral, eps
    real(8)    :: t1, t2, time
    integer(8) :: nint
    real(8), parameter :: pi = 4.0*atan(1.0)
    real(8), external  :: f
    
! Integration limits and accuraccy.
    a = 0.0
    b = pi
    eps = 1.0e-8

    write(*,99)

    call cpu_time(t1)
    call gauss2(f,a,b,eps,integral,nint)
    call cpu_time(t2)
    time = t2 - t1
    
    write(*,103) nint, integral
    write(*,104) time
     99 format ('Integration of a function using Gauss 8 and 16 points with doubling', /, &
                'number of intervals till  error = |I_16 - I_8| < eps.',/)
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
    real(8) :: f, x

    f = sin(x)
!    f = x*cos(10.0*x**2)/(x**2 + 1.0)
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
    real(8)    :: a, b, eps, integral, gauss8, gauss16
    real(8)    :: s1, s2, h, ax, bx
    integer(8) :: nint, n, i
    integer(8), parameter :: nmax = 16384   ! max number of intervals
    real(8), external     :: f

! loop over number of intervals (starting from 1 interval)
    n = 1
    write(*,100)
    100 format ('   i           ax              bx              s1              s2')
    do while (n <= nmax)
        s1 = 0.0
        s2 = 0.0
        h  = (b-a)/dfloat(n)
        do i = 1, n
            ax =  a + h*dfloat(i-1)
            bx = ax + h
            s1 = s1 + gauss8( f,ax,bx)
            s2 = s2 + gauss16(f,ax,bx)
        end do
        write(*,101) i, ax, bx, s1, s2
        101 format (i5, 4f15.10)
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

function gauss8(f,a,b)
!==========================================================
! Integration of f(x) on [a,b]
! Method: Gauss 8 points  
! written by: Alex Godunov (October 2009)
!----------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a	  - Lower limit of integration
! b	  - Upper limit of integration
! OUT:
! gauss8 - Result of integration
!==========================================================
    implicit none
    integer(8), parameter :: n = 4
    real(8)    :: gauss8, f, a, b
    real(8)    :: ti(n), ci(n)
    data ti/0.1834346424, 0.5255324099, 0.7966664774, 0.9602898564/
    data ci/0.3626837833, 0.3137066458, 0.2223810344, 0.1012285362/ 
    real(8)    :: r, m, c
    integer(8) :: i
    
    r = 0.0;
    m = (b-a)/2.0;
    c = (b+a)/2.0;
    
    do i = 1, n 
        r = r + ci(i)*(f(m*(-1.0)*ti(i) + c) + f(m*ti(i) + c))
    end do
    gauss8 = r*m
    return
end function gauss8

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
    integer(8), parameter :: n = 8
    real(8)    :: gauss16, f, a, b
    real(8)    :: ti(n), ci(n)
    data ti/0.0950125098, 0.2816035507, 0.4580167776, 0.6178762444, &  
            0.7554044083, 0.8656312023, 0.9445750230, 0.9894009349/ 
    data ci/0.1894506104, 0.1826034150, 0.1691565193, 0.1495959888, &
            0.1246289712, 0.0951585116, 0.0622535239, 0.0271524594/ 
    real(8)    :: r, m, c
    integer(8) :: i
    
    r = 0.0;
    m = (b-a)/2.0;
    c = (b+a)/2.0;
    
    do i = 1, n 
        r = r + ci(i)*(f(m*(-1.0)*ti(i) + c) + f(m*ti(i) + c))
    end do
    gauss16 = r*m
    return
end function gauss16
    