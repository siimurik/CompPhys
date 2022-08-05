program main
!=============================================
! Integration of a function using Simpson rule
!=============================================
    implicit none
    real(8)    :: a, b, integral
    real(8)    :: t1, t2, time
    integer(8) :: n, i
    real(8), parameter :: pi = 4.0*atan(1.0)
    real(8), external  :: f
    
! Integration limits.
    a = 0.0
    b = pi

! Initial number of bars under the integral.
    n = 2

    write(*,99)
    write(*,100)

! Calling the simpson subroutine and formatting the answer.
    call cpu_time(t1)
    do i = 1, 16
        call simpson(f,a,b,integral,n)
        write (*,101) i, n, integral
        n = n*2
    end do
    call cpu_time(t2)
    time = t2 - t1

! Number of loops, number of intervals, result
    write(*,102) time
     99 format ('   Integration of a function using Simpson rule',/,/)
    100 format('        i    nint    simpson')
    !101 format((i9, i9, 1pe15.6))
    101 format((i9, i9, 1pe26.17))
    102 format (/'Calculation time is ', 1pe10.4, ' seconds.')
end program main

function f(x)
!----------------------------------------
! Function for integration
!----------------------------------------
    implicit none
    real(8) :: f, x

!    f = sin(x)
    f = x*cos(10.0*x**2)/(x**2 + 1.0)
    return

end function f

subroutine simpson(f,a,b,integral,n)
!==========================================================
! Integration of f(x) on [a,b]
! Method: Simpson rule for n intervals  
! written by: Alex Godunov (October 2009)
!----------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a	  - Lower limit of integration
! b	  - Upper limit of integration
! n   - number of intervals
! OUT:
! integral - Result of integration
!==========================================================

    implicit none 
    real(8)    :: f, a, b, integral, s
    real(8)    :: h, x
    integer(8) :: n, i

! if n is odd we add +1 to make it even
    if (n/2*2.ne.n) n = n + 1

! loop over n (number of intervals)
    s = 0.0
    h = (b-a)/dfloat(n)
    do i = 2, n-2, 2
        x = a + dfloat(i)*h
        s = s + 2.0*f(x) + 4.0*f(x+h)
    end do
    integral = (s + f(a) + f(b) + 4.0*f(a+h))*h/3.0
    return
end subroutine simpson