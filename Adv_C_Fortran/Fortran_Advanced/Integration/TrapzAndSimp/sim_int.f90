! Compile and execute with:
!   $ gfortran sim_int.f90 -o int
!   $ ./int
program Simpsons
    implicit none

    ! Declare variables
    real(8) :: a, b, h, x, s, result
    real(8) :: t1, t2, time
    integer :: n, i, j

    ! Set the lower and upper limits of integration
    a = -10.D0
    b =  10.D0

    ! Set the number of subintervals to use
    n = 2
    
    ! Perform the integration using the Simpson's method
    call cpu_time(t1)
    do i = 1, 16
        s = 0.0
        h = (b-a)/dfloat(n) ! Calculate the step size
        do j = 2, n-2, 2
            x = a + dfloat(j)*h
            s = s + 2.D0*func(x) + 4.D0*func(x+h)
        end do
        result = (s + func(a) + 4.D0*func(a+h) + func(b))*h/3.0
        write (*,101) i, n, result
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

contains

    ! Define the function to integrate
    real(8) function func(x)
    real(8), intent(in) :: x

        ! Insert the expression for the function here
        !func = x**2.D0
        func = exp(-x**2)

    end function func

end program
