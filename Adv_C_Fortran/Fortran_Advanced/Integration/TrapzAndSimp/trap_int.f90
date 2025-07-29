!   $ gfortran trap_int.f90 -o trap
!   $ ./trap

module mymath

    contains

    function func(x) result(y)
        implicit none
        double precision, intent(in) :: x
        double precision             :: y

        y = dexp(-x**2.D0)

    end function func

    pure function integrate(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    double precision, intent(in)  :: x(:)         !! Variable x
    double precision, intent(in)  :: y(size(x))   !! Function y(x)
    double precision              :: r            !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
        r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2.D0
    end associate
  end function
    
end module mymath

program main
    use mymath
    implicit none
    integer, parameter  :: DIM = 1000000
    double precision    :: x(DIM), y(DIM), result
    double precision    :: a, b, nx, hx
    integer             :: i

    ! Define integration bounds
    a = -2.0D0  ! left side
    b =  2.0D0  ! right side

    ! Set appropiate step size based on DIM
    nx = DBLE(DIM)
    hx = abs(b-a)/(nx-1.D0)

    ! Defing the length of x and y vectors
    do i = 1, DIM
        x(i) = a + (i-1.0D0)*hx
        y(i) = func(x(i))
    end do

    !do i = 1, DIM
    !    WRITE (*,*) x(i), ',',y(i)
    !end do

    ! Define values for x (in this example, an array of 10 elements)
    !do i = 1, DIM
    !    x(i) = 0.1D0 * i
    !end do

    ! Evaluate the function at each x value
    !do i = 1, DIM
    !    y(i) = func(x(i))
    !end do

    ! Call the integration function
    result = integrate(x, y)

    ! Print the result
    print *, "Result of the integration: ", result
end program main
