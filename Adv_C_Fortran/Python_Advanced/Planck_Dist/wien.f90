!=====================================================
! Compile and execute with
!   $ gfortran wien.f90 -o wien
!   $ ./wien
!=====================================================
! Output:
!   x(0) = 5.00000000
!   x(1) = 4.96513570     |x(1) - x(0)| = 0.3486E-01
!   x(2) = 4.96511423     |x(2) - x(1)| = 0.2146E-04
!   x(3) = 4.96511423     |x(3) - x(2)| = 0.8327E-11
!   
!   Wien's displacement constant b = lambda_max * T
!   b = 2.897752E-03 m * T
!=====================================================
program wien
    implicit none

    double precision :: h, c, k, b
    double precision :: x, xvana, limit, f, df, vahed
    integer :: count

    ! Starting guess
    x = 5.0
    xvana = 0.0
    vahed = abs(x-xvana)
    count = 0
    write(*,9) count, x

    ! Define necessary constants
    h = 6.626068e-34 ! m^2 kg / s
    c = 2.99792458e8 ! m / s
    k = 1.380658e-23 ! J / K Boltzmann constant

    ! Desired accuracy
    limit = 1e-11

    ! Apply Newton's root-finding algorithm
    do while ( abs(x-xvana) .gt. limit )
        !print*, 'x_', count,'=', x, '|x - xvana| =', vahed
        xvana = x
        f  = 5.0 - 5.0*exp(-x) - x
        df = 5.0*exp(-x) - 1.0
        x = x - f/df
        vahed = abs(x-xvana)
        count = count + 1
        write(*,10) count, x, count, count-1, vahed
    enddo

    ! Printing out iterative values
9   format('x('i1') = ', 1F10.8)
10  format( 'x('i1') = ', 1F10.8, '     |x('i1') - x('i1')| = ', E10.4)

    ! Calculating Wien's diplacement costant
    b = h*c/(k*x)
    write(*, 11) b
11  format( /"Wien's displacement constant b = lambda_max * T",/, "b = ", ES12.6, " m * T" )

end program wien
