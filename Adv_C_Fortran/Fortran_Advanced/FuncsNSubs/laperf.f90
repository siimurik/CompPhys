PROGRAM laperf
    implicit none
    integer :: num, i
    double precision, parameter :: pi = 4.D0*atan(1.D0)
    double precision :: val, x, phi
    print *, "Input x: "
    read *, x

    num = 100
    val = 0.0D0
    
    !do i = 0, num
    !    val = val + (-x**2.D0/2.D0)**i/((2.D0*i+1.D0)*gamma(dble(i)+1.D0))
    !end do

    do i = 0, num
        val = val + (-1.D0)**i * x**(2.D0*i+1.D0) / (gamma(dble(i)+1.D0) * (2.D0*i+1.D0))
    end do
    
    !phi = x/sqrt(2.D0*pi) * val
    phi = 2.D0/sqrt(pi) * val
    print *, "Estimated value of Phi: ", phi

END PROGRAM laperf