!==============================================
! Compile and execute with:
!   $ gfortran romb_int.f90 -o romb
!   $ ./romb
!==============================================
program main
    implicit none
    double precision           :: a, b
    double precision, external :: func
    double precision           :: s
    integer                    :: start_time, end_time, elapsed_time, rate
    double precision           :: elapsed_seconds

    ! Define integration bounds
    a =  0.0D0  ! left side
    b =  2.0D0  ! right side

    ! Simple example for using trapzd() on its own
    !m = 20
    !do j = 1, m
    !    call trapzd(func,a,b,s,j)
    !enddo 

    ! Get starting time
    call SYSTEM_CLOCK(count=start_time, count_rate=rate)

    ! Romberg Integration Method
    call qromb(func,a,b,s)

    ! Get end time and calculate the difference
    call SYSTEM_CLOCK(count=end_time)
    elapsed_time = end_time - start_time
    elapsed_seconds = dble(elapsed_time)/dble(rate)

    write(*,*) "Computations completed."
    write(*,*) ""
    write(*,*) 'answer = ', s
    write(*,*) ""
    write (*,'(a, 1PE12.3, a)') 'Elapsed time:', elapsed_seconds, ' seconds.'

end program main

FUNCTION func(x)
    !----------------------------------------
    ! Function for integration
    !----------------------------------------
    implicit none
    double precision, parameter :: pi = 4.D0*atan(1.D0)
    double precision :: func, x

    func = 1.D0/4.D0 * pi * x*x*x*x * &
            cos(1.D0/4.D0 * pi * x)
    return
    
END FUNCTION func

SUBROUTINE trapzd(func,a,b,s,n)
    !==============================================================================================
    ! This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
    ! input as the name of the function to be integrated between limits a and b, also input. When
    ! called with n=1, the routine returns as s the crudest estimate of ∫f(x)dx from a to b. 
    ! Subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy of s 
    ! by adding 2n-2 additional interior points. s should not be modified between sequential calls.
    !==============================================================================================
    IMPLICIT NONE
    INTEGER n
    DOUBLE PRECISION           :: a,b,s
    DOUBLE PRECISION, EXTERNAL :: func
    DOUBLE PRECISION           :: del,sum,tnm,x
    INTEGER                    :: it, j
    if (n.eq.1) then
        s = 0.5*(b-a)*(func(a)+func(b))
    else
        it  = 2**(n-2)
        tnm = it
        del = (b-a)/tnm ! This is the spacing of the points to be added.
        x   = a+0.5*del
        sum = 0.0
        do j=1, it
            sum = sum+func(x)
            x   = x+del
        end do 
        s = 0.5*(s+(b-a)*sum/tnm) ! This replaces s by its refined value.
    end if
    return
END SUBROUTINE trapzd

SUBROUTINE polint(xa,ya,n,x,y,dy)
    !======================================================================================
    ! Neville’s algorithm:
    ! Given arrays xa and ya, each of length n, and given a value x, this routine returns a
    ! value y, and an error estimate dy. If P(x) is the polynomial of degree N-1 such that
    ! P(xai) = yai, i = 1, ..., n, then the returned value y = P(x).
    !======================================================================================
    IMPLICIT NONE
    INTEGER            :: n
    INTEGER, PARAMETER :: NMAX=10 ! Largest anticipated value of n.
    DOUBLE PRECISION   :: dy,x,y,xa(n),ya(n)
    INTEGER          :: i, m, ns
    DOUBLE PRECISION :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns   = 1
    dif  = abs(x-xa(1))
    do i = 1, n ! Here we find the index ns of the closest table entry,
        dift = abs(x-xa(i))
        if (dift.lt.dif) then
            ns  = i
            dif = dift
        end if
        c(i) = ya(i)    ! and initialize the tableau of c's and d's.
        d(i) = ya(i)
    end do
    y  = ya(ns)         ! This is the initial approximation to y.
    ns = ns-1
    do m = 1, n-1       ! For each column of the tableau,
        do i = 1,n-m    ! we loop over the current c's and d's and update them.
            ho  = xa(i)-x
            hp  = xa(i+m)-x
            w   = c(i+1)-d(i)
            den = ho-hp
            if (den .eq. 0.0) stop 'failure in polint'
            !This error can occur only if two input xa’s are (to within roundoff) identical.
            den  = w/den
            d(i) = hp*den       ! Here the c's and d's are updated.
            c(i) = ho*den
        end do
        if (2*ns .lt. n-m) then ! After each column in the tableau is completed, we decide
            dy = c(ns+1)        ! which correction, c or d, we want to add to our accu-
        else                    ! mulating value of y, i.e., which path to take through
            dy = d(ns)          ! the tableau—forking up or down. We do this in such a
            ns = ns-1           ! way as to take the most “straight line” route through the
        end if                  ! tableau to its apex, updating ns accordingly to keep track
        y=y+dy                  ! of where we are. This route keeps the partial approxima-
    end do                      ! tions centered (insofar as possible) on the target x. The
    return                      ! last dy added is thus the error indication.
END SUBROUTINE polint

SUBROUTINE qromb(func,a,b,ss)
    !========================================================================================
    ! USES polint, trapzd
    !Returns as ss the integral of the function func from a to b. Integration is performed by
    !Romberg’s method of order 2K, where, e.g., K=2 is Simpson’s rule.
    !Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation
    !error estimate; JMAX limits the total number of steps; K is the number of points used in
    !the extrapolation.
    !========================================================================================
    IMPLICIT NONE
    INTEGER, PARAMETER          :: JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1
    DOUBLE PRECISION            :: a,b,ss
    DOUBLE PRECISION, EXTERNAL  :: func
    DOUBLE PRECISION, PARAMETER :: EPS=1.D-14
    INTEGER          :: j
    DOUBLE PRECISION :: dss,h(JMAXP),s(JMAXP) ! These store the successive trapezoidal approximations
    !and their relative stepsizes.
    h(1) = 1.0
    do j = 1, JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
            call polint(h(j-KM), s(j-KM), K, 0.D0, ss, dss)
            if (abs(dss) .le. EPS*abs(ss)) return
        end if
        s(j+1) = s(j)
        h(j+1) = 0.25*h(j)          ! This is a key step: The factor is 0.25 even though
    end do                          ! the stepsize is decreased by only 0.5. This makes
    stop 'too many steps in qromb'  ! the extrapolation a polynomial in h2 as allowed
END SUBROUTINE qromb                ! by equation (4.2.1), not just a polynomial in h.