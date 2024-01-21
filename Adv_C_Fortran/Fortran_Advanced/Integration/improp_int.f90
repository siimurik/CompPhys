!==============================================
! Compile and execute with:
!   $ gfortran improp_int.f90 -o im
!   $ ./im 
!==============================================
program main
    implicit none
!    integer, parameter          :: DIM = 12
!    double precision    :: x(DIM), y(DIM)
    double precision            :: a, b
    double precision, external  :: func
    character(6)                :: midpnt = 'midpnt', midinf = 'midinf' 
    double precision            :: s1, s2, answer
    integer                     :: start_time, end_time, elapsed_time, rate
    double precision            :: elapsed_seconds

    ! Define integration bounds
    a = -2.0D0  ! left side
    b =  2.0D0  ! right side

    ! Get starting time
    call SYSTEM_CLOCK(count=start_time, count_rate=rate)

    ! Integrate func from a to b using midpoint rule
    call qromo(func,-2.D0, 2.0D0, s1, midpnt)

    ! Integrate func from b to infinity using mid-infinty rule
    ! Be cautious with the choice of bounds for midinf
    !call qromo(func, 1.D-4, 1.D30, s2, midinf)
    ! Sum up the results from both integrals
    !answer = s1 + s2

    call SYSTEM_CLOCK(count=end_time)
    elapsed_time = end_time - start_time
    elapsed_seconds = dble(elapsed_time) / dble(rate)

    write(*,*) "Computations completed."
    write(*,*) ""
    write(*,*) 'answer = ', s1
    write(*,*) ""
    write (*,'(a, 1PE11.4, a)') 'Elapsed time:', elapsed_seconds, ' seconds.'
    
end program main

FUNCTION func(x)
    !----------------------------------------
    ! Function for integration
    !----------------------------------------
    implicit none
    double precision :: func, x

    func = exp(-x**2.D0)
    return
    
END FUNCTION func

SUBROUTINE polint(xa,ya,n,x,y,dy)
    IMPLICIT NONE
    INTEGER             :: n
    DOUBLE PRECISION    :: dy,x,y,xa(n),ya(n)
    INTEGER, PARAMETER  :: NMAX=10 ! Largest anticipated value of n.
    !Given arrays xa and ya, each of length n, and given a value x, this routine returns a
    !value y, and an error estimate dy. If P (x) is the polynomial of degree N − 1 such that
    !P (xai ) = yai, i = 1, . . . , n, then the returned value y = P (x).
    INTEGER             :: i,m,ns
    DOUBLE PRECISION    :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns=1
    dif=abs(x-xa(1))
    do i=1, n       ! Here we find the index ns of the closest table entry,
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
            ns=i
            dif=dift
        endif
        c(i)=ya(i)  ! and initialize the tableau of c's and d's.
        d(i)=ya(i)
    enddo
    ! y=ya(ns) This is the initial approximation to y.
    ns=ns-1
    do m=1,n-1      ! For each column of the tableau,
        do i=1,n-m  ! we loop over the current c's and d's and update them.
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if(den.eq.0.) STOP !'failure in polint'
            ! This error can occur only if two input xa's are (to within roundoff) identical.
            den=w/den
            d(i)=hp*den ! Here the c's and d's are updated.
            c(i)=ho*den
        enddo
        if (2*ns.lt.n-m)then    ! After each column in the tableau is completed, we decide
            dy=c(ns+1)          ! which correction, c or d, we want to add to our accu-
        else                    ! mulating value of y, i.e., which path to take through
            dy=d(ns)            ! the tableau—forking up or down. We do this in such a    
            ns=ns-1             ! way as to take the most “straight line” route through the
        endif                   ! tableau to its apex, updating ns accordingly to keep track
    y=y+dy                      ! of where we are. This route keeps the partial approxima-
    enddo                       ! tions centered (insofar as possible) on the target x. The
    return                      ! last dy added is thus the error indication.
END SUBROUTINE polint

SUBROUTINE midpnt(func,a,b,s,n)
    IMPLICIT NONE
    INTEGER                     :: n
    DOUBLE PRECISION            :: a,b,s
    DOUBLE PRECISION, EXTERNAL  :: func
    ! This routine computes the nth stage of refinement of an extended midpoint rule. func is
    ! input as the name of the function to be integrated between limits a and b, also input. When
    ! called with n=1, the routine returns as s the crudest estimate of ∫ b a f (x)dx. Subsequent 
    ! calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 
    ! (2/3) × 3n-1 additional interior points. s should not be modified between sequential calls.
    INTEGER          :: it,j
    DOUBLE PRECISION :: ddel,del,sum,tnm,x
    if (n.eq.1) then
        s = (b-a)*func(0.5D0*(a+b))
    else
        it   = 3**(n-2)
        tnm  = it
        del  = (b-a)/(3.0D0*tnm)
        ddel = del+del !The added points alternate in spacing between del and ddel.
        x    = a + 0.5D0*del
        sum  = 0.D0
        do j = 1, it
            sum = sum+func(x)
            x   = x+ddel
            sum = sum+func(x)
            x   = x+del
        end do
        s = (s+(b-a)*sum/tnm)/3.D0  !The new sum is combined with the old integral to give a
                                    !refined integral.
    end if
    return
END SUBROUTINE

SUBROUTINE midinf(funk,aa,bb,s,n)
    INTEGER          :: n
    DOUBLE PRECISION :: aa,bb,s
    DOUBLE PRECISION, EXTERNAL :: funk
    ! This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
    ! of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
    ! points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
    ! as the computer allows, or the lower limit aa to be as large and negative, but not both.
    ! aa and bb must have the same sign.
    INTEGER          :: it,j
    DOUBLE PRECISION :: a,b,ddel,del,sum,tnm,func,x
    func(x) = funk(1./x)/x**2.D0  ! This statement function effects the change of variable.
    b = 1.D0/aa                 ! These two statements change the limits of integration ac-
    a = 1.D0/bb                 ! cordingly.
    if (n.eq.1) then            ! From this point on, the routine is exactly identical to midpnt.
        s=(b-a)*func(0.5D0*(a+b))
    else
        it   = 3**(n-2)
        tnm  = it
        del  = (b-a)/(3.D0*tnm)
        ddel = del + del
        x    = a + 0.5D0*del
        sum  = 0.D0
        do j = 1, it
            sum = sum+func(x)
            x   = x+ddel
            sum = sum+func(x)
            x   = x+del
        enddo
        s = (s+(b-a)*sum/tnm)/3.D0
    endif
    return
END SUBROUTINE midinf

SUBROUTINE qromo(func,a,b,ss,choose)
    IMPLICIT NONE
    INTEGER, PARAMETER          :: JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1
    DOUBLE PRECISION            :: a,b,ss
    DOUBLE PRECISION, EXTERNAL  :: func
    CHARACTER(6)                :: choose
    DOUBLE PRECISION, PARAMETER :: EPS=1.e-6
    !C USES polint
    !Romberg integration on an open interval. Returns as ss the integral of the function func
    !from a to b, using any specified integrating subroutine choose and Romberg's method.
    !Normally choose will be an open formula, not evaluating the function at the endpoints. It
    !is assumed that choose triples the number of steps on each call, and that its error series
    !contains only even powers of the number of steps. The routines midpnt, midinf, midsql,
    !midsqu, are possible choices for choose. The parameters have the same meaning as in
    !qromb.
    INTEGER             :: j
    DOUBLE PRECISION    :: dss,h(JMAXP),s(JMAXP)
    h(1)=1.
    do j=1,JMAX
        if (choose .eq. 'midpnt') then
            call midpnt(func, a, b, s(j), j)
        else if (choose .eq. 'midinf') then
            call midinf(func, a, b, s(j), j)
        else
            ! Handle other choices if needed
            write(*,*) 'Invalid choice!'
            return
        endif
        if (j.ge.K) then
            call polint(h(j-KM), s(j-KM), K, 0.D0, ss, dss)
            if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1) = s(j)
        h(j+1) = h(j)/9.D0  ! This is where the assumption of step tripling and an even
                            ! error series is used.
    enddo 
    STOP 'too many steps in qromo'
END SUBROUTINE

