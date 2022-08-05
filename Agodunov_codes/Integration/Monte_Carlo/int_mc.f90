program rinteg
!================================================================
! 1D intergration using Monte-Carlo method
! random numbers are generated either by rand() or urand
!
! written by Alex Godunov (October 2006)
!================================================================
    implicit none
    real(4)     :: rint, x, xint, a, b
    real(4)     :: rand
    integer(4)  :: i, nmax, key, ix
    real(4), external   :: fint, urand
    real(4), parameter  :: pi = 4.0*atan(1.0)

    write(*,*) ' enter number of random points'
    read(*,*) nmax

!*** initial values (rng, limits of integration)
    key = 2
    a   = 0.0
    b   = pi
    rint = 0.0
! initialization for random number generators
    call srand(1000)
    ix = 1000
!***
    do i = 1, nmax
        if(key .eq. 1) call random_number(x)
        if(key .eq. 2) x = urand(ix)
        xint = a + (b-a)*x
        rint = rint + fint(xint)
    end do
    rint = rint*(b-a)/dfloat(nmax)
    write(*,*) nmax, rint, x

    stop
end program

function fint(x)
!----------------------------------------
! Function for integration
!----------------------------------------
    implicit none
    real(4) :: fint, x
    fint = sin(x)
!      fint = 0.2/((x-6.0)**2 + 0.02)
!      fint = x*cos(10.0*x**2)/(x**2 + 1.0)
    return
end

function urand(iy)
!  urand is a uniform random number generator based  on  theory  and
!  suggestions  given  in  d.e. knuth (1969),  vol  2.   the integer  iy
!  should be initialized to an arbitrary integer prior to the first call
!  to urand.  the calling program should  not  alter  the  value  of  iy
!  between  subsequent calls to urand.  values of urand will be returned
!  in the interval (0,1).
    implicit none
    integer(4)  :: iy
    integer(4)  :: ia, ic, itwo, m2, m, mic
    real(4)     :: halfm, s, urand
    data m2 / 0 / , itwo / 2 /

    if (m2 .ne. 0) go to 20

!  if first entry, compute machine integer word length

        m = 1
    10  m2 = m
        m = itwo*m2
        if (m .gt. m2) go to 10
        halfm = m2

!  compute multiplier and increment for linear congruential method

    ia = 8*int(halfm*atan(1.D0)/8.D0) + 5
    ic = 2*int(halfm*(0.5d0-sqrt(3.D0)/6.D0)) + 1
    mic = (m2 - ic) + m2

!  s is the scale factor for converting to floating point

    s = 0.5/halfm

!  compute next random number

    20 iy = iy*ia

!  the following statement is for computers which do not allow
!  integer overflow on addition

    if (iy .gt. mic) iy = (iy - m2) - m2
!
    iy = iy + ic

!  the following statement is for computers where the
!  word length for addition is greater than for multiplication

    if (iy/2 .gt. m2) iy = (iy - m2) - m2

!  the following statement is for computers where integer
!  overflow affects the sign bit

    if (iy .lt. 0) iy = (iy + m2) + m2
    urand = float(iy)*s
    return

end function