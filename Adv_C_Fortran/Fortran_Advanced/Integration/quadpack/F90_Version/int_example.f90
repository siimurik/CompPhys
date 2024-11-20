!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Doing Integration in FORTRAN usind QUADPACK functions with
! double precision 
!-------------------------------------------------------------------
! These functions are written by complete MANIACS because compiling 
! this simple code with
!
!       gfortran int_example.f90 dquad.f90 -o int
!
! Gives a peculiar warning:
!
!   "/usr/bin/ld: warning: /tmp/ccKX5W65.o: requires executable &
!   stack (because the .note.GNU-stack section is executable)"
!
! According to some research this "indicates that some part of this 
! code requires an executable stack. This is not necessarily a 
! critical issue, but it can be a SECURITY CONCERN, as having an 
! executable stack can make the program more vulnerable to certain 
! types of exploits."
!
! My brother in Christ... I am just trying to solve an integral here.
!
! Also, the supposed error is somehow on the scale of 10E-21...
! How is this even possible when I am using double precision????
!-------------------------------------------------------------------
! Compile and execute with a more "cleaner" version:
!       gfortran -Wall -Wextra -o int int_example.f90 dquadSAFER.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program integrate_example
    implicit none
    double precision :: a, b, epsabs, epsrel, result, abserr
    integer :: key, limit, neval, ier
    double precision, dimension(100) :: alist, blist, rlist, elist
    integer, dimension(100) :: iord
    integer :: last

    ! Define the integration limits
    a = 0.0d0
    b = 1.0d0

    ! Set the desired absolute and relative accuracy
    epsabs = 1.0d-6
    epsrel = 1.0d-6

    ! Choose the integration rule (1 to 6)
    key = 2  ! Using 10-21 point Gauss-Kronrod pair

    ! Set the limit for the number of subintervals
    limit = 100

    ! Call the dqage subroutine
    call dqage(f, a, b, epsabs, epsrel, key, limit, result, abserr, &
                neval, ier, alist, blist, rlist, elist, iord, last)

    ! Check for errors and print results
    if (ier == 0) then
        print *, 'Integral result: ', result
        print *, 'Estimated absolute error: ', abserr
        print *, 'Number of function evaluations: ', neval
    else
        print *, 'Error in integration: ier = ', ier
    end if

contains

    ! Define the integrand function
    double precision function f(x)
        double precision :: x
        f = exp(-x**2)  ! Example: f(x) = e^(-x^2)
    end function f

end program integrate_example

! dqage subroutine should be included here or in a separate file.