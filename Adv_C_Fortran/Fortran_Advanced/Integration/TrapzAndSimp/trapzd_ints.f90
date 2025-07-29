!==============================================
! Compile and execute with:
!   $ gfortran trapzd_ints.f90 -o im
!   $ ./im 
!==============================================
program main
    implicit none
    double precision    :: a, b
    double precision, external  :: func
    double precision            :: s
    integer                     :: start_time, end_time, elapsed_time, rate
    double precision            :: elapsed_seconds

    ! Define integration bounds
    a = -2.0D0  ! left side
    b =  2.0D0  ! right side

    ! Simple example for using trapzd() on its own
    !m = 20
    !do j = 1, m
    !    call trapzd(func,a,b,s,j)
    !enddo 

    ! Get starting time
    call SYSTEM_CLOCK(count=start_time, count_rate=rate)

    ! Choose a solver: trapezoid or Simpson
    !call qtrap(func,a,b,s)
    call qsimp(func,a,b,s)  ! 100x faster than qtrap()

    ! Get end time and calculate the difference
    call SYSTEM_CLOCK(count=end_time)
    elapsed_time = end_time - start_time
    elapsed_seconds = dble(elapsed_time) / dble(rate)

    write(*,*) "Computations completed."
    write(*,*) ""
    write(*,*) 'answer = ', s
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
    DOUBLE PRECISION :: del,sum,tnm,x
    INTEGER :: it, j
    if (n.eq.1) then
        s = 0.5D0*(b-a)*(func(a)+func(b))
    else
        it  = 2**(n-2)
        tnm = it
        del = (b-a)/tnm ! This is the spacing of the points to be added.
        x   = a+0.5D0*del
        sum = 0.0D0
        do j=1, it
            sum = sum+func(x)
            x   = x+del
        end do 
        s = 0.5D0*(s+(b-a)*sum/tnm) ! This replaces s by its refined value.
    end if
    return
END SUBROUTINE trapzd

SUBROUTINE qtrap(func,a,b,s)
    !==========================================================================================
    ! USES trapzd()
    ! Returns as s the integral of the function func from a to b. The parameters EPS can be set
    ! to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
    ! allowed number of steps. Integration is performed by the trapezoidal rule.
    !==========================================================================================
    IMPLICIT NONE
    INTEGER, PARAMETER          :: JMAX = 26
    DOUBLE PRECISION            :: a,b,s
    DOUBLE PRECISION, EXTERNAL  :: func
    DOUBLE PRECISION, PARAMETER :: EPS=1.0D-14
    DOUBLE PRECISION :: olds
    INTEGER          :: j
    olds = -1.e30   ! Any number that is unlikely to be the average of the function
    do j = 1, JMAX  ! at its endpoints will do here.
        call trapzd(func,a,b,s,j)
        if (j.gt.5) then    ! Avoid spurious early convergence.
            if (abs(s-olds).lt.EPS*abs(olds).or.(s.eq.0..and.olds.eq.0.)) return
        end if
        olds = s
    end do
    stop 'too many steps in qtrap (j > JMAX = 26)'
END SUBROUTINE qtrap

SUBROUTINE qsimp(func,a,b,s)
    !=========================================================================================
    ! USES trapzd
    !Returns as s the integral of the function func from a to b. The parameters EPS can be set
    !to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
    !allowed number of steps. Integration is performed by Simpson's rule.
    !=========================================================================================
    IMPLICIT NONE
    INTEGER, PARAMETER          :: JMAX = 26
    DOUBLE PRECISION            :: a,b,s
    DOUBLE PRECISION, EXTERNAL  :: func
    DOUBLE PRECISION, PARAMETER :: EPS=1.0D-14
    INTEGER          :: j
    DOUBLE PRECISION :: os,ost,st
    ost  = -1.e30
    os   = -1.e30
    do j = 1,JMAX
        call trapzd(func,a,b,st,j)
        s = (4.*st-ost)/3.  ! Compare equation (4.2.4), above.
        if (j.gt.5) then    ! Avoid spurious early convergence.
            if (abs(s-os).lt.EPS*abs(os).or.(s.eq.0..and.os.eq.0.)) return
        end if
        os  = s
        ost = st
    end do
    stop 'too many steps in qsimp'
END SUBROUTINE  qsimp