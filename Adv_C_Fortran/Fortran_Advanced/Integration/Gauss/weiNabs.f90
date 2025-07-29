!==============================================
! Compile and execute with:
!   $ gfortran weiNabs.f90 -o wa
!   $ ./wa
!==============================================
PROGRAM main
    INTEGER, PARAMETER :: n = 10  ! or any desired number of points
    DOUBLE PRECISION :: x(n), w(n)
    DOUBLE PRECISION :: x1, x2
    x1 = -1.0d0   ! Set lower limit
    x2 =  1.0d0   ! Set upper limit

    ! Subroutine which generates weights and abscissae 
    ! for use in performing Gauss-Legendre quadrature 
    ! integral approximation
    CALL gauleg(x1, x2, x, w, n)

    PRINT *, "Weights and Abscissae:"
    DO i = 1, n
        WRITE(*, '(a,i3,a,F19.16,a,i3,a,F19.16)') "w(", i, ") = ", w(i), "   x(", i, ") = ", x(i)
    END DO

END PROGRAM main

SUBROUTINE gauleg(x1,x2,x,w,n)
    INTEGER                     :: n
    DOUBLE PRECISION            :: x1,x2,x(n),w(n)
    DOUBLE PRECISION, PARAMETER :: EPS = 3.d-14, pi = 4.0D0*atan(1.0D0) !EPS is the relative precision.
    !Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
    !arrays x(1:n) and w(1:n) of length n, containing the abscissas and weights of the Gauss-
    !Legendre n-point quadrature formula.
    INTEGER             :: i,j,m
    DOUBLE PRECISION    :: p1,p2,p3,pp,xl,xm,z,z1
    !High precision is a good idea for this routine.
    m  = (n+1)/2        !The roots are symmetric in the interval, so we
    xm = 0.5d0*(x2+x1)  !only have to find half of them.
    xl = 0.5d0*(x2-x1)
    do i=1, m           !Loop over the desired roots.
        z = cos(pi*(i-0.25d0)/(n+0.5d0))
        ! Starting with the above approximation to the ith root, we enter the main loop of re-
        ! finement by Newton's method.
        1 continue
        p1 = 1.d0
        p2 = 0.d0
        do j=1,n        ! Loop up the recurrence relation to get the Leg-
            p3=p2       ! endre polynomial evaluated at z.
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        end do
        !p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by
        !a standard relation involving also p2, the polynomial of one lower order.
        pp = n*(z*p1-p2)/(z*z-1.d0)
        z1 = z
        z  = z1 - p1/pp ! Newton's method.
        if(abs(z-z1).gt.EPS) goto 1
        x(i)     = xm-xl*z                      !Scale the root to the desired interval,
        x(n+1-i) = xm+xl*z                      !and put in its symmetric counterpart.
        w(i)     = 2.d0*xl/((1.d0-z*z)*pp*pp)   !Compute the weight
        w(n+1-i) = w(i)                         ! and its symmetric counterpart.
    end do
    return
END SUBROUTINE 