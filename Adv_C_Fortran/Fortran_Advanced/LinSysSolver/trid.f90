!=============================================================================
!   $ gfortran trid.f90 -o trid
!   $ ./trid
! In GNU Plot:
!   $ plot 'data.csv' with lines
!=============================================================================
! This code solves the general N-by-N system A x = b
! where A is symmetric tridiagonal (N >= 2). The off-
! diagonal vector du and dl must be one element shorter 
! than the diagonal vector diag. The form of A for the 
! 4-by-4 case is shown below
!        [ b1   c1   0     0     0 ]   [u1]     [r1]  
!        [ a2  b2   c2     0     0 ]   [u2]     [r2] 
!        [ 0    a3  b3    c3     0 ] * [u3]  =  [r3] 
!        [ 0    0   a4    b4    c4 ]   [u4]     [r4] 
!        [ 0    0    0    a5    b5 ]   [u5]     [r5]  
! Note how lower diagonal 'a' starts at 2 and ends at N, 
! while the upper diagonal 'b' starts at 1, but ends at 
! N-1. The solution is storet in vector 'u'.
!-----------------------------------------------------------------------------
! Suprisingly, using the tridag subroutine to solve this 
! problem is actually faster than LAPACK's dgtsv solver.
!=============================================================================

PROGRAM main
    INTEGER, PARAMETER :: N = 5001
    DOUBLE PRECISION   :: a(N), b(N), c(N), r(N), u(N), x(N)
    DOUBLE PRECISION   :: h

    ! Printing out the dimensions
    write (*,1) N-1, N-1
1   format ("Dimension of tridiagonal system:", i5," x", i5, ".",/)

    ! Initialize tridiagonal system coefficients and right-hand side
    a = 0.D0
    b = 0.D0
    c = 0.D0
    r = 0.D0

    ! Initializing 
    h = 2.0D0/dble(N-1.D0)
    do i = 1, N
        x(i) = (i-1.D0)*h
    end do

    ! Boundary conditions
    b(1) =   1.D0/h + h/2.D0 * (4.D0-x(1))
    r(1) = h/2.D0*(x(1)+5.D0) - 1.D0

    ! Elements of the main diagonal b and 
    ! righthand side vector r
    do i = 2, N-1
        b(i) =  2.D0/h +h*(4.D0-x(i))
        r(i) =  h*(x(i) + 5.D0)
    end do

    ! Elements of off-diagonal elements
    ! Upper diagonal
    do i = 1, N-1
        c(i) = -1.D0/h
    end do
    ! Lower diagonal
    do i = 2, N
        a(i) = -1.D0/h
    end do

    ! Final values
    b(N) =   1.D0/h + h/2.D0 * (4.D0-x(N))
    r(N) = h/2.D0*(x(N)+5.D0) - 1.D0

    ! record the starting time
    call cpu_time(start_time)

    ! Call the tridag subroutine to solve the tridiagonal system
    CALL tridag(a, b, c, r, u, N)

    ! record the ending time
    call cpu_time(end_time)

    ! Print the solution vector u
    ! printing out the first and last 3 elements
    write(*,*)'DGTSV Program Results'
    write(*,*) '	x	    u'
    do i = 1, 3
        write(*,'(6f12.6)') x(i), u(i)
    end do
    write(*,*) '            ...'
    do i = N-2, N
        write(*,'(6f12.6)') x(i), u(i)
    end do

    ! Print the solution into a separate file
    open(unit = 11, file="data.csv")
    do i = 1, N
        write(11,*) x(i), u(i)
    end do
    close(11)

    ! compute the elapsed time
    elapsed_time = end_time - start_time

    ! print the elapsed time
    write(*,9) elapsed_time
9   format (/'Elapsed time: ', E10.4, ' seconds.')
!9   format (/'Elapsed time: ', F8.6 ' seconds.')

END PROGRAM main


SUBROUTINE tridag(a,b,c,r,u,n)
    !INTEGER, PARAMETER :: NMAX=n
    INTEGER          :: n
    DOUBLE PRECISION :: a(n),b(n),c(n),r(n),u(n)
    !Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1).
    !a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are not modified.
    !Parameter: NMAX is the maximum expected value of n.
    INTEGER          :: j
    DOUBLE PRECISION :: bet, gam(n) !One vector of workspace, gam is needed.
    if (b(1) .eq. 0.D0) stop 'tridag: rewrite equations'
    !If this happens then you should rewrite your equations as a set of order N − 1, with u2
    !trivially eliminated.
    bet  = b(1)
    u(1) = r(1)/bet
    do j = 2, n ! Decomposition and forward substitution.
        gam(j) = c(j-1)/bet
        bet    = b(j) - a(j)*gam(j)
        if (bet .eq. 0.0D0) stop 'tridag failed' ! Algorithm fails; see below.
        u(j)   = (r(j)-a(j)*u(j-1))/bet
    end do
    do j = n-1, 1, -1 ! Backsubstitution.
        u(j) = u(j) - gam(j+1)*u(j+1)
    end do
    return
END SUBROUTINE tridag

! There is no pivoting in tridag. It is for this reason that tridag can fail
! (pause) even when the underlying matrix is nonsingular: A zero pivot can be
! encountered even for a nonsingular matrix. In practice, this is not something to lose
! sleep about. The kinds of problems that lead to tridiagonal linear sets usually have
! additional properties which guarantee that the algorithm in tridag will succeed.
! For example, if
!           |bj | > |aj | + |cj | j = 1, . . . , N (2.4.2)
! (called diagonal dominance) then it can be shown that the algorithm cannot encounter
! a zero pivot.
! It is possible to construct special examples in which the lack of pivoting in the
! algorithm causes numerical instability. In practice, however, such instability is almost
! never encountered — unlike the general matrix problem where pivoting is essential.
! The tridiagonal algorithm is the rare case of an algorithm that, in practice, is
! more robust than theory says it should be. Of course, should you ever encounter a
! problem for which tridag fails, you can instead use the more general method for
! band diagonal systems (routines bandec and banbks).