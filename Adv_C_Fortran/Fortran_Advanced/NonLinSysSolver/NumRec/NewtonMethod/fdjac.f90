SUBROUTINE fdjac(n,x,fvec,np,df)
    INTEGER n,np,NMAX
    REAL df(np,np),fvec(n),x(n),EPS
    PARAMETER (NMAX=40, EPS=1.e-4)
    ! USES funcv
    ! Computes forward-difference approximation to Jacobian. On input, x(1:n) is the point
    ! at which the Jacobian is to be evaluated, fvec(1:n) is the vector of function values at
    ! the point, and np is the physical dimension of the Jacobian array df(1:n,1:n) which is
    ! output. subroutine funcv(n,x,f) is a fixed-name, user-supplied routine that returns
    ! the vector of functions at x.
    ! Parameters: NMAX is the maximum value of n; EPS is the approximate square root of the
    ! machine precision.
    INTEGER i,j
    REAL h,temp,f(NMAX)
    do j = 1, n
        temp = x(j)
        h = EPS*abs(temp)
        if(h.eq.0.) h = EPS
        x(j) = temp+h ! Trick to reduce finite precision error.
        h    = x(j)-temp
        call funcv(n,x,f)
        x(j) = temp
        do i = 1, n ! Forward difference formula.
            df(i,j) = (f(i)-fvec(i))/h
        enddo
    enddo
    return
END SUBROUTINE