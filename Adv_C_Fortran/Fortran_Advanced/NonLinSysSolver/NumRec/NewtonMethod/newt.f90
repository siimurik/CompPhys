SUBROUTINE newt(x, n, check, eval_count)
    IMPLICIT NONE
    INTEGER n,nn,NP,MAXITS,eval_count
    LOGICAL check
    REAL x(n),fvec,TOLF,TOLMIN,TOLX,STPMX
    PARAMETER (NP=40,MAXITS=200,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7,STPMX=100.)
    COMMON /newtv/ fvec(NP),nn
    SAVE /newtv/
    ! USES fdjac,fmin,lnsrch,lubksb,ludcmp
    ! Given an initial guess x(1:n) for a root in n dimensions, find the root by a globally
    ! convergent Newton's method. The vector of functions to be zeroed, called fvec(1:n)
    ! in the routine below, is returned by a user-supplied subroutine that must be called funcv
    ! and have the declaration subroutine funcv(n,x,fvec). The output quantity check
    ! is false on a normal return and true if the routine has converged to a local minimum of the
    ! function fmin defined below. In this case try restarting from a different initial guess.
    ! Parameters: NP is the maximum expected value of n; MAXITS is the maximum number of
    ! iterations; TOLF sets the convergence criterion on function values; TOLMIN sets the criterion
    ! for deciding whether spurious convergence to a minimum of fmin has occurred; TOLX is
    ! the convergence criterion on δx; STPMX is the scaled maximum step length allowed in line
    ! searches.
    INTEGER i,its,j,indx(NP)
    REAL d,den,f,fold,stpmax,sum,temp,test
    REAL fjac(NP,NP),g(NP),p(NP),xold(NP),fmin
    EXTERNAL fmin
    nn=n
    f=fmin(x)  ! The vector fvec is also computed by this call.
    test=0.    ! Test for initial guess being a root. Use more stringent test than simply TOLF.
    do i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
    enddo
    if(test.lt.0.01*TOLF)then
        check=.false.
        return
    endif
    sum=0.     ! Calculate stpmax for line searches.
    do i=1,n
        sum=sum+x(i)**2
    enddo
    stpmax=STPMX*max(sqrt(sum),float(n))
    eval_count = 0
    do its=1,MAXITS  ! Start of iteration loop.
        call fdjac(n,x,fvec,NP,fjac)
        ! If analytic Jacobian is available, you can replace the routine fdjac below with your own routine.
        do i=1,n     ! Compute ∇f for the line search.
            sum=0.
            do j=1,n
                sum=sum+fjac(j,i)*fvec(j)
            enddo
            g(i)=sum
        enddo
        do i=1,n     ! Store x,
            xold(i)=x(i)
        enddo
        fold=f          ! and f.
        do i=1,n     ! Right-hand side for linear equations.
            p(i)=-fvec(i)
        enddo
        call ludcmp(fjac,n,NP,indx,d)  ! Solve linear equations by LU decomposition.
        call lubksb(fjac,n,NP,indx,p)
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
        ! lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.
        test=0.  ! Test for convergence on function values.
        do i=1,n
            if(abs(fvec(i)).gt.test)test=abs(fvec(i))
        enddo
        if(test.lt.TOLF)then
            check=.false.
            return
        endif
        if(check)then   ! Check for gradient of f zero, i.e., spurious convergence.
            test=0.
            den=max(f,.5*n)
            do i=1,n
                temp=abs(g(i))*max(abs(x(i)),1.)/den
                if(temp.gt.test)test=temp
            enddo
            if(test.lt.TOLMIN)then
                check=.true.
            else
                check=.false.
            endif
            return
        endif
        test=0.  ! Test for convergence on δx.
        do i=1,n
            temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
            if(temp.gt.test)test=temp
        enddo
        if(test.lt.TOLX)return
        eval_count = eval_count + 1
    enddo
    PRINT *, 'MAXITS exceeded in newt'
END SUBROUTINE