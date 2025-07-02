SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
    INTEGER n
    LOGICAL check
    REAL f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
    PARAMETER (ALF=1.e-4,TOLX=1.e-7)
    EXTERNAL func
    ! USES func
    ! Given an n-dimensional point xold(1:n), the value of the function and gradient there,
    ! fold and g(1:n), and a direction p(1:n), finds a new point x(1:n) along the direction
    ! p from xold where the function func has decreased “sufficiently.” The new function value
    ! is returned in f. stpmax is an input quantity that limits the length of the steps so that you
    ! do not try to evaluate the function in regions where it is undefined or subject to overflow.
    ! p is usually the Newton direction. The output quantity check is false on a normal exit.
    ! It is true when x is too close to xold. In a minimization algorithm, this usually signals
    ! convergence and can be ignored. However, in a zero-finding algorithm the calling program
    ! should check whether the convergence is spurious.
    ! Parameters: ALF ensures sufficient decrease in function value; TOLX is the convergence
    ! criterion on ∆x.
    INTEGER i
    REAL a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam
    check=.false.
    sum=0.
    do i=1,n
        sum=sum+p(i)*p(i)
    enddo
    sum=sqrt(sum)
    if(sum.gt.stpmax) then !! Scale if attempted step is too big.
        do i=1,n
            p(i)=p(i)*stpmax/sum
        enddo
    endif
    slope=0.
    do i=1,n
        slope=slope+g(i)*p(i)
    enddo
    if(slope.ge.0.) stop 'roundoff problem in lnsrch'
    test=0. !! Compute λmin.
    do i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.)
        if (temp.gt.test) test=temp
    enddo
    alamin=TOLX/test
    alam=1.  !! Always try full Newton step first.
1   continue !! Start of iteration loop.
        do i=1,n
            x(i)=xold(i)+alam*p(i)
        enddo
        f=func(x)
        if(alam.lt.alamin) then !! Convergence on ∆x. For zero finding,
            do i=1,n !! the calling program should verify the convergence.
                x(i)=xold(i)
            enddo
            check=.true.
            return
        else if(f.le.fold+ALF*alam*slope) then !! Sufficient function decrease.
            return
        else !! Backtrack.
            if(alam.eq.1.)then !! First time.
                tmplam=-slope/(2.*(f-fold-slope))
            else !! Subsequent backtracks.
            rhs1=f-fold-alam*slope
            rhs2=f2-fold-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.)then
                tmplam=-slope/(2.*b)
            else
                disc=b*b-3.*a*slope
                if(disc.lt.0.)then
                    tmplam=.5*alam
                else if(b.le.0.)then
                    tmplam=(-b+sqrt(disc))/(3.*a)
                else
                    tmplam=-slope/(b+sqrt(disc))
                endif
            endif
            if(tmplam.gt..5*alam)tmplam=.5*alam !! λ ≤ 0.5λ1.
        endif
    endif
    alam2=alam
    f2=f
    alam=max(tmplam,.1*alam) !! λ ≥ 0.1λ1.
    goto 1 !! Try again.
END SUBROUTINE