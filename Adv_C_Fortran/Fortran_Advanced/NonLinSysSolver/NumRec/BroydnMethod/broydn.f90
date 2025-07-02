SUBROUTINE broydn(x, n, check, eval_count)
    IMPLICIT NONE
    INTEGER n,nn,NP,MAXITS, eval_count
    REAL x(n),fvec,EPS,TOLF,TOLMIN,TOLX,STPMX
    LOGICAL check
    PARAMETER (NP=40,MAXITS=200,EPS=1.e-7,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=EPS,STPMX=100.)
    COMMON /newtv/ fvec(NP),nn ! Communicates with fmin.
    SAVE /newtv/
    ! USES fdjac,fmin,lnsrch,qrdcmp,qrupdt,rsolv
    ! Given an initial guess x(1:n) for a root in n dimensions, find the root by Broyden's method
    ! embedded in a globally convergent strategy. The vector of functions to be zeroed, called
    ! fvec(1:n) in the routine below, is returned by a user-supplied subroutine that must be
    ! called funcv and have the declaration subroutine funcv(n,x,fvec). The subroutine
    ! fdjac and the function fmin from newt are used. The output quantity check is false on
    ! a normal return and true if the routine has converged to a local minimum of the function
    ! fmin or if Broyden's method can make no further progress. In this case try restarting from
    ! a different initial guess.
    ! Parameters: NP is the maximum expected value of n; MAXITS is the maximum number of
    ! iterations; EPS is close to the machine precision; TOLF sets the convergence criterion on
    ! function values; TOLMIN sets the criterion for deciding whether spurious convergence to a
    ! minimum of fmin has occurred; TOLX is the convergence criterion on δx; STPMX is the
    ! scaled maximum step length allowed in line searches.
    INTEGER i,its,j,k
    REAL den,f,fold,stpmax,sum,temp,test,c(NP),d(NP),fvcold(NP), &
                g(NP),p(NP),qt(NP,NP),r(NP,NP),s(NP),t(NP),w(NP),&
                xold(NP),fmin
    LOGICAL restrt,sing,skip
    EXTERNAL fmin
    nn = n
    f = fmin(x)   !! The vector fvec is also computed by this call.
    test = 0.     !! Test for initial guess being a root. Use more strin-
    do i = 1,n    !! gent test than simply TOLF.
        if(abs(fvec(i)).gt.test) test=abs(fvec(i))
    enddo
    if (test .lt. .01*TOLF) then
        check=.false.
        return
    endif
    sum = 0.  !! Calculate stpmax for line searches.
    do i = 1,n
        sum = sum+x(i)**2
    enddo
    stpmax = STPMX*max(sqrt(sum),float(n))
    restrt=.true.       !! Ensure initial Jacobian gets computed.
    eval_count = 0
    do its = 1, MAXITS  !! Start of iteration loop.
        if(restrt)then
            call fdjac(n,x,fvec,NP,r)    !! Initialize or reinitialize Jacobian in r.
            call qrdcmp(r,n,NP,c,d,sing) !! QR decomposition of Jacobian.
            if(sing) stop 'singular Jacobian in broydn'
                do i = 1, n                  !! Form QT explicitly.
                    do j = 1, n
                        qt(i,j) = 0.
                    enddo
                    qt(i,i) = 1.
                enddo
                do k = 1, n-1
                    if(c(k) .ne. 0.) then
                        do j=1,n
                            sum=0.
                            do  i=k,n
                                sum=sum+r(i,k)*qt(i,j)
                            enddo
                            sum=sum/c(k)
                            do i = k, n
                                qt(i,j)=qt(i,j)-sum*r(i,k)
                            enddo
                        enddo
                    endif
                enddo
                do i = 1, n !! Form R explicitly.
                    r(i,i) = d(i)
                    do j = 1, i-1
                        r(i,j) = 0.
                    enddo
                enddo
            else        !! Carry out Broyden update.
            do i = 1, n !! s = δx.
                s(i) = x(i)-xold(i)
            enddo
            do i = 1, n  !! t = R · s.
                sum = 0.
                do j = i, n
                    sum = sum + r(i,j)*s(j)
                enddo
                t(i) = sum
            enddo
            skip = .true.
            do i = 1, n !! w = δF − B · s.
                sum = 0.
                do j = 1, n
                    sum = sum + qt(j,i)*t(j)
                enddo
                w(i) = fvec(i)-fvcold(i)-sum
                if (abs(w(i)).ge.EPS*(abs(fvec(i))+abs(fvcold(i)))) then
                    skip = .false. !! Don't update with noisy components of w.
                else
                    w(i) = 0.
                endif
            enddo
            if(.not.skip)then
                do i = 1, n !! t = QT · w.
                    sum=0.
                    do j = 1, n
                        sum = sum + qt(i,j)*w(j)
                    enddo
                    t(i)=sum
                enddo
                den = 0.
                do i = 1, n
                    den = den + s(i)**2
                enddo
                do i = 1, n !! Store s/(s · s) in s.
                    s(i) = s(i)/den
                enddo
                call qrupdt(r,qt,n,NP,t,s) !! Update R and QT.
                do i = 1, n
                    if(r(i,i).eq.0.) stop 'r singular in broydn'
                    d(i)=r(i,i) !! Diagonal of R stored in d.
                enddo
            endif
        endif
        do i = 1, n !! Compute ∇f ≈ (Q·R)T·F for the line search.
            sum = 0.
            do j = 1, n
                sum = sum + qt(i,j)*fvec(j)
            enddo
            g(i) = sum
        enddo
        do i = n, 1, -1
            sum=0.
            do j = 1, i
                sum = sum + r(j,i)*g(j)
            enddo
            g(i)=sum
        enddo
        do i = 1, n !! Store x and F.
            xold(i) = x(i)
            fvcold(i) = fvec(i)
        enddo
        fold = f !! Store f .
        do i = 1, n !! Right-hand side for linear equations is −QT · F.
            sum = 0.
            do j = 1, n
                sum = sum + qt(i,j)*fvec(j)
            enddo
            p(i)=-sum
        enddo
        call rsolv(r,n,NP,d,p) !! Solve linear equations.
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)   !! lnsrch returns new x and f. 
                                !! It also calculates fvec at the new x when it calls fmin.
        test = 0. !! Test for convergence on function values.
        do i = 1, n
            if(abs(fvec(i)) .gt. test) test = abs(fvec(i))
        enddo
        if (test .lt. TOLF) then
            check = .false.
        return
        endif
        if (check) then  !! True if line search failed to find a new x.
            if (restrt) then !! Failure; already tried reinitializing the Jacobian.
                return
            else          !! Check for gradient of f zero, i.e., spurious con-
                test = 0. !! vergence.
                den = max(f,.5*n)
                do i = 1, n
                    temp = abs(g(i))*max(abs(x(i)),1.)/den
                    if (temp .gt. test) test = temp
                enddo
                if(test.lt.TOLMIN)then
                        return
                else !! Try reinitializing the Jacobian.
                    restrt=.true.
                endif
            endif
        else !! Successful step; will use Broyden update for next
            restrt=.false. !! step.
            test=0. !! Test for convergence on δx.
            do i = 1, n
                temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
                if(temp.gt.test) test = temp
            enddo
            if (test.lt.TOLX) return
        endif
        eval_count = eval_count + 1
    enddo
    stop 'MAXITS exceeded in broydn'
END SUBROUTINE
