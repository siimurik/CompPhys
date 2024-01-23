PROGRAM test_erf
    REAL, EXTERNAL :: erf
    REAL :: x, y

    x = 0.5
    y = erf(x)

    print *, y

END PROGRAM test_erf

FUNCTION erf(x)
    REAL erf,x
    ! USES gammp
    !Returns the error function erf(x).
    REAL gammp
    if(x.lt.0.)then
    erf=-gammp(.5, x**2)
    else
    erf=gammp(.5, x**2)
    endif
    return
END FUNCTION erf

FUNCTION gammln(xx)
    REAL gammln,xx
    !Returns the value ln[Γ(xx)] for xx > 0.
    INTEGER j
    DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
    !Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
    !accuracy is good enough.
    SAVE cof,stp
    DATA cof,stp/   76.18009172947146d0,-86.50532032941677d0,   &
                    24.01409824083091d0,-1.231739572450155d0,   &
                    .1208650973866179d-2, -.5395239384953d-5,   &
                    2.5066282746310005d0/
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
    end do
    gammln=tmp+log(stp*ser/x)
    return
END FUNCTION gammln

SUBROUTINE gser(gamser,a,x,gln)
    INTEGER, PARAMETER :: ITMAX=100
    REAL a,gamser,gln,x
    real, PARAMETER :: EPS=3.e-7
    !USES gammln
    !Returns the incomplete gamma function P (a, x) evaluated by its series representation as
    !gamser. Also returns ln Γ(a) as gln.
    INTEGER n
    REAL ap,del,sum,gammln
    gln=gammln(a)
    if(x.le.0.)then
        if(x.lt.0.)stop !'x < 0 in gser'
        gamser=0.
        return
    end if
    ap=a
    sum=1./a
    del=sum
    do n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
    end do
    stop 'a too large, ITMAX too small in gser'
1   gamser=sum*exp(-x+a*log(x)-gln)
    return
END SUBROUTINE gser

SUBROUTINE gcf(gammcf,a,x,gln)
    INTEGER, PARAMETER :: ITMAX=100
    REAL a,gammcf,gln,x
    REAL,PARAMETER :: EPS=3.e-7,FPMIN=1.e-30
    !USES gammln
    !Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction repre-
    !sentation as gammcf. Also returns ln Γ(a) as gln.
    !Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accu-
    !racy; FPMIN is a number near the smallest representable floating-point number.
    INTEGER i
    REAL an,b,c,d,del,h,gammln
    gln=gammln(a)
    b=x+1.-a
    !Set up for evaluating continued fraction by modified
    !Lentz's method (§5.2) with b0 = 0.
    c=1./FPMIN
    d=1./b
    h=d
    do i=1,ITMAX
        !Iterate to convergence.
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS) goto 1
    end do
    stop 'a too large, ITMAX too small in gcf'
1   gammcf=exp(-x+a*log(x)-gln)*h
    !Put factors in front.
    return
END SUBROUTINE gcf

FUNCTION gammp(a,x)
    REAL a,gammp,x
    !USES gcf,gser
    !Returns the incomplete gamma function P (a, x).
    REAL gammcf,gamser,gln
    if(x.lt.0..or.a.le.0.)stop 'bad arguments in gammp'
    if(x.lt.a+1.)then
        !Use the series representation.
        call gser(gamser,a,x,gln)
        gammp=gamser
        else
        !Use the continued fraction representation
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
        !and take its complement.
    end if
    return
END FUNCTION gammp