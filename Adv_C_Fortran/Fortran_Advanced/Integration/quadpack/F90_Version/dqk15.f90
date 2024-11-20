    subroutine dqk15(f,a,b,result,abserr,resabs,resasc)
!***begin prologue  dqk15
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  15-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description

!           integration rules
!           standard fortran subroutine
!           double precision version

!           parameters
!            on entry
!              f      - double precision
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.

!              a      - double precision
!                       lower limit of integration

!              b      - double precision
!                       upper limit of integration

!            on return
!              result - double precision
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the7-point gauss rule(resg).

!              abserr - double precision
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)

!              resabs - double precision
!                       approximation to the integral j

!              resasc - double precision
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)

!***references  (none)
!***routines called  d1mach
!***end prologue  dqk15

    double precision :: a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1, &
    d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
    integer :: j,jtw,jtwm1
    external f

    dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)

!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.

!           xgk    - abscissae of the 15-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point gauss rule

!           wgk    - weights of the 15-point kronrod rule

!           wg     - weights of the 7-point gauss rule


! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    data wg  (  1) / 0.129484966168869693270611432679082d0 /
    data wg  (  2) / 0.279705391489276667901467771423780d0 /
    data wg  (  3) / 0.381830050505118944950369775488975d0 /
    data wg  (  4) / 0.417959183673469387755102040816327d0 /

    data xgk (  1) / 0.991455371120812639206854697526329d0 /
    data xgk (  2) / 0.949107912342758524526189684047851d0 /
    data xgk (  3) / 0.864864423359769072789712788640926d0 /
    data xgk (  4) / 0.741531185599394439863864773280788d0 /
    data xgk (  5) / 0.586087235467691130294144838258730d0 /
    data xgk (  6) / 0.405845151377397166906606412076961d0 /
    data xgk (  7) / 0.207784955007898467600689403773245d0 /
    data xgk (  8) / 0.000000000000000000000000000000000d0 /

    data wgk (  1) / 0.022935322010529224963732008058970d0 /
    data wgk (  2) / 0.063092092629978553290700663189204d0 /
    data wgk (  3) / 0.104790010322250183839876322541518d0 /
    data wgk (  4) / 0.140653259715525918745189590510238d0 /
    data wgk (  5) / 0.169004726639267902826583426598550d0 /
    data wgk (  6) / 0.190350578064785409913256402421014d0 /
    data wgk (  7) / 0.204432940075298892414161999234649d0 /
    data wgk (  8) / 0.209482141084727828012999174891714d0 /


!           list of major variables
!           -----------------------

!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 7-point gauss formula
!           resk   - result of the 15-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)

!           machine dependent constants
!           ---------------------------

!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.

!***first executable statement  dqk15
    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 15-point kronrod approximation to
!           the integral, and estimate the absolute error.

    fc = f(centr)
    resg = fc*wg(4)
    resk = fc*wgk(8)
    resabs = dabs(resk)
    do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    10 END DO
    do 15 j = 1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    15 END DO
    reskh = resk*0.5d+00
    resasc = wgk(8)*dabs(fc-reskh)
    do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    20 END DO
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc /= 0.0d+00 .AND. abserr /= 0.0d+00) &
    abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
    if(resabs > uflow/(0.5d+02*epmach)) abserr = dmax1 &
    ((epmach*0.5d+02)*resabs,abserr)
    return
    end subroutine dqk15
