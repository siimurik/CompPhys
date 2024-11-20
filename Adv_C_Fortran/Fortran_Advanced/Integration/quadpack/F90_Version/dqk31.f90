    subroutine dqk31(f,a,b,result,abserr,resabs,resasc)
!***begin prologue  dqk31
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  31-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b) with error
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
!                       result is computed by applying the 31-point
!                       gauss-kronrod rule (resk), obtained by optimal
!                       addition of abscissae to the 15-point gauss
!                       rule (resg).

!              abserr - double precison
!                       estimate of the modulus of the modulus,
!                       which should not exceed abs(i-result)

!              resabs - double precision
!                       approximation to the integral j

!              resasc - double precision
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)

!***references  (none)
!***routines called  d1mach
!***end prologue  dqk31
    double precision :: a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1, &
    d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
    integer :: j,jtw,jtwm1
    external f

    dimension fv1(15),fv2(15),xgk(16),wgk(16),wg(8)

!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.

!           xgk    - abscissae of the 31-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 15-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 15-point gauss rule

!           wgk    - weights of the 31-point kronrod rule

!           wg     - weights of the 15-point gauss rule


! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    data wg  (  1) / 0.030753241996117268354628393577204d0 /
    data wg  (  2) / 0.070366047488108124709267416450667d0 /
    data wg  (  3) / 0.107159220467171935011869546685869d0 /
    data wg  (  4) / 0.139570677926154314447804794511028d0 /
    data wg  (  5) / 0.166269205816993933553200860481209d0 /
    data wg  (  6) / 0.186161000015562211026800561866423d0 /
    data wg  (  7) / 0.198431485327111576456118326443839d0 /
    data wg  (  8) / 0.202578241925561272880620199967519d0 /

    data xgk (  1) / 0.998002298693397060285172840152271d0 /
    data xgk (  2) / 0.987992518020485428489565718586613d0 /
    data xgk (  3) / 0.967739075679139134257347978784337d0 /
    data xgk (  4) / 0.937273392400705904307758947710209d0 /
    data xgk (  5) / 0.897264532344081900882509656454496d0 /
    data xgk (  6) / 0.848206583410427216200648320774217d0 /
    data xgk (  7) / 0.790418501442465932967649294817947d0 /
    data xgk (  8) / 0.724417731360170047416186054613938d0 /
    data xgk (  9) / 0.650996741297416970533735895313275d0 /
    data xgk ( 10) / 0.570972172608538847537226737253911d0 /
    data xgk ( 11) / 0.485081863640239680693655740232351d0 /
    data xgk ( 12) / 0.394151347077563369897207370981045d0 /
    data xgk ( 13) / 0.299180007153168812166780024266389d0 /
    data xgk ( 14) / 0.201194093997434522300628303394596d0 /
    data xgk ( 15) / 0.101142066918717499027074231447392d0 /
    data xgk ( 16) / 0.000000000000000000000000000000000d0 /

    data wgk (  1) / 0.005377479872923348987792051430128d0 /
    data wgk (  2) / 0.015007947329316122538374763075807d0 /
    data wgk (  3) / 0.025460847326715320186874001019653d0 /
    data wgk (  4) / 0.035346360791375846222037948478360d0 /
    data wgk (  5) / 0.044589751324764876608227299373280d0 /
    data wgk (  6) / 0.053481524690928087265343147239430d0 /
    data wgk (  7) / 0.062009567800670640285139230960803d0 /
    data wgk (  8) / 0.069854121318728258709520077099147d0 /
    data wgk (  9) / 0.076849680757720378894432777482659d0 /
    data wgk ( 10) / 0.083080502823133021038289247286104d0 /
    data wgk ( 11) / 0.088564443056211770647275443693774d0 /
    data wgk ( 12) / 0.093126598170825321225486872747346d0 /
    data wgk ( 13) / 0.096642726983623678505179907627589d0 /
    data wgk ( 14) / 0.099173598721791959332393173484603d0 /
    data wgk ( 15) / 0.100769845523875595044946662617570d0 /
    data wgk ( 16) / 0.101330007014791549017374792767493d0 /


!           list of major variables
!           -----------------------
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 15-point gauss formula
!           resk   - result of the 31-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)

!           machine dependent constants
!           ---------------------------
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!***first executable statement  dqk31
    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 31-point kronrod approximation to
!           the integral, and estimate the absolute error.

    fc = f(centr)
    resg = wg(8)*fc
    resk = wgk(16)*fc
    resabs = dabs(resk)
    do 10 j=1,7
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
    do 15 j = 1,8
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
    resasc = wgk(16)*dabs(fc-reskh)
    do 20 j=1,15
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
    end subroutine dqk31
