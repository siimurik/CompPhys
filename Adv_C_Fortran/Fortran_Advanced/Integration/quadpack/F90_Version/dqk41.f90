    subroutine dqk41(f,a,b,result,abserr,resabs,resasc)
!***begin prologue  dqk41
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  41-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
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
!                       result is computed by applying the 41-point
!                       gauss-kronrod rule (resk) obtained by optimal
!                       addition of abscissae to the 20-point gauss
!                       rule (resg).

!              abserr - double precision
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)

!              resabs - double precision
!                       approximation to the integral j

!              resasc - double precision
!                       approximation to the integal of abs(f-i/(b-a))
!                       over (a,b)

!***references  (none)
!***routines called  d1mach
!***end prologue  dqk41

    double precision :: a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1, &
    d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
    integer :: j,jtw,jtwm1
    external f

    dimension fv1(20),fv2(20),xgk(21),wgk(21),wg(10)

!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.

!           xgk    - abscissae of the 41-point gauss-kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 20-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 20-point gauss rule

!           wgk    - weights of the 41-point gauss-kronrod rule

!           wg     - weights of the 20-point gauss rule


! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    data wg  (  1) / 0.017614007139152118311861962351853d0 /
    data wg  (  2) / 0.040601429800386941331039952274932d0 /
    data wg  (  3) / 0.062672048334109063569506535187042d0 /
    data wg  (  4) / 0.083276741576704748724758143222046d0 /
    data wg  (  5) / 0.101930119817240435036750135480350d0 /
    data wg  (  6) / 0.118194531961518417312377377711382d0 /
    data wg  (  7) / 0.131688638449176626898494499748163d0 /
    data wg  (  8) / 0.142096109318382051329298325067165d0 /
    data wg  (  9) / 0.149172986472603746787828737001969d0 /
    data wg  ( 10) / 0.152753387130725850698084331955098d0 /

    data xgk (  1) / 0.998859031588277663838315576545863d0 /
    data xgk (  2) / 0.993128599185094924786122388471320d0 /
    data xgk (  3) / 0.981507877450250259193342994720217d0 /
    data xgk (  4) / 0.963971927277913791267666131197277d0 /
    data xgk (  5) / 0.940822633831754753519982722212443d0 /
    data xgk (  6) / 0.912234428251325905867752441203298d0 /
    data xgk (  7) / 0.878276811252281976077442995113078d0 /
    data xgk (  8) / 0.839116971822218823394529061701521d0 /
    data xgk (  9) / 0.795041428837551198350638833272788d0 /
    data xgk ( 10) / 0.746331906460150792614305070355642d0 /
    data xgk ( 11) / 0.693237656334751384805490711845932d0 /
    data xgk ( 12) / 0.636053680726515025452836696226286d0 /
    data xgk ( 13) / 0.575140446819710315342946036586425d0 /
    data xgk ( 14) / 0.510867001950827098004364050955251d0 /
    data xgk ( 15) / 0.443593175238725103199992213492640d0 /
    data xgk ( 16) / 0.373706088715419560672548177024927d0 /
    data xgk ( 17) / 0.301627868114913004320555356858592d0 /
    data xgk ( 18) / 0.227785851141645078080496195368575d0 /
    data xgk ( 19) / 0.152605465240922675505220241022678d0 /
    data xgk ( 20) / 0.076526521133497333754640409398838d0 /
    data xgk ( 21) / 0.000000000000000000000000000000000d0 /

    data wgk (  1) / 0.003073583718520531501218293246031d0 /
    data wgk (  2) / 0.008600269855642942198661787950102d0 /
    data wgk (  3) / 0.014626169256971252983787960308868d0 /
    data wgk (  4) / 0.020388373461266523598010231432755d0 /
    data wgk (  5) / 0.025882133604951158834505067096153d0 /
    data wgk (  6) / 0.031287306777032798958543119323801d0 /
    data wgk (  7) / 0.036600169758200798030557240707211d0 /
    data wgk (  8) / 0.041668873327973686263788305936895d0 /
    data wgk (  9) / 0.046434821867497674720231880926108d0 /
    data wgk ( 10) / 0.050944573923728691932707670050345d0 /
    data wgk ( 11) / 0.055195105348285994744832372419777d0 /
    data wgk ( 12) / 0.059111400880639572374967220648594d0 /
    data wgk ( 13) / 0.062653237554781168025870122174255d0 /
    data wgk ( 14) / 0.065834597133618422111563556969398d0 /
    data wgk ( 15) / 0.068648672928521619345623411885368d0 /
    data wgk ( 16) / 0.071054423553444068305790361723210d0 /
    data wgk ( 17) / 0.073030690332786667495189417658913d0 /
    data wgk ( 18) / 0.074582875400499188986581418362488d0 /
    data wgk ( 19) / 0.075704497684556674659542775376617d0 /
    data wgk ( 20) / 0.076377867672080736705502835038061d0 /
    data wgk ( 21) / 0.076600711917999656445049901530102d0 /


!           list of major variables
!           -----------------------

!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 20-point gauss formula
!           resk   - result of the 41-point kronrod formula
!           reskh  - approximation to mean value of f over (a,b), i.e.
!                    to i/(b-a)

!           machine dependent constants
!           ---------------------------

!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.

!***first executable statement  dqk41
    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 41-point gauss-kronrod approximation to
!           the integral, and estimate the absolute error.

    resg = 0.0d+00
    fc = f(centr)
    resk = wgk(21)*fc
    resabs = dabs(resk)
    do 10 j=1,10
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
    do 15 j = 1,10
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
    resasc = wgk(21)*dabs(fc-reskh)
    do 20 j=1,20
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    20 END DO
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc /= 0.0d+00 .AND. abserr /= 0.d+00) &
    abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
    if(resabs > uflow/(0.5d+02*epmach)) abserr = dmax1 &
    ((epmach*0.5d+02)*resabs,abserr)
    return
    end subroutine dqk41
