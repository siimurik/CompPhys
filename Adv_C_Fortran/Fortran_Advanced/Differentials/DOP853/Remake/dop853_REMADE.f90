module dop853_mod
    implicit none
    private
    public :: dop853, contd8
    
    ! Module-level constants
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    ! Common block replacement
    type :: condo8_type
        real(dp) :: xold, hout
    end type
    type(condo8_type), save :: condo8_data
    
    ! Butcher tableau coefficients as module parameters
    real(dp), parameter :: &
        c2  = 0.526001519587677318785587544488e-01_dp, &
        c3  = 0.789002279381515978178381316732e-01_dp, &
        c4  = 0.118350341907227396726757197510e+00_dp, &
        c5  = 0.281649658092772603273242802490e+00_dp, &
        c6  = 0.333333333333333333333333333333e+00_dp, &
        c7  = 0.25e+00_dp, &
        c8  = 0.307692307692307692307692307692e+00_dp, &
        c9  = 0.651282051282051282051282051282e+00_dp, &
        c10 = 0.6e+00_dp, &
        c11 = 0.857142857142857142857142857142e+00_dp, &
        c14 = 0.1e+00_dp, &
        c15 = 0.2e+00_dp, &
        c16 = 0.777777777777777777777777777778e+00_dp
    
    real(dp), parameter :: &
        b1 =   5.42937341165687622380535766363e-2_dp, &
        b6 =   4.45031289275240888144113950566e0_dp, &
        b7 =   1.89151789931450038304281599044e0_dp, &
        b8 =  -5.8012039600105847814672114227e0_dp, &
        b9 =   3.1116436695781989440891606237e-1_dp, &
        b10 = -1.52160949662516078556178806805e-1_dp, &
        b11 =  2.01365400804030348374776537501e-1_dp, &
        b12 =  4.47106157277725905176885569043e-2_dp
    
    real(dp), parameter :: &
        bhh1 = 0.244094488188976377952755905512e+00_dp, &
        bhh2 = 0.733846688281611857341361741547e+00_dp, &
        bhh3 = 0.220588235294117647058823529412e-01_dp
    
    real(dp), parameter :: &
        er1 =  0.1312004499419488073250102996e-01_dp, &
        er6 = -0.1225156446376204440720569753e+01_dp, &
        er7 = -0.4957589496572501915214079952e+00_dp, &
        er8 =  0.1664377182454986536961530415e+01_dp, &
        er9 = -0.3503288487499736816886487290e+00_dp, &
        er10 =  0.3341791187130174790297318841e+00_dp, &
        er11 =  0.8192320648511571246570742613e-01_dp, &
        er12 = -0.2235530786388629525884427845e-01_dp
    
    ! A-matrix coefficients
    real(dp), parameter :: &
        a21 =    5.26001519587677318785587544488e-2_dp, &
        a31 =    1.97250569845378994544595329183e-2_dp, &
        a32 =    5.91751709536136983633785987549e-2_dp, &
        a41 =    2.95875854768068491816892993775e-2_dp, &
        a43 =    8.87627564304205475450678981324e-2_dp, &
        a51 =    2.41365134159266685502369798665e-1_dp, &
        a53 =   -8.84549479328286085344864962717e-1_dp, &
        a54 =    9.24834003261792003115737966543e-1_dp, &
        a61 =    3.7037037037037037037037037037e-2_dp, &
        a64 =    1.70828608729473871279604482173e-1_dp, &
        a65 =    1.25467687566822425016691814123e-1_dp, &
        a71 =    3.7109375e-2_dp, &
        a74 =    1.70252211019544039314978060272e-1_dp, &
        a75 =    6.02165389804559606850219397283e-2_dp, &
        a76 =   -1.7578125e-2_dp
    
    real(dp), parameter :: &
        a81 =    3.70920001185047927108779319836e-2_dp, &
        a84 =    1.70383925712239993810214054705e-1_dp, &
        a85 =    1.07262030446373284651809199168e-1_dp, &
        a86 =   -1.53194377486244017527936158236e-2_dp, &
        a87 =    8.27378916381402288758473766002e-3_dp, &
        a91 =    6.24110958716075717114429577812e-1_dp, &
        a94 =   -3.36089262944694129406857109825e0_dp, &
        a95 =   -8.68219346841726006818189891453e-1_dp, &
        a96 =    2.75920996994467083049415600797e1_dp, &
        a97 =    2.01540675504778934086186788979e1_dp, &
        a98 =   -4.34898841810699588477366255144e1_dp, &
        a101 =   4.77662536438264365890433908527e-1_dp, &
        a104 =  -2.48811461997166764192642586468e0_dp, &
        a105 =  -5.90290826836842996371446475743e-1_dp, &
        a106 =   2.12300514481811942347288949897e1_dp, &
        a107 =   1.52792336328824235832596922938e1_dp, &
        a108 =  -3.32882109689848629194453265587e1_dp, &
        a109 =  -2.03312017085086261358222928593e-2_dp
    
    real(dp), parameter :: &
        a111 =  -9.3714243008598732571704021658e-1_dp, &
        a114 =   5.18637242884406370830023853209e0_dp, &
        a115 =   1.09143734899672957818500254654e0_dp, &
        a116 =  -8.14978701074692612513997267357e0_dp, &
        a117 =  -1.85200656599969598641566180701e1_dp, &
        a118 =   2.27394870993505042818970056734e1_dp, &
        a119 =   2.49360555267965238987089396762e0_dp, &
        a1110 = -3.0467644718982195003823669022e0_dp, &
        a121 =   2.27331014751653820792359768449e0_dp, &
        a124 =  -1.05344954667372501984066689879e1_dp, &
        a125 =  -2.00087205822486249909675718444e0_dp, &
        a126 =  -1.79589318631187989172765950534e1_dp, &
        a127 =   2.79488845294199600508499808837e1_dp, &
        a128 =  -2.85899827713502369474065508674e0_dp, &
        a129 =  -8.87285693353062954433549289258e0_dp, &
        a1210 =  1.23605671757943030647266201528e1_dp, &
        a1211 =  6.43392746015763530355970484046e-1_dp
    
    real(dp), parameter :: &
        a141 =  5.61675022830479523392909219681e-2_dp, &
        a147 =  2.53500210216624811088794765333e-1_dp, &
        a148 = -2.46239037470802489917441475441e-1_dp, &
        a149 = -1.24191423263816360469010140626e-1_dp, &
        a1410 =  1.5329179827876569731206322685e-1_dp, &
        a1411 =  8.20105229563468988491666602057e-3_dp, &
        a1412 =  7.56789766054569976138603589584e-3_dp, &
        a1413 = -8.298e-3_dp
    
    real(dp), parameter :: &
        a151 =  3.18346481635021405060768473261e-2_dp, &
        a156 =  2.83009096723667755288322961402e-2_dp, &
        a157 =  5.35419883074385676223797384372e-2_dp, &
        a158 = -5.49237485713909884646569340306e-2_dp, &
        a1511 = -1.08347328697249322858509316994e-4_dp, &
        a1512 =  3.82571090835658412954920192323e-4_dp, &
        a1513 = -3.40465008687404560802977114492e-4_dp, &
        a1514 =  1.41312443674632500278074618366e-1_dp, &
        a161 = -4.28896301583791923408573538692e-1_dp, &
        a166 = -4.69762141536116384314449447206e0_dp, &
        a167 =  7.68342119606259904184240953878e0_dp, &
        a168 =  4.06898981839711007970213554331e0_dp, &
        a169 =  3.56727187455281109270669543021e-1_dp, &
        a1613 = -1.39902416515901462129418009734e-3_dp, &
        a1614 =  2.9475147891527723389556272149e0_dp, &
        a1615 = -9.15095847217987001081870187138e0_dp
    
    ! Dense output coefficients
    real(dp), parameter :: &
        d41  = -0.84289382761090128651353491142e+01_dp, &
        d46  =  0.56671495351937776962531783590e+00_dp, &
        d47  = -0.30689499459498916912797304727e+01_dp, &
        d48  =  0.23846676565120698287728149680e+01_dp, &
        d49  =  0.21170345824450282767155149946e+01_dp, &
        d410 = -0.87139158377797299206789907490e+00_dp, &
        d411 =  0.22404374302607882758541771650e+01_dp, &
        d412 =  0.63157877876946881815570249290e+00_dp, &
        d413 = -0.88990336451333310820698117400e-01_dp, &
        d414 =  0.18148505520854727256656404962e+02_dp, &
        d415 = -0.91946323924783554000451984436e+01_dp, &
        d416 = -0.44360363875948939664310572000e+01_dp
    
    real(dp), parameter :: &
        d51  =  0.10427508642579134603413151009e+02_dp, &
        d56  =  0.24228349177525818288430175319e+03_dp, &
        d57  =  0.16520045171727028198505394887e+03_dp, &
        d58  = -0.37454675472269020279518312152e+03_dp, &
        d59  = -0.22113666853125306036270938578e+02_dp, &
        d510 =  0.77334326684722638389603898808e+01_dp, &
        d511 = -0.30674084731089398182061213626e+02_dp, &
        d512 = -0.93321305264302278729567221706e+01_dp, &
        d513 =  0.15697238121770843886131091075e+02_dp, &
        d514 = -0.31139403219565177677282850411e+02_dp, &
        d515 = -0.93529243588444783865713862664e+01_dp, &
        d516 =  0.35816841486394083752465898540e+02_dp
    
    real(dp), parameter :: &
        d61 =  0.19985053242002433820987653617e+02_dp, &
        d66 = -0.38703730874935176555105901742e+03_dp, &
        d67 = -0.18917813819516756882830838328e+03_dp, &
        d68 =  0.52780815920542364900561016686e+03_dp, &
        d69 = -0.11573902539959630126141871134e+02_dp, &
        d610 =  0.68812326946963000169666922661e+01_dp, &
        d611 = -0.10006050966910838403183860980e+01_dp, &
        d612 =  0.77771377980534432092869265740e+00_dp, &
        d613 = -0.27782057523535084065932004339e+01_dp, &
        d614 = -0.60196695231264120758267380846e+02_dp, &
        d615 =  0.84320405506677161018159903784e+02_dp, &
        d616 =  0.11992291136182789328035130030e+02_dp
    
    real(dp), parameter :: &
        d71  = -0.25693933462703749003312586129e+02_dp, &
        d76  = -0.15418974869023643374053993627e+03_dp, &
        d77  = -0.23152937917604549567536039109e+03_dp, &
        d78  =  0.35763911791061412378285349910e+03_dp, &
        d79  =  0.93405324183624310003907691704e+02_dp, &
        d710 = -0.37458323136451633156875139351e+02_dp, &
        d711 =  0.10409964950896230045147246184e+03_dp, &
        d712 =  0.29840293426660503123344363579e+02_dp, &
        d713 = -0.43533456590011143754432175058e+02_dp, &
        d714 =  0.96324553959188282948394950600e+02_dp, &
        d715 = -0.39177261675615439165231486172e+02_dp, &
        d716 = -0.14972683625798562581422125276e+03_dp

contains

    subroutine dop853(n, fcn, x, y, xend, rtol, atol, itol, &
                      solout, iout, work, lwork, iwork, liwork, &
                      rpar, ipar, idid)
        ! ----------------------------------------------------------
        ! NUMERICAL SOLUTION OF A SYSTEM OF FIRST ORDER
        ! ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y).
        ! THIS IS AN EXPLICIT RUNGE-KUTTA METHOD OF ORDER 8(5,3)
        ! DUE TO DORMAND & PRINCE (WITH STEPSIZE CONTROL AND
        ! DENSE OUTPUT)
        !
        ! AUTHORS: E. HAIRER AND G. WANNER
        !          UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
        !          CH-1211 GENEVE 24, SWITZERLAND
        !          E-MAIL:  Ernst.Hairer@unige.ch
        !                   Gerhard.Wanner@unige.ch
        !
        ! MODERNIZED TO FORTRAN 90 STYLE
        ! ----------------------------------------------------------
        
        implicit none
        
        ! Arguments
        integer, intent(in) :: n, itol, iout, lwork, liwork
        real(dp), intent(inout) :: x
        real(dp), intent(in) :: xend
        real(dp), intent(inout) :: y(n)
        real(dp), intent(in) :: rtol(*), atol(*)
        real(dp), intent(inout) :: work(lwork)
        integer, intent(inout) :: iwork(liwork)
        real(dp), intent(inout) :: rpar(*)
        integer, intent(inout) :: ipar(*)
        integer, intent(out) :: idid
        
        ! External procedures
        external :: fcn, solout
        
        ! Local variables
        integer :: nfcn, nstep, naccpt, nrejct
        integer :: iprint, nmax, meth, nstiff, nrdens
        integer :: iek1, iek2, iek3, iek4, iek5, iek6, iek7, iek8, iek9, iek10
        integer :: iey1, ieco, icomp, istore, i
        real(dp) :: uround, safe, fac1, fac2, beta, hmax, h
        logical :: arret
        
        ! Initialize counters
        nfcn = 0
        nstep = 0
        naccpt = 0
        nrejct = 0
        arret = .false.
        
        ! Set print unit
        if (iwork(3) == 0) then
            iprint = 6
        else
            iprint = iwork(3)
        end if
        
        ! Maximum number of steps
        if (iwork(1) == 0) then
            nmax = 100000
        else
            nmax = iwork(1)
            if (nmax <= 0) then
                if (iprint > 0) write(iprint, *) 'WRONG INPUT IWORK(1)=', iwork(1)
                arret = .true.
            end if
        end if
        
        ! Method coefficients
        if (iwork(2) == 0) then
            meth = 1
        else
            meth = iwork(2)
            if (meth <= 0 .or. meth >= 4) then
                if (iprint > 0) write(iprint, *) 'CURIOUS INPUT IWORK(2)=', iwork(2)
                arret = .true.
            end if
        end if
        
        ! Stiffness detection parameter
        nstiff = iwork(4)
        if (nstiff == 0) nstiff = 1000
        if (nstiff < 0) nstiff = nmax + 10
        
        ! Number of dense output components
        nrdens = iwork(5)
        if (nrdens < 0 .or. nrdens > n) then
            if (iprint > 0) write(iprint, *) 'CURIOUS INPUT IWORK(5)=', iwork(5)
            arret = .true.
        else
            if (nrdens > 0 .and. iout < 2) then
                if (iprint > 0) write(iprint, *) 'WARNING: PUT IOUT=2 OR IOUT=3 FOR DENSE OUTPUT'
            end if
            if (nrdens == n) then
                do i = 1, nrdens
                    iwork(i + 20) = i
                end do
            end if
        end if
        
        ! Unit roundoff
        if (work(1) == 0.0_dp) then
            uround = 2.3e-16_dp
        else
            uround = work(1)
            if (uround <= 1.0e-35_dp .or. uround >= 1.0_dp) then
                if (iprint > 0) write(iprint, *) 'WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:', work(1)
                arret = .true.
            end if
        end if
        
        ! Safety factor
        if (work(2) == 0.0_dp) then
            safe = 0.9_dp
        else
            safe = work(2)
            if (safe >= 1.0_dp .or. safe <= 1.0e-4_dp) then
                if (iprint > 0) write(iprint, *) 'CURIOUS INPUT FOR SAFETY FACTOR WORK(2)=', work(2)
                arret = .true.
            end if
        end if
        
        ! Step size selection parameters
        if (work(3) == 0.0_dp) then
            fac1 = 0.333_dp
        else
            fac1 = work(3)
        end if
        
        if (work(4) == 0.0_dp) then
            fac2 = 6.0_dp
        else
            fac2 = work(4)
        end if
        
        ! Beta for step control stabilization
        if (work(5) == 0.0_dp) then
            beta = 0.0_dp
        else
            if (work(5) < 0.0_dp) then
                beta = 0.0_dp
            else
                beta = work(5)
                if (beta > 0.2_dp) then
                    if (iprint > 0) write(iprint, *) 'CURIOUS INPUT FOR BETA: WORK(5)=', work(5)
                    arret = .true.
                end if
            end if
        end if
        
        ! Maximum step size
        if (work(6) == 0.0_dp) then
            hmax = xend - x
        else
            hmax = work(6)
        end if
        
        ! Initial step size
        h = work(7)
        
        ! Prepare entry points for arrays in work
        iek1 = 21
        iek2 = iek1 + n
        iek3 = iek2 + n
        iek4 = iek3 + n
        iek5 = iek4 + n
        iek6 = iek5 + n
        iek7 = iek6 + n
        iek8 = iek7 + n
        iek9 = iek8 + n
        iek10 = iek9 + n
        iey1 = iek10 + n
        ieco = iey1 + n
        
        ! Check storage requirement
        istore = ieco + 8 * nrdens - 1
        if (istore > lwork) then
            if (iprint > 0) write(iprint, *) 'INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=', istore
            arret = .true.
        end if
        
        icomp = 21
        istore = icomp + nrdens - 1
        if (istore > liwork) then
            if (iprint > 0) write(iprint, *) 'INSUFFICIENT STORAGE FOR IWORK, MIN. LIWORK=', istore
            arret = .true.
        end if
        
        ! Return with error if problems occurred
        if (arret) then
            idid = -1
            return
        end if
        
        ! Call core integrator
        call dp86co(n, fcn, x, y, xend, hmax, h, rtol, atol, itol, iprint, &
                    solout, iout, idid, nmax, uround, meth, nstiff, safe, beta, fac1, fac2, &
                    work(iek1), work(iek2), work(iek3), work(iek4), work(iek5), &
                    work(iek6), work(iek7), work(iek8), work(iek9), work(iek10), &
                    work(iey1), work(ieco), iwork(icomp), nrdens, rpar, ipar, &
                    nfcn, nstep, naccpt, nrejct)
        
        ! Store results
        work(7) = h
        iwork(17) = nfcn
        iwork(18) = nstep
        iwork(19) = naccpt
        iwork(20) = nrejct
        
    end subroutine dop853

    subroutine dp86co(n, fcn, x, y, xend, hmax, h, rtol, atol, itol, iprint, &
                      solout, iout, idid, nmax, uround, meth, nstiff, safe, beta, fac1, fac2, &
                      k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, y1, cont, icomp, nrd, rpar, ipar, &
                      nfcn, nstep, naccpt, nrejct)
        ! ----------------------------------------------------------
        ! CORE INTEGRATOR FOR DOP853
        ! PARAMETERS SAME AS IN DOP853 WITH WORKSPACE ADDED
        ! FIXED: Initialized NONSTI variable to prevent uninitialized access
        ! ----------------------------------------------------------
        
        implicit none
        
        ! Arguments
        integer, intent(in) :: n, itol, iprint, iout, nmax, meth, nstiff, nrd
        real(dp), intent(inout) :: x, h, hmax
        real(dp), intent(in) :: xend, uround, safe, beta, fac1, fac2
        real(dp), intent(inout) :: y(n), y1(n)
        real(dp), intent(inout) :: k1(n), k2(n), k3(n), k4(n), k5(n), k6(n)
        real(dp), intent(inout) :: k7(n), k8(n), k9(n), k10(n)
        real(dp), intent(in) :: atol(*), rtol(*)
        real(dp), intent(inout) :: cont(8*nrd)
        integer, intent(in) :: icomp(nrd)
        real(dp), intent(inout) :: rpar(*)
        integer, intent(inout) :: ipar(*)
        integer, intent(inout) :: nfcn, nstep, naccpt, nrejct
        integer, intent(out) :: idid
        
        ! External procedures
        external :: fcn, solout
        
        ! Local variables
        real(dp) :: facold, expo1, facc1, facc2, posneg
        real(dp) :: atoli, rtoli, hlamb, xold, hout, xout
        real(dp) :: err, err2, deno, fac11, fac, hnew, xph
        real(dp) :: sk, erri, stnum, stden, ydiff, bspl
        integer :: iord, irtrn, iasti, nonsti  ! FIX: Added nonsti declaration
        integer :: i, j
        logical :: reject, last, event
        
        ! Initialize variables (FIX: Initialize nonsti to prevent uninitialized access)
        facold = 1.0e-4_dp
        expo1 = 1.0_dp/8.0_dp - beta * 0.2_dp
        facc1 = 1.0_dp / fac1
        facc2 = 1.0_dp / fac2
        posneg = sign(1.0_dp, xend - x)
        
        ! Initial preparations
        atoli = atol(1)
        rtoli = rtol(1)
        last = .false.
        hlamb = 0.0_dp
        iasti = 0
        nonsti = 0  ! FIX: Initialize nonsti to prevent uninitialized variable access
        
        call fcn(n, x, y, k1, rpar, ipar)
        hmax = abs(hmax)
        iord = 8
        
        if (h == 0.0_dp) then
            h = hinit(n, fcn, x, y, xend, posneg, k1, k2, k3, iord, &
                      hmax, atol, rtol, itol, rpar, ipar)
        end if
        
        nfcn = nfcn + 2
        reject = .false.
        xold = x
        
        if (iout /= 0) then
            irtrn = 1
            hout = 1.0_dp
            condo8_data%xold = xold
            condo8_data%hout = hout
            call solout(naccpt + 1, xold, x, y, n, cont, icomp, nrd, &
                        rpar, ipar, irtrn, xout)
            if (irtrn < 0) goto 79
        end if
        
        ! Basic integration step
        integration_loop: do
            if (nstep > nmax) goto 78
            if (0.1_dp * abs(h) <= abs(x) * uround) goto 77
            
            if ((x + 1.01_dp * h - xend) * posneg > 0.0_dp) then
                h = xend - x
                last = .true.
            end if
            
            nstep = nstep + 1
            
            ! The twelve stages
            if (irtrn >= 2) then
                call fcn(n, x, y, k1, rpar, ipar)
            end if
            
            ! Stage 2
            do i = 1, n
                y1(i) = y(i) + h * a21 * k1(i)
            end do
            call fcn(n, x + c2 * h, y1, k2, rpar, ipar)
            
            ! Stage 3
            do i = 1, n
                y1(i) = y(i) + h * (a31 * k1(i) + a32 * k2(i))
            end do
            call fcn(n, x + c3 * h, y1, k3, rpar, ipar)
            
            ! Stage 4
            do i = 1, n
                y1(i) = y(i) + h * (a41 * k1(i) + a43 * k3(i))
            end do
            call fcn(n, x + c4 * h, y1, k4, rpar, ipar)
            
            ! Stage 5
            do i = 1, n
                y1(i) = y(i) + h * (a51 * k1(i) + a53 * k3(i) + a54 * k4(i))
            end do
            call fcn(n, x + c5 * h, y1, k5, rpar, ipar)
            
            ! Stage 6
            do i = 1, n
                y1(i) = y(i) + h * (a61 * k1(i) + a64 * k4(i) + a65 * k5(i))
            end do
            call fcn(n, x + c6 * h, y1, k6, rpar, ipar)
            
            ! Stage 7
            do i = 1, n
                y1(i) = y(i) + h * (a71 * k1(i) + a74 * k4(i) + a75 * k5(i) + a76 * k6(i))
            end do
            call fcn(n, x + c7 * h, y1, k7, rpar, ipar)
            
            ! Stage 8
            do i = 1, n
                y1(i) = y(i) + h * (a81 * k1(i) + a84 * k4(i) + a85 * k5(i) + a86 * k6(i) + a87 * k7(i))
            end do
            call fcn(n, x + c8 * h, y1, k8, rpar, ipar)
            
            ! Stage 9
            do i = 1, n
                y1(i) = y(i) + h * (a91 * k1(i) + a94 * k4(i) + a95 * k5(i) + a96 * k6(i) + &
                                   a97 * k7(i) + a98 * k8(i))
            end do
            call fcn(n, x + c9 * h, y1, k9, rpar, ipar)
            
            ! Stage 10
            do i = 1, n
                y1(i) = y(i) + h * (a101 * k1(i) + a104 * k4(i) + a105 * k5(i) + a106 * k6(i) + &
                                   a107 * k7(i) + a108 * k8(i) + a109 * k9(i))
            end do
            call fcn(n, x + c10 * h, y1, k10, rpar, ipar)
            
            ! Stage 11
            do i = 1, n
                y1(i) = y(i) + h * (a111 * k1(i) + a114 * k4(i) + a115 * k5(i) + a116 * k6(i) + &
                                   a117 * k7(i) + a118 * k8(i) + a119 * k9(i) + a1110 * k10(i))
            end do
            call fcn(n, x + c11 * h, y1, k2, rpar, ipar)
            
            ! Stage 12
            xph = x + h
            do i = 1, n
                y1(i) = y(i) + h * (a121 * k1(i) + a124 * k4(i) + a125 * k5(i) + a126 * k6(i) + &
                                   a127 * k7(i) + a128 * k8(i) + a129 * k9(i) + a1210 * k10(i) + &
                                   a1211 * k2(i))
            end do
            call fcn(n, xph, y1, k3, rpar, ipar)
            nfcn = nfcn + 11
            
            ! Compute solution and embedded solution
            do i = 1, n
                k4(i) = b1 * k1(i) + b6 * k6(i) + b7 * k7(i) + b8 * k8(i) + b9 * k9(i) + &
                        b10 * k10(i) + b11 * k2(i) + b12 * k3(i)
                k5(i) = y(i) + h * k4(i)
            end do
            
            ! Error estimation
            err = 0.0_dp
            err2 = 0.0_dp
            
            if (itol == 0) then
                do i = 1, n
                    sk = atoli + rtoli * max(abs(y(i)), abs(k5(i)))
                    erri = k4(i) - bhh1 * k1(i) - bhh2 * k9(i) - bhh3 * k3(i)
                    err2 = err2 + (erri / sk)**2
                    erri = er1 * k1(i) + er6 * k6(i) + er7 * k7(i) + er8 * k8(i) + er9 * k9(i) + &
                           er10 * k10(i) + er11 * k2(i) + er12 * k3(i)
                    err = err + (erri / sk)**2
                end do
            else
                do i = 1, n
                    sk = atol(i) + rtol(i) * max(abs(y(i)), abs(k5(i)))
                    erri = k4(i) - bhh1 * k1(i) - bhh2 * k9(i) - bhh3 * k3(i)
                    err2 = err2 + (erri / sk)**2
                    erri = er1 * k1(i) + er6 * k6(i) + er7 * k7(i) + er8 * k8(i) + er9 * k9(i) + &
                           er10 * k10(i) + er11 * k2(i) + er12 * k3(i)
                    err = err + (erri / sk)**2
                end do
            end if
            
            deno = err + 0.01_dp * err2
            if (deno <= 0.0_dp) deno = 1.0_dp
            err = abs(h) * err * sqrt(1.0_dp / (n * deno))
            
            ! Computation of new step size
            fac11 = err**expo1
            
            ! Lund-stabilization
            fac = fac11 / facold**beta
            
            ! We require fac1 <= hnew/h <= fac2
            fac = max(facc2, min(facc1, fac / safe))
            hnew = h / fac
            
            if (err <= 1.0_dp) then
                ! Step is accepted
                facold = max(err, 1.0e-4_dp)
                naccpt = naccpt + 1
                call fcn(n, xph, k5, k4, rpar, ipar)
                nfcn = nfcn + 1
                
                ! Stiffness detection
                if (mod(naccpt, nstiff) == 0 .or. iasti > 0) then
                    stnum = 0.0_dp
                    stden = 0.0_dp
                    do i = 1, n
                        stnum = stnum + (k4(i) - k3(i))**2
                        stden = stden + (k5(i) - y1(i))**2
                    end do
                    if (stden > 0.0_dp) hlamb = abs(h) * sqrt(stnum / stden)
                    if (hlamb > 6.1_dp) then
                        nonsti = 0
                        iasti = iasti + 1
                        if (iasti == 15) then
                            if (iprint > 0) write(iprint, *) 'THE PROBLEM SEEMS TO BECOME STIFF AT X = ', x
                            if (iprint <= 0) goto 76
                        end if
                    else
                        nonsti = nonsti + 1
                        if (nonsti == 6) iasti = 0
                    end if
                end if
                
                ! Final preparation for dense output
                event = (iout == 3) .and. (xout <= xph)
                
                if (iout == 2 .or. event) then
                    ! Save the first function evaluations
                    do j = 1, nrd
                        i = icomp(j)
                        cont(j) = y(i)
                        ydiff = k5(i) - y(i)
                        cont(j + nrd) = ydiff
                        bspl = h * k1(i) - ydiff
                        cont(j + nrd * 2) = bspl
                        cont(j + nrd * 3) = ydiff - h * k4(i) - bspl
                        cont(j + nrd * 4) = d41 * k1(i) + d46 * k6(i) + d47 * k7(i) + d48 * k8(i) + &
                                            d49 * k9(i) + d410 * k10(i) + d411 * k2(i) + d412 * k3(i)
                        cont(j + nrd * 5) = d51 * k1(i) + d56 * k6(i) + d57 * k7(i) + d58 * k8(i) + &
                                            d59 * k9(i) + d510 * k10(i) + d511 * k2(i) + d512 * k3(i)
                        cont(j + nrd * 6) = d61 * k1(i) + d66 * k6(i) + d67 * k7(i) + d68 * k8(i) + &
                                            d69 * k9(i) + d610 * k10(i) + d611 * k2(i) + d612 * k3(i)
                        cont(j + nrd * 7) = d71 * k1(i) + d76 * k6(i) + d77 * k7(i) + d78 * k8(i) + &
                                            d79 * k9(i) + d710 * k10(i) + d711 * k2(i) + d712 * k3(i)
                    end do
                    
                    ! The next three function evaluations
                    do i = 1, n
                        y1(i) = y(i) + h * (a141 * k1(i) + a147 * k7(i) + a148 * k8(i) + &
                                           a149 * k9(i) + a1410 * k10(i) + a1411 * k2(i) + &
                                           a1412 * k3(i) + a1413 * k4(i))
                    end do
                    call fcn(n, x + c14 * h, y1, k10, rpar, ipar)
                    
                    do i = 1, n
                        y1(i) = y(i) + h * (a151 * k1(i) + a156 * k6(i) + a157 * k7(i) + &
                                           a158 * k8(i) + a1511 * k2(i) + a1512 * k3(i) + &
                                           a1513 * k4(i) + a1514 * k10(i))
                    end do
                    call fcn(n, x + c15 * h, y1, k2, rpar, ipar)
                    
                    do i = 1, n
                        y1(i) = y(i) + h * (a161 * k1(i) + a166 * k6(i) + a167 * k7(i) + &
                                           a168 * k8(i) + a169 * k9(i) + a1613 * k4(i) + &
                                           a1614 * k10(i) + a1615 * k2(i))
                    end do
                    call fcn(n, x + c16 * h, y1, k3, rpar, ipar)
                    nfcn = nfcn + 3
                    
                    ! Final preparation
                    do j = 1, nrd
                        i = icomp(j)
                        cont(j + nrd * 4) = h * (cont(j + nrd * 4) + d413 * k4(i) + d414 * k10(i) + &
                                                d415 * k2(i) + d416 * k3(i))
                        cont(j + nrd * 5) = h * (cont(j + nrd * 5) + d513 * k4(i) + d514 * k10(i) + &
                                                d515 * k2(i) + d516 * k3(i))
                        cont(j + nrd * 6) = h * (cont(j + nrd * 6) + d613 * k4(i) + d614 * k10(i) + &
                                                d615 * k2(i) + d616 * k3(i))
                        cont(j + nrd * 7) = h * (cont(j + nrd * 7) + d713 * k4(i) + d714 * k10(i) + &
                                                d715 * k2(i) + d716 * k3(i))
                    end do
                    hout = h
                    condo8_data%hout = hout
                end if
                
                do i = 1, n
                    k1(i) = k4(i)
                    y(i) = k5(i)
                end do
                xold = x
                x = xph
                condo8_data%xold = xold
                
                if (iout == 1 .or. iout == 2 .or. event) then
                    call solout(naccpt + 1, xold, x, y, n, cont, icomp, nrd, &
                                rpar, ipar, irtrn, xout)
                    if (irtrn < 0) goto 79
                end if
                
                ! Normal exit
                if (last) then
                    h = hnew
                    idid = 1
                    return
                end if
                
                if (abs(hnew) > hmax) hnew = posneg * hmax
                if (reject) hnew = posneg * min(abs(hnew), abs(h))
                reject = .false.
            else
                ! Step is rejected
                hnew = h / min(facc1, fac11 / safe)
                reject = .true.
                if (naccpt >= 1) nrejct = nrejct + 1
                last = .false.
            end if
            
            h = hnew
        end do integration_loop
        
        ! Error exits
        76 continue
        idid = -4
        return
        
        77 continue
        if (iprint > 0) write(iprint, 979) x
        if (iprint > 0) write(iprint, *) ' STEP SIZE TOO SMALL, H=', h
        idid = -3
        return
        
        78 continue
        if (iprint > 0) write(iprint, 979) x
        if (iprint > 0) write(iprint, *) 'MORE THAN NMAX =', nmax, 'STEPS ARE NEEDED'
        idid = -2
        return
        
        79 continue
        if (iprint > 0) write(iprint, 979) x
        979 format(' EXIT OF DOP853 AT X=', e18.4)
        idid = 2
        return
        
    end subroutine dp86co

    function hinit(n, fcn, x, y, xend, posneg, f0, f1, y1, iord, &
                   hmax, atol, rtol, itol, rpar, ipar) result(hinit_result)
        ! ----------------------------------------------------------
        ! COMPUTATION OF AN INITIAL STEP SIZE GUESS
        ! ----------------------------------------------------------
        
        implicit none
        
        ! Arguments
        integer, intent(in) :: n, iord, itol
        real(dp), intent(in) :: x, xend, posneg, hmax
        real(dp), intent(in) :: y(n), f0(n), atol(*), rtol(*)
        real(dp), intent(inout) :: f1(n), y1(n)
        real(dp), intent(inout) :: rpar(*)
        integer, intent(inout) :: ipar(*)
        real(dp) :: hinit_result
        
        ! External procedure
        external :: fcn
        
        ! Local variables
        real(dp) :: dnf, dny, atoli, rtoli, sk, h, der2, der12, h1
        integer :: i
        
        ! Compute a first guess for explicit Euler as
        ! H = 0.01 * NORM(Y0) / NORM(F0)
        dnf = 0.0_dp
        dny = 0.0_dp
        atoli = atol(1)
        rtoli = rtol(1)
        
        if (itol == 0) then
            do i = 1, n
                sk = atoli + rtoli * abs(y(i))
                dnf = dnf + (f0(i) / sk)**2
                dny = dny + (y(i) / sk)**2
            end do
        else
            do i = 1, n
                sk = atol(i) + rtol(i) * abs(y(i))
                dnf = dnf + (f0(i) / sk)**2
                dny = dny + (y(i) / sk)**2
            end do
        end if
        
        if (dnf <= 1.0e-10_dp .or. dny <= 1.0e-10_dp) then
            h = 1.0e-6_dp
        else
            h = sqrt(dny / dnf) * 0.01_dp
        end if
        
        h = min(h, hmax)
        h = sign(h, posneg)
        
        ! Perform an explicit Euler step
        do i = 1, n
            y1(i) = y(i) + h * f0(i)
        end do
        call fcn(n, x + h, y1, f1, rpar, ipar)
        
        ! Estimate the second derivative of the solution
        der2 = 0.0_dp
        
        if (itol == 0) then
            do i = 1, n
                sk = atoli + rtoli * abs(y(i))
                der2 = der2 + ((f1(i) - f0(i)) / sk)**2
            end do
        else
            do i = 1, n
                sk = atol(i) + rtol(i) * abs(y(i))
                der2 = der2 + ((f1(i) - f0(i)) / sk)**2
            end do
        end if
        
        der2 = sqrt(der2) / h
        
        ! Step size is computed such that
        ! H**IORD * MAX(NORM(F0), NORM(DER2)) = 0.01
        der12 = max(abs(der2), sqrt(dnf))
        
        if (der12 <= 1.0e-15_dp) then
            h1 = max(1.0e-6_dp, abs(h) * 1.0e-3_dp)
        else
            h1 = (0.01_dp / der12)**(1.0_dp / iord)
        end if
        
        h = min(100.0_dp * abs(h), h1, hmax)
        hinit_result = sign(h, posneg)
        
    end function hinit

    function contd8(ii, x, con, icomp, nd) result(contd8_result)
        ! ----------------------------------------------------------
        ! THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION
        ! WITH THE OUTPUT-SUBROUTINE FOR DOP853. IT PROVIDES AN
        ! APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X.
        ! ----------------------------------------------------------
        
        implicit none
        
        ! Arguments
        integer, intent(in) :: ii, nd
        real(dp), intent(in) :: x
        real(dp), intent(in) :: con(8*nd)
        integer, intent(in) :: icomp(nd)
        real(dp) :: contd8_result
        
        ! Local variables
        integer :: i, j
        real(dp) :: s, s1, conpar
        
        ! Compute place of ii-th component
        i = 0
        do j = 1, nd
            if (icomp(j) == ii) i = j
        end do
        
        if (i == 0) then
            write(6, *) ' NO DENSE OUTPUT AVAILABLE FOR COMP.', ii
            contd8_result = 0.0_dp
            return
        end if
        
        s = (x - condo8_data%xold) / condo8_data%hout
        s1 = 1.0_dp - s
        conpar = con(i + nd * 4) + s * (con(i + nd * 5) + s1 * (con(i + nd * 6) + s * con(i + nd * 7)))
        contd8_result = con(i) + s * (con(i + nd) + s1 * (con(i + nd * 2) + s * (con(i + nd * 3) + s1 * conpar)))
        
    end function contd8

end module dop853_mod