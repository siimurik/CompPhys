SUBROUTINE ODEX(N,FCN,X,Y,XEND,H, &
    RTOL,ATOL,ITOL, &
    SOLOUT,IOUT, &
    WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
! ----------------------------------------------------------
!     NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER
!     ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y).
!     THIS IS AN EXTRAPOLATION-ALGORITHM (GBS), BASED ON THE
!     EXPLICIT MIDPOINT RULE (WITH STEPSIZE CONTROL,
!     ORDER SELECTION AND DENSE OUTPUT).

!     AUTHORS: E. HAIRER AND G. WANNER
!              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
!              CH-1211 GENEVE 24, SWITZERLAND
!              E-MAIL:  Ernst.Hairer@unige.ch
!                       Gerhard.Wanner@unige.ch
!              DENSE OUTPUT WRITTEN BY E. HAIRER AND A. OSTERMANN

!     THIS CODE IS DESCRIBED IN SECTION II.9 OF THE BOOK:
!         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
!         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
!         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!         SPRINGER-VERLAG (1993)

!     BASED ON VERSION SEPTEMBER 30, 1995
!         SMALL CORRECTIONS ON OCTOBER 11, 2009
!         SMALL CHANGE FOR FEATURE IDID=2, OCTOBER 4, 2015
!               INTERRUPTED BY SOLOUT (BY C. LUDWIG)

!     INPUT PARAMETERS
!     ----------------
!     N           DIMENSION OF THE SYSTEM

!     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
!                 VALUE OF F(X,Y):
!                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
!                    DOUBLE PRECISION X,Y(N),F(N)
!                    F(1)=...   ETC.

!     X           INITIAL X-VALUE

!     Y(N)        INITIAL VALUES FOR Y

!     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)

!     H           INITIAL STEP SIZE GUESS;
!                 H=1.D0/(NORM OF F'), USUALLY 1.D-1 OR 1.D-3, IS GOOD.
!                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
!                 ADAPTS ITS STEP SIZE. WHEN YOU ARE NOT SURE, THEN
!                 STUDY THE CHOSEN VALUES FOR A FEW
!                 STEPS IN SUBROUTINE "SOLOUT".
!                 (IF H=0.D0, THE CODE PUTS H=1.D-4).

!     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
!                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.

!     ITOL        SWITCH FOR RTOL AND ATOL:
!                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
!                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
!                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
!                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
!                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
!                     RTOL(I)*ABS(Y(I))+ATOL(I).

!     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
!                 NUMERICAL SOLUTION DURING INTEGRATION.
!                 IF IOUT.GE.1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
!                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
!                 IT MUST HAVE THE FORM
!                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,NCON,ICOMP,ND,
!                                       RPAR,IPAR,IRTRN)
!                    DIMENSION X,Y(N),CON(NCON),ICOMP(ND)
!                    ....
!                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
!                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
!                    THE FIRST GRID-POINT).
!                 "XOLD" IS THE PRECEEDING GRID-POINT.
!                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
!                    IS SET <0, ODEX WILL RETURN TO THE CALLING PROGRAM.

!          -----  CONTINUOUS OUTPUT (IF IOUT=2): -----
!                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
!                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
!                 THE DOUBLE PRECISION FUNCTION
!                    >>>   CONTEX(I,S,CON,NCON,ICOMP,ND)   <<<
!                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
!                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
!                 S SHOULD LIE IN THE INTERVAL [XOLD,X].

!     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
!                    IOUT=0: SUBROUTINE IS NEVER CALLED
!                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
!                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT

!     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
!                 SERVES AS WORKING SPACE FOR ALL VECTORS.
!                 "LWORK" MUST BE AT LEAST
!                    N*(KM+5)+5*KM+20+(2*KM*(KM+2)+5)*NRDENS
!                 WHERE NRDENS=IWORK(8) (SEE BELOW) AND
!                        KM=9                IF IWORK(2)=0
!                        KM=IWORK(2)         IF IWORK(2).GT.0
!                 WORK(1),...,WORK(20) SERVE AS PARAMETERS
!                 FOR THE CODE. FOR STANDARD USE, SET THESE
!                 PARAMETERS TO ZERO BEFORE CALLING.

!     LWORK       DECLARED LENGTH OF ARRAY "WORK".

!     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
!                 "LIWORK" MUST BE AT LEAST
!                               2*KM+21+NRDENS
!                 IWORK(1),...,IWORK(20) SERVE AS PARAMETERS
!                 FOR THE CODE. FOR STANDARD USE, SET THESE
!                 PARAMETERS TO ZERO BEFORE CALLING.

!     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".

!     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
!                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
!                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES.

!-----------------------------------------------------------------------

!     SOPHISTICATED SETTING OF PARAMETERS
!     -----------------------------------
!              SEVERAL PARAMETERS (WORK(1),...,IWORK(1),...) ALLOW
!              TO ADAPT THE CODE TO THE PROBLEM AND TO THE NEEDS OF
!              THE USER. FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES.

!    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 2.3D-16.

!    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X.

!    WORK(3)   STEP SIZE IS REDUCED BY FACTOR WORK(3), IF THE
!              STABILITY CHECK IS NEGATIVE, DEFAULT 0.5.

!    WORK(4), WORK(5)   PARAMETERS FOR STEP SIZE SELECTION
!              THE NEW STEP SIZE FOR THE J-TH DIAGONAL ENTRY IS
!              CHOSEN SUBJECT TO THE RESTRICTION
!                 FACMIN/WORK(5) <= HNEW(J)/HOLD <= 1/FACMIN
!              WHERE FACMIN=WORK(4)**(1/(2*J-1))
!              DEFAULT VALUES: WORK(4)=0.02D0, WORK(5)=4.D0

!    WORK(6), WORK(7)   PARAMETERS FOR THE ORDER SELECTION
!              STEP SIZE IS DECREASED IF    W(K-1) <= W(K)*WORK(6)
!              STEP SIZE IS INCREASED IF    W(K) <= W(K-1)*WORK(7)
!              DEFAULT VALUES: WORK(6)=0.8D0, WORK(7)=0.9D0

!    WORK(8), WORK(9)   SAFETY FACTORS FOR STEP CONTROL ALGORITHM
!             HNEW=H*WORK(9)*(WORK(8)*TOL/ERR)**(1/(J-1))
!             DEFAULT VALUES: WORK(8)=0.65D0,
!                        WORK(9)=0.94D0  IF "HOPE FOR CONVERGENCE"
!                        WORK(9)=0.90D0  IF "NO HOPE FOR CONVERGENCE"

!    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
!              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 10000.

!    IWORK(2)  THE MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION
!              TABLE. THE DEFAULT VALUE (FOR IWORK(2)=0) IS 9.
!              IF IWORK(2).NE.0 THEN IWORK(2) SHOULD BE .GE.3.

!    IWORK(3)  SWITCH FOR THE STEP SIZE SEQUENCE (EVEN NUMBERS ONLY)
!              IF IWORK(3).EQ.1 THEN 2,4,6,8,10,12,14,16,...
!              IF IWORK(3).EQ.2 THEN 2,4,8,12,16,20,24,28,...
!              IF IWORK(3).EQ.3 THEN 2,4,6,8,12,16,24,32,...
!              IF IWORK(3).EQ.4 THEN 2,6,10,14,18,22,26,30,...
!              IF IWORK(3).EQ.5 THEN 4,8,12,16,20,24,28,32,...
!              THE DEFAULT VALUE IS IWORK(3)=1 IF IOUT.LE.1;
!              THE DEFAULT VALUE IS IWORK(3)=4 IF IOUT.GE.2.

!    IWORK(4)  STABILITY CHECK IS ACTIVATED AT MOST IWORK(4) TIMES IN
!              ONE LINE OF THE EXTRAP. TABLE, DEFAULT IWORK(4)=1.

!    IWORK(5)  STABILITY CHECK IS ACTIVATED ONLY IN THE LINES
!              1 TO IWORK(5) OF THE EXTRAP. TABLE, DEFAULT IWORK(5)=1.

!    IWORK(6)  IF  IWORK(6)=0  ERROR ESTIMATOR IN THE DENSE
!              OUTPUT FORMULA IS ACTIVATED. IT CAN BE SUPPRESSED
!              BY PUTTING IWORK(6)=1.
!              DEFAULT IWORK(6)=0  (IF IOUT.GE.2).

!    IWORK(7)  DETERMINES THE DEGREE OF INTERPOLATION FORMULA
!              MU = 2 * KAPPA - IWORK(7) + 1
!              IWORK(7) SHOULD LIE BETWEEN 1 AND 6
!              DEFAULT IWORK(7)=4  (IF IWORK(7)=0).

!    IWORK(8)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT
!              IS REQUIRED

!    IWORK(21),...,IWORK(NRDENS+20) INDICATE THE COMPONENTS, FOR WHICH
!              DENSE OUTPUT IS REQUIRED

!----------------------------------------------------------------------
!     OUTPUT PARAMETERS
!     -----------------
!     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
!                 (AFTER SUCCESSFUL RETURN X=XEND).

!     Y(N)        NUMERICAL SOLUTION AT X

!     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP

!     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
!                   IDID=1  COMPUTATION SUCCESSFUL,
!                   IDID=2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
!                   IDID=-1 COMPUTATION UNSUCCESSFUL.

!   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS
!   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS
!   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS
!   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
!                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
!-----------------------------------------------------------------------
! *** *** *** *** *** *** *** *** *** *** *** *** ***
!          DECLARATIONS
! *** *** *** *** *** *** *** *** *** *** *** *** ***
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, ITOL, IOUT, LWORK, LIWORK
    DOUBLE PRECISION, INTENT(INOUT) :: X, H
    DOUBLE PRECISION, INTENT(IN) :: XEND
    DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
    DOUBLE PRECISION, INTENT(IN) :: ATOL(*), RTOL(*)
    DOUBLE PRECISION, INTENT(INOUT) :: WORK(LWORK)
    INTEGER, INTENT(INOUT) :: IWORK(LIWORK)
    DOUBLE PRECISION, INTENT(INOUT) :: RPAR(*)
    INTEGER, INTENT(INOUT) :: IPAR(*)
    INTEGER, INTENT(OUT) :: IDID
    
    ! Local variables
    INTEGER :: NFCN, NSTEP, NACCPT, NREJCT
    INTEGER :: NMAX, KM, NSEQU, MSTAB, JSTAB, IDERR, MUDIF, NRDENS
    INTEGER :: LFSAFE, IEDY, IEYH1, IEYH2, IEDZ, IESCAL, IET, IEFS
    INTEGER :: IEYS, IEHH, IEW, IEA, IEFAC, IECO, ISTORE
    INTEGER :: ICOM, IENJ, IEIP, NRD, NCOM
    INTEGER :: I
    DOUBLE PRECISION :: UROUND, HMAX, SAFE3, FAC1, FAC2, FAC3, FAC4
    DOUBLE PRECISION :: SAFE1, SAFE2
    LOGICAL :: ARRET
    
    EXTERNAL FCN, SOLOUT
    
! *** *** *** *** *** *** ***
!        SETTING THE PARAMETERS
! *** *** *** *** *** *** ***
    NFCN = 0
    NSTEP = 0
    NACCPT = 0
    NREJCT = 0
    ARRET = .FALSE.
    
! -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
    IF (IWORK(1) == 0) THEN
        NMAX = 10000
    ELSE
        NMAX = IWORK(1)
        IF (NMAX <= 0) THEN
            WRITE(6,*) ' WRONG INPUT IWORK(1)=', IWORK(1)
            ARRET = .TRUE.
        END IF
    END IF
    
! -------- KM     MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION
    IF (IWORK(2) == 0) THEN
        KM = 9
    ELSE
        KM = IWORK(2)
        IF (KM <= 2) THEN
            WRITE(6,*) ' CURIOUS INPUT IWORK(2)=', IWORK(2)
            ARRET = .TRUE.
        END IF
    END IF
    
! -------- NSEQU     CHOICE OF STEP SIZE SEQUENCE
    NSEQU = IWORK(3)
    IF (IWORK(3) == 0 .AND. IOUT <= 1) NSEQU = 1
    IF (IWORK(3) == 0 .AND. IOUT >= 2) NSEQU = 4
    IF (NSEQU <= 0 .OR. NSEQU >= 6) THEN
        WRITE(6,*) ' CURIOUS INPUT IWORK(3)=', IWORK(3)
        ARRET = .TRUE.
    END IF
    IF (NSEQU <= 3 .AND. IOUT >= 2) THEN
        WRITE(6,*) ' IWORK(3) NOT COMPATIBLE WITH IOUT'
        ARRET = .TRUE.
    END IF
    
! -------- MSTAB     PARAMETER FOR STABILITY CHECK
    IF (IWORK(4) == 0) THEN
        MSTAB = 1
    ELSE
        MSTAB = IWORK(4)
    END IF
    
! -------- JSTAB     PARAMETER FOR STABILITY CHECK
    IF (IWORK(5) == 0) THEN
        JSTAB = 2
    ELSE
        JSTAB = IWORK(5)
    END IF
    
! -------- IDERR  PARAMETER FOR ERROR ESTIMATION IN DENSE OUTPUT
    IF (IWORK(6) == 0) THEN
        IF (IOUT <= 1) IDERR = 1
        IF (IOUT >= 2) IDERR = 0
    ELSE
        IDERR = IWORK(6)
        IF (IOUT <= 1) THEN
            WRITE(6,*) ' ERROR ESTIMATION IN DENSE OUTPUT', &
            ' NOT POSSIBLE, WRONG IWORK(6)=', IWORK(6)
            ARRET = .TRUE.
        END IF
    END IF
    
! -------- MUDIF
    IF (IWORK(7) == 0) THEN
        MUDIF = 4
    ELSE
        MUDIF = IWORK(7)
        IF (MUDIF <= 0 .OR. MUDIF >= 7) THEN
            WRITE(6,*) ' WRONG INPUT IWORK(7)=', IWORK(7)
            ARRET = .TRUE.
        END IF
    END IF
    
! -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS
    NRDENS = IWORK(8)
    IF (NRDENS < 0 .OR. NRDENS > N) THEN
        WRITE(6,*) ' CURIOUS INPUT IWORK(8)=', IWORK(8)
        ARRET = .TRUE.
    END IF
    IF (NRDENS == N) THEN
        DO I = 1, NRDENS
            IWORK(20+I) = I
        END DO
    END IF
    
! -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0
    IF (WORK(1) == 0.D0) THEN
        UROUND = 2.3D-16
    ELSE
        UROUND = WORK(1)
        IF (UROUND <= 1.D-35 .OR. UROUND >= 1.D0) THEN
            WRITE(6,*) ' WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:', &
            WORK(1)
            ARRET = .TRUE.
        END IF
    END IF
    
! -------- MAXIMAL STEP SIZE
    IF (WORK(2) == 0.D0) THEN
        HMAX = XEND - X
    ELSE
        HMAX = ABS(WORK(2))
    END IF
    
! -------- STEP SIZE REDUCTION FACTOR
    IF (WORK(3) == 0.D0) THEN
        SAFE3 = 0.5D0
    ELSE
        SAFE3 = WORK(3)
        IF (SAFE3 <= UROUND .OR. SAFE3 >= 1.D0) THEN
            WRITE(6,*) ' CURIOUS INPUT WORK(3)=', WORK(3)
            ARRET = .TRUE.
        END IF
    END IF
    
! -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
    IF (WORK(4) == 0.D0) THEN
        FAC1 = 0.02D0
    ELSE
        FAC1 = WORK(4)
    END IF
    IF (WORK(5) == 0.D0) THEN
        FAC2 = 4.0D0
    ELSE
        FAC2 = WORK(5)
    END IF
    
! -------  FAC3, FAC4   PARAMETERS FOR THE ORDER SELECTION
    IF (WORK(6) == 0.D0) THEN
        FAC3 = 0.8D0
    ELSE
        FAC3 = WORK(6)
    END IF
    IF (WORK(7) == 0.D0) THEN
        FAC4 = 0.9D0
    ELSE
        FAC4 = WORK(7)
    END IF
    
! ------- SAFE1, SAFE2 SAFETY FACTORS FOR STEP SIZE PREDICTION
    IF (WORK(8) == 0.D0) THEN
        SAFE1 = 0.65D0
    ELSE
        SAFE1 = WORK(8)
    END IF
    IF (WORK(9) == 0.D0) THEN
        SAFE2 = 0.94D0
    ELSE
        SAFE2 = WORK(9)
    END IF
    
! ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
    LFSAFE = 2*KM*KM + KM
    IEDY = 21
    IEYH1 = IEDY + N
    IEYH2 = IEYH1 + N
    IEDZ = IEYH2 + N
    IESCAL = IEDZ + N
    IET = IESCAL + N
    IEFS = IET + KM*N
    IEYS = IEFS + LFSAFE*NRDENS
    IEHH = IEYS + KM*NRDENS
    IEW = IEHH + KM
    IEA = IEW + KM
    IEFAC = IEA + KM
    
! ------ TOTAL STORAGE REQUIREMENT -----------
    IECO = IEFAC + 2*KM
    ISTORE = IECO + (2*KM+5)*NRDENS - 1
    IF (ISTORE > LWORK) THEN
        WRITE(6,*) ' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=', ISTORE
        ARRET = .TRUE.
    END IF
    
! ------- ENTRY POINTS FOR INTEGER WORKSPACE -----
    ICOM = 21
    IENJ = ICOM + NRDENS
    
! --------- TOTAL REQUIREMENT ---------------
    IEIP = IENJ + KM
    ISTORE = IEIP + KM + 1 - 1
    IF (ISTORE > LIWORK) THEN
        WRITE(6,*) ' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=', ISTORE
        ARRET = .TRUE.
    END IF
    
! ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
    IF (ARRET) THEN
        IDID = -1
        RETURN
    END IF
    
! -------- CALL TO CORE INTEGRATOR ------------
    NRD = MAX(1, NRDENS)
    NCOM = MAX(1, (2*KM+5)*NRDENS)
    CALL ODXCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,KM, &
    SOLOUT,IOUT,IDID,NMAX,UROUND,WORK(IEDY),WORK(IEYH1), &
    WORK(IEYH2),WORK(IEDZ),WORK(IESCAL),WORK(IEFS), &
    WORK(IEYS),WORK(IET),WORK(IEHH),WORK(IEW),WORK(IEA), &
    WORK(IECO),NCOM,IWORK(ICOM), &
    IWORK(IENJ),IWORK(IEIP),NSEQU,MSTAB,JSTAB,LFSAFE, &
    SAFE1,SAFE2,SAFE3,FAC1,FAC2,FAC3,FAC4,IDERR,WORK(IEFAC), &
    MUDIF,NRD,RPAR,IPAR,NFCN,NSTEP,NACCPT,NREJCT)
    
    IWORK(17) = NFCN
    IWORK(18) = NSTEP
    IWORK(19) = NACCPT
    IWORK(20) = NREJCT
    
! ----------- RETURN -----------
    RETURN
    END SUBROUTINE ODEX



!  ----- ... AND HERE IS THE CORE INTEGRATOR  ----------

    SUBROUTINE ODXCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,KM, &
    SOLOUT,IOUT,IDID,NMAX,UROUND,DY,YH1,YH2,DZ,SCAL,FSAFE, &
    YSAFE,T,HH,W,A,DENS,NCOM,ICOMP,NJ,IPOINT,NSEQU,MSTAB,JSTAB, &
    LFSAFE,SAFE1,SAFE2,SAFE3,FAC1,FAC2,FAC3,FAC4,IDERR,ERRFAC, &
    MUDIF,NRD,RPAR,IPAR,NFCN,NSTEP,NACCPT,NREJCT)
! ----------------------------------------------------------
!     CORE INTEGRATOR FOR ODEX
!     PARAMETERS SAME AS IN ODEX WITH WORKSPACE ADDED
! ----------------------------------------------------------
!         DECLARATIONS
! ----------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, ITOL, KM, IOUT, NMAX, NSEQU, MSTAB, JSTAB
    INTEGER, INTENT(IN) :: LFSAFE, IDERR, MUDIF, NRD, NCOM
    INTEGER, INTENT(OUT) :: IDID
    INTEGER, INTENT(INOUT) :: NFCN, NSTEP, NACCPT, NREJCT
    DOUBLE PRECISION, INTENT(INOUT) :: X, H
    DOUBLE PRECISION, INTENT(IN) :: XEND, HMAX, UROUND
    DOUBLE PRECISION, INTENT(IN) :: SAFE1, SAFE2, SAFE3
    DOUBLE PRECISION, INTENT(IN) :: FAC1, FAC2, FAC3, FAC4
    DOUBLE PRECISION, INTENT(INOUT) :: Y(N), DY(N), YH1(N), YH2(N), DZ(N), SCAL(N)
    DOUBLE PRECISION, INTENT(IN) :: RTOL(*), ATOL(*)
    DOUBLE PRECISION, INTENT(INOUT) :: T(KM,N), HH(KM), W(KM), A(KM)
    DOUBLE PRECISION, INTENT(INOUT) :: FSAFE(LFSAFE,NRD), YSAFE(KM,NRD)
    DOUBLE PRECISION, INTENT(INOUT) :: ERRFAC(2*KM), DENS(NCOM)
    INTEGER, INTENT(INOUT) :: NJ(KM), IPOINT(KM+1), ICOMP(NRD)
    DOUBLE PRECISION, INTENT(INOUT) :: RPAR(*)
    INTEGER, INTENT(INOUT) :: IPAR(*)
    
    ! Local variables
    LOGICAL :: REJECT, LAST, ATOV
    INTEGER :: K, KC, J, I, IPT, IRTRN, MU, NJADD, NCON
    INTEGER :: KMI, KBEG, KK, L, LBEG, LEND, KOPT, KRN, KLN, KMIT
    DOUBLE PRECISION :: POSNEG, ERR, ERROLD, HOPTDE, XOLD, FACNJ
    DOUBLE PRECISION :: FAC, ERRX, PROD, DBLENJ, FACTOR, ERRINT
    
    EXTERNAL FCN, SOLOUT
    COMMON /CONODX/XOLDD,HHH,KMIT
    DOUBLE PRECISION :: XOLDD, HHH
    
! --- DEFINE THE STEP SIZE SEQUENCE
    IF (NSEQU == 1) THEN
        DO I = 1, KM
            NJ(I) = 2*I
        END DO
    END IF
    IF (NSEQU == 2) THEN
        NJ(1) = 2
        DO I = 2, KM
            NJ(I) = 4*I - 4
        END DO
    END IF
    IF (NSEQU == 3) THEN
        NJ(1) = 2
        NJ(2) = 4
        NJ(3) = 6
        DO I = 4, KM
            NJ(I) = 2*NJ(I-2)
        END DO
    END IF
    IF (NSEQU == 4) THEN
        DO I = 1, KM
            NJ(I) = 4*I - 2
        END DO
    END IF
    IF (NSEQU == 5) THEN
        DO I = 1, KM
            NJ(I) = 4*I
        END DO
    END IF
    
! --- DEFINE THE A(I) FOR ORDER SELECTION
    A(1) = 1.D0 + NJ(1)
    DO I = 2, KM
        A(I) = A(I-1) + NJ(I)
    END DO
    
! --- INITIAL SCALING
    DO I = 1, N
        IF (ITOL == 0) THEN
            SCAL(I) = ATOL(1) + RTOL(1)*ABS(Y(I))
        ELSE
            SCAL(I) = ATOL(I) + RTOL(I)*ABS(Y(I))
        END IF
    END DO
    
! --- INITIAL PREPARATIONS
    POSNEG = SIGN(1.D0, XEND-X)
    K = MAX(2, MIN(KM-1, INT(-LOG10(RTOL(1)+1.0D-40)*0.6D0+1.5D0)))
    H = MAX(ABS(H), 1.D-4)
    H = POSNEG*MIN(H, HMAX, ABS(XEND-X)/2.D0)
    
    IF (IOUT >= 1) THEN
        IF (IOUT >= 2) THEN
            IPOINT(1) = 0
            DO I = 1, KM
                NJADD = 4*I - 2
                IF (NJ(I) > NJADD) NJADD = NJADD + 1
                IPOINT(I+1) = IPOINT(I) + NJADD
            END DO
            DO MU = 1, KM*2
                ERRX = SQRT(REAL(MU, KIND=KIND(1.0D0))/(MU+4.D0))*0.5D0
                PROD = 1.D0/(MU+4.D0)**2
                DO J = 1, MU
                    PROD = PROD*ERRX/J
                END DO
                ERRFAC(MU) = PROD
            END DO
            IPT = 0
        END IF
        IRTRN = 0
        XOLD = X
        NCON = NCOM
        CALL SOLOUT(NACCPT+1, XOLD, X, Y, N, DENS, NCON, ICOMP, NRD, &
        RPAR, IPAR, IRTRN)
        IF (IRTRN < 0) GOTO 130
    END IF
    
    ERR = 0.D0
    ERROLD = 1.D10
    HOPTDE = POSNEG*HMAX
    W(1) = 0.D0
    REJECT = .FALSE.
    LAST = .FALSE.
    
10  ATOV = .FALSE.
! --- IS XEND REACHED IN THE NEXT STEP?
    IF (0.1D0*ABS(XEND-X) <= ABS(X)*UROUND) GOTO 110
    H = POSNEG*MIN(ABS(H), ABS(XEND-X), HMAX, ABS(HOPTDE))
    IF ((X+1.01D0*H-XEND)*POSNEG > 0.D0) THEN
        H = XEND - X
        LAST = .TRUE.
    END IF
    IF (NSTEP == 0 .OR. IOUT /= 2) CALL FCN(N, X, Y, DZ, RPAR, IPAR)
    NFCN = NFCN + 1
    
! --- THE FIRST AND LAST STEP
    IF (NSTEP == 0 .OR. LAST) THEN
        IPT = 0
        NSTEP = NSTEP + 1
        DO J = 1, K
            KC = J
            CALL MIDEX(J,X,Y,H,HMAX,N,FCN,DY,YH1,YH2,DZ,T,NJ,HH,W, &
            ERR,FAC,A,SAFE1,UROUND,FAC1,FAC2,SAFE2,SCAL,ATOV,SAFE3, &
            REJECT,KM,RTOL,ATOL,ITOL,MSTAB,JSTAB,ERROLD,FSAFE,LFSAFE, &
            IOUT,IPT,YSAFE,ICOMP,NRD,RPAR,IPAR,NFCN)
            IF (ATOV) GO TO 10
            IF (J > 1 .AND. ERR <= 1.D0) GO TO 60
        END DO
        GO TO 55
    END IF
    
! --- BASIC INTEGRATION STEP
30  CONTINUE
    IPT = 0
    NSTEP = NSTEP + 1
    IF (NSTEP >= NMAX) GO TO 120
    KC = K - 1
    DO J = 1, KC
        CALL MIDEX(J,X,Y,H,HMAX,N,FCN,DY,YH1,YH2,DZ,T,NJ,HH,W, &
        ERR,FAC,A,SAFE1,UROUND,FAC1,FAC2,SAFE2,SCAL,ATOV,SAFE3, &
        REJECT,KM,RTOL,ATOL,ITOL,MSTAB,JSTAB,ERROLD,FSAFE,LFSAFE, &
        IOUT,IPT,YSAFE,ICOMP,NRD,RPAR,IPAR,NFCN)
        IF (ATOV) GO TO 10
    END DO
    
! --- CONVERGENCE MONITOR
    IF (K == 2 .OR. REJECT) GO TO 50
    IF (ERR <= 1.D0) GO TO 60
    IF (ERR > ((NJ(K+1)*NJ(K))/4.D0)**2) GO TO 100
    
50  CONTINUE
    CALL MIDEX(K,X,Y,H,HMAX,N,FCN,DY,YH1,YH2,DZ,T,NJ,HH,W, &
    ERR,FAC,A,SAFE1,UROUND,FAC1,FAC2,SAFE2,SCAL,ATOV,SAFE3, &
    REJECT,KM,RTOL,ATOL,ITOL,MSTAB,JSTAB,ERROLD,FSAFE,LFSAFE, &
    IOUT,IPT,YSAFE,ICOMP,NRD,RPAR,IPAR,NFCN)
    IF (ATOV) GO TO 10
    KC = K
    IF (ERR <= 1.D0) GO TO 60
    
! --- HOPE FOR CONVERGENCE IN LINE K+1
55  CONTINUE
    IF (ERR > (NJ(K+1)/2.D0)**2) GO TO 100
    KC = K + 1
    CALL MIDEX(KC,X,Y,H,HMAX,N,FCN,DY,YH1,YH2,DZ,T,NJ,HH,W, &
    ERR,FAC,A,SAFE1,UROUND,FAC1,FAC2,SAFE2,SCAL,ATOV,SAFE3, &
    REJECT,KM,RTOL,ATOL,ITOL,MSTAB,JSTAB,ERROLD,FSAFE,LFSAFE, &
    IOUT,IPT,YSAFE,ICOMP,NRD,RPAR,IPAR,NFCN)
    IF (ATOV) GO TO 10
    IF (ERR > 1.D0) GO TO 100
    
! --- STEP IS ACCEPTED
60  XOLD = X
    X = X + H
    IF (IOUT >= 2) THEN
    ! ---  KMIT = MU OF THE PAPER
        KMIT = 2*KC - MUDIF + 1
        DO I = 1, NRD
            DENS(I) = Y(ICOMP(I))
        END DO
        XOLDD = XOLD
        HHH = H
        DO I = 1, NRD
            DENS(NRD+I) = H*DZ(ICOMP(I))
        END DO
        KLN = 2*NRD
        DO I = 1, NRD
            DENS(KLN+I) = T(1, ICOMP(I))
        END DO
        
    ! --- COMPUTE SOLUTION AT MID-POINT ----
        DO J = 2, KC
            DBLENJ = NJ(J)
            DO L = J, 2, -1
                FACTOR = (DBLENJ/NJ(L-1))**2 - 1.D0
                DO I = 1, NRD
                    YSAFE(L-1,I) = YSAFE(L,I) + (YSAFE(L,I)-YSAFE(L-1,I))/FACTOR
                END DO
            END DO
        END DO
        KRN = 4*NRD
        DO I = 1, NRD
            DENS(KRN+I) = YSAFE(1, I)
        END DO
        
    ! --- COMPUTE FIRST DERIVATIVE AT RIGHT END ----
        DO I = 1, N
            YH1(I) = T(1, I)
        END DO
        CALL FCN(N, X, YH1, YH2, RPAR, IPAR)
        KRN = 3*NRD
        DO I = 1, NRD
            DENS(KRN+I) = YH2(ICOMP(I))*H
        END DO
        
    ! --- THE LOOP ---
        DO KMI = 1, KMIT
        ! --- COMPUTE KMI-TH DERIVATIVE AT MID-POINT ----
            KBEG = (KMI+1)/2
            DO KK = KBEG, KC
                FACNJ = (NJ(KK)/2.D0)**(KMI-1)
                IPT = IPOINT(KK+1) - 2*KK + KMI
                DO I = 1, NRD
                    YSAFE(KK,I) = FSAFE(IPT,I)*FACNJ
                END DO
            END DO
            DO J = KBEG+1, KC
                DBLENJ = NJ(J)
                DO L = J, KBEG+1, -1
                    FACTOR = (DBLENJ/NJ(L-1))**2 - 1.D0
                    DO I = 1, NRD
                        YSAFE(L-1,I) = YSAFE(L,I) + (YSAFE(L,I)-YSAFE(L-1,I))/FACTOR
                    END DO
                END DO
            END DO
            KRN = (KMI+4)*NRD
            DO I = 1, NRD
                DENS(KRN+I) = YSAFE(KBEG,I)*H
            END DO
            IF (KMI == KMIT) GOTO 180
            
        ! --- COMPUTE DIFFERENCES
            DO KK = (KMI+2)/2, KC
                LBEG = IPOINT(KK+1)
                LEND = IPOINT(KK) + KMI + 1
                IF (KMI == 1 .AND. NSEQU == 4) LEND = LEND + 2
                DO L = LBEG, LEND, -2
                    DO I = 1, NRD
                        FSAFE(L,I) = FSAFE(L,I) - FSAFE(L-2,I)
                    END DO
                END DO
                IF (KMI == 1 .AND. NSEQU == 4) THEN
                    L = LEND - 2
                    DO I = 1, NRD
                        FSAFE(L,I) = FSAFE(L,I) - DZ(ICOMP(I))
                    END DO
                END IF
            END DO
            
        ! --- COMPUTE DIFFERENCES
            DO KK = (KMI+2)/2, KC
                LBEG = IPOINT(KK+1) - 1
                LEND = IPOINT(KK) + KMI + 2
                DO L = LBEG, LEND, -2
                    DO I = 1, NRD
                        FSAFE(L,I) = FSAFE(L,I) - FSAFE(L-2,I)
                    END DO
                END DO
            END DO
180     END DO
        CALL INTERP(NRD, DENS, KMIT)
        
    ! --- ESTIMATION OF INTERPOLATION ERROR
        IF (IDERR == 0 .AND. KMIT >= 1) THEN
            ERRINT = 0.D0
            DO I = 1, NRD
                ERRINT = ERRINT + (DENS((KMIT+4)*NRD+I)/SCAL(ICOMP(I)))**2
            END DO
            ERRINT = SQRT(ERRINT/NRD)*ERRFAC(KMIT)
            HOPTDE = H/MAX((ERRINT)**(1.D0/(KMIT+4)), 0.01D0)
            IF (ERRINT > 10.D0) THEN
                H = HOPTDE
                X = XOLD
                NREJCT = NREJCT + 1
                REJECT = .TRUE.
                GOTO 10
            END IF
        END IF
        DO I = 1, N
            DZ(I) = YH2(I)
        END DO
    END IF
    
    DO I = 1, N
        Y(I) = T(1, I)
    END DO
    NACCPT = NACCPT + 1
    IF (IOUT >= 1) THEN
        NCON = NCOM
        CALL SOLOUT(NACCPT+1, XOLD, X, Y, N, DENS, NCON, ICOMP, NRD, &
        RPAR, IPAR, IRTRN)
        IF (IRTRN < 0) GOTO 130
    END IF
    
! --- COMPUTE OPTIMAL ORDER
    IF (KC == 2) THEN
        KOPT = MIN(3, KM-1)
        IF (REJECT) KOPT = 2
        GO TO 80
    END IF
    IF (KC <= K) THEN
        KOPT = KC
        IF (W(KC-1) < W(KC)*FAC3) KOPT = KC - 1
        IF (W(KC) < W(KC-1)*FAC4) KOPT = MIN(KC+1, KM-1)
    ELSE
        KOPT = KC - 1
        IF (KC > 3 .AND. W(KC-2) < W(KC-1)*FAC3) KOPT = KC - 2
        IF (W(KC) < W(KOPT)*FAC4) KOPT = MIN(KC, KM-1)
    END IF
    
! --- AFTER A REJECTED STEP
80  IF (REJECT) THEN
        K = MIN(KOPT, KC)
        H = POSNEG*MIN(ABS(H), ABS(HH(K)))
        REJECT = .FALSE.
        GO TO 10
    END IF
    
! --- COMPUTE STEPSIZE FOR NEXT STEP
    IF (KOPT <= KC) THEN
        H = HH(KOPT)
    ELSE
        IF (KC < K .AND. W(KC) < W(KC-1)*FAC4) THEN
            H = HH(KC)*A(KOPT+1)/A(KC)
        ELSE
            H = HH(KC)*A(KOPT)/A(KC)
        END IF
    END IF
    K = KOPT
    H = POSNEG*ABS(H)
    GO TO 10
    
! --- STEP IS REJECTED
100 CONTINUE
    K = MIN(K, KC, KM-1)
    IF (K > 2 .AND. W(K-1) < W(K)*FAC3) K = K - 1
    NREJCT = NREJCT + 1
    H = POSNEG*HH(K)
    REJECT = .TRUE.
    GO TO 30
    
! --- SOLUTION EXIT
110 CONTINUE
    IDID = 1
    RETURN
    
! --- INTERRUPTED BY SOLOUT
130 CONTINUE
    IDID = 2
    RETURN
    
! --- FAIL EXIT
120 WRITE(6, 979) X, H
979 FORMAT(' EXIT OF ODEX AT X=', D14.7, '   H=', D14.7)
    IDID = -1
    RETURN
    END SUBROUTINE ODXCOR

    SUBROUTINE MIDEX(J,X,Y,H,HMAX,N,FCN,DY,YH1,YH2,DZ,T,NJ,HH,W, &
    ERR,FAC,A,SAFE1,UROUND,FAC1,FAC2,SAFE2,SCAL,ATOV,SAFE3, &
    REJECT,KM,RTOL,ATOL,ITOL,MSTAB,JSTAB,ERROLD,FSAFE,LFSAFE, &
    IOUT,IPT,YSAFE,ICOMP,NRD,RPAR,IPAR,NFCN)
! --- THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE
! --- EXTRAPOLATION TABLE AND PROVIDES AN ESTIMATION
! --- OF THE OPTIMAL STEPSIZE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: J, N, KM, ITOL, MSTAB, JSTAB, LFSAFE, IOUT, NRD
    INTEGER, INTENT(INOUT) :: NFCN, IPT
    DOUBLE PRECISION, INTENT(IN) :: X, HMAX, SAFE1, UROUND, FAC1, FAC2
    DOUBLE PRECISION, INTENT(IN) :: SAFE2, SAFE3
    DOUBLE PRECISION, INTENT(INOUT) :: Y(N), DY(N), YH1(N), YH2(N), DZ(N), SCAL(N)
    DOUBLE PRECISION, INTENT(INOUT) :: T(KM,N), HH(KM), W(KM), A(KM)
    DOUBLE PRECISION, INTENT(IN) :: RTOL(*), ATOL(*)
    DOUBLE PRECISION, INTENT(INOUT) :: ERR, FAC, ERROLD, H
    DOUBLE PRECISION, INTENT(INOUT) :: FSAFE(LFSAFE,NRD), YSAFE(KM,NRD)
    INTEGER, INTENT(IN) :: NJ(KM), ICOMP(NRD)
    DOUBLE PRECISION, INTENT(INOUT) :: RPAR(*)
    INTEGER, INTENT(INOUT) :: IPAR(*)
    LOGICAL, INTENT(INOUT) :: REJECT, ATOV
    
    ! Local variables
    INTEGER :: I, M, MM, NJMID, L
    DOUBLE PRECISION :: HJ, YS, DEL1, DEL2, QUOT, DBLENJ, T1I
    DOUBLE PRECISION :: EXPO, FACMIN
    
    EXTERNAL FCN
    
    HJ = H/NJ(J)
    
! --- EULER STARTING STEP
    DO I = 1, N
        YH1(I) = Y(I)
        YH2(I) = Y(I) + HJ*DZ(I)
    END DO
    
! --- EXPLICIT MIDPOINT RULE
    M = NJ(J) - 1
    NJMID = NJ(J)/2
    DO MM = 1, M
        IF (IOUT >= 2 .AND. MM == NJMID) THEN
            DO I = 1, NRD
                YSAFE(J,I) = YH2(ICOMP(I))
            END DO
        END IF
        CALL FCN(N, X+HJ*MM, YH2, DY, RPAR, IPAR)
        IF (IOUT >= 2 .AND. ABS(MM-NJMID) <= 2*J-1) THEN
            IPT = IPT + 1
            DO I = 1, NRD
                FSAFE(IPT,I) = DY(ICOMP(I))
            END DO
        END IF
        DO I = 1, N
            YS = YH1(I)
            YH1(I) = YH2(I)
            YH2(I) = YS + 2.D0*HJ*DY(I)
        END DO
        IF (MM <= MSTAB .AND. J <= JSTAB) THEN
        ! --- STABILITY CHECK
            DEL1 = 0.D0
            DO I = 1, N
                DEL1 = DEL1 + (DZ(I)/SCAL(I))**2
            END DO
            DEL2 = 0.D0
            DO I = 1, N
                DEL2 = DEL2 + ((DY(I)-DZ(I))/SCAL(I))**2
            END DO
            QUOT = DEL2/MAX(UROUND, DEL1)
            IF (QUOT > 4.D0) THEN
                NFCN = NFCN + 1
                GOTO 79
            END IF
        END IF
    END DO
    
! --- FINAL SMOOTHING STEP
    CALL FCN(N, X+H, YH2, DY, RPAR, IPAR)
    IF (IOUT >= 2 .AND. NJMID <= 2*J-1) THEN
        IPT = IPT + 1
        DO I = 1, NRD
            FSAFE(IPT,I) = DY(ICOMP(I))
        END DO
    END IF
    DO I = 1, N
        T(J,I) = (YH1(I) + YH2(I) + HJ*DY(I))/2.D0
    END DO
    NFCN = NFCN + NJ(J)
    
! --- POLYNOMIAL EXTRAPOLATION
    IF (J == 1) RETURN
    DBLENJ = NJ(J)
    DO L = J, 2, -1
        FAC = (DBLENJ/NJ(L-1))**2 - 1.D0
        DO I = 1, N
            T(L-1,I) = T(L,I) + (T(L,I)-T(L-1,I))/FAC
        END DO
    END DO
    ERR = 0.D0
    
! --- SCALING
    DO I = 1, N
        T1I = MAX(ABS(Y(I)), ABS(T(1,I)))
        IF (ITOL == 0) THEN
            SCAL(I) = ATOL(1) + RTOL(1)*T1I
        ELSE
            SCAL(I) = ATOL(I) + RTOL(I)*T1I
        END IF
        ERR = ERR + ((T(1,I)-T(2,I))/SCAL(I))**2
    END DO
    ERR = SQRT(ERR/N)
    IF (ERR*UROUND >= 1.D0) GOTO 79
    IF (J > 2 .AND. ERR >= ERROLD) GOTO 79
    ERROLD = MAX(4*ERR, 1.D0)
    
! --- COMPUTE OPTIMAL STEPSIZES
    EXPO = 1.D0/(2*J-1)
    FACMIN = FAC1**EXPO
    FAC = MIN(FAC2/FACMIN, MAX(FACMIN, (ERR/SAFE1)**EXPO/SAFE2))
    FAC = 1.D0/FAC
    HH(J) = MIN(ABS(H)*FAC, HMAX)
    W(J) = A(J)/HH(J)
    RETURN
    
79  ATOV = .TRUE.
    H = H*SAFE3
    REJECT = .TRUE.
    RETURN
    END SUBROUTINE MIDEX

    SUBROUTINE INTERP(N, Y, IMIT)
! --- COMPUTES THE COEFFICIENTS OF THE INTERPOLATION FORMULA
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, IMIT
    DOUBLE PRECISION, INTENT(INOUT) :: Y(N*(IMIT+5))
    
    ! Local variables
    INTEGER :: I, IM
    DOUBLE PRECISION :: Y0, Y1, YP0, YP1, YDIFF, ASPL, BSPL
    DOUBLE PRECISION :: PH0, PH1, PH2, PH3, FAC1, FAC2
    DOUBLE PRECISION :: A(0:30)
    
! --- BEGIN WITH HERMITE INTERPOLATION
    DO I = 1, N
        Y0 = Y(I)
        Y1 = Y(2*N+I)
        YP0 = Y(N+I)
        YP1 = Y(3*N+I)
        YDIFF = Y1 - Y0
        ASPL = -YP1 + YDIFF
        BSPL = YP0 - YDIFF
        Y(N+I) = YDIFF
        Y(2*N+I) = ASPL
        Y(3*N+I) = BSPL
        IF (IMIT < 0) GOTO 100
        
    ! --- COMPUTE THE DERIVATIVES OF HERMITE AT MIDPOINT
        PH0 = (Y0+Y1)*0.5D0 + 0.125D0*(ASPL+BSPL)
        PH1 = YDIFF + (ASPL-BSPL)*0.25D0
        PH2 = -(YP0-YP1)
        PH3 = 6.D0*(BSPL-ASPL)
        
    ! --- COMPUTE THE FURTHER COEFFICIENTS
        IF (IMIT < 1) GOTO 20
        A(1) = 16.D0*(Y(5*N+I)-PH1)
        IF (IMIT < 3) GOTO 20
        A(3) = 16.D0*(Y(7*N+I)-PH3+3*A(1))
        IF (IMIT < 5) GOTO 20
        DO IM = 5, IMIT, 2
            FAC1 = IM*(IM-1)/2.D0
            FAC2 = FAC1*(IM-2)*(IM-3)*2.D0
            A(IM) = 16.D0*(Y((IM+4)*N+I)+FAC1*A(IM-2)-FAC2*A(IM-4))
        END DO
        
20      CONTINUE
        A(0) = (Y(4*N+I)-PH0)*16.D0
        IF (IMIT < 2) GOTO 60
        A(2) = (Y(N*6+I)-PH2+A(0))*16.D0
        IF (IMIT < 4) GOTO 60
        DO IM = 4, IMIT, 2
            FAC1 = IM*(IM-1)/2.D0
            FAC2 = IM*(IM-1)*(IM-2)*(IM-3)
            A(IM) = (Y(N*(IM+4)+I)+A(IM-2)*FAC1-A(IM-4)*FAC2)*16.D0
        END DO
        
60      CONTINUE
        DO IM = 0, IMIT
            Y(N*(IM+4)+I) = A(IM)
        END DO
100 END DO
    RETURN
    END SUBROUTINE INTERP

    FUNCTION CONTEX(II, X, Y, NCON, ICOMP, N)
! ----------------------------------------------------------
!     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT IN CONECTION
!     WITH THE OUTPUT-SUBROUTINE FOR ODEX. IT PROVIDES AN
!     APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X.
! ----------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: II, NCON, N
    DOUBLE PRECISION, INTENT(IN) :: X
    DOUBLE PRECISION, INTENT(IN) :: Y(NCON)
    INTEGER, INTENT(IN) :: ICOMP(N)
    DOUBLE PRECISION :: CONTEX
    
    ! Local variables
    INTEGER :: I, J, IM
    DOUBLE PRECISION :: THETA, THETA1, PHTHET, THETAH
    
    COMMON /CONODX/XOLD, H, IMIT
    DOUBLE PRECISION :: XOLD, H
    INTEGER :: IMIT
    
! ----- COMPUTE PLACE OF II-TH COMPONENT
    I = 0
    DO J = 1, N
        IF (ICOMP(J) == II) I = J
    END DO
    IF (I == 0) THEN
        WRITE(6, *) ' NO DENSE OUTPUT AVAILABLE FOR COMP.', II
        CONTEX = 0.D0
        RETURN
    END IF
    
! ----- COMPUTE THE INTERPOLATED VALUE
    THETA = (X-XOLD)/H
    THETA1 = 1.D0 - THETA
    PHTHET = Y(I) + THETA*(Y(N+I) + THETA1*(Y(2*N+I)*THETA + Y(3*N+I)*THETA1))
    IF (IMIT < 0) THEN
        CONTEX = PHTHET
        RETURN
    END IF
    THETAH = THETA - 0.5D0
    CONTEX = Y(N*(IMIT+4)+I)
    DO IM = IMIT, 1, -1
        CONTEX = Y(N*(IM+3)+I) + CONTEX*THETAH/IM
    END DO
    CONTEX = PHTHET + (THETA*THETA1)**2*CONTEX
    RETURN
    END FUNCTION CONTEX