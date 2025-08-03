*DECK DVJAC
      SUBROUTINE dvjac (Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM, F, JAC,
     1                 ierpj, rpar, ipar)
      EXTERNAL f, jac
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM, RPAR
      INTEGER LDYH, IWM, IERPJ, IPAR
      dimension y(*), yh(ldyh,*), ewt(*), ftem(*), savf(*),
     1   wm(*), iwm(*), rpar(*), ipar(*)
C-----------------------------------------------------------------------
C Call sequence input -- Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM,
C                        F, JAC, RPAR, IPAR
C Call sequence output -- WM, IWM, IERPJ
C COMMON block variables accessed:
C     /DVOD01/  CCMXJ, DRC, H, RL1, TN, UROUND, ICF, JCUR, LOCJS,
C               MITER, MSBJ, N, NSLJ
C     /DVOD02/  NFE, NST, NJE, NLU
C
C Subroutines called by DVJAC: F, JAC, DACOPY, DCOPY, DGBFA, DGEFA,
C                              DSCAL
C Function routines called by DVJAC: DVNORM
C-----------------------------------------------------------------------
C DVJAC is called by DVNLSD to compute and process the matrix
C P = I - h*rl1*J , where J is an approximation to the Jacobian.
C Here J is computed by the user-supplied routine JAC if
C MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
C If MITER = 3, a diagonal approximation to J is used.
C If JSV = -1, J is computed from scratch in all cases.
C If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is
C considered acceptable, then P is constructed from the saved J.
C J is stored in wm and replaced by P.  If MITER .ne. 3, P is then
C subjected to LU decomposition in preparation for later solution
C of linear systems with P as coefficient matrix. This is done
C by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
C
C Communication with DVJAC is done with the following variables.  (For
C more details, please see the comments in the driver subroutine.)
C Y          = Vector containing predicted values on entry.
C YH         = The Nordsieck array, an LDYH by LMAX array, input.
C LDYH       = A constant .ge. N, the first dimension of YH, input.
C EWT        = An error weight vector of length N.
C SAVF       = Array containing f evaluated at predicted y, input.
C WM         = Real work space for matrices.  In the output, it containS
C              the inverse diagonal matrix if MITER = 3 and the LU
C              decomposition of P if MITER is 1, 2 , 4, or 5.
C              Storage of matrix elements starts at WM(3).
C              Storage of the saved Jacobian starts at WM(LOCJS).
C              WM also contains the following matrix-related data:
C              WM(1) = SQRT(UROUND), used in numerical Jacobian step.
C              WM(2) = H*RL1, saved for later use if MITER = 3.
C IWM        = Integer work space containing pivot information,
C              starting at IWM(31), if MITER is 1, 2, 4, or 5.
C              IWM also contains band parameters ML = IWM(1) and
C              MU = IWM(2) if MITER is 4 or 5.
C F          = Dummy name for the user supplied subroutine for f.
C JAC        = Dummy name for the user supplied Jacobian subroutine.
C RPAR, IPAR = Dummy names for user's real and integer work arrays.
C RL1        = 1/EL(2) (input).
C IERPJ      = Output error flag,  = 0 if no trouble, 1 if the P
C              matrix is found to be singular.
C JCUR       = Output flag to indicate whether the Jacobian matrix
C              (or approximation) is now current.
C              JCUR = 0 means J is not current.
C              JCUR = 1 means J is current.
C-----------------------------------------------------------------------
C
      include "vode.H"
C
C Type declarations for local variables --------------------------------
C
      DOUBLE PRECISION CON, DI, FAC, HRL1, ONE, PT1, R, R0, SRUR, THOU,
     1     yi, yj, yjj, zero
      INTEGER I, I1, I2, IER, II, J, J1, JJ, JOK, LENP, MBA, MBAND,
     1        meb1, meband, ml, ml3, mu, np1
C
C Type declaration for function subroutines called ---------------------
C
      DOUBLE PRECISION DVNORM
      parameter(one = 1.0d0, thou = 1000.0d0, zero = 0.0d0, pt1 = 0.1d0)
C
      ierpj = 0
      hrl1 = h*rl1
C See whether J should be evaluated (JOK = -1) or not (JOK = 1). -------
      jok = jsv
      IF (jsv .EQ. 1) THEN
        IF (nst .EQ. 0 .OR. nst .GT. nslj+msbj) jok = -1
        IF (icf .EQ. 1 .AND. drc .LT. ccmxj) jok = -1
        IF (icf .EQ. 2) jok = -1
      ENDIF
C End of setting JOK. --------------------------------------------------
C
      IF (jok .EQ. -1 .AND. miter .EQ. 1) THEN
C If JOK = -1 and MITER = 1, call JAC to evaluate Jacobian. ------------
      nje = nje + 1
      nslj = nst
      jcur = 1
      lenp = n*n
      DO 110 i = 1,lenp
 110    wm(i+2) = zero
      CALL jac (n, tn, y, 0, 0, wm(3), n, rpar, ipar)
      IF (jsv .EQ. 1) CALL dcopy (lenp, wm(3), 1, wm(locjs), 1)
      ENDIF
C
      IF (jok .EQ. -1 .AND. miter .EQ. 2) THEN
C If MITER = 2, make N calls to F to approximate the Jacobian. ---------
      nje = nje + 1
      nslj = nst
      jcur = 1
      fac = dvnorm(n, savf, ewt)
      r0 = thou*abs(h)*uround*REAL(n)*FAC
      IF (r0 .EQ. zero) r0 = one
      srur = wm(1)
      j1 = 2
      DO 230 j = 1,n
        yj = y(j)
        r = max(srur*abs(yj),r0/ewt(j))
        y(j) = y(j) + r
        fac = one/r
        CALL f (n, tn, y, ftem, rpar, ipar)
        DO 220 i = 1,n
 220      wm(i+j1) = (ftem(i) - savf(i))*fac
        y(j) = yj
        j1 = j1 + n
 230    CONTINUE
      nfe = nfe + n
      lenp = n*n
      IF (jsv .EQ. 1) CALL dcopy (lenp, wm(3), 1, wm(locjs), 1)
      ENDIF
C
      IF (jok .EQ. 1 .AND. (miter .EQ. 1 .OR. miter .EQ. 2)) THEN
      jcur = 0
      lenp = n*n
      CALL dcopy (lenp, wm(locjs), 1, wm(3), 1)
      ENDIF
C
      IF (miter .EQ. 1 .OR. miter .EQ. 2) THEN
C Multiply Jacobian by scalar, add identity, and do LU decomposition. --
      con = -hrl1
      CALL dscal (lenp, con, wm(3), 1)
      j = 3
      np1 = n + 1
      DO 250 i = 1,n
        wm(j) = wm(j) + one
 250    j = j + np1
      nlu = nlu + 1
      CALL dgefa (wm(3), n, n, iwm(31), ier)
      IF (ier .NE. 0) ierpj = 1
      RETURN
      ENDIF
C End of code block for MITER = 1 or 2. --------------------------------
C
      IF (miter .EQ. 3) THEN
C If MITER = 3, construct a diagonal approximation to J and P. ---------
      nje = nje + 1
      jcur = 1
      wm(2) = hrl1
      r = rl1*pt1
      DO 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      CALL f (n, tn, y, wm(3), rpar, ipar)
      nfe = nfe + 1
      DO 320 i = 1,n
        r0 = h*savf(i) - yh(i,2)
        di = pt1*r0 - h*(wm(i+2) - savf(i))
        wm(i+2) = one
        IF (abs(r0) .LT. uround/ewt(i)) GO TO 320
        IF (abs(di) .EQ. zero) GO TO 330
        wm(i+2) = pt1*r0/di
 320    CONTINUE
      RETURN
 330  ierpj = 1
      RETURN
      ENDIF
C End of code block for MITER = 3. -------------------------------------
C
C Set constants for MITER = 4 or 5. ------------------------------------
      ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
C
      IF (jok .EQ. -1 .AND. miter .EQ. 4) THEN
C If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian. ------------
      nje = nje + 1
      nslj = nst
      jcur = 1
      DO 410 i = 1,lenp
 410    wm(i+2) = zero
      CALL jac (n, tn, y, ml, mu, wm(ml3), meband, rpar, ipar)
      IF (jsv .EQ. 1)
     1   CALL dacopy (mband, n, wm(ml3), meband, wm(locjs), mband)
      ENDIF
C
      IF (jok .EQ. -1 .AND. miter .EQ. 5) THEN
C If MITER = 5, make ML+MU+1 calls to F to approximate the Jacobian. ---
      nje = nje + 1
      nslj = nst
      jcur = 1
      mba = min(mband,n)
      meb1 = meband - 1
      srur = wm(1)
      fac = dvnorm(n, savf, ewt)
      r0 = thou*abs(h)*uround*REAL(n)*FAC
      IF (r0 .EQ. zero) r0 = one
      DO 560 j = 1,mba
        DO 530 i = j,n,mband
          yi = y(i)
          r = max(srur*abs(yi),r0/ewt(i))
 530      y(i) = y(i) + r
        CALL f (n, tn, y, ftem, rpar, ipar)
        DO 550 jj = j,n,mband
          y(jj) = yh(jj,1)
          yjj = y(jj)
          r = max(srur*abs(yjj),r0/ewt(jj))
          fac = one/r
          i1 = max(jj-mu,1)
          i2 = min(jj+ml,n)
          ii = jj*meb1 - ml + 2
          DO 540 i = i1,i2
 540        wm(ii+i) = (ftem(i) - savf(i))*fac
 550      CONTINUE
 560    CONTINUE
      nfe = nfe + mba
      IF (jsv .EQ. 1)
     1   CALL dacopy (mband, n, wm(ml3), meband, wm(locjs), mband)
      ENDIF
C
      IF (jok .EQ. 1) THEN
      jcur = 0
      CALL dacopy (mband, n, wm(locjs), mband, wm(ml3), meband)
      ENDIF
C
C Multiply Jacobian by scalar, add identity, and do LU decomposition.
      con = -hrl1
      CALL dscal (lenp, con, wm(3), 1 )
      ii = mband + 2
      DO 580 i = 1,n
        wm(ii) = wm(ii) + one
 580    ii = ii + meband
      nlu = nlu + 1
      CALL dgbfa (wm(3), meband, n, ml, mu, iwm(31), ier)
      IF (ier .NE. 0) ierpj = 1
      RETURN
C End of code block for MITER = 4 or 5. --------------------------------
C
C----------------------- End of Subroutine DVJAC -----------------------
      END
