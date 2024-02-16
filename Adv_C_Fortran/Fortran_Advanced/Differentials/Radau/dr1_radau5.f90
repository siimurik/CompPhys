! Compile and execute with:
!    $ gfortran dr1_radau5.f90 radau5.f radaua.f -o vdp
!    $ ./vdp
! * * * * * * * * * * * * * * * * * * * * * * * * *
! --- DRIVER FOR RADAU5 AT VAN DER POL'S EQUATION
! * * * * * * * * * * * * * * * * * * * * * * * * *
! USES radau5.f and radaua.f
PROGRAM main
    !REAL(kind=8) :: A-H,O-Z
    IMPLICIT REAL(kind=8) (A-H,O-Z)
    ! --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
    PARAMETER  (ND=2,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
    DIMENSION  Y(ND),WORK(LWORK),IWORK(LIWORK)
    EXTERNAL   FVPOL,JVPOL,SOLOUT
    ! --- PARAMETER IN THE DIFFERENTIAL EQUATION
    RPAR  = 1.0D-6
    ! --- DIMENSION OF THE SYSTEM
    N     = 2
    ! --- COMPUTE THE JACOBIAN ANALYTICALLY
    IJAC  = 1
    ! --- JACOBIAN IS A FULL MATRIX
    MLJAC = N
    ! --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
    IMAS  = 0
    ! --- OUTPUT ROUTINE IS USED DURING INTEGRATION
    IOUT  = 1
    ! --- INITIAL VALUES
    X    = 0.0D0        ! t
    Y(1) = 2.0D0        ! x
    Y(2) = -0.66D0      ! y
    ! --- ENDPOINT OF INTEGRATION
    XEND = 2.0D0
    ! --- REQUIRED TOLERANCE
    RTOL = 1.0D-14      ! max 1.0D-14
    ATOL = 1.0D0*RTOL
    ITOL = 0
    ! --- INITIAL STEP SIZE
    H    = 1.0D-6 
    ! --- SET DEFAULT VALUES 
    DO I = 1, 20
        IWORK(I) = 0
        WORK(I)  = 0.D0
    END DO
    ! --- CALL OF THE SUBROUTINE RADAU5
    CALL RADAU5(N,FVPOL,X,Y,XEND,H,             &
                RTOL,ATOL,ITOL,                 &
                JVPOL,IJAC,MLJAC,MUJAC,         &
                FVPOL,IMAS,MLMAS,MUMAS,         &
                SOLOUT,IOUT,                    &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
    ! --- PRINT FINAL SOLUTION
    WRITE (6,99) X,Y(1),Y(2)
99  FORMAT(1X,'X =',F6.3,'    Y =',2E18.10)
    ! --- PRINT STATISTICS
    WRITE (6,90) RTOL
90  FORMAT('       rtol=',D8.2)
    WRITE (6,91) (IWORK(J),J=14,20)
91  FORMAT(' fcn=',I5,' jac=',I5,' step=',I5,' accpt=',I5, &
           ' rejct=',I3,' dec=',I5,' sol=',I5)
    STOP
END PROGRAM main
!
!
SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,DT,LRC,N,RPAR,IPAR,IRTRN)
    ! --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS BY USING "CONTR5"
    IMPLICIT REAL(kind=8) (A-H,O-Z)
    DIMENSION Y(N),CONT(LRC)
    COMMON /INTERN/XOUT
    IF (NR.EQ.1) THEN
        WRITE (6,99) X,Y(1),Y(2),NR-1
        DT   = 1.D-03    ! Step size
        XOUT = DT
    ELSE
10  CONTINUE
        IF (X.GE.XOUT) THEN
        ! --- CONTINUOUS OUTPUT FOR RADAU5
            WRITE (6,99) XOUT,CONTR5(1,XOUT,CONT,LRC), &
                    CONTR5(2,XOUT,CONT,LRC),NR-1
            XOUT=XOUT+DT
            GOTO 10
        END IF
    END IF
99  FORMAT(1X,'X =',F6.3,'    Y =',2E18.10,'    NSTEP =',I5)
    RETURN
END SUBROUTINE
!
!
SUBROUTINE FVPOL(N,X,Y,F,RPAR,IPAR)
    ! --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
    IMPLICIT REAL(kind=8) (A-H,O-Z)
    DIMENSION Y(N),F(N)
    F(1)=Y(2)
    F(2)=((1-Y(1)**2)*Y(2)-Y(1))/RPAR
    RETURN
END SUBROUTINE
!
!
SUBROUTINE JVPOL(N,X,Y,DFY,LDFY,RPAR,IPAR)
    ! --- JACOBIAN OF VAN DER POL'S EQUATION
    IMPLICIT REAL(kind=8) (A-H,O-Z)
    DIMENSION Y(N),DFY(LDFY,N)
    DFY(1,1)=0.0D0
    DFY(1,2)=1.0D0
    DFY(2,1)=(-2.0D0*Y(1)*Y(2)-1.0D0)/RPAR
    DFY(2,2)=(1.0D0-Y(1)**2)/RPAR
    RETURN
END SUBROUTINE

