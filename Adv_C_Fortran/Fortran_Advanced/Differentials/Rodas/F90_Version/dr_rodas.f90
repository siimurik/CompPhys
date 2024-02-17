! * * * * * * * * * * * * * * * * * * * * * * * * *
! --- DRIVER FOR ROSENBROCK CODE RODAS ON VAN DER POL
! * * * * * * * * * * * * * * * * * * * * * * * * *
! link dr_rodas rodas decsol dc_decsol
! link dr_rodas rodas lapack lapackc dc_lapack
! Compilation commands:
!       $ gfortran -O2 dr_rodas.f90 rodas.f90 decsol.f90 dc_decsol.f90 -o decrod
!       $ gfortran -O2 dr_rodas.f90 rodas.f90 lapack.f90 lapackc.f90 dc_lapack.f90 -o laprod
! NOTE!
! Even with second level compiler optimization, the 'decsol'
! SUBROUTINES are faster than using 'lapack' SUBROUTINES.
    PROGRAM main
    IMPLICIT REAL(KIND=8) (A-H,O-Z)
! --- PARAMETERS FOR RODAS (FULL JACOBIAN)
    PARAMETER (ND=2,LWORK=2*ND*ND+14*ND+20,LIWORK=ND+20)
! --- DECLARATIONS
    DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
    EXTERNAL FVPOL,JVPOL,SOLOUT
! --- DIMENSION OF THE SYSTEM
    RPAR=1.D-6
    N=2
! --- PROBLEM IS AUTONOMOUS
    IFCN=0
! --- COMPUTE THE JACOBIAN ANALYTICALLY
    IJAC=1
! --- JACOBIAN IS A FULL MATRIX
    MLJAC=N
! --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
    IMAS=0
! --- OUTPUT ROUTINE IS USED DURING INTEGRATION
    IOUT=1
! --- INITIAL VALUES
    X=0.0D0
    Y(1)=2.0D0
    Y(2)=-0.66D0
! --- ENDPOINT OF INTEGRATION
    XEND=2.0D0
! --- REQUIRED TOLERANCE
! MAX: RTOL = 1.0D-13, increase NMAX in rodas.f90 if necessary 
    RTOL=1.0D-13    
    ATOL=1.0D-6*RTOL
    ITOL=0
! --- INITIAL STEP SIZE
    H=1.0D-6
! --- SET DEFAULT VALUES
    DO I=1,20
        IWORK(I)=0
        WORK(I)=0.D0
    END DO
! --- CALL OF THE SUBROUTINE RODAS
    CALL RODAS(N,FVPOL,IFCN,X,Y,XEND,H, &
    RTOL,ATOL,ITOL, &
    JVPOL,IJAC,MLJAC,MUJAC,FVPOL,IDFX, &
    FVPOL,IMAS,MLMAS,MUMAS, &
    SOLOUT,IOUT, &
    WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
! --- PRINT FINAL SOLUTION
    WRITE (6,99) X,Y(1),Y(2)
    99 FORMAT(1X,'X =',F5.2,'    Y =',2E18.10)
! --- PRINT STATISTICS
    WRITE (6,90) RTOL
    90 FORMAT('       rtol=',D8.2)
    WRITE (6,91) (IWORK(J),J=14,20)
    91 FORMAT(' fcn=',  I7,' jac=',  I6,' step=',I6, &
              ' accpt=',I6,' rejct=',I4,' dec=', I6, &
              ' sol=',  I7)
    STOP
    END PROGRAM


    SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
! --- PRINTS SOLUTION
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION Y(N),CONT(LRC)
    COMMON /INTERN/XOUT
    IF (NR == 1) THEN
        WRITE (6,99) X,Y(1),Y(2),NR-1
        XOUT=0.2D0
    ELSE
        IF (X >= XOUT) THEN
            Y1=CONTRO(1,XOUT,CONT,LRC)
            Y2=CONTRO(2,XOUT,CONT,LRC)
            WRITE (6,99) XOUT,Y1,Y2,NR-1
            XOUT=XOUT+0.2D0
        END IF
    END IF
    99 FORMAT(1X,'X =',F5.2,'    Y =',2E18.10,'    NSTEP =',I6)
    RETURN
    END SUBROUTINE SOLOUT


    SUBROUTINE FVPOL(N,X,Y,F,RPAR,IPAR)
! --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION Y(N),F(N)
    F(1)=Y(2)
    F(2)=((1-Y(1)**2)*Y(2)-Y(1))/RPAR
    RETURN
    END SUBROUTINE FVPOL


    SUBROUTINE JVPOL(N,X,Y,DFY,LDFY,RPAR,IPAR)
! --- JACOBIAN OF VAN DER POL'S EQUATION
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION Y(N),DFY(LDFY,N)
    DFY(1,1)=0.0D0
    DFY(1,2)=1.0D0
    DFY(2,1)=(-2.0D0*Y(1)*Y(2)-1.0D0)/RPAR
    DFY(2,2)=(1.0D0-Y(1)**2)/RPAR
    RETURN
    END SUBROUTINE JVPOL
