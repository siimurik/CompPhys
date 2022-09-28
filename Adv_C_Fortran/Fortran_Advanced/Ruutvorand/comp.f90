!Ax^2+Bx+C=0 
PROGRAM delta1 
    IMPLICIT NONE 
    REAL(8) :: A,B,C 
    REAL(8) :: DELTA,X1,X2 
    COMPLEX(8) :: Z1, Z2, z
    complex, parameter :: i = (0, 1)   ! sqrt(-1)
    PRINT *,"Enter the a,b,and c" 
    READ *,A,B,C 
    delta=b**2.D0-4.D0*a*c 
   

    IF (delta>0)THEN 
        PRINT *,"The roots are real" 
        X1=(-B+SQRT(DELTA))/(2.D0*A) 
        X2=(-B-SQRT(DELTA))/(2.D0*A) 
        PRINT *,"X1 = ",X1,"X2 = ",X2 
    ELSE IF(delta==0)THEN 
        PRINT *,"There is one real root" 
        X1 = -B/(2.D0*A) 
        PRINT *,"X = ",X1 
    ELSE

        PRINT *,"The roots are complex" 
        Z1=((-B+SQRT(DELTA))/(2.D0*A), 0.)
        Z2=((-B-SQRT(DELTA))/(2.D0*A), 0.) 
        PRINT *,"X1 = ",z1,"X2 = ",z2
    END IF 
END PROGRAM delta1

