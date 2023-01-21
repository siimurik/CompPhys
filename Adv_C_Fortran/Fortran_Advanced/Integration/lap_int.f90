program LapackIntegration
    USE QDAGS_INT
    USE UMACH_INT
    IMPLICIT NONE
    INTEGER    NOUT
    REAL       A, ABS, B, ERRABS, ERREST, ERROR, ERRREL, EXACT, F, RESULT

    INTRINSIC  ABS
    EXTERNAL   F
!                                 Get output unit number
    CALL UMACH (2, NOUT)
!                                 Set limits of integration
    A = 0.0
    B = 1.0
!                                 Set error tolerances

    ERRABS = 0.0
    CALL QDAGS (F, A, B, RESULT, ERRABS=ERRABS, ERREST=ERREST)
!                                 Print results
    EXACT = -4.0
    ERROR = ABS(RESULT-EXACT)
    WRITE (NOUT,99999) RESULT, EXACT, ERREST, ERROR

99999   FORMAT (' Computed =', F8.3, 13X, ' Exact =', F8.3, /, /, &
                ' Error estimate =', 1PE10.3, 6X, 'Error =', 1PE10.3)
END program LapackIntegration

!

    REAL FUNCTION F (X)
    REAL       X
    REAL       ALOG, SQRT

    INTRINSIC  ALOG, SQRT

    F = ALOG(X)/SQRT(X)

    RETURN

    END