MODULE LinspaceModule

    CONTAINS

    FUNCTION linspace(start, stop, array) RESULT(arrayFinal)
        DOUBLE PRECISION, INTENT(IN)  :: start, stop, array(:)
        DOUBLE PRECISION, allocatable :: arrayFinal(:)

        INTEGER :: num, i
        DOUBLE PRECISION :: step

        num = SIZE(array)
        IF (num <= 1) THEN
            PRINT *, "Error: Size of array must be greater than 1 in linspace subroutine."
            RETURN
        END IF

        step = (stop - start) / DBLE(num - 1)
        allocate(arrayFinal(num))
        DO i = 1, num
            arrayFinal(i) = start + (i - 1) * step
        END DO
    END FUNCTION linspace
END MODULE LinspaceModule

PROGRAM TestLinspace
    USE LinspaceModule
    DOUBLE PRECISION :: a, b, c(10), x(10)
    a = 0.0D0
    b = 1.0D0
    call random_number(c)
    x =  linspace(a, b, c)
    PRINT *, x
END PROGRAM TestLinspace

!I understand the confusion. The difference lies in the scoping rules of Fortran.
!
!In the case where linspace is within a module (LinspaceModule), it automatically has an explicit interface, and any program or subroutine that USES the module can access that interface. This is Fortran's way of encapsulating procedures and providing a clean way to share them among different parts of your code.
!
!In the second case, where linspace is outside any module, Fortran relies on explicit interfaces for procedures. When a subroutine is called, the compiler needs to know the number and types of arguments it takes. Without an explicit interface, the compiler cannot check if the arguments passed in the call match the subroutine's expectations.
!
!To make it work without a module, you can provide an explicit interface for linspace before the PROGRAM statement in the main program. This informs the compiler about the subroutine's interface. Here's the modified version: