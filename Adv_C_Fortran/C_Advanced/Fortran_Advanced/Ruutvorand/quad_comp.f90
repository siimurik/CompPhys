! Quadratic formula a*x**2 + b*x + c = 0 solver 
! for REAL and COMPLEX roots.
!-------------------------------------------------
! Compile with:
!   > gfortran quad_comp.f90 -o comp
!   >./comp
!-------------------------------------------------
PROGRAM quadratic
    IMPLICIT NONE 
    REAL(8)     :: a, b, c, x1_real, x2_real, delta
    COMPLEX(8)  :: x1_num, x2_num, x_denum
    COMPLEX(8)  :: z1, z2, D
    complex, parameter :: i = (0, 1)   ! sqrt(-1)

    PRINT *,"Enter the a,b,and c" 
    READ *, a, b, c

    if (a .eq. 0.0) then
        print *, "Boy that ain't a quadratic anymore..."
        print *, "anyway, the asnwer is x = ", -c/b
        stop
    end if

    delta = b**2 - 4.0 * a * c

    if (delta < 0) then
        PRINT *,"The roots are complex"
        call comproots(a,b,c,z1,z2)
        PRINT *, "Z1 = ", z1
        print *, "Z2 = ", z2
    else if (delta > 0) then
        PRINT *,"The roots are real"
        call comproots(a,b,c,z1,z2)
        x1_real = z1
        x2_real = z2
        PRINT *, "X1 = ", x1_real
        print *, "X2 = ", x2_real
    else
        x1_real = -b/(2.D0*a) 
        PRINT *, 'Double Root:		x1=', x1_real
    end if
END PROGRAM quadratic

subroutine comproots(a,b,c,z1,z2)
    implicit none
    REAL(8)     :: a,b,c, delta
    COMPLEX(8)  :: x1_num, x2_num, x_denum
    COMPLEX(8)  :: z1, z2, D

    D = b**2 - 4.0 * a * c
    x1_num = (-b + SQRT(D))
    x2_num = (-b - SQRT(D))
    x_denum = 2.D0*A
    z1 = x1_num/x_denum
    z2 = x2_num/x_denum

end subroutine comproots



