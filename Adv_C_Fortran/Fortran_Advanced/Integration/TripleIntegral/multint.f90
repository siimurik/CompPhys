!==============================================
! Compile and execute with:
!   $ gfortran multint.f90 -o mi
!   $ ./mi
!==============================================
PROGRAM multint
    REAL :: x1, x2, result
    ! Set integration limits
    x1 = -1.0
    x2 =  1.0

    ! Call quad3d to perform the integration
    CALL quad3d(x1, x2, result)

    ! Display the result
    WRITE(*,*) 'Result of 3D integration:', result
END PROGRAM multint

SUBROUTINE qgaus(func,a,b,ss)
    REAL :: a,b,ss
    REAL, EXTERNAL :: func
    !Returns as ss the integral of the function func between a and b, by 64-point Gauss-
    !Legendre integration: the function is evaluated exactly 64 times at interior points in the
    !range of integration.
    INTEGER, PARAMETER  :: NDIM = 32
    INTEGER             :: j
    REAL    :: dx,xm,xr,w(NDIM),x(NDIM) ! The weights and abscissae
    SAVE w, x
    DATA x/ 0.0243502926634244,    0.0729931217877990,    0.1214628192961206,    0.1696444204239928, &
            0.2174236437400071,    0.2646871622087674,    0.3113228719902110,    0.3572201583376681, &
            0.4022701579639916,    0.4463660172534641,    0.4894031457070530,    0.5312794640198946, &
            0.5718956462026340,    0.6111553551723933,    0.6489654712546573,    0.6852363130542333, &
            0.7198818501716109,    0.7528199072605319,    0.7839723589433414,    0.8132653151227975, &
            0.8406292962525803,    0.8659993981540928,    0.8893154459951141,    0.9105221370785028, &
            0.9295691721319396,    0.9464113748584028,    0.9610087996520538,    0.9733268277899110, &
            0.9833362538846260,    0.9910133714767443,    0.9963401167719553,    0.9993050417357722  /
    DATA w/ 0.0486909570091397,    0.0485754674415034,    0.0483447622348030,    0.0479993885964583, &
            0.0475401657148303,    0.0469681828162100,    0.0462847965813144,    0.0454916279274181, &
            0.0445905581637566,    0.0435837245293235,    0.0424735151236536,    0.0412625632426235, &
            0.0399537411327203,    0.0385501531786156,    0.0370551285402400,    0.0354722132568824, &
            0.0338051618371416,    0.0320579283548516,    0.0302346570724025,    0.0283396726142595, &
            0.0263774697150547,    0.0243527025687109,    0.0222701738083833,    0.0201348231535302, &
            0.0179517157756973,    0.0157260304760247,    0.0134630478967186,    0.0111681394601311, &
            0.0088467598263639,    0.0065044579689784,    0.0041470332605625,    0.0017832807216964/
    xm = 0.5*(b+a)
    xr = 0.5*(b-a)
    ss = 0.0          ! Will be twice the average value of the function, since the 64
    do j = 1, NDIM  ! weights (five numbers above each used twice) sum to 2.
        dx = xr*x(j)
        ss = ss+w(j)*(func(xm+dx)+func(xm-dx))
    end do
    ss = xr*ss !Scale the answer to the range of integration.
    return
END SUBROUTINE

SUBROUTINE quad3d(x1,x2,ss)
    REAL           :: ss,x1,x2
    REAL, EXTERNAL :: h
    !C USES h,qgausx
    !Returns as ss the integral of a user-supplied function func over a three-dimensional region
    !specified by the limits x1, x2, and by the user-supplied functions y1, y2, z1, and z2, as
    !defined in (4.6.2).
    call qgaus(h,x1,x2,ss)
    return
END SUBROUTINE quad3d
    

FUNCTION func(x,y,z)
    REAL :: func, x, y, z

    func = x*x + y*y - z*z
    return

END FUNCTION func

FUNCTION y1(x)
    REAL :: y1, x
    ! Define the lower limit of integration along the y-axis
    y1 = -sqrt(1.0 - x*x) ! Adjust this based on your problem
    RETURN
END FUNCTION y1

FUNCTION y2(x)
    REAL :: y2, x
    ! Define the upper limit of integration along the y-axis
    y2 = sqrt(1.0 - x*x)  ! Adjust this based on your problem
    RETURN
END FUNCTION y2

FUNCTION z1(x, y)
    REAL :: z1, x, y
    ! Define the lower limit of integration along the z-axis
    z1 = -sqrt(1.0 - x*x - y*y)  ! Adjust this based on your problem
    RETURN
END FUNCTION z1

FUNCTION z2(x, y)
    REAL :: z2, x, y
    ! Define the upper limit of integration along the z-axis
    z2 = sqrt(1.0 - x*x - y*y)  ! Adjust this based on your problem
    RETURN
END FUNCTION z2

FUNCTION f(zz)
    REAL f,zz,func,x,y,z
    COMMON /xyz/ x,y,z
    ! USES func
    !Called by qgausz. Calls func.
    z=zz
    f=func(x,y,z)
    return
END FUNCTION f
    
FUNCTION g(yy)
    REAL :: g,yy,z1,z2,x,y,z
    REAL, EXTERNAL :: f
    COMMON /xyz/ x,y,z
    ! USES f,qgausz,z1,z2
    !Called by qgausy. Calls qgausz.
    REAL :: ss
    y=yy
    call qgaus(f,z1(x,y),z2(x,y),ss)
    g=ss
    return
END FUNCTION g

FUNCTION h(xx)
    REAL :: h,xx,y1,y2,x,y,z
    REAL, EXTERNAL :: g
    COMMON /xyz/ x,y,z
    ! USES g,qgausy,y1,y2
    !Called by qgausx. Calls qgausy.
    REAL :: ss
    x=xx
    call qgaus(g,y1(x),y2(x),ss)
    h=ss
    return
END FUNCTION h