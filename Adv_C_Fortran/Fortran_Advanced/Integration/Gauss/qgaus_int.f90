!==============================================
! Compile and execute with:
!   $ gfortran qgaus_int.f90 -o qg
!   $ ./qg
!==============================================
PROGRAM main
    DOUBLE PRECISION           :: a, b, s
    DOUBLE PRECISION, EXTERNAL :: func

    a = -2.0
    b =  2.0
    
    call qgaus(func, a, b, s)

    print *, 'Result: ', s

END PROGRAM main

FUNCTION func(x)
    DOUBLE PRECISION func, x

    func = exp(-x**2)
    return
END FUNCTION func

SUBROUTINE qgaus(func,a,b,ss)
    DOUBLE PRECISION :: a,b,ss
    DOUBLE PRECISION, EXTERNAL :: func
    !Returns as ss the integral of the function func between a and b, by 64-point Gauss-
    !Legendre integration: the function is evaluated exactly 64 times at interior points in the
    !range of integration.
    INTEGER, PARAMETER  :: NDIM = 32
    INTEGER             :: j
    DOUBLE PRECISION    :: dx,xm,xr,w(NDIM),x(NDIM) ! The weights and abscissae
    SAVE w, x
    DATA x/ 0.0243502926634244D0,    0.0729931217877990D0,    0.1214628192961206D0,    0.1696444204239928D0, &
            0.2174236437400071D0,    0.2646871622087674D0,    0.3113228719902110D0,    0.3572201583376681D0, &
            0.4022701579639916D0,    0.4463660172534641D0,    0.4894031457070530D0,    0.5312794640198946D0, &
            0.5718956462026340D0,    0.6111553551723933D0,    0.6489654712546573D0,    0.6852363130542333D0, &
            0.7198818501716109D0,    0.7528199072605319D0,    0.7839723589433414D0,    0.8132653151227975D0, &
            0.8406292962525803D0,    0.8659993981540928D0,    0.8893154459951141D0,    0.9105221370785028D0, &
            0.9295691721319396D0,    0.9464113748584028D0,    0.9610087996520538D0,    0.9733268277899110D0, &
            0.9833362538846260D0,    0.9910133714767443D0,    0.9963401167719553D0,    0.9993050417357722D0  /
    DATA w/ 0.0486909570091397D0,    0.0485754674415034D0,    0.0483447622348030D0,    0.0479993885964583D0, &
            0.0475401657148303D0,    0.0469681828162100D0,    0.0462847965813144D0,    0.0454916279274181D0, &
            0.0445905581637566D0,    0.0435837245293235D0,    0.0424735151236536D0,    0.0412625632426235D0, &
            0.0399537411327203D0,    0.0385501531786156D0,    0.0370551285402400D0,    0.0354722132568824D0, &
            0.0338051618371416D0,    0.0320579283548516D0,    0.0302346570724025D0,    0.0283396726142595D0, &
            0.0263774697150547D0,    0.0243527025687109D0,    0.0222701738083833D0,    0.0201348231535302D0, &
            0.0179517157756973D0,    0.0157260304760247D0,    0.0134630478967186D0,    0.0111681394601311D0, &
            0.0088467598263639D0,    0.0065044579689784D0,    0.0041470332605625D0,    0.0017832807216964D0/
    xm = 0.5*(b+a)
    xr = 0.5*(b-a)
    ss = 0          ! Will be twice the average value of the function, since the 64
    do j = 1, NDIM  ! weights (five numbers above each used twice) sum to 2.
        dx = xr*x(j)
        ss = ss+w(j)*(func(xm+dx)+func(xm-dx))
    end do
    ss = xr*ss !Scale the answer to the range of integration.
    return
END SUBROUTINE