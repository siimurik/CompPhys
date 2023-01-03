program main
implicit none

    integer :: i,j,info1,info2
    integer :: neqn, n ! number of equations
    double precision ,dimension(201,201) :: coeff
    double precision ,dimension (201,1)  :: lhs
    double precision ,dimension (201,1)  :: ipiv
    double precision                     :: h, sum
    double precision, dimension(201,201) :: A
    double precision, dimension(201)     :: y, x, u

    h = 2.0/200.0
    do i = 1, n+1
        x(i) = (i-1)*h
    end do
    A = 0.0
    y = 0.0
    n = 200

    A(1,1) = 1.0/h + h/2*(4.0-x(1))
    A(1,2) = -1.0/h
    y(1)   = h/2*(x(1)+5.0) - 1.0

    do i = 2, n
        A(i,i-1) = -1.0/h
        A(i,i)   =  2.0/h + h*(4.0-x(i))
        A(i,i+1) = -1.0/h
        y(i)     = h*(x(i)+5.0)
    end do
    A(n+1,n)   = -1/h
    A(n+1,n+1) = 1/h + h/2*(4-x(n+1))
    y(n+1)     = h/2*(x(n+1)+5) - 1

    neqn=201

!coeff = reshape( (/4,1,3,-2,3,1,3,-4,2/),(/3,3/))
!lhs = reshape ( (/1,-7,5/),(/3,1/) )

    do i = 1, 201
        lhs(i,1) = y(i)
    enddo

    call SGETRF (neqn,201,A,neqn,ipiv,infO1)
        if (info1==0) then
            call SGETRS ('N',neqn,1,A,neqn,ipiv,lhs,neqn,info2) !Error
        else
        end if

    write (*,*) 'Answer: '
        do j=1,5,1
            write (*,100) lhs(j,1)
            100 format (F12.5,' ')
        end do

    write (*,100) (lhs)

end program