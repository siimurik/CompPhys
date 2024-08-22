!=========================================================
!   $ gfortran -o invmat inverse_matrix.f90 -llapack -lblas
!   $ ./invmat > output.csv
! In GNU Plot:
!   $ plot 'output.csv' with lines
!---------------------------------------------------------
!    A*u = y
!   u = A^-1 * y
!=========================================================
program main
    implicit none
    integer, parameter                   :: n = 200
    integer                              :: i
    double precision                     :: h
    double precision, dimension(n+1,n+1) :: A, C
    double precision, dimension(n+1)     :: y, x, u
    double precision :: start_time, end_time, elapsed_time

    h = 2.0/n
    do i = 1, n+1
        x(i) = (i-1)*h
    end do
    A = 0.0
    y = 0.0

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

    ! record the starting time
    call cpu_time(start_time)

    !call inverse(A, C, n+1)! inverse(a,c,n) - slower, uses LU factorization 
    C = inv(A)              ! inv(A) - much faster; uses LAPACK
    ! C(n+1, n+1) = A(n+1,n+1)^-1
    ! u(i)  = C(n+1,n+1) * y(i)

    ! perform matrix-vector multiplication
    ! TODO: write a matrix-vector multiplcation 
    ! algorithm using LAPACK

    ! Basic and slow method of matrix-vector multiplication
    !do i = 1, n+1
    !    sum = 0.0
    !    do j = 1, n+1
    !        sum = sum + C(i,j)*y(j)
    !    end do
    !    u(i) = sum
    !end do

    ! Perform matrix-vector multiplication
    u = matmul(C, y)

    ! record the ending time
    call cpu_time(end_time)

    ! compute the elapsed time
    elapsed_time = end_time - start_time

    WRITE(*,*)'DGETRF Inverse Matrix Program Results'
    do i = 1, 3
        write(*,'(6f12.5)') x(i), u(i)
    end do
    WRITE(*,*) '            ...'
    do i = n-1, n+1
        write(*,'(6f12.5)') x(i), u(i)
    end do

    ! print the elapsed time
    write(*,9) elapsed_time
9   format ( 'Elapsed time: ', F8.6, ' seconds.')
    ! print the inverse matrix C = A^{-1} 
    !write (*,202)
    !do i = 1,n+1
    !    write (*,201)  (c(i,j),j=1,n+1)
    !end do
    !200 format (' Computing Inverse matrix ',/,/, &
    !            ' Matrix A')
    !201 format (6f12.6)
    !202 format (/,' Inverse matrix A^{-1}')
!========================================================================
    ! print the solution
    open(unit=11, file="output.csv")
    do i = 1, n+1
        write(11,*) x(i), u(i)
    end do
    close(11)
!========================================================================
    contains
! -- Returns the inverse of a general squared matrix A
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
    function inv(Amat) result(Ainv)
        double precision, dimension(:,:), intent(in) :: Amat
        double precision, dimension(size(Amat,1),size(Amat,2)) :: Ainv
        double precision, dimension(size(Amat,1)) :: work  ! work array for LAPACK
        integer, dimension(size(Amat,1)) :: ipiv   ! pivot indices
        integer :: m, info

        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI

        ! Store Amat in Ainv to prevent it from being overwritten by LAPACK
        Ainv = Amat
        m = size(Amat,1)

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(m, m, Ainv, m, ipiv, info)

        if (info /= 0) then
            stop 'Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(m, Ainv, m, ipiv, work, m, info)

        if (info /= 0) then
            stop 'Matrix inversion failed!'
        end if
    end function inv
end program main
!========================================================================
subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
    do i=k+1,n
        coeff=a(i,k)/a(k,k)
        L(i,k) = coeff
        do j=k+1,n
            a(i,j) = a(i,j)-coeff*a(k,j)
        end do
    end do
    end do

    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
    L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
    do i=1,j
        U(i,j) = a(i,j)
    end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
    b(k)=1.0
    d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
    do i=2,n
        d(i)=b(i)
        do j=1,i-1
        d(i) = d(i) - L(i,j)*d(j)
        end do
    end do
    ! Step 3b: Solve Ux=d using the back substitution
    x(n)=d(n)/U(n,n)
    do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
        x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
    end do
    ! Step 3c: fill the solutions x(n) into column k of C
    do i=1,n
        c(i,k) = x(i)
    end do
    b(k)=0.0
    end do
end subroutine inverse

