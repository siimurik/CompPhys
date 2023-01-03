!==================================================================
!   $ gfortran dgemv_solver.f90 -o matvec -llapack -lblas
!   $ ./matvec
!==================================================================
! DGEMV  performs one of the matrix-vector operations
!
!    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!
! where alpha and beta are scalars, x and y are vectors and A is an
! m by n matrix.
!==================================================================
program matvec
    implicit none
    integer          N
    parameter        ( N = 4 )
    integer          LDA, INCX, INCY, i, j
    parameter        ( LDA = N )
    double precision A( LDA, N )
    double precision ALPHA, BETA
    double precision X(N), Y(N)
    double precision start_time, end_time, elapsed_time

    ALPHA = 1.0
    BETA  = 1.0
    INCX  = 1
    INCY  = INCX

    ! initialize the matrix A and vectors x and y
    do i = 1, N
        do j = 1, N
        A(i,j) = (i-1)*N - j
        end do
        x(i) = 1.0*i
        y(i) = 0.0
    end do

    ! print the matrix A
    write (*,10)
    do i = 1, N
        write (*,11)  (A(i,j),j = 1, N)
    end do
    10 format (' Matrix A')
    11 format (6f12.6)

    ! print vector x
    write(*,*)'DGEMV Program Results'
    do i = 1, N
        write(*,'(6g12.5)') X(i)
    end do
    !write(*,*) '    ...'
    !do i = N-2, N
    !    write(*,'(6g12.5)') X(i)
    !end do

    ! record the starting time
    call cpu_time(start_time)

    ! perform the matrix-vector multiplication y = A*x
    call dgemv( 'N', N, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )

    ! record the ending time
    call cpu_time(end_time)

    ! printing out the first and last 3 elements
    write(*,*)'DGEMV Program Results'
    do i = 1, N
        write(*,'(6g12.5)') Y(i)
    end do
    !write(*,*) '    ...'
    !do i = N-2, N
    !    write(*,'(6g12.5)') Y(i)
    !end do

    ! compute the elapsed time
    elapsed_time = end_time - start_time

    ! print the elapsed time
    write(*,20) elapsed_time
20  format ( 'Elapsed time: ', F8.6, ' seconds.')
end program matvec
