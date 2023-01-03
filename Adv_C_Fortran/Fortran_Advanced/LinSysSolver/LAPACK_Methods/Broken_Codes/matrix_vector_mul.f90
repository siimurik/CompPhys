!   $ gfortran matrix_vector_mul.f90 -o mvm -llapack -lblas
program matrix_vector_mult
    implicit none

    integer, parameter :: n = 3, m = 3  ! Dimensions of matrix A
    real*8, dimension(n,m) :: A          ! Matrix A
    real*8, dimension(m) :: x            ! Vector x
    real*8, dimension(n) :: y            ! Vector y

    ! Initialize matrix A and vector x
    A = reshape((/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 /), shape(A))
    x = (/ 1.0, 2.0, 3.0 /)

    ! Perform matrix-vector multiplication using dgemv from LAPACK
    y = matmul(A,x)

    ! Print result
    write(*,'(3F8.2)') y

end program matrix_vector_mult

