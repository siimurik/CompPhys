!======================================================
! Compile and execute with
!   $ gfortran lu_solver.f90 -o lu_solver -llapack
!   $ ./lu_solver <N> <NRHS>
!
!<N>    - The number of linear equations, i.e., 
!         the order of the matrix A.
!<NRHS> - The number of right hand sides, i.e., 
!         the number of columns of the matrix B. 
!------------------------------------------------------
! Example:
! $ ./lu_solver 4 1
! Enter the values for matrix A (4 rows expected):
! 4, 2, 0, 1,
! 3, 0, 0, 2,
! 0, 1, 1, 1,
! 0, 2, 1, 0
! Enter the values for vector B (4 values expected):
! -1, -0.5, -1, 2
!  Solution vector X:
!   12.50000
!  -16.00000
!   34.00000
!  -19.00000
!======================================================
program lu_solve
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    implicit none
  
    ! Declarations
    character(len=100) :: arg1, arg2
    integer(int32)     :: n, nrhs, lda, ldb, info
    integer(int32), allocatable :: ipiv(:)
    real(real64),   allocatable :: A(:,:), B(:)
    
    ! Command-line input for matrix dimensions
    call get_command_argument(1, arg1)
    call get_command_argument(2, arg2)

    read(arg1, *) n     ! The number of linear equations, i.e., the order of the matrix A.
    read(arg2, *) nrhs  ! The number of right hand sides, i.e., the number of columns of the matrix B.
    
    if (n <= 0 .or. nrhs <= 0) then
        print *, "Matrix dimensions must be positive integers."
        stop
    endif

    lda = n
    ldb = n

    ! Allocate arrays based on input dimensions
    allocate(A(lda, n))
    allocate(B(n))
    allocate(ipiv(n))

    ! Read matrix and vector from user input
    call read_matrix(A, n, lda)
    call read_vector(B, n)

    ! Solve the system of equations A * X = B
    call solve_system(A, B, n, nrhs, ipiv, info)

    ! Check if the solution was successful
    if (info == 0) then
        print *, "Solution vector X:"
        call print_vector(B, n)
    else if (info > 0) then
        print *, "Matrix is singular; solution could not be computed."
    else
        print *, "An error occurred; check the input values."
    endif

    ! Deallocate memory
    deallocate(A)
    deallocate(B)
    deallocate(ipiv)

end program lu_solve

! Subroutine to read matrix values from user input
subroutine read_matrix(A, n, lda)
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    integer(int32), intent(in)  :: n, lda
    real(real64), intent(inout) :: A(lda, n)
    integer(int32) :: i, j

    ! Read matrix A row by row
    write (*, '(a, i1, a)') 'Enter the values for matrix A (', lda, ' rows expected):'
        do i = 1, lda
        read(*,*) (A(i,j), j = 1, lda)  ! Read each row of the matrix
    end do
end subroutine read_matrix

! Subroutine to read vector values from user input
subroutine read_vector(B, n)
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    integer(int32), intent(in)  :: n
    real(real64), intent(inout) :: B(n)
    integer(int32) :: i

    write (*, '(a, i1, a)') 'Enter the values for vector B (', n, ' values expected):'
    read(*,*) (B(i), i = 1, n)
end subroutine read_vector

! Subroutine to solve the system using LAPACK's dgesv
subroutine solve_system(A, B, n, nrhs, ipiv, info)
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    integer(int32), intent(in)  :: n, nrhs
    real(real64), intent(inout) :: A(n, n), B(n)
    integer(int32), intent(out) :: ipiv(n)
    integer(int32), intent(out) :: info

    ! Call LAPACK's dgesv to perform LU decomposition and solve
    ! https://netlib.org/lapack/explore-html-3.6.1/d7/d3b/group__double_g_esolve_ga5ee879032a8365897c3ba91e3dc8d512.html
    call dgesv(n, nrhs, A, n, ipiv, B, n, info)
end subroutine solve_system

! Subroutine to print a vector
subroutine print_vector(B, n)
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    integer(int32), intent(in) :: n
    real(real64), intent(in)   :: B(n)
    integer(int32) :: i

    do i = 1, n
        print '(F10.5)', B(i)
    end do
end subroutine print_vector