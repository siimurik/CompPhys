!========================================================
! Compile and execute with:
!   $ gfortran -o lapack_bicgstab lapack_bicgstab.f90 -llapack
!   $ ./bicgstab_solver <matrix dim> <tolerance> <max iterations>
!--------------------------------------------------------
! For testing:
! ./bicgstab_solver 4 1.0E-6 1000
! Matrix A:
!   4, 2, 0, 1,
!   3, 0, 0, 2,
!   0, 1, 1, 1,
!   0, 2, 1, 0
! Vector b:
!   -1, -0.5, -1, 2
! Solution x:
! x(  1) =  12.50000000
! x(  2) = -16.00000000
! x(  3) =  34.00000000
! x(  4) = -19.00000000
!========================================================
program BiCGSTABSolver
    implicit none
    logical :: custom_input
    integer :: n, max_iter
    real(8) :: tol
    real(8), allocatable :: A(:,:), b(:), x(:), r(:), r0(:), p(:), v(:), s(:), t(:)
    character(len=100)   :: arg1, arg2, arg3, input_choice
    
    ! Command-line inputs: matrix size, tolerance, and maximum iterations
    call get_command_argument(1, arg1)
    call get_command_argument(2, arg2)
    call get_command_argument(3, arg3)
    
    read(arg1, *) n
    read(arg2, *) tol
    read(arg3, *) max_iter
    
    ! Allocate memory for the matrix and vectors
    allocate(A(n,n), b(n), x(n), r(n), r0(n), p(n), v(n), s(n), t(n))
    
    ! Ask user if they want to input their own matrix and vector
    print *, 'Do you want to input your own matrix A and vector b? (y/n)'
    read(*, '(A)') input_choice
    custom_input = trim(adjustl(input_choice)) == 'y'
    
    if (custom_input) then
        call ReadMatrixAndVector(A, b, n)
    else
        ! Initialize the matrix A and vector b (example initialization)
        call InitializeMatrixAndVector(A, b, n)
        call PrintMatrix(A, n)  ! (optional)
        call PrintVector(b, n)  ! (optional)
    end if

    ! Solve the linear system Ax = b
    call BiCGSTAB(A, b, x, n, tol, max_iter)    ! Alternative ver: BiCGSTAB_v2()
                                                ! A standalone and self-contained version
    ! Deallocate memory
    deallocate(A, b, x, r, r0, p, v, s, t)

    CONTAINS
    
    ! The reason why 'matvec' and 'BICGSTAB' are here is because
    ! 'matvec' requires explicit implementation and 'BICGSTAB' is
    ! the only other subroutine which uses it. 
    ! There is also an option to create a module, e.g. 
    ! module bicgstab_module
    !    contains
    !        ! Function to perform matrix-vector multiplication
    !        function matvec(matA, vecX, NDIM) result(vecAx)
    !            implicit none
    !            integer, intent(in) :: NDIM
    !        ...
    !        end function
    ! end module

    ! BLAS boosted version to perform matrix-vector multiplication
    function matvec(matA, vecX, NDIM) result(vecAx)
        !use, intrinsic :: iso_c_binding, only: c_double
        implicit none
        integer, intent(in) :: NDIM
        real(8), intent(in) :: matA(NDIM,NDIM), vecX(NDIM)
        real(8) :: vecAx(NDIM)
        integer :: one
        real(8) :: alpha, beta
        
        ! Parameters for DGEMV
        one = 1
        alpha = 1.0d0
        beta = 0.0d0
    
        ! Call to DGEMV (y = alpha*A*x + beta*y)
        call dgemv('N', NDIM, NDIM, alpha, matA, NDIM, vecX, one, beta, vecAx, one)
    end function matvec    

    ! Subroutine which uses BIConjugate Gradient STABilized iteration to solve Ax = b
    subroutine BiCGSTAB(A_inout, b_inout, x_inout, m, tolerence, maxIter)
        implicit none
        integer, intent(in) :: m, maxIter
        real(8), intent(in) :: tolerence
        real(8), intent(inout) :: A_inout(m,m), b_inout(m), x_inout(m)
        real(8) :: r_og(m), r0_og(m), p_og(m), v_og(m), s_og(m), t_og(m)
        real(8) :: alpha, beta, omega, rho, rho_old, res_norm
        integer :: l
        logical :: converged
        
        ! Initial guess for the solution
        x_inout = 0.0
        
        ! Start the BiCGSTAB algorithm
        r_og = b_inout - matvec(A_inout, x_inout, m)
        r0_og = r_og
        rho_old = 1.0
        alpha = 1.0
        omega = 1.0
        p_og = 0.0
        converged = .false.
        
        do l = 1, maxIter
            rho = dot_product(r0_og, r_og)
            
            if (rho == 0.0) then
                print *, 'Method failed: rho = 0.'
                exit
            end if
            
            if (l == 1) then
                p_og = r_og
            else
                beta = (rho / rho_old) * (alpha / omega)
                p_og = r_og + beta * (p_og - omega * v_og)
            end if
            
            v_og = matvec(A_inout, p_og, m)
            alpha = rho / dot_product(r0_og, v_og)
            
            s_og = r_og - alpha * v_og
            t_og = matvec(A_inout, s_og, m)
            omega = dot_product(t_og, s_og) / dot_product(t_og, t_og)
            
            x_inout = x_inout + alpha * p_og + omega * s_og
            r_og = s_og - omega * t_og
            
            res_norm = norm2(r_og)
            print *, 'Iteration:', l, ' Residual norm:', res_norm
            
            if (res_norm < tolerence) then
                converged = .true.
                exit
            end if
            
            rho_old = rho
        end do
        
        if (converged) then
            print *, 'Converged after', l, 'iterations.'
            print *, 'Result:'
            do l = 1, m
                write(*, '(A, I3, A, F12.8)') 'x(', l, ') = ', x_inout(l)
            end do
        else
            print *, 'Failed to converge within', maxIter, 'iterations.'
        end if
    end subroutine BiCGSTAB    
 
end program BiCGSTABSolver
 
! Subroutine to initialize the matrix A and vector b
subroutine InitializeMatrixAndVector(A_out, b_out, dim)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(out) :: A_out(dim,dim), b_out(dim)
    integer :: i, j
    
    ! Example: A is a diagonally dominant matrix
    A_out = 0.0
    do i = 1, dim
        do j = 1, dim
            if (i == j) then
                A_out(i,j) = 2.0
            else
                A_out(i,j) = -1.0 / (abs(i-j) + 1)
            end if
        end do
        b_out(i) = 1.0
    end do
end subroutine InitializeMatrixAndVector

! Subroutine to read the matrix A and vector b from user input
subroutine ReadMatrixAndVector(A_out, b_out, dim)
    implicit none
    integer, intent(in)  :: dim
    real(8), intent(out) :: A_out(dim,dim), b_out(dim)
    integer :: i, j
    
    ! Read matrix A row by row
    print *, 'Enter the values for matrix A (', dim, 'rows expected):'
    do i = 1, dim
        read(*,*) (A_out(i,j), j = 1, dim)  ! Read each row of the matrix
    end do
    
    ! Read vector b in one line
    print *, 'Enter the values for vector b (', dim, 'values expected):'
    read(*,*) (b_out(i), i = 1, dim)
end subroutine ReadMatrixAndVector

! Subroutine to print matrix A
subroutine PrintMatrix(A, n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: A(n,n)
    integer :: i, j
    integer :: print_size

    ! Determine the size to print
    if (n > 5) then
        print_size = 5
    else
        print_size = n
    end if

    print *, 'Matrix A:'
    do i = 1, print_size
        do j = 1, print_size
            write(*, '(F8.4)', advance='no') A(i,j)
            if (j < print_size) then
                write(*, '(A)', advance='no') '  '
            end if
        end do
        write(*, *) ''  ! Move to the next line after each row
    end do

    ! If matrix is larger than 5x5, indicate that the matrix is truncated
    if (n > 5) then
        print *, '... (Matrix is larger than 5x5, showing only the top-left 5x5 corner) '
        print *, ''
    end if
end subroutine PrintMatrix

! Subroutine to print vector v (b)
subroutine PrintVector(v, n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: v(n)
    integer :: i
    integer :: print_size

    ! Determine the number of elements to print
    if (n > 5) then
        print_size = 5
    else
        print_size = n
    end if

    print *, 'Vector b:'
    do i = 1, print_size
        write(*, '(A, I2, A, F12.8)') 'b(', i, ') = ', v(i)
    end do

    ! If vector has more than 5 elements, indicate that the vector is truncated
    if (n > 5) then
        print *, '... (Vector is larger than 5 elements, showing only the first 5 elements)'
        print *, ''
    end if
end subroutine PrintVector
