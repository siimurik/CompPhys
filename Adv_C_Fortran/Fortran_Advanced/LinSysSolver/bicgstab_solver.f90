! Compile and execute with:
!  $ gfortran -o bicgstab_solver bicgstab_solver.f90
!  $ ./bicgstab_solver 10 1.0e-6 1000
program BiCGSTABSolver
    implicit none
    integer :: n, max_iter, k
    real(8), allocatable :: A(:,:), b(:), x(:), r(:), r0(:), p(:), v(:), s(:), t(:)
    real(8) :: tol, alpha, beta, omega, rho, rho_old, res_norm
    logical :: converged
    character(len=100) :: arg1, arg2, arg3, input_choice
    logical :: custom_input

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
        ! Custom input for matrix A
        print *, 'Enter the values for matrix A, row by row (', n*n, 'values expected):'
        call ReadMatrix(A, n)
  
        ! Custom input for vector b
        print *, 'Enter the values for vector b (', n, 'values expected):'
        call ReadVector(b, n)
    else
        ! Initialize the matrix A and vector b (example initialization)
        call InitializeMatrixAndVector(A, b, n)
    end if
 
    ! Initial guess for the solution
    x = 0.0
 
    ! Start the BiCGSTAB algorithm
    r = b - matvec(A, x, n)
    r0 = r
    rho_old = 1.0
    alpha = 1.0
    omega = 1.0
    p = 0.0
    converged = .false.
 
    do k = 1, max_iter
       rho = dot_product(r0, r)
 
       if (rho == 0.0) then
          print *, 'Method failed: rho = 0.'
          exit
       end if
 
       if (k == 1) then
          p = r
       else
          beta = (rho / rho_old) * (alpha / omega)
          p = r + beta * (p - omega * v)
       end if
 
       v = matvec(A, p, n)
       alpha = rho / dot_product(r0, v)
 
       s = r - alpha * v
       t = matvec(A, s, n)
       omega = dot_product(t, s) / dot_product(t, t)
 
       x = x + alpha * p + omega * s
       r = s - omega * t
 
       res_norm = norm2(r)
       print *, 'Iteration:', k, ' Residual norm:', res_norm
 
       if (res_norm < tol) then
          converged = .true.
          exit
       end if
 
       rho_old = rho
    end do
 
    if (converged) then
        print *, 'Converged after', k, 'iterations.'
        print *, 'Result:'
        do k = 1, n
        write(*, '(A, I3, A, F12.8)') 'x(', k, ') = ', x(k)
        end do
    else
        print *, 'Failed to converge within', max_iter, 'iterations.'
    end if
 
    ! Deallocate memory
    deallocate(A, b, x, r, r0, p, v, s, t)
 
CONTAINS
 
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
 
! Subroutine to read matrix A from user input
subroutine ReadMatrix(A_out, dim)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(out) :: A_out(dim,dim)
    integer :: i, j
 
    do i = 1, dim
       print *, 'Row ', i, ':'
       do j = 1, dim
          read(*, *) A_out(i,j)
       end do
    end do
end subroutine ReadMatrix
 
! Subroutine to read vector b from user input
subroutine ReadVector(b_out, dim)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(out) :: b_out(dim)
    integer :: i
 
    print *, 'Enter ', dim, ' values for vector b:'
    do i = 1, dim
       read(*, *) b_out(i)
    end do
end subroutine ReadVector
 
 ! Function to perform matrix-vector multiplication
function matvec(matA, x_vec, NDIM) result(vecAx)
    implicit none
    integer, intent(in) :: NDIM
    real(8), intent(in) :: matA(NDIM,NDIM), x_vec(NDIM)
    real(8) :: vecAx(NDIM)
    integer :: i, j
 
    vecAx = 0.0
    do i = 1, NDIM
       do j = 1, NDIM
          vecAx(i) = vecAx(i) + matA(i,j) * x_vec(j)
       end do
    end do
end function matvec
 
 ! Function to calculate the 2-norm of a vector
 !function norm2(v) result(norm)
 !   implicit none
 !   real(8), intent(in) :: v(:)
 !   real(8) :: norm
 !
 !   norm = sqrt(sum(v**2))
 !end function norm2
 
end program BiCGSTABSolver
 