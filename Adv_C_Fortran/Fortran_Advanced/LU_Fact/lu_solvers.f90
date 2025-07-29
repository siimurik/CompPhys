module lu_solvers
    !------------------ LU solvers --------------------
    contains

        subroutine lu_doolittle(A, L, U, success)
            !=============================================================
            ! Key features:
            !   * L has 1's on the diagonal.
            !   * U has non-zero diagonal.
            !   * No pivoting (use only for diagonally dominant matrices).
            !=============================================================
            implicit none
            double precision, intent(in) ::  A(:, :)
            double precision, intent(out) :: L(:, :), U(:, :)
            logical, intent(out) :: success
            integer :: i, j, k, n
            double precision :: sum_lu

            n = size(A, 1)

            success = .true.
            L = 0.0d0
            U = 0.0d0

            do k = 1, n
                ! Compute U(k, k:n)
                do j = k, n
                    sum_lu = 0.0d0
                    do i = 1, k-1
                        sum_lu = sum_lu + L(k, i) * U(i, j)
                    end do
                    U(k, j) = A(k, j) - sum_lu
                end do

                ! Check for zero pivot
                if (abs(U(k, k)) < tiny(1.0d0)) then
                    success = .false.
                    return
                end if

                ! Compute L(k+1:n, k)
                do i = k+1, n
                    sum_lu = 0.0d0
                    do j = 1, k-1
                        sum_lu = sum_lu + L(i, j) * U(j, k)
                    end do
                    L(i, k) = (A(i, k) - sum_lu) / U(k, k)
                end do
                L(k, k) = 1.0d0  ! Doolittle's L has 1's on diagonal
            end do
        end subroutine lu_doolittle


        subroutine lu_crout(A, L, U, success)
            !=============================================================
            ! Key features:
            !   * L has non-zero diagonal.
            !   * U has 1's on the diagonal.
            !   * More stable for tridiagonal/banded matrices.
            !=============================================================
            implicit none
            double precision, intent(in) :: A(:, :)
            double precision, intent(out) :: L(:, :), U(:, :)
            logical, intent(out) :: success
            integer :: i, j, k, n
            double precision :: sum_lu

            n = size(A, 1)

            success = .true.
            L = 0.0d0
            U = 0.0d0

            ! Crout's Algorithm
            do j = 1, n
                ! Compute L(i,j) for i >= j
                do i = j, n
                    sum_lu = 0.0d0
                    do k = 1, j-1
                        sum_lu = sum_lu + L(i, k) * U(k, j)
                    end do
                    L(i, j) = A(i, j) - sum_lu
                end do

                ! Check for zero pivot
                if (abs(L(j, j)) < tiny(1.0d0)) then
                    success = .false.
                    return
                end if

                ! Compute U(j,i) for i > j
                do i = j+1, n
                    sum_lu = 0.0d0
                    do k = 1, j-1
                        sum_lu = sum_lu + L(j, k) * U(k, i)
                    end do
                    U(j, i) = (A(j, i) - sum_lu) / L(j, j)
                end do

                ! Set U(j,j) = 1
                U(j, j) = 1.0d0
            end do
        end subroutine lu_crout


        subroutine cholesky(A, L, success)
            ! Cholesky decomposition of a symmetric positive definite matrix A
            ! Computes the lower triangular matrix L such that A = L * L^T
            !
            ! Parameters:
            ! A(n,n)  - Input: symmetric positive definite matrix
            ! L(n,n)  - Output: lower triangular Cholesky factor
            ! n       - Input: dimension of the matrices
            ! success - Output: logical flag indicating success
            
            double precision, intent(in) ::  A(:,:)
            double precision, intent(out) :: L(:,:)
            logical, intent(out) :: success
            integer :: i, j, k, n
            double precision :: sum_val, temp
            
            n = size(A,1)

            ! Initialize output matrix L to zero
            L = 0.0_8
            success = .true.
            
            ! Perform Cholesky decomposition
            do i = 1, n
                ! Compute diagonal element L(i,i)
                sum_val = 0.0_8
                do k = 1, i-1
                    sum_val = sum_val + L(i,k) * L(i,k)
                end do
                
                temp = A(i,i) - sum_val
                
                ! Check for positive definiteness
                if (temp <= 0.0_8) then
                    success = .false.
                    return
                end if
                
                L(i,i) = sqrt(temp)
                
                ! Compute off-diagonal elements L(j,i) for j > i
                do j = i+1, n
                    sum_val = 0.0_8
                    do k = 1, i-1
                        sum_val = sum_val + L(j,k) * L(i,k)
                    end do
                    L(j,i) = (A(j,i) - sum_val) / L(i,i)
                end do
            end do
                
        end subroutine cholesky
 

        subroutine forward_substitution(L, b, y)
            !--------------------------------------------------------------
            ! Forward substitution (solve Ly = b where L is lower-triangular)
            !--------------------------------------------------------------
            implicit none
            double precision, intent(in)  :: L(:,:), b(:)
            double precision, intent(out) :: y(:)
            integer :: i, j, n

            n = size(L,1)
            y(1) = b(1) / L(1,1)
            do i = 2, n
                y(i) = b(i)
                do j = 1, i-1
                    y(i) = y(i) - L(i,j) * y(j)
                end do
                y(i) = y(i) / L(i,i)
            end do
        end subroutine forward_substitution

        !--------------------------------------------------------------
        ! Backward substitution (solve Ux = y where U is upper-triangular)
        !--------------------------------------------------------------
        subroutine backward_substitution(U, y, x)
            implicit none
            double precision, intent(in)  :: U(:,:), y(:)
            double precision, intent(out) :: x(:)
            integer :: i, j, n

            n = size(U,1)

            x(n) = y(n) / U(n,n)
            do i = n-1, 1, -1
                x(i) = y(i)
                do j = i+1, n
                    x(i) = x(i) - U(i,j) * x(j)
                end do
                x(i) = x(i) / U(i,i)
            end do
        end subroutine backward_substitution

        subroutine solve_lu(A, b, x, method, success)
            implicit none
            ! Input/Output arrays with automatic dimension detection
            double precision, intent(in)  :: A(:,:)   ! Matrix (square)
            double precision, intent(in)  :: b(:)     ! RHS vector
            double precision, intent(out) :: x(:)     ! Solution vector
            character(*), intent(in)      :: method   ! "doolittle", "crout", "cholesky"
            logical, intent(out)          :: success
            
            ! Local variables with automatic sizing
            integer :: n
            double precision, allocatable :: L(:,:), U(:,:), y(:)

            ! Get system size from input
            n = size(A, 1)

            ! Validate dimensions
            if (size(A,2) /= n)       error stop "Matrix A must be square"
            if (size(b) /= n)         error stop "Vector b has wrong size"
            if (size(x) /= n)         error stop "Vector x has wrong size"

            ! Allocate workspace
            allocate(L(n,n), U(n,n), y(n))

            ! Select decomposition method
            select case (method)
            case ("doolittle")
                call lu_doolittle(A, L, U, success)
            case ("crout")
                call lu_crout(A, L, U, success)
            case ("cholesky")
                call cholesky(A, L, success)
                U = transpose(L)  ! A = L·Lᵀ
            case default
                success = .false.
                return
            end select


            ! Solve Ly = b, then Ux = y
            call forward_substitution(L, b, y)
            call backward_substitution(U, y, x)

            deallocate(L, U, y)
        end subroutine solve_lu

end module lu_solvers