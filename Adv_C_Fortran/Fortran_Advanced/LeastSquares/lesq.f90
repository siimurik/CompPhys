!   gfortran lesq.f90 -o lesq
program main
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), dimension(9) :: x, y, w, wy
    real(dp), dimension(:, :), allocatable :: V, Q
    real(dp), dimension(:, :), allocatable :: R
    real(dp), dimension(:), allocatable :: coefs, qty
    integer :: rows, cols, n, i
    real(dp), dimension(:,:), allocatable :: Q_house
    real(dp), dimension(:,:), allocatable :: R_house
    real(dp), dimension(:)  , allocatable :: coefs_house, qty_house

    ! Initialize data
    x = [1.2D0, 1.3D0, 1.4D0, 1.5D0, 1.6D0, 1.7D0, 1.8D0, 1.9D0, 2.0D0]
    y = [9.08D0, 10.43D0, 11.9D0, 13.48D0, 15.19D0, 17.03D0, 19.01D0, 21.13D0, 23.39D0]
    w = [1.0D0, 1.0D0, 2.0D0, 5.0D0, 1.0D0, 4.0D0, 2.0D0, 2.0D0, 1.0D0]

    rows = size(x)
    n = 3
    cols = n + 1

    allocate(V(rows,cols), Q(rows,cols), R(cols,cols))
    allocate(coefs(cols), qty(cols))
    allocate(Q_house(rows,rows), R_house(cols,cols))
    allocate(coefs_house(cols),qty_house(cols))

    call vander(x, y, w, n, rows, V, wy)
    call print_matrix(V, "Weighted Vandermonde matrix V:")

    write(*, '(a)') ""
    write (*,'(A)') "Gram-Schmidt method:"
    call qr_gramschmidt(rows, cols, V, Q, R)
    call print_matrix(Q, "Orthogonal matrix Q:")
    call print_matrix(R, "Upper triangular matrix R:")

    ! Step 1: Compute Qᵀ·wy
    qty = matmul(transpose(Q), wy)
    call print_vector(qty, "Qᵀ·wy")


    ! V·x = wy --> Q·R·x = wy --> Qᵀ | Q·R·x = wy
    ! Step 2: Solve R·X = Qᵀ·wy using back substitution
    call backsub(R, qty, coefs, cols) 

    ! OPTIONAL
    ! Step 3: Flip coefs to match polynomial order [c0, c1, c2, c3]
    !coefs = coefs(cols:1:-1)
    !print *, "Polynomial coefs (highest power first):"
    
    print *, "Polynomial coefs (lowest power first):"
    print "(4(f14.6,1x))", coefs

    write(*, '(a)') ""
    write (*,'(A)') "Householder method:"
    call qr_householder(rows, cols, V, Q_house, R_house)
    call print_matrix(Q_house, "Orthogonal matrix Q:")
    call print_matrix(R_house, "Upper triangular matrix R:")
    qty_house = matmul(transpose(Q_house(:,1:cols)), wy)  ! Economy Qᵀy
    call print_vector(qty_house, "Qᵀ·wy")
    call backsub(R_house, qty_house, coefs_house, cols)
    print *, "Polynomial coefs (lowest power first):"
    print "(4(f14.6,1x))", coefs_house

    deallocate(V, Q, R, coefs, qty)
    deallocate(Q_house, R_house, coefs_house, qty_house)

    contains

        subroutine print_matrix(matrix, title)
            implicit none
            integer :: j
            real(dp), dimension(:,:), intent(in) :: matrix
            character(len=*), intent(in) :: title   
            character(len=15) :: mformat = '(100(F14.6,1x))'
        
            print *, title
            do i = 1, size(matrix, 1)
                write(*, mformat) (matrix(i, j), j = 1, size(matrix, 2))
            end do
        end subroutine print_matrix

        ! Helper subroutine to print vectors
        subroutine print_vector(vec, description)
            implicit none
            character(len=*), intent(in) :: description
            double precision, intent(in) :: vec(:)
            character(len=15) :: mformat = '(100(F14.6,1x))'
            
            write(*,*)
            write(*,*) description
            write(*, mformat) vec
        end subroutine print_vector
end program main


subroutine vander(x, y, w, n, m, V, wy)
    implicit none
    integer, parameter :: dp = kind(1.0d0)  ! Double precision kind
    integer, intent(in) :: n, m
    real(dp), dimension(m), intent(in) :: x, y, w
    real(dp), dimension(m, n+1), intent(out) :: V
    real(dp), dimension(m), intent(out) :: wy
    
    integer :: i, j
    
    ! Construct the weighted Vandermonde matrix
    do i = 1, m
        ! Set the last column as the weight
        V(i, n+1) = w(i)  ! Last column is the weight (Fortran is 1-based)
        
        ! Fill the Vandermonde row starting from the highest power to the lowest
        do j = n, 1, -1
            V(i, j) = x(i) * V(i, j+1)
        end do
        
        ! Alternative implementation commented out:
        ! Fill the Vandermonde row starting from the lowest power to the highest
        !do j = 0, n
        !    if (j == 0) then
        !        V(i, j+1) = 1.0d0  ! First element is 1 (adjust for 1-based index)
        !    else
        !        V(i, j+1) = x(i) * V(i, j)
        !    end if
        !end do
        
        ! Calculate weighted y values
        wy(i) = w(i) * y(i)
    end do
end subroutine vander

subroutine qr_gramschmidt(m, n, a, q, r)
    !**********************************************************************************
    !This subroutine performs QR decomposition using Gram-Schmidt orthogonalization
    ![input]   m,n is the dimension of [m * n] input 'a' matrix
    !          a is the input matrix
    ![output]  q is the orthogonal matrix
    !          r is the upper triangular matrix
    !**********************************************************************************
    implicit none
    integer, intent(in) :: m, n
    integer, parameter :: dp = kind(1.0d0)  ! Double precision kind
    real(dp), dimension(m, n), intent(in) :: a
    real(dp), dimension(m, n), intent(out) :: q
    real(dp), dimension(n, n), intent(out) :: r
    real(dp), dimension(m, n) :: u
    integer :: k, i

    q = 0.0_dp
    r = 0.0_dp
    u = 0.0_dp

    ! Gram-Schmidt process
    do k = 1, n
        u(:, k) = a(:, k)
        do i = 1, k-1
            r(i, k) = dot_product(q(:, i), a(:, k))
            u(:, k) = u(:, k) - r(i, k) * q(:, i)
        end do
        r(k, k) = norm2(u(:, k))
        if (r(k, k) > epsilon(1.0_dp)) then
            q(:, k) = u(:, k) / r(k, k)
        else
            q(:, k) = u(:, k)
        end if
    end do
end subroutine qr_gramschmidt

subroutine qr_householder(m, n, a, q, r)
    !**********************************************************************************
    ! This subroutine performs QR decomposition using Householder reflections
    ! [input]   m,n is the dimension of [m × n] input 'a' matrix (m ≥ n)
    !           a is the input matrix
    ! [output]  q is the orthogonal matrix (m × m)
    !           r is the upper triangular matrix (n × n)
    !**********************************************************************************
    implicit none
    integer, intent(in) :: m, n
    integer, parameter :: dp = kind(1.0d0)
    real(dp), dimension(m,n), intent(in) :: a
    real(dp), dimension(m,m), intent(out) :: q
    real(dp), dimension(n,n), intent(out) :: r
    real(dp), dimension(m,n) :: a_copy
    real(dp), dimension(m) :: x, v
    real(dp) :: alpha, norm_x
    integer :: k, j

    a_copy = a
    q = 0.0_dp
    ! Initialize Q as identity matrix
    do k = 1, m
        q(k,k) = 1.0_dp
    end do

    ! Householder triangularization
    do k = 1, n
        x = 0.0_dp
        x(k:m) = a_copy(k:m,k)
        norm_x = norm2(x(k:m))
        
        ! Compute Householder vector
        alpha = -sign(1.0_dp, x(k)) * norm_x
        v = x
        v(k) = v(k) - alpha
        v = v / norm2(v)
        
        ! Apply reflection to remaining columns
        do j = k, n
            a_copy(:,j) = a_copy(:,j) - 2.0_dp * v * dot_product(v, a_copy(:,j))
        end do
        
        ! Accumulate reflections into Q
        do j = 1, m
            !q(:,j) = q(:,j) - 2.0_dp * v * dot_product(v, q(:,j))
            q(j,:) = q(j,:) - 2.0_dp * v * dot_product(v, q(j,:))
        end do
    end do

    ! Extract R from upper triangle of transformed a_copy
    r = 0.0_dp
    do k = 1, n
        r(1:k,k) = a_copy(1:k,k)
    end do
end subroutine qr_householder

subroutine backsub(r, b, x, n)
    !**********************************************************************************
    ! Solves Rx = b where R is upper triangular
    ! [input]   r - n×n upper triangular matrix
    !           b - right-hand side vector (length n)
    !           n - system size
    ! [output]  x - solution vector
    !**********************************************************************************
    implicit none
    integer, intent(in) :: n
    integer, parameter :: dp = kind(1.0d0)
    real(dp), dimension(n,n), intent(in) :: r
    real(dp), dimension(n), intent(in) :: b
    real(dp), dimension(n), intent(out) :: x
    integer :: i, j

    x = 0.0_dp
    do i = n, 1, -1
        x(i) = b(i)
        do j = i+1, n
            x(i) = x(i) - r(i,j) * x(j)
        end do
        if (abs(r(i,i)) > epsilon(1.0_dp)) then
            x(i) = x(i) / r(i,i)
        else
            x(i) = 0.0_dp  ! Handle zero diagonal (rank-deficient case)
        end if
    end do
end subroutine backsub