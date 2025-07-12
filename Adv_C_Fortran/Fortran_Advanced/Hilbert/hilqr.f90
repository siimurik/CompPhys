!   gfortran hilqr.f90 -o hilqr -llapack
program qr_failure_demo
    implicit none
    integer, parameter :: dp = kind(1.0d0), dim = 5
    real(dp), dimension(dim,dim) :: Q_gs, Q_hh, R_gs, R_hh, eye5, H
    real(dp) :: cond_num
    integer :: k
    
    call hilb(dim, H)
    call print_matrix(H, "Hilbert Matrix H:")
    
    cond_num = norm2(H) * norm2(inv(H))
    print *, "Condition number:", cond_num
    
    call qr_gramschmidt(dim, dim, H, Q_gs, R_gs)
    call qr_householder(dim, dim, H, Q_hh, R_hh)
    
    eye5 = eye(dim)
    print *, "Gram-Schmidt ‖QᵀQ - I‖:", norm2(matmul(transpose(Q_gs), Q_gs) - eye5)
    print *, "Householder ‖QᵀQ - I‖:", norm2(matmul(transpose(Q_hh), Q_hh) - eye5)
    
    ! Reconstruction error
    print *, "Gram-Schmidt ‖QR - H‖:", norm2(matmul(Q_gs, R_gs) - H)
    print *, "Householder ‖QR - H‖:", norm2(matmul(Q_hh(:,1:dim), R_hh) - H)
    
    ! Check diagonal of R for stability
    print *, "Gram-Schmidt R diag:", [(R_gs(k,k), k=1,dim)]
    print *, "Householder R diag:", [(R_hh (k,k), k=1,dim)]

    ! Check first column orthogonality
    print *, "GS q1·q2:", dot_product(Q_gs(:,1), Q_gs(:,2))
    print *, "HH q1·q2:", dot_product(Q_hh(:,1), Q_hh(:,2))

    contains

        subroutine print_matrix(matrix, title)
            implicit none
            integer :: i, j
            real(dp), dimension(:,:), intent(in) :: matrix
            character(len=*), intent(in) :: title   
            character(len=15) :: mformat = '(100(F14.6,1x))'
        
            print *, title
            do i = 1, size(matrix, 1)
                write(*, mformat) (matrix(i, j), j = 1, size(matrix, 2))
            end do
        end subroutine print_matrix

        ! Returns the inverse of a matrix calculated by finding the LU
        ! decomposition.  Depends on LAPACK.
        function inv(A) result(Ainv)
            !integer, parameter :: dp = kind(1.0d0)  ! Double precision kind
            real(dp), dimension(:,:), intent(in) :: A
            real(dp), dimension(size(A,1),size(A,2)) :: Ainv

            real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
            integer, dimension(size(A,1)) :: ipiv   ! pivot indices
            integer :: n, info

            ! External procedures defined in LAPACK
            external DGETRF
            external DGETRI

            ! Store A in Ainv to prevent it from being overwritten by LAPACK
            Ainv = A
            n = size(A,1)

            ! DGETRF computes an LU factorization of a general M-by-N matrix A
            ! using partial pivoting with row interchanges.
            call DGETRF(n, n, Ainv, n, ipiv, info)

            if (info /= 0) then
                stop 'Matrix is numerically singular!'
            end if

            ! DGETRI computes the inverse of a matrix using the LU factorization
            ! computed by DGETRF.
            call DGETRI(n, Ainv, n, ipiv, work, n, info)

            if (info /= 0) then
                stop 'Matrix inversion failed!'
            end if
        end function inv

        function eye(n) result(mat_i)
            implicit none
            integer, intent(in) :: n
            real(dp), dimension(n,n) :: mat_i
            integer :: i

            mat_i = 0.0_dp

            do i = 1, n
                mat_i(i,i) = 1.0_dp
            end do
            
        end function eye

end program

subroutine hilb(n, H)
    !=========================================================== 
    ! Evaluate the values of the Hilbert matrix following the
    ! formula:
    ! H(i,j) = 1 / (i + j - 1)
    ! 
    ! Input:
    ! ========
    ! n - dimension of the matrix
    !
    ! Output:
    ! ========
    ! H - Hilbert matrix with dimension n x n
    !
    ! Example for n = 3:
    ! =======
    !
    ! Hilbert matrix with dimension  3
    !      1.0000      0.5000      0.3333
    !      0.5000      0.3333      0.2500
    !      0.3333      0.2500      0.2000
    !
    !=========================================================== 
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), dimension(n,n), intent(out) :: H
    integer, intent(in) :: n
    integer :: i, j

    do i = 1, n
        do j = 1, n
            H(i,j) = 1.0_dp / real(i + j - 1)
        end do
    end do
    return

end subroutine hilb

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