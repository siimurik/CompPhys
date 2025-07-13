program main
    use lu_solvers
    implicit none
    integer :: m, n
    double precision, dimension(:,:), allocatable :: A, L, U, expmat
    double precision, dimension(:)  , allocatable :: y, ux
    double precision, dimension(:)  , allocatable :: x, f, b
    double precision :: k0, k1, kn, h
    !double precision :: val_main, val_sub, val_super
    integer :: i, j
    logical :: succ

    !------------------ intro text-----------------------------------
    write (*,'(A)') "======================================================="
    write (*,'(A)') "                  ██████╗ ███████╗ ██╗"
    write (*,'(A)') "                 ██╔════╝ ██╔════╝ ██║"
    write (*,'(A)') "                 ██║      ███████╗ ██║"
    write (*,'(A)') "                 ██║      ╚════██║ ██║"
    write (*,'(A)') "                 ╚██████╗ ███████║ ██║"
    write (*,'(A)') "                  ╚═════╝ ╚══════╝ ╚═╝"
    write (*,'(A)') "-------------------------------------------------------"
    write (*,'(A)') "  C U B I C   S P L I N E   I N T E R P O L A T I O N  "
    write (*,'(A)') "-------------------------------------------------------"
    write (*,'(A)') "Computes smooth curves through (x,y) points  "
    write (*,'(A)') "using LU decomposition (Doolittle/Crout/Cholesky)  "
    write (*,'(A)') "======================================================="
    write (*,*) ""
    !------------------ intro text-----------------------------------

    ! Size of initial dataset arrays x and f
    m = 5   ! CHANGE THIS ACCORDING TO DATASET SIZE
    n = m-1 ! 

    allocate(x(m), f(m), b(m))
    allocate(A(n-1,n-1), L(n-1,n-1), U(n-1,n-1), y(n-1), ux(n-1))
    allocate(expmat(m-1, 4))
    
    x = [-2.0D0, -1.0D0,  0.0D0,  1.0D0, 2.0D0]
    f = [ 0.0D0,  0.0D0,  1.0D0,  0.0D0, 0.0D0]
    k0 = 0.0D0
    kn = 0.0D0

    !f = [1.0, 9.0, 41.0, 41.0]
    !x = [0.0, 2.0,  4.0,  6.0]
    !k0 =   0.0
    !kn = -12.0

    if (size(x,1) .ne. m) then
        error stop "Mismatch between input dimension and dataset size."
    end if

    ! Step size based on x
    h = abs(x(2) - x(1))
    
    if (n .gt. 2) then
        A = 0.0D0
        b = 0.0D0
        y = 0.0D0

        call fill_diagonals(A, val_main=4.0d0)

        ! Fill sub- and superdiagonal
        if (n .gt. 3) then
            call fill_diagonals(A, val_sub=1.0d0, val_super=1.0d0)
        end if

        call print_matrix(A, "Matrix A")

        ! Right-hand side of Ax = y
        call init_rightside(h, f, y)
        !call print_vector(y, "Vector y")

        ! Adjust first and last equations for boundary conditions
        y(1) = y(1) - k0
        y(n-1) = y(n-1) - kn

        call print_vector(y, "Vector y")
        write (*,*) ""

        ! Solve the system using LU-decomp
        L = 0.0D0
        U = 0.0D0
                                                            !  | - how many el in ux
        ! Solve with all methods                            !  |
        call solve_lu(A, y, ux, "doolittle", succ)          ! \/ 
        if (succ) print '(a, 3e14.4)', "Doolittle x =", ux  ! 3e14.4               
        !if (succ) write(*,'(a, 3e12.4)') "Doolittle x = ", (ux(i), i=1,n)
        
        call solve_lu(A, y, ux, "crout", succ)
        if (succ) print '(a, 3e14.4)', "Crout     x =", ux

        call solve_lu(A, y, ux, "cholesky", succ)
        if (succ) print '(a, 3e14.4)', "Cholesky  x =", ux
        

        ! Adjust first and last equations for boundary conditions
        b(1) = k0
        do i = 1, n-1
            b(i+1) = ux(i)
        end do
        b(m) = kn

    else
        ! Special case for only 3 points (1 interior point)
        k1 = (3.0D0/h) * (f(3) - f(1)) - k0 - kn
        b = [k0, k1, kn]
    end if

    call print_vector(b, "Vector b")
    write (*,*) ""

    call calc_spline_coefs(x, f, b, h, expmat)
    write (*,*) ""
    call print_matrix(expmat, "Expanded coefs (highest to lowest):")

    deallocate(x, f, b)
    deallocate(A, L, U, y, ux)
    deallocate(expmat)

    !------------- end of main code -------------!

    contains

        subroutine fill_diagonals(mat, val_main, val_sub, val_super)
            implicit none
            double precision, dimension(:,:), intent(inout) :: mat
            double precision, intent(in), optional :: val_main, val_sub, val_super
            integer :: ndim
            
            ndim = size(mat, 1)
            
            ! Check if matrix is square and at least 2x2 for sub/super diagonals
            if (size(mat, 1) .ne. size(mat, 2)) then
                print *, "Matrix must be square."
                stop
            end if
            
            if (ndim < 2 .and. (present(val_sub) .or. present(val_super))) then
                print *, "Matrix must be at least 2x2 for sub/super diagonals."
                stop
            end if
            
            ! Fill main diagonal if value is provided
            if (present(val_main)) then
                do i = 1, ndim
                    mat(i,i) = val_main
                end do
            end if
            
            ! Fill subdiagonal if value is provided
            if (present(val_sub)) then
                do i = 1, ndim-1
                    mat(i+1,i) = val_sub
                end do
            end if
            
            ! Fill superdiagonal if value is provided
            if (present(val_super)) then
                do i = 1, ndim-1
                    mat(i,i+1) = val_super
                end do
            end if
        end subroutine fill_diagonals

        subroutine init_rightside(step_size, f_data, y_vec)
            implicit none
            double precision, intent(in) :: step_size
            double precision, intent(in) :: f_data(:)
            double precision, intent(inout) :: y_vec(:)
            integer :: mdim

            mdim = size(f_data,1)

            ! Right-hand side
            do j = 2, mdim-1
                y_vec(j-1) = (3.0D0/step_size) * (f_data(j+1) - f_data(j-1))
            end do

        end subroutine init_rightside

        subroutine print_matrix(matrix, title)
            implicit none
            double precision, dimension(:,:)       :: matrix
            character(len=*), intent(in), optional :: title   
            ! intent(in) - information for the programmer and the 
            ! compiler that this is an INPUT; optional - not ne-
            ! cessary for the user to input it to work correctly.
            integer :: irows, jcols
            character(len=15) :: mformat = '(100(F14.6,1x))'
        
            write(*, *) title
            do irows = 1, min(5, size(matrix, 1))
                write(*, mformat) (matrix(irows, jcols), jcols = 1, min(5, size(matrix, 2)))
            end do
        end subroutine print_matrix


        ! Helper subroutine to print vectors
        subroutine print_vector(vec, description)
            implicit none
            character(len=*), intent(in), optional :: description
            double precision, intent(in) :: vec(:)
            character(len=15) :: mformat = '(100(F14.6,1x))'
            
            write(*,*)
            write(*,*) description
            write(*, mformat) vec
        end subroutine print_vector

        function expand_poly(a0, a1, a2, a3, xj) result(expanded)
            double precision, intent(in) :: a0, a1, a2, a3, xj
            double precision :: expanded(4)
            double precision :: xj_sq, xj_cb
            
            xj_sq = xj * xj
            xj_cb = xj * xj_sq
            
            expanded(1) = a3
            expanded(2) = a2 - 3.0d0*a3*xj
            expanded(3) = a1 - 2.0d0*a2*xj + 3.0d0*a3*xj_sq
            expanded(4) = a0 - a1*xj + a2*xj_sq - a3*xj_cb
        end function expand_poly

        subroutine calc_spline_coefs(x_data, f_data, b_vec, step_size, expanded_mat)
            double precision, intent(in) :: x_data(:), f_data(:), b_vec(I), step_size
            double precision, dimension(:,:), allocatable, intent(inout) :: expanded_mat
            double precision :: a0, a1, a2, a3
            integer :: ndim
            
            ndim = size(x_data)-1
            
            ! Print header
            write (*, '(a)') "---------------------------------------------"
            write (*, '(a)') "|   j   |   aj0  |   aj1  |   aj2  |   aj3  |"
            write (*, '(a)') "---------------------------------------------"
            
            do j = 1, ndim
                ! Calculate coefficients
                a0 = f_data(j)
                a1 = b_vec(j)
                a2 = (3.0d0/(step_size**2)) * (f_data(j+1) - f_data(j  )) - (1.0d0/step_size) * (b_vec(j+1) + 2.0d0*b_vec(j))
                a3 = (2.0d0/(step_size**3)) * (f_data(j  ) - f_data(j+1)) + (1.0d0/(step_size**2)) * (b_vec(j+1) +  b_vec(j))
                
                ! Print current row
                print '(a,i3,a,f5.2,a,f5.2,a,f5.2,a,f5.2,a)', &
                    "| ", j, "   | ", a0, "  | ", a1, "  | ", a2, "  | ", a3, "  |"
                
                expanded_mat(j,:) = expand_poly(a0, a1, a2, a3, x(j))
            end do
            
            write (*, '(a)') "---------------------------------------------"

        end subroutine calc_spline_coefs

end program main

! Unused
subroutine thomas_algorithm(a, b, c, d, x, n)
    implicit none
    !------------------------------------------------------------------
    ! Solves a tridiagonal system Ax = d using the Thomas Algorithm.
    !
    ! Args:
    !   a(n) : Sub-diagonal (a(1) is unused, a(2:n) are the entries).
    !   b(n) : Main diagonal.
    !   c(n) : Super-diagonal (c(n) is unused, c(1:n-1) are the entries).
    !   d(n) : Right-hand side vector.
    !   n    : Size of the system.
    ! Returns:
    !   x(n) : Solution vector.
    !------------------------------------------------------------------
    integer, intent(in) :: n
    double precision, intent(in) :: a(n), b(n), c(n), d(n)
    double precision, intent(out) :: x(n)
    double precision :: cp(n), dp(n)  ! Modified coefficients
    integer :: i


    !-------------------------------------------
    ! Forward Elimination (Compute c' and d')
    !-------------------------------------------
    ! First row (i=1)
    cp(1) = c(1) / b(1)
    dp(1) = d(1) / b(1)

    ! Rows i=2 to i=n-1
    do i = 2, n-1
        cp(i) = c(i) / (b(i) - a(i) * cp(i-1))
        dp(i) = (d(i) - a(i) * dp(i-1)) / (b(i) - a(i) * cp(i-1))
    end do

    ! Last row (i=n)
    dp(n) = (d(n) - a(n) * dp(n-1)) / (b(n) - a(n) * cp(n-1))

    !-------------------------------------------
    ! Back Substitution (Compute x)
    !-------------------------------------------
    x(n) = dp(n)
    do i = n-1, 1, -1
        x(i) = dp(i) - cp(i) * x(i+1)
    end do
end subroutine thomas_algorithm