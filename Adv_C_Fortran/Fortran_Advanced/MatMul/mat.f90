!=============================================
! Compile and execute with:
!   $ gfortran -O3 mat.f90 -o mat -lblas -fopenmp
!   $ ./mat
!=============================================
program main
    implicit none
    external :: DGEMM  ! Explicit interface for BLAS
    double precision, parameter :: alpha = 1.0D0, beta = 0.0D0
    integer, parameter          :: dim = 5000
    integer                     :: m, k, n, j
    double precision, allocatable, dimension(:,:) :: a, b, c
    double precision :: elapsed_seconds
    integer          :: start_time, end_time, elapsed_time, rate

    ! Declare symmetric matrix
    m = dim; k = m; n = m

    ! Dynamically allocate space for matrices
    allocate(a(m,k))
    allocate(b(k,n))
    allocate(c(m,n))

    ! Parallel initialization of matrices
    print *, "Initializing matrices with random numbers..."
    call system_clock(count=start_time, count_rate=rate)
    
    !$OMP PARALLEL DO PRIVATE(j)
    do j = 1, k
        call random_number(a(:,j))
        call random_number(b(:,j))
    end do
    !$OMP END PARALLEL DO
    
    c = 0.0D0  ! Initialize matrix C as a zero matrix
    
    call system_clock(count=end_time)
    elapsed_time = end_time - start_time
    elapsed_seconds = real(elapsed_time) / real(rate)
    write(*, 100) elapsed_seconds
    write(*, 110)
100 format('Matrix initialization completed in ', E9.3, ' seconds')
110 format('')

    call print_matrix(a, 'Matrix A')
    call print_matrix(b, 'Matrix B')

    ! Get starting time for multiplication
    write(*, 120) 
    write(*, 130) m, k, k, n, m, n
120 format(/,'Computing matrix product using BLAS DGEMM subroutine')
130 format('Matrix dimensions: A(',I0,'×',I0,') × B(',I0,'×',I0,') → C(',I0,'×',I0,')')
    
    call system_clock(count=start_time, count_rate=rate)

    ! Perform matrix multiplication with a subroutine
    call matrix_multiply(a, b, c)
    
    ! Get end time and calculate elapsed time in seconds.
    call system_clock(count=end_time)
    elapsed_time = end_time - start_time
    elapsed_seconds = real(elapsed_time) / real(rate)

    ! Print out computation time
    print *, "Matrix multiplication completed."
    write (*,'(/"Calculation time is ", F6.3, " seconds.")') elapsed_seconds
    print *, ""

    ! Print out the result
    call print_matrix(c, 'Matrix C')

    ! Free used memory
    deallocate(a)
    deallocate(b)
    deallocate(c)

contains

    subroutine print_matrix(matrix, title)
        implicit none
        double precision, dimension(:,:)       :: matrix
        character(len=*), intent(in), optional :: title   
        integer :: i, j
        character(len=15) :: mformat = '(100(F14.6,1x))'
    
        if (present(title)) write(*, *) title
        do i = 1, min(5, size(matrix, 1))
            write(*, mformat) (matrix(i, j), j = 1, min(5, size(matrix, 2)))
        end do
    end subroutine print_matrix

    subroutine matrix_multiply(matA, matB, matC)
        implicit none
        double precision, dimension(:,:), intent(in)  :: matA, matB
        double precision, dimension(:,:), intent(out) :: matC
        integer :: rows, cols, cdim

        ! Get dimensions of matrices
        rows = size(matA, 1)
        cols = size(matA, 2)
        cdim = size(matB, 2)

        ! Check dimensions
        if (size(matB, 1) /= cols) then
            write(*,*) 'Error: Matrix dimensions do not match for multiplication.'
            return
        end if

        ! Perform matrix multiplication
        call dgemm('N', 'N', rows, cdim, cols, alpha, matA, rows, &
                   matB, cols, beta, matC, rows)
    end subroutine matrix_multiply
    
end program main