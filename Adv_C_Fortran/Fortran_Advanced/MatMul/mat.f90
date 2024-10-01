!==============================================
! Compile and execute with:
!   $ gfortran -O3 mat.f90 -o mat -lblas -fopenmp
!   $ ./mat
!==============================================
program main
    implicit none
    !integer :: i, j
    double precision, parameter :: alpha = 1.0D0, beta = 0.0D0
    integer, parameter          :: dim = 5000
    integer                     :: m, k, n
    double precision, allocatable, dimension(:,:) :: a, b, c
    double precision :: elapsed_seconds
    integer          :: start_time, end_time, elapsed_time, rate

    ! Declare symmetric matrix
    m = dim; k = m; n = m

    ! Dynamically allocate space for matrices
    allocate(a(m,k))
    allocate(b(k,n))
    allocate(c(m,n))

    ! Generate random matrices with real values between 0 and 1.
    call random_number(a)
    call random_number(b)
    c = 0.0D0;  ! Decalare matrix C as a zero matrix.

    call print_matrix(a, 'Matrix A')
    call print_matrix(b, 'Matrix B')

    ! Get starting time
    PRINT *, "Computing matrix product using BLAS DGEMM subroutine"
    CALL SYSTEM_CLOCK(count=START_TIME, count_rate=RATE)

    ! Perform matrix multplication with a subroutine
    call matrix_multiply(A, B, C)
    
    ! Get end time and calculate elapsed time in seconds.
    CALL SYSTEM_CLOCK(count=END_TIME)
    ELAPSED_TIME = END_TIME - START_TIME
    ELAPSED_SECONDS = REAL(ELAPSED_TIME) / REAL(RATE)

    ! Print out computation time
    PRINT *, "Computations completed."
    WRITE (*,15) ELAPSED_SECONDS
    PRINT *, ""
15  FORMAT(/'Calculation time is ', F6.3, ' seconds.')

    ! Print out the result
    call print_matrix(c, 'Matrix C')

    ! Free used memory
    deallocate(a)
    deallocate(b)
    deallocate(c)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains

!----------------------------------------------------------------------------
    
    subroutine print_matrix(matrix, title)
        implicit none
        double precision, dimension(:,:)       :: matrix
        character(len=*), intent(in), optional :: title   
        ! intent(in) - information for the programmer and the 
        ! compiler that this is an INPUT; optional - not ne-
        ! cessary for the user to input it to work correctly.
        integer :: i, j
        character(len=15) :: mformat = '(100(F14.6,1x))'
    
        write(*, *) title
        do i = 1, min(5, size(matrix, 1))
            write(*, mformat) (matrix(i, j), j = 1, min(5, size(matrix, 2)))
        end do
    end subroutine print_matrix

!----------------------------------------------------------------------------

    subroutine matrix_multiply(matA, matB, matC)
        implicit none

        double precision, dimension(:,:), intent(in)  :: matA, matB
        double precision, allocatable, dimension(:,:), intent(inout) :: matC
        integer :: rows, cols, cdim
        double precision, parameter :: pALPHA = 1.0D0, pBETA = 0.0D0;

        ! Get dimensions of matrices
        rows = size(matA, 1)
        cols = size(matA, 2)
        cdim = size(matB, 2)

        ! Check dimensions
        if (size(matB, 1) .ne. cols) then
            write(*,*) 'Error: Matrix dimensions do not match for multiplication.'
            return
        end if

        ! Perform matrix multiplication
        ! 'N' or 'n',  op( A ) = A.
        call DGEMM('N', 'N', rows, cdim, cols, pALPHA, matA, rows, &
                                    matB, cols, pBETA, matC, rows)

        ! You can use the resulting matrix C as needed
    end subroutine matrix_multiply
    
end program main
