!==============================================
! Compile and execute with:
!   $ gfortran matmult.f90 -o mm -llapack -lblas
!   $ ./mm
!==============================================
module matrix_module
    implicit none

    contains

    subroutine matrix_multiply(A, B, C)
        implicit none

        double precision, dimension(:,:), intent(in)  :: A, B
        double precision, allocatable, dimension(:,:) :: C
        integer :: M, N, K
        double precision :: ALPHA, BETA

        ! Get dimensions of matrices
        M = size(A, 1)
        K = size(A, 2)
        N = size(B, 2)

        ! Set default values
        ALPHA = 1.0D0
        BETA  = 0.0D0

        ! Check dimensions
        if (size(B, 1) .ne. K) then
            write(*,*) 'Error: Matrix dimensions do not match for multiplication.'
            return
        end if

        ! Allocate C within the subroutine
        allocate(C(M, N))

        ! Perform matrix multiplication
        call DGEMM('N', 'N', M, N, K, ALPHA, A, M, B, K, BETA, C, M)

        ! You can use the resulting matrix C as needed
    end subroutine matrix_multiply
!====================================================================
    function matfunc(A, B) result(C)
        double precision, dimension(:,:), intent(in)  :: A, B    ! intent(in) - information for the programmer and the compiler that this is an INPUT
        double precision, allocatable, dimension(:,:) :: C
        integer :: M, N, K
        double precision :: ALPHA, BETA

        ! Get dimensions of matrices
        M = size(A, 1)
        K = size(A, 2)
        N = size(B, 2)

        ! Set default values
        ALPHA = 1.0D0
        BETA  = 0.0D0

        ! Check dimensions
        if (size(B, 1) /= K) then
            write(*,*) 'Error: Matrix dimensions do not match for multiplication.'
            C = 0.0D0
            return
        end if

        ! Perform matrix multiplication
        allocate(C(M, N))
        call DGEMM('N', 'N', M, N, K, ALPHA, A, M, B, K, BETA, C, M)

        ! You can use the resulting matrix C as needed
    end function matfunc
!====================================================================
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
    
end module matrix_module

program main
    use matrix_module
    implicit none
    !integer             :: i, j
    integer, parameter  :: M=5000, K=5000, N=5000
    double precision, dimension(M,K) :: A
    double precision, dimension(K,N) :: B
    double precision, allocatable, dimension(:,:) :: C1!, C2
    integer :: start_time, end_time, elapsed_time, rate
    real    :: elapsed_seconds
    !character(len=15) :: mformat='(100(F14.6,1x))'

    print *, "Initializing matrix data."
    print *, ""
    call random_number(A)
    call random_number(B)

    ! Allocate C before calling the subroutine
    !-----------------------------------------------------------------
    ! If C is already allocated in the subroutine or function then
    ! this section along with 'deallocate' at the end won't have much
    ! point. This is here just in case to show an alternative approach
    ! and provide an option for more control in memory management.
    !-----------------------------------------------------------------
    !allocate(C1(M, N))
    !allocate(C2(M, N))

    call SYSTEM_CLOCK(count=start_time, count_rate=rate)

    ! Two approaches to do matrix multplication
    call matrix_multiply(A, B, C1)  ! subroutine
    !C2 = matfunc(A, B)              ! function

    call SYSTEM_CLOCK(count=end_time)
    elapsed_time = end_time - start_time
    elapsed_seconds = real(elapsed_time) / real(rate)

    write(*,*) "Computations completed."
    write (*,15) elapsed_seconds
    write(*,*) ""

    !write(*,*) 'Matrix C values:'
    !do i = 1, 5
    !    write(*,mformat) (C1(i,j), j = 1, 5)
    !end do
    call print_matrix(C1, 'Matrix C1 values:')
    !call print_matrix(C2)

15  format(/'Calculation time is', F6.3, ' seconds.')

    ! Deallocate matrix C
    !deallocate(C1)
    !deallocate(C2)

end program main
