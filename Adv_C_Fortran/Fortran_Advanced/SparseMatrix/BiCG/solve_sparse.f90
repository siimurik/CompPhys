!=========================================================
!   gfortran solve_sparse.f90 linbcg.f90 -o solsparse
!=========================================================
PROGRAM main
    IMPLICIT NONE

    ! -------------------------------
    ! PARAMETERS AND CONSTANTS
    ! -------------------------------
    integer, parameter :: NMAX = 1000  !! Make sure this matches
    integer, parameter :: n = 3 !! Dimension of the system of equations
    
    ! -------------------------------
    ! VARIABLE DECLARATIONS
    ! -------------------------------
    ! Matrix and vectors
    double precision, allocatable, dimension(:,:) :: a
    double precision, allocatable, dimension(:)   :: b
    double precision, allocatable, dimension(:)   :: x
    
    ! Sparse storage
    double precision, dimension(NMAX) :: sa  = 0.0D0
    integer,          dimension(NMAX) :: ija = 0
    common /mat/ sa, ija    !! Matrix A data in sparse format, which gets shared
                            !! in subroutines atimes and solve. Check linbcg().
    ! Solver parameters
    integer :: itol = 1     !! convergence test options: 1, 2, 3, 4
    integer :: iter = 0     !! Initialized to zero to avoid silent errors.
    integer :: np   = n     !! EXTREMELY IMPORTANT. DO NOT TOUCH. 
    integer :: itmax = 1000
    double precision :: thresh = 0.D0   !! EXTREMELY IMPORTANT. DO NOT TOUCH. 
    double precision :: err = 0.0D-6    !! Initialized to zero to avoid silent errors.
    double precision :: tol = 1.0D-6
    
    ! Utility variables
    integer :: i, k, m_ija, idebug
    character(len=15) :: mformat = '(100(F14.6,1x))'

    ! Status variables
    integer :: alloc_stat

    ! Timing variables
    real :: start_time, end_time

    ! -------------------------------
    ! INITIALIZATION
    ! -------------------------------
    ! Start timer
    call cpu_time(start_time)

    ! Allocate arrays with error checking
    allocate(a(n,n), x(n), b(n), stat=alloc_stat)
    if (alloc_stat .ne. 0) then
        write(*,*) "Error: Failed to allocate memory!"
        stop 1
    end if

    ! Initialize matrix - ROW-MAJOR ORDER
    a = reshape([4.0d0, -1.0d0,  0.0d0,  &  ! First ROW
                -1.0d0,  4.0d0, -1.0d0,  &  ! Second ROW
                 0.0d0, -1.0d0,  4.0d0], &  ! Third ROW
                [n,n])

    ! Initialize vectors
    b = [2.d0, 4.d0, 10.d0]      ! Right-hand side
    x = [0.d0, 0.d0, 0.d0]       ! Initial guess

    ! -------------------------------
    ! PRINT SYSTEM INFORMATION
    ! -------------------------------
    write(*,*) "========================================"
    write(*,*) " ITERATIVE BICONJUGATE GRADIENT METHOD"
    write(*,*) "----------------------------------------"
    write(*,*) " SPARSE MATRIX SOLVER"
    write(*,'(A,I2,A,I2)') "  System size: ", n, " x", n
    write(*,*) " Max iterations: ", NMAX
    write(*,*) " Tolerance: ", TOL
    write(*,*) "========================================"
    write(*,*)

    ! TODO: Add subroutine for printing matrix A
    ! do i = 1, min(5, size(a, 1))
    !     write(*, mformat) (a(i, j), j = 1, min(5, size(a, 2)))
    ! end do

    call print_matrix("Original dense matrix A:", a)
    call print_vector("Right-hand side vector b:", b)
    call print_vector("Initial guess x:", x)

    ! -------------------------------
    ! CONVERT TO SPARSE FORMAT
    ! -------------------------------
    ! A super sneaky step, which sends 'sa' and 'ija' to be used in
    ! subroutines atimes() and asolve(), which themselves are used  
    ! in the middle of the biconjugate gradient solver linbcg().
    call dsprsin(a, n, np, thresh, nmax, sa, ija)
    
    idebug = 0
    if (idebug .ne. 0) then
        m_ija = ija(ija(1)-1)-1 !! 11
        WRITE(*,'(A,100(I5))')   'index k:', (k, k = 1, m_ija)
        WRITE(*,'(A,100(I5))')   'ija(k) :', (ija(k), k = 1, m_ija)
        WRITE(*,'(A,100(F5.1))') 'sa(k)  :', (sa(k), k = 1, m_ija)
        call print_coo(n, sa, ija)
    end if

    ! -------------------------------
    ! SOLVE LINEAR SYSTEM
    ! -------------------------------
    write(*,*)
    write(*,*) "Solving system Ax = b..."
    
    ! Set callback for progress reporting
    !call set_progress_callback(progress_report)

    ! -------------------------------
    ! Call the solver
    ! -------------------------------
    ! linbcg() uses asolve() and atimes() which stores 
    ! matrix A into a new vector x thanks to the COMMON
    ! operator, which shares sparse matrix values 'sa'
    ! and index location data 'ija' between different 
    ! subroutines. This can be a little confusing since
    ! we don't see exactly how A gets stored into x in
    ! the main section of the program.
    call linbcg(n, b, x, itol, tol, itmax, iter, err)

    ! -------------------------------
    ! OUTPUT RESULTS
    ! -------------------------------
    call cpu_time(end_time)
    
    ! Print results
    write(*,*)
    write(*,*) "========================================"
    write(*,*) " SOLVER RESULTS"
    write(*,*) "========================================"
    write(*,'(A,I5)') " Iterations: ", iter
    write(*,'(A,E12.4)') " Final error: ", err
    write(*,'(A,E10.4,A)') " Solution time: ", end_time-start_time, " seconds"
    write(*,*)

    call print_vector("Solution vector x:", x)

    ! -------------------------------
    ! CLEAN UP
    ! -------------------------------
    deallocate(a, x, b, stat=alloc_stat)
    if (alloc_stat /= 0) then
        write(*,*) "Warning: Memory deallocation error!"
    end if
    
    write(*,*)
    write(*,*) "Program completed successfully."

    CONTAINS

        ! Helper subroutine to print matrices
        subroutine print_matrix(description, mat)
            character(len=*), intent(in) :: description
            double precision, intent(in) :: mat(:,:)
            
            write(*,*)
            write(*,*) description
            do i = 1, size(mat,1)
                write(*, mformat) mat(i,:)
            end do
        end subroutine print_matrix

        ! Helper subroutine to print vectors
        subroutine print_vector(description, vec)
            character(len=*), intent(in) :: description
            double precision, intent(in) :: vec(:)
            
            write(*,*)
            write(*,*) description
            write(*, mformat) vec
        end subroutine print_vector

END PROGRAM main

SUBROUTINE dsprsin(a, n, np, thresh, nmax, sa, ija)
    IMPLICIT NONE
    ! Converts a square matrix a(1:n,1:n) with physical dimension np into row-indexed sparse
    ! storage mode. Only elements of a with magnitude ≥thresh are retained. Output is in
    ! two linear arrays with physical dimension nmax (an input parameter): sa(1:) contains
    ! array values, indexed by ija(1:). The logical sizes of sa and ija on output are both
    ! ija(ija(1)-1)-1 (see text).
    
    ! Input parameters
    INTEGER, INTENT(IN) :: n, np, nmax
    DOUBLE PRECISION, INTENT(IN) :: thresh, a(np,np)
    
    ! Output arrays
    DOUBLE PRECISION, INTENT(OUT) :: sa(nmax)
    INTEGER, INTENT(OUT) :: ija(nmax)
    
    ! Local variables
    INTEGER :: i, j, k

    do j = 1, n         !! Store diagonal elements.
        sa(j) = a(j,j)
    end do
    ija(1) = n+2        !! Index to 1st row off-diagonal element, if any.
    k = n+1
    do i = 1, n         !! Loop over rows.
        do j = 1, n     !! Loop over columns.
            if (abs(a(i,j)) .ge. thresh) then
                if(i.ne.j)then  !! Store off-diagonal elements and their columns.
                    k=k+1
                    if (k .gt. nmax) stop 'nmax too small in sprsin'
                    sa (k) = a(i,j)
                    ija(k) = j
                endif
            endif
        end do
        ija(i+1)=k+1    !! As each row is completed, store index to next.
    end do
    return
END SUBROUTINE dsprsin

SUBROUTINE print_coo(n, sa, ija)
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: n
    INTEGER, INTENT(IN)          :: ija(*)
    DOUBLE PRECISION, INTENT(IN) :: sa(*)
    
    INTEGER :: i, k_start, k_end, k

    ! Print only nonzero diagonals
    DO i = 1, n
        IF (sa(i) .NE. 0.0) THEN
            WRITE(*,'(A,I0,A,I0,A,F8.2)') '(', i, ',', i, ') ', sa(i)
        END IF
    END DO

    ! Print off-diagonal elements
    DO i = 1, n
        k_start = ija(i)
        k_end = ija(i+1) - 1
        DO k = k_start, k_end
            WRITE(*,'(A,I0,A,I0,A,F8.2)') '(', i, ',', ija(k), ') ', sa(k)
        END DO
    END DO
END SUBROUTINE print_coo