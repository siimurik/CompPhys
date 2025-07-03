PROGRAM test_matvec
    IMPLICIT NONE
    integer, parameter :: n = 3
    double precision, allocatable, dimension(:,:) :: a
    double precision, allocatable, dimension(:) :: x, b
    integer :: i, j
    character(len=15) :: mformat = '(100(F14.6,1x))'
    integer ::  np, k, m_ija
    double precision    :: thresh
    integer, parameter :: nmax = 50    ! 
    INTEGER, allocatable :: ija(:)
    double precision, allocatable :: sa(:)

    allocate(a(n,n))
    allocate(x(n))
    allocate(b(n))
    allocate(sa(nmax))
    allocate(ija(nmax))
    sa = 0.d0
    ija = 0

    ! Initialize matrix - ROW-MAJOR ORDER
    a = reshape([4.0d0, -1.0d0,  0.0d0, &   ! First ROW
                -1.0d0,  4.0d0, -1.0d0, &   ! Second ROW
                 0.0d0, -1.0d0,  4.0d0], &  ! Third ROW
                [n,n])

    do i = 1, min(5, size(a, 1))
        write(*, mformat) (a(i, j), j = 1, min(5, size(a, 2)))
    end do

    x = [1.0d0, 2.0d0, 3.0d0]
    np = n          !! Two most important values
    thresh = 0.d0   !! Two most important values
    call dsprsin(a, n, np, thresh, nmax, sa, ija)
    
    m_ija = ija(ija(1)-1)-1 !! 11
    WRITE(*,'(A,100(I5))')   'index k:', (k, k = 1, m_ija)
    WRITE(*,'(A,100(I5))')   'ija(k) :', (ija(k), k = 1, m_ija)
    WRITE(*,'(A,100(F5.1))') 'sa(k)  :', (sa(k), k = 1, m_ija)
    call print_coo(n, sa, ija)

    b = [0.d0, 0.d0, 0.d0]
    call dsprsax(sa, ija, x, b, n)
    do i = 1, n
        print *, b(i)
    end do

    deallocate(a,x,b,ija,sa)

END PROGRAM test_matvec

SUBROUTINE dsprsin(a, n, np, thresh, nmax, sa, ija)
    INTEGER n, nmax, np, ija(nmax)
    double precision thresh, a(np,np), sa(nmax)
    ! Converts a square matrix a(1:n,1:n) with physical dimension np into row-indexed sparse
    ! storage mode. Only elements of a with magnitude ≥thresh are retained. Output is in
    ! two linear arrays with physical dimension nmax (an input parameter): sa(1:) contains
    ! array values, indexed by ija(1:). The logical sizes of sa and ija on output are both
    ! ija(ija(1)-1)-1 (see text).
    INTEGER i, j, k
    do j = 1, n     !! Store diagonal elements.
        sa(j) = a(j,j)
    end do
    ija(1) = n+2      !! Index to 1st row off-diagonal element, if any.
    k = n+1
    do i = 1, n     !! Loop over rows.
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
    INTEGER, INTENT(IN) :: n
    double precision,    INTENT(IN) :: sa(*)
    INTEGER, INTENT(IN) :: ija(*)
    
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

SUBROUTINE dsprsax(sa,ija,x,b,n)
    INTEGER n,ija(*)
    DOUBLE PRECISION b(n),sa(*),x(n)
    ! Multiply a matrix in row-index sparse storage arrays sa and ija by a vector x(1:n), giving
    ! a vector b(1:n).
    INTEGER i,k
    if (ija(1).ne.n+2) stop 'mismatched vector and matrix in dsprsax'
    do i = 1, n
        b(i) = sa(i)*x(i)           !! Start with diagonal term.
        do k = ija(i), ija(i+1)-1   !! Loop over off-diagonal terms.
            b(i) = b(i) + sa(k)*x(ija(k))
        end do
    end do
    return
END SUBROUTINE dsprsax
