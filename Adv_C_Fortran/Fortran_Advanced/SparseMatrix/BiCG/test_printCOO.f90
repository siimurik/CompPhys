PROGRAM main
    IMPLICIT NONE
    INTEGER :: n, nmax, np, k, m_ija
    REAL    :: thresh
    INTEGER, ALLOCATABLE :: ija(:)
    REAL, ALLOCATABLE :: a(:,:), sa(:)

    n = 5
    np = n                ! Because a is n×n
    thresh = 1.0          ! Set your threshold (example: 1.0)
    nmax = 50             ! Set max size for sparse arrays (enough room)

    allocate(a(n,n))
    allocate(sa(nmax))
    allocate(ija(nmax))

    a = reshape(            &
        (/                  &
        3., 0., 1., 0., 0., &
        0., 4., 0., 0., 0., &
        0., 7., 5., 9., 0., &
        0., 0., 0., 0., 2., &
        0., 0., 0., 6., 5.  &
        /), [n,n]           &
    ); a = transpose(a)     !! ✨just Fortran things✨
    

    CALL sprsin(a, n, np, thresh, nmax, sa, ija)

    ! Print sparse matrix content
    !WRITE(*,*) "Sparse representation (SA):"
    !DO i = 1, ija(1)-1
    !    WRITE(*,*) "sa(", i, ") = ", sa(i)
    !END DO
!
    !WRITE(*,*) "Index array (IJA):"
    !DO i = 1, ija(1)
    !    WRITE(*,*) "ija(", i, ") = ", ija(i)
    !END DO

    ! Print k
    m_ija = ija(ija(1)-1)-1 !! 11
    WRITE(*,'(A,100(I5))')   'index k:', (k, k = 1, m_ija)

    ! Print ija
    WRITE(*,'(A,100(I5))')   'ija(k) :', (ija(k), k = 1, m_ija)

    ! Print sa
    WRITE(*,'(A,100(F5.1))') 'sa(k)  :', (sa(k), k = 1, m_ija)

    call print_coo(n, sa, ija)

    deallocate(a, sa, ija)
END PROGRAM main

SUBROUTINE print_coo(n, sa, ija)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL,    INTENT(IN) :: sa(*)
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


!        ija:
!        ↓      ↓      ↓      ↓      ↓      ↓
!      [ 7,     8,     8,    10,    11,    12,  ... ]
!        │      │      │     │      │      │
!        │      │      │     │      │      └─→ end of row 5's off-diags
!        │      │      │     │      └────────→ row 4 starts at index 11
!        │      │      │     └──────────────→ row 3 starts at index 10
!        │      │      └────────────────────→ row 2 has no off-diagonals
!        │      └──────────────────────────→ row 1 starts at index 8
!        └─────────────────────────────────→ off-diagonals begin at index 7
