SUBROUTINE sprsin(a, n, np, thresh, nmax, sa, ija)
    INTEGER n, nmax, np, ija(nmax)
    REAL thresh, a(np,np), sa(nmax)
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
END SUBROUTINE sprsin