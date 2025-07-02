SUBROUTINE ludcmp(a,n,np,indx,d)
    INTEGER             :: n, np, indx(n)
    REAL                :: d, a(np,np)
    INTEGER, PARAMETER  :: NMAX = 500
    REAL, PARAMETER     :: TINY = 1.0e-20 ! Largest expected n, and a small number.
    ! Given a matrix a(1:n,1:n), with physical dimension 'np' by 'np', this routine replaces it by
    ! the LU decomposition of a rowwise permutation of itself. 'a' and 'n' are input. 'a' is output,
    ! arranged as in equation (2.3.14) above; 'indx(1:n)' is an output vector that records the
    ! row permutation effected by the partial pivoting; 'd' is output as ±1 depending on whether
    ! the number of row interchanges was even or odd, respectively. This routine is used in
    ! combination with lubksb to solve linear equations or invert a matrix.
    INTEGER :: i,imax,j,k
    REAL    :: aamax,dum,sum,vv(NMAX) !vv stores the implicit scaling of each row.
    
    d = 1.0     ! No row interchanges yet.
    do i = 1, n ! Loop over rows to get the implicit scaling informa-tion.
        aamax = 0.0
        do j = 1, n
            if (abs(a(i,j)) .gt. aamax) aamax=abs(a(i,j))
        end do
        if (aamax .eq. 0.0) stop 'singular matrix in ludcmp'! No nonzero largest element.
        vv(i) = 1.0/aamax ! Save the scaling.
    end do
    do j = 1, n ! This is the loop over columns of Crout’s method.
        do i = 1, j-1 ! This is equation (2.3.12) except for i = j.
            sum = a(i,j)
            do k = 1,i-1
                sum = sum-a(i,k)*a(k,j)
            end do
            a(i,j) = sum
        end do
        aamax = 0.0 ! Initialize for the search for largest pivot element.
        do i = j, n ! This is i = j of equation (2.3.12) and i = j + 1 ... N
            sum = a(i,j)  ! of equation (2.3.13).
            do k = 1, j-1
                sum = sum - a(i,k)*a(k,j)
            end do
            a(i,j) = sum
            dum = vv(i)*abs(sum) ! Figure of merit for the pivot.
            if (dum .ge. aamax) then ! Is it better than the best so far?
                imax  = i
                aamax = dum
            end if
        end do
        if (j .ne. imax) then ! Do we need to interchange rows?
            do k = 1, n ! Yes, do so...
                dum       = a(imax,k)
                a(imax,k) = a(j,k)
                a(j,k)    = dum
            end do
            d = -d ! ...and change the parity of d.
            vv(imax) = vv(j) ! Also interchange the scale factor.
        end if
        indx(j) = imax
        if (a(j,j) .eq. 0.0) a(j,j) = TINY
        ! If the pivot element is zero the matrix is singular (at least to the precision of the al-
        ! gorithm). For some applications on singular matrices, it is desirable to substitute TINY
        ! for zero.
        if (j .ne. n) then ! Now, finally, divide by the pivot element.
            dum = 1./a(j,j)
            do i = j+1, n
                a(i,j) = a(i,j)*dum
            end do
        endif
    end do ! Go back for the next column in the reduction.
    return
END SUBROUTINE ludcmp