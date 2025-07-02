SUBROUTINE lubksb(a,n,np,indx,b)
    INTEGER n,np,indx(n)
    REAL a(np,np),b(n)
    ! Solves the set of n linear equations A · X = B. Here a is input, not as the matrix A but
    ! rather as its LU decomposition, determined by the routine ludcmp. indx is input as the
    ! permutation vector returned by ludcmp. b(1:n) is input as the right-hand side vector B,
    ! and returns with the solution vector X. a, n, np, and indx are not modified by this routine
    ! and can be left in place for successive calls with different right-hand sides b. This routine
    ! takes into account the possibility that b will begin with many zero elements, so it is efficient
    ! for use in matrix inversion.
    INTEGER i,ii,j,ll
    REAL sum
    ii=0            ! When ii is set to a positive value, it will become the in-
    do i=1,n        ! dex of the first nonvanishing element of b. We now do
        ll=indx(i)  ! the forward substitution, equation (2.3.6). The only new
        sum=b(ll)   ! wrinkle is to unscramble the permutation as we go.
        b(ll)=b(i)
        if (ii.ne.0)then
            do j=ii,i-1
                sum=sum-a(i,j)*b(j)
            end do
        else if (sum.ne.0.) then
            ii=i ! A nonzero element was encountered, so from now on we will
        endif    ! have to do the sums in the loop above.
        b(i)=sum
    end do
    do i=n,1,-1 ! Now we do the backsubstitution, equation (2.3.7).
        sum=b(i)
        do j=i+1,n
            sum=sum-a(i,j)*b(j)
        end do
        b(i)=sum/a(i,i) ! Store a component of the solution vector X.
    end do
    return ! All done!
END SUBROUTINE lubksb