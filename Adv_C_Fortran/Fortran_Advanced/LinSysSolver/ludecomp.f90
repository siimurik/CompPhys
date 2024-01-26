!====================================
! Compile and execute with:
!   $ gfortran ludecomp.f90 -o lu
!   $ ./lu
!====================================
PROGRAM LU_Decomposition_Example
    INTEGER, PARAMETER :: N = 3
    REAL :: A(N, N), B(N), d
    INTEGER :: i, j, indx(N)

    ! Initialize matrix A
    ! NB! 'DATA' command actually transposes the initial matrix.
!    DATA A /  2.0,  3.0,  1.0, &
!              3.0,  5.0, -1.0, &
!             -1.0,  2.0,  3.0 /

    ! To get around the column-major order vs. row-major order 
    ! problem, matrix A data is written into vector, which is 
    ! convenient to read and easy to understand. It is then 
    ! reshaped, which messes up the ordering. This happens with
    ! DATA as well as reshape. To fix the matrix and portray it #
    ! in the form that it was originally written, it is necessary 
    ! to use the transpose command right after.
    A = reshape((/ 2.0,  3.0, -1.0, &
                   3.0,  5.0,  2.0, &
                   1.0, -1.0,  3.0 /),(/n,n/))
    ! Transpose the matrix
    A = transpose(A)

    WRITE(*,*) 'Initial matrix:'
    DO i = 1, n
        WRITE(*, '(4(F8.4, 2X))') (A(i, j), j = 1, n)
    END DO


    ! Initialize vector B
    DATA B / 1.0, 1.0, 1.0 /  ! Adjust the values as needed
  
    call ludcmp(A, N, N, indx, d)
    call lubksb(A, N, N, indx, B)

    WRITE(*,*) 'Solution Vector:'
    do i = 1, n
        WRITE(*,'(F12.6)') b(i)
    end do
  
END PROGRAM LU_Decomposition_Example
  
SUBROUTINE ludcmp(a,n,np,indx,d)
    INTEGER             :: n,np,indx(n)
    REAL                :: d,a(np,np)
    INTEGER, PARAMETER  :: NMAX=500
    REAL, PARAMETER     :: TINY=1.0e-20 ! Largest expected n, and a small number.
    !Given a matrix a(1:n,1:n), with physical dimension 'np' by 'np', this routine replaces it by
    !the LU decomposition of a rowwise permutation of itself. 'a' and 'n' are input. 'a' is output,
    !arranged as in equation (2.3.14) above; 'indx(1:n)' is an output vector that records the
    !row permutation effected by the partial pivoting; 'd' is output as ±1 depending on whether
    !the number of row interchanges was even or odd, respectively. This routine is used in
    !combination with lubksb to solve linear equations or invert a matrix.
    INTEGER i,imax,j,k
    REAL aamax,dum,sum,vv(NMAX) !vv stores the implicit scaling of each row.
    
    d=1.0     !No row interchanges yet.
    do i=1,n !Loop over rows to get the implicit scaling informa-tion.
        aamax=0.0
        do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        end do
        if (aamax.eq.0.) stop  'singular matrix in ludcmp'! No nonzero largest element.
        vv(i)=1./aamax ! Save the scaling.
    end do
    do j=1,n !This is the loop over columns of Crout’s method.
        do i=1,j-1 !This is equation (2.3.12) except for i = j.
            sum=a(i,j)
            do k=1,i-1
                sum=sum-a(i,k)*a(k,j)
            end do
            a(i,j)=sum
        end do
        aamax=0. !Initialize for the search for largest pivot element.
        do i=j,n !This is i = j of equation (2.3.12) and i = j + 1 . . . N
            sum=a(i,j)  !of equation (2.3.13).
            do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
            end do
            a(i,j)=sum
            dum=vv(i)*abs(sum) !Figure of merit for the pivot.
            if (dum.ge.aamax) then !Is it better than the best so far?
            imax=i
            aamax=dum
            endif
        end do
        if (j.ne.imax)then !Do we need to interchange rows?
            do k=1,n !Yes, do so...
                dum=a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k)=dum
            end do
            d=-d !...and change the parity of d.
            vv(imax)=vv(j) !Also interchange the scale factor.
        end if
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        !If the pivot element is zero the matrix is singular (at least to the precision of the al-
        !gorithm). For some applications on singular matrices, it is desirable to substitute TINY
        !for zero.
        if(j.ne.n)then !Now, finally, divide by the pivot element.
            dum=1./a(j,j)
            do i=j+1, n
                a(i,j)=a(i,j)*dum
            end do
        endif
    end do !Go back for the next column in the reduction.
    return
END SUBROUTINE ludcmp

SUBROUTINE lubksb(a,n,np,indx,b)
    INTEGER n,np,indx(n)
    REAL a(np,np),b(n)
    !Solves the set of n linear equations A · X = B. Here a is input, not as the matrix A but
    !rather as its LU decomposition, determined by the routine ludcmp. indx is input as the
    !permutation vector returned by ludcmp. b(1:n) is input as the right-hand side vector B,
    !and returns with the solution vector X. a, n, np, and indx are not modified by this routine
    !and can be left in place for successive calls with different right-hand sides b. This routine
    !takes into account the possibility that b will begin with many zero elements, so it is efficient
    !for use in matrix inversion.
    INTEGER i,ii,j,ll
    REAL sum
    ii=0 !When ii is set to a positive value, it will become the in-
    !dex of the first nonvanishing element of b. We now do
    !the forward substitution, equation (2.3.6). The only new
    !wrinkle is to unscramble the permutation as we go.
    do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
            do j=ii,i-1
                sum=sum-a(i,j)*b(j)
            end do
        else if (sum.ne.0.) then
            ii=i !A nonzero element was encountered, so from now on we will
            !have to do the sums in the loop above.
        endif
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

SUBROUTINE PRINT_VECTOR( DESC, N, A )
    CHARACTER*(*)    DESC
    INTEGER          N
    REAL          A( N )
    INTEGER          I

    WRITE(*,*)
    WRITE(*,*) DESC
    WRITE(*,9999) ( A( I ), I = 1, N )
9999 FORMAT( 11(:,1X,I6) )
    RETURN
END