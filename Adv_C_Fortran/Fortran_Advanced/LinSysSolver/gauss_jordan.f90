!====================================
! Compile and execute with:
!   $ gfortran gauss_jordan.f90 -o gj
!   $ ./gj
!====================================
PROGRAM GaussJordanExample
    !use MOD4GJ
    implicit none
    INTEGER, PARAMETER :: N = 3
    REAL               :: A(N, N), B(N)
    INTEGER            :: i

    ! Initialize matrix A
    DATA A / 2.0,  3.0, -1.0, &
             3.0,  5.0,  2.0, &
             1.0, -1.0,  3.0 /
    A = transpose(A)    ! ✨just Fortran things✨
                        ! Or in other words: welome to column-major ordering
    ! Initialize vector B
    DATA B / 1.0, 1.0, 1.0 /  ! Adjust the values as needed

    ! Call the gaussj subroutine to solve the system of linear equations
    CALL gaussj(A, N, N, B, 1, N)

    ! Output the inverse matrix
    call print_matrix(A, 'Inverse Matrix:')

    ! Output the solution vector
    WRITE(*,*) 'Solution Vector:'
    DO i = 1, N
        WRITE(*,'(F12.6)') B(i)
    END DO

    contains 
    ! Needs to be contained, because the subroutine 
    ! requires explicit interface.
    subroutine print_matrix(matrix, title)
        implicit none
        REAL, dimension(:,:)                   :: matrix
        character(len=*), intent(in), optional :: title   
        ! intent(in) - information for the programmer and the 
        ! compiler that this is an INPUT; optional - not ne-
        ! cessary for the user to input it to work correctly.
        integer :: j, k
        character(len=15) :: mformat = '(100(F14.6,1x))'

        write(*, *) title
        do j = 1, min(5, size(matrix, 1))
            write(*, mformat) (matrix(j, k), k = 1, min(5, size(matrix, 2)))
        end do
    end subroutine print_matrix

END PROGRAM GaussJordanExample

SUBROUTINE gaussj(a,n,np,b,m,mp)
    INTEGER :: m,mp,n,np
    REAL    :: a(np,np),b(np,mp)
    INTEGER, PARAMETER :: NMAX=50
    ! Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a(1:n,1:n)
    ! is an input matrix stored in an array of physical dimensions np by np. b(1:n,1:m) is an in-
    ! put matrix containing the 'm' right-hand side vectors, stored in an array of physical dimensions
    ! 'np' by 'mp'. On output, a(1:n,1:n) is replaced by its matrix inverse, and b(1:n,1:m) is
    ! replaced by the corresponding set of solution vectors.
    ! Parameter: NMAX is the largest anticipated value of n.
    INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX), ipiv(NMAX) !The integer arrays ipiv, indxr, and indxc are used
    !for bookkeeping on the pivoting.
    REAL big,dum,pivinv
    do j=1,n
        ipiv(j)=0
    end do
    do i=1,n !This is the main loop over the columns to be reduced.
        big=0.
        do j=1,n !This is the outer loop of the search for a pivot element.
            if(ipiv(j).ne.1)then
                do k=1,n
                    if (ipiv(k).eq.0) then
                        if (abs(a(j,k)).ge.big)then
                            big=abs(a(j,k))
                            irow=j
                            icol=k
                        end if
                    else if (ipiv(k).gt.1) then
                        stop !’singular matrix in gaussj’
                    end if
                end do
            end if
        end do
        ipiv(icol)=ipiv(icol)+1
        ! We now have the pivot element, so we interchange rows, if needed, to put the pivot
        ! element on the diagonal. The columns are not physically interchanged, only relabeled:
        ! indxc(i), the column of the ith pivot element, is the ith column that is reduced, while
        ! indxr(i) is the row in which that pivot element was originally located. If indxr(i) 6 =
        ! indxc(i) there is an implied column interchange. With this form of bookkeeping, the
        ! solution b’s will end up in the correct order, and the inverse matrix will be scrambled
        ! by columns.
        if (irow.ne.icol) then
            do l=1,n
                dum=a(irow,l)
                a(irow,l)=a(icol,l)
                a(icol,l)=dum
            end do
            do l=1,m
                dum=b(irow,l)
                b(irow,l)=b(icol,l)
                b(icol,l)=dum
            end do
        end if
        indxr(i)=irow   ! We are now ready to divide the pivot row by the pivot
        indxc(i)=icol   ! element, located at irow and icol.
        if (a(icol,icol).eq.0.) stop ! ’singular matrix in gaussj’
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do l=1,n
            a(icol,l)=a(icol,l)*pivinv
        end do
        do l=1,m
            b(icol,l)=b(icol,l)*pivinv
        end do
        do ll=1,n ! Next, we reduce the rows...
            if(ll.ne.icol)then ! ...except for the pivot one, of course.
                dum=a(ll,icol)
                a(ll,icol)=0.
                do l=1,n
                    a(ll,l)=a(ll,l)-a(icol,l)*dum
                end do
                do l=1,m
                    b(ll,l)=b(ll,l)-b(icol,l)*dum
                end do
            end if
        end do
    enddo   ! This is the end of the main loop over columns of the reduction.
    do l=n,1,-1  ! It only remains to unscramble the solution in view
        if(indxr(l).ne.indxc(l))then    ! of the column interchanges. We do this by in-
            do k=1,n                    ! terchanging pairs of columns in the reverse order
                dum=a(k,indxr(l))       ! that the permutation was built up.
                a(k, indxr(l))=a(k,indxc(l)) 
                a(k, indxc(l))=dum           
            end do
        end if
    end do
    return !And we are done.
END SUBROUTINE gaussj
