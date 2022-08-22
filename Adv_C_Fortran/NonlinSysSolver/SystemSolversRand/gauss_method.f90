program main
!====================================================================
!
!====================================================================
    implicit none
    integer(4), parameter :: n = 1500
    integer(4) :: i, j
    real(8)    :: a(n,n), b(n), x(n)
    real(8)    :: t1, t2, time

    do i = 1, n-1
        do j = 1, n-1
            a(i, j) = 4.D0*(dfloat(i) + 9.D0 + dfloat(j))
            b(j)    = 4.D0*(11.D0-dfloat(i)+30.D0*dfloat(j))
        enddo
    enddo

!    do i = 1, n
!        print *, (a(i,j),j=1,n), b(i)
!    end do

    ! Gauss elination and time measurement
    call cpu_time(t1)
    call gauss_2(a,b,x,n)
    call cpu_time(t2)
    time = t2 - t1

    ! print matrix A and vector b after the elimination 
!    write (*,202)
!        do i = 1,n
!        write (*,201)  (a(i,j),j=1,n), b(i)
!    end do

    ! print solutions
!    write (*,200)
!    write (*,201) (x(i),i=1,n)
!    write (*,203)
    write (*,204) time
!    200 format (' Gauss elimination with scaling and pivoting ',/,/,/, &
!                ' Matrix A and vector b')
!    201 format (6f12.6)
!    202 format (/,' Matrix A and vector b after elimination')
!    203 format (/,' Solutions x(n)')
    204 format (/'Calculation time is ', e9.3, ' seconds.')
!    204 format (/'Calculation time is ', f10.8, ' seconds.')

end program main
    
subroutine gauss_2(a,b,x,n)
!============================================================
! Solutions to a system of linear equations A*x=b
! Method: the basic elimination (simple Gauss elimination)
! Alex G. November 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! b(n)   - vector of the right hand coefficients b
! n      - number of equations (size of matrix A)
! output ...
! x(n)   - solutions
! comments ...
! the original arrays a(n,n) and b(n) will be destroyed 
! during the calculation
!===========================================================
    implicit none
    integer(4)  :: n
    real(8)     :: a(n,n), b(n), x(n), s(n)
    real(8)     :: c, pivot, store
    integer(4)  :: i, j, k, l

! step 1: begin forward elimination
    do k = 1, n-1

! step 2: "scaling"
! s(i) will have the largest element from row 1
        do i = k, n             ! loop over rows
            s(i) = 0.0
            do j = k, n         ! loop over elements of row i
                s(i) = max(s(i), abs(a(i,j)))
            end do
        end do

! step 3: "pivoting 1"
! find a row with the largest pivoting element
        pivot = abs(a(k,k)/s(k))
        l = k
        do j = k+1, n
            if(abs(a(j,k)/s(j)) > pivot) then 
                pivot = abs(a(j,k)/s(j))
                l = j
            end if
        end do

! Check if the system has a sigular matrix
        if (pivot == 0.0) then
            write(*,*) ' The matrix is sigular '
            return
        end if

! step 4: "pivoting 2" interchange rows k and l (if needed)
        if (l /= k) then
            do j = k, n
                store  = a(k,j)
                a(k,j) = a(l,j)
                a(l,j) = store
            end do
            store = b(k)
            b(k)  = b(l)
            b(l)  = store
        end if
    
! step 5: the elimination (after scalinga and pivoting)
        do i = k+1, n
            c      = a(i,k)/a(k,k)
            a(i,k) = 0.0
            b(i)   = b(i) - c*b(k)
            do j = k+1, n
                a(i,j) = a(i,j) - c*a(k,j)
            end do
        end do
    end do      ! end of forward elimination
    
! step 6: back substitution
    x(n) = b(n)/a(n,n)
    do i = n-1, 1, -1
        c = 0.0
        do j = i+1, n
            c = c + a(i,j)*x(j)
        end do
        x(i) = (b(i) - c)/a(i,i)
    end do

end subroutine gauss_2