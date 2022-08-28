!========================================================
! Compile and execute with
!       $ gfortran hilbert_matrix.f90 -o hilbert_matrix
!       $ ./hilbert_matrix
!========================================================

program hilbert_matrix

    implicit none
    integer(4), parameter      :: n = 4
    integer(4)                 :: i, j
    real(8), dimension(n,n)    :: H

    call hilb(n, H)

    write (*,9) n
    do i = 1, n 
        write (*,10) ( H(i,j), j = 1, n ) 
    end do
    write (*,11) 

9   format (/' Hilbert matrix with dimension ', i2)
10  format ( 6f12.4 )
11  format (/)

end program hilbert_matrix

subroutine hilb(n, H)
!=========================================================== 
! Evaluate the values of the Hilbert matrix following the
! formula:
! H(i,j) = 1 / (i + j - 1)
! 
! Input:
! ========
! n - dimension of the matrix
!
! Output:
! ========
! H - Hilbert matrix with dimension n x n
!
! Example for n = 3:
! =======
!
! Hilbert matrix with dimension  3
!      1.0000      0.5000      0.3333
!      0.5000      0.3333      0.2500
!      0.3333      0.2500      0.2000
!
!=========================================================== 
    implicit none
    real(8), dimension(n,n) :: H
    integer(4) :: i, j, n

    do i = 1, n
        do j = 1, n
            H(i,j) = 1.0 / real(i + j - 1)
        end do
    end do
    return

end subroutine hilb