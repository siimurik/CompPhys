program addmats
        implicit none

        integer, parameter :: dimmat = 3
        real, dimension(dimmat, dimmat) :: a, b, c
        integer :: i, j

        ! Important to initialize the matrices.
        ! If this is not done, values are assumed
        ! as empty, which can result in weird 
        ! values in the matrix elements.
        a = 0.0; b = a; c = a
        
        ! This creates the matrices
        a(1,2) = 2.0
        do i = 2, dimmat-1
                a(i,i+1) = 2.0
                b(i,i-1) = 1.0
        enddo

        b(dimmat, dimmat-1) = 1.0                 
        
        !call print_matrix(a, dimmat)
        
        ! This adds the matrices a and b
        do i = 1, dimmat
                do j = 1, dimmat
                        c(i,j) = a(i,j) + b(i,j)
                enddo
        enddo
        
        ! This prints c
        do i = 1, dimmat
                write(*,*) ( c(i,j), j = 1, dimmat )
        enddo

end program addmats

subroutine print_matrix(matrix, dimmat)
        implicit none
        integer :: dimmat 
        real, dimension(dimmat,dimmat) :: matrix
        integer :: i, j

        do i = 1, size(matrix, 1)
                write (*, '(999F12.6)') (matrix(i, j), j = 1, size(matrix, 2))
        end do
end subroutine print_matrix
