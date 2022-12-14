! Computes the solution to system of linear equations A * X = B for GE matrices
! by useing DGESV function from the LAPACK libraries.
! Compile and execute using these commands:
!       $ gfortran -o lapack_ex lapack_ex.f90 -L/usr/local/lib -llapack -lblas
!       $ ./lapack_ex 
program solve

        implicit none
        integer(4), parameter   :: n = 3
        real(8), dimension(n)   :: x, b
        real(8), dimension(n,n) :: a
        integer(4)              :: i, info, lda, ldb, nrhs
        integer, dimension(n)   :: ipiv

        a = reshape((/   3.0, 2.0,-1.0,&
                         2.0,-2.0, 0.5,&
                        -1.0, 4.0,-1.0/),&
                        (/n,n/))
        b = (/1.0, -2.0, 0.0/)
        x = b

        nrhs = 1        ! number of right hand sides in b
        lda  = n        ! leading dimension of a
        ldb  = n        ! leading dimension of b

        call dgesv(n, nrhs, a, lda, ipiv, x, ldb, info)

        print *, 'The solution using the LAPACK subroutine is:'

        do i = 1, n
                print '("X",i1," is:",f16.6)', i, x(i)
        enddo

end program solve
