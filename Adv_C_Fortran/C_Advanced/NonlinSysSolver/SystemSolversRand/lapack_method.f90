! Computes the solution to system of linear equations A * X = B for GE matrices
! by useing DGESV function from the LAPACK libraries.
! Compile and execute using these commands:
!       $ gfortran -o lapack_ex lapack_ex.f90 -L/usr/local/lib -llapack -lblas
!       $ ./lapack_ex 
program solve

        implicit none
        integer(4), parameter    :: n = 5000
        real(8), dimension(n)    :: x, b
        real(8), dimension(n,n)  :: a
        integer(4)               :: i, j, info, lda, ldb, nrhs
        integer(4), dimension(n) :: ipiv
        real(8)                  :: t1, t2, time

        do i = 1, n-1
            do j = 1, n-1
                a(i, j) = 4.D0*(dfloat(i) + 9.D0 + dfloat(j))
                b(j)    = 4.D0*(11.D0-dfloat(i)+30.D0*dfloat(j))
            enddo
        enddo
        x = b

        nrhs = 1        ! number of right hand sides in b
        lda  = n        ! leading dimension of a
        ldb  = n        ! leading dimension of b

        call cpu_time(t1)
        call dgesv(n, nrhs, a, lda, ipiv, x, ldb, info)
        call cpu_time(t2)
        time = t2 - t1

        write (*,100) time
        100 format (/'Calculation time is ', e9.3, ' seconds.')

!        do i = 1, n
!                print '("X",i1," is:",f16.6)', i, x(i)
!        enddo

end program solve
