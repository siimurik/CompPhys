! Solves a system of linear equations using LAPACK.
! Compile and execute with:
!       $ gfortran -o example example.f90 -L/usr/local/lib -llapack -lblas
!       $ ./example
program LinearEquations
        
        implicit none
        integer(4), parameter   :: n = 3
        real(8), dimension(n,n) :: A(n,n) 
        real(8), dimension(n)   :: b(n), x(n)
        integer(4)              :: i, info, ldA, ldb, nrhs
        integer(4), dimension(n):: ipiv
        
        ! Matrix A
        A(1,:) = (/3.0, 1.0, 3.0/)
        A(2,:) = (/1.0, 5.0, 9.0/)
        A(3,:) = (/2.0, 6.0, 5.0/)
        
        ! Vector b
        b = (/-1.0, 3.0, -3.0/)
        
        ! Initialize vector x as zeros won't work
        ! for some odd reason. Setting x == b will
        ! do the trick though.
        x = b

        nrhs = 1        ! number of right hand sides in b
        ldA  = n        ! leading dimension of A
        ldb  = n        ! leading dimension of b

        ! Find the solution using the LAPACK routine DGESV
        call DGESV(n, nrhs, A, ldA, ipiv, x, ldb, info)
        
        print *, 'The solution to the system of equations:'
        do i = 1, n
                write(*,9) i, x(i)
        end do

9       format('x[',i1, ']= ', f5.2)
end program LinearEquations
