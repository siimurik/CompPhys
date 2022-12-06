! Compile and execute with
!       $ gfortran read_matrix.f90 -o read_matrix
!       $ ./read_matrix

!======================================================================
! This code's main purpose is to demonstrate how to read
! data values from a separate file into a matrix and use
! that symmetrical matrix to find its eigenvalues and
! eigenvectors.
!
! The "data.dat" files contains the following values
!    1.96  -6.49  -0.47  -7.20  -0.65
!   -6.49   3.80  -6.39   1.50  -6.34
!   -0.47  -6.39   4.17  -1.51   2.67
!   -7.20   1.50  -1.51   5.70   1.80
!   -0.65  -6.34   2.67   1.80  -7.10
!
! Description.
! ============
!
! The routine computes all eigenvalues and, optionally, eigenvectors of an
! n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies
!
! A*v(j) = lambda(j)*v(j)
!
! where lambda(j) is its eigenvalue. The computed eigenvectors are NOT
! orthonormal!

!  Example Program Results.
!  ========================
!
! Jacobi(...) Example Program Results
!
! Eigenvalues
!    8.865457   -6.228747    0.864028   16.094837  -11.065575
!
! Eigenvectors
!    0.374481    0.607513    0.402620    0.489637   -0.298067
!    0.357169    0.287968   -0.406586   -0.605255   -0.507798
!   -0.500764    0.384320   -0.659965    0.399148   -0.081606
!   -0.620365    0.446730    0.455290   -0.456375   -0.003589
!   -0.310768   -0.448032    0.172458    0.162248   -0.804130
!======================================================================

program read_mat

        implicit none
        
        integer, parameter       :: m = 5, n = 5   ! the file where we read from has 5x5 values.
        real(8), parameter       :: abserr = 1.0e-9
        real(8), dimension(m, n) :: A, x
        integer :: i, j, count

! opening the file for reading
        open (2, file = 'data.dat', status = 'old')
        do i = 1, m
                read(2,*) ( A( i, j ), j = 1, n )
        end do
        close(2)

        write (*,3)
! Format the values from the file like a matrix.
        do i = 1, m
                write(*,4) ( A( i, j ), j = 1, n )
        end do
! Printing out the initial values
3       format (/'Values from the data.dat file')
4       format ( (6f12.6) )

        call Jacobi(A, x, abserr, n, count)
        ! print solutions
        write (*,200) count
        write (*,201) 
        write (*,202) (a(i,i),i=1,n)
        write (*,203)
        do i = 1,n 
                write (*,204) (x(i,j),j=1,n) 
        end do
        write (*, 205)

200     format (/'Eigenvalues and eigenvectors (Jacobi method) ',/,&
                'Jacobi method took',i2,' loops.')
201     format (/,' Eigenvalues')
202     format (6f12.6) 
203     format (/,' Eigenvectors')
204     format (6f12.6)
205     format ('NB! Eigenvalues and vectors are correct, but poorly organised.',/,& 
                'Hence the orthonormality thing.')

end program read_mat

subroutine Jacobi(a,x,abserr,n,count) 
!=========================================================== 
! Evaluate eigenvalues and eigenvectors 
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009) 
!----------------------------------------------------------- 
! input ... 
! a(n,n) - array of coefficients for matrix A 
! n - number of equations 
! abserr - abs tolerance [sum of (off-diagonal elements)^2] 
! output ... 
! a(i,i) - eigenvalues 
! x(i,j) - eigenvectors 
! comments ... 
!=========================================================== 
        implicit none 
        integer i, j, k, n, count
        real(8) a(n,n),x(n,n) 
        real(8) abserr, b2, bar 
        real(8) beta, coeff, c, s, cs, sc 
        
! initialize x(i,j)=0, x(i,i)=1 
! *** the array operation x=0.0 is specific for Fortran 90/95 
        x = 0.0 
        do i = 1, n 
                x(i,i) = 1.0 
        end do 

! find the sum of all off-diagonal elements (squared) 
        b2 = 0.0 
        do i = 1, n 
                do j = 1, n 
                        if (i.ne.j) b2 = b2 + a(i,j)**2 
                end do 
        end do 
        
        if (b2 <= abserr) return 
        
! average for off-diagonal elements /2 
        bar = 0.5*b2/float(n*n) 
        
        count = 0
        do while (b2.gt.abserr) 
                count = count + 1
                do i = 1, n-1 
                        do j=i+1,n 
                                if (a(j,i)**2 <= bar) cycle ! do not touch small elements
                                b2 = b2 - 2.0*a(j,i)**2 
                                bar = 0.5*b2/float(n*n) 
                                
                        ! calculate coefficient c and s for Givens matrix 
                                beta = (a(j,j)-a(i,i))/(2.0*a(j,i)) 
                                coeff = 0.5*beta/sqrt(1.0+beta**2) 
                                s = sqrt(max(0.5+coeff,0.0)) 
                                c = sqrt(max(0.5-coeff,0.0)) 
                        ! recalculate rows i and j 
                                do k = 1, n 
                                        cs = c*a(i,k)+s*a(j,k) 
                                        sc = -s*a(i,k)+c*a(j,k) 
                                        a(i,k) = cs 
                                        a(j,k) = sc 
                                end do 
                        ! new matrix a_{k+1} from a_{k}, and eigenvectors 
                                do k = 1, n 
                                        cs = c*a(k,i)+s*a(k,j) 
                                        sc = -s*a(k,i)+c*a(k,j) 
                                        a(k,i) = cs 
                                        a(k,j) = sc 
                                        cs = c*x(k,i)+s*x(k,j) 
                                        sc = -s*x(k,i)+c*x(k,j) 
                                        x(k,i) = cs 
                                        x(k,j) = sc 
                                end do 
                        end do 
                end do 
        end do 
        return 
end



