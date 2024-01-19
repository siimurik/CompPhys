!==================================================
! Compile and execute with:
!   $ gfortran test_evs.f90 -o test -llapack -lblas
!   $ ./test
!--------------------------------------------------
! DSYEV documentation:
!   $ man dsyev
    !--------------------------------------------------
! Dowload LAPACK with:
!   $ sudo apt-get install liblapack3 liblapack-doc 
!     liblapack-dev.
!==================================================
program test_evs
    implicit none
    integer, parameter :: P = 100 ! P = LDA
    integer, parameter :: LWORK = 3*P-1
    double precision :: A(P,P), W(P), WORK(LWORK)
    integer :: N ! DSYEV diagonalizes A(N,N)
    integer :: i, j
    integer :: LDA, INFO
    character(1) :: JOBZ, UPLO
    character(len=15) :: mformat='(100(F14.6,1x))'
    
    ! Initialize values     
    N = 4       ! Matrix dimension
    A = 0.0     
    W = 0.0
    WORK = 0.0    

    ! Define the **symmetric** matrix to be diagonalized
    ! The subroutine uses the upper triangularpart (UPLO= 'U')
    ! therefore the lower triangular part needs not to be defined
    A(1,1) = -7.7;      A(1,2) =  2.1;      A(1,3) = -3.7;      A(1,4) =  4.4;
                        A(2,2) =  8.3;      A(2,3) = -16.0;     A(2,4) =  4.6;
                                            A(3,3) = -12.0;     A(3,4) = -1.04;
                                                                A(4,4) = -3.7;

    ! We print marix A before calling DSYEV since it is 
    ! destroyed adter the call.

    !do i = 1, N
    !    do j = 1, N
    !        print *, 'A( ',i,' , ', j,' )= ', A(i,j)
    !    enddo
    !enddo
    write(*,*) 'Matrix A values:'
    do i = 1, N
        write(*,mformat) (A(i,j), j = 1, N)
    end do

    ! We ask for eigenvalues AND eigenvectors (JOBZ='V')
    JOBZ = 'V'; UPLO = 'U'
    print *, 'COMPUTING WITH DSYEV:'
    LDA = P     ! notice that LDA -> P>N !!

    ! DSYEV computes the eigenvalues and eigenvectors    
    call DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
    
    print *, 'DSYEV: DONE. CHECKING NOW:'

    ! If INFO is nonzero, then there is an error:
    if (INFO .ne. 0) then
        print *, 'DSYEV FAILED. INFO = ', INFO
        stop
    end if

    ! Print results: W(I) has the eigenvalues:
    print *, 'DSYEV: DONE.:'
    print *, 'EIGENVALUES OF MATRIX:'
    do i = 1, size(W)
        if (W(i) .ne. 0.0) then 
            !print *,  'LAMBDA(',i,')=', W(i)
            write (*,'(a,i3,a,F21.16)') 'LAMBDA(', i, ') =', W(i)
        end if
    end do
    
    ! Eigenvectors are stored in coluns of A:
    print *, 'EIGENVECTORS OF MATRIX:'
    do j = 1, N
        !print *, 'EIGENVECTOR ', j, 'FOR EIGEVALUE', W(j)
        write (*,'(a,i2,a,F21.16)') 'EIGENVECTOR', j, ' FOR EIGENVALUE', W(j)
        do i = 1, N
            !print *, 'V_',j,'(',i,')=', A(i,j)
            write (*, '(a,i1,a,i1,a,F20.16)') 'V_',j,'(',i,') =', A(i,j)
        end do
    end do

end program test_evs
