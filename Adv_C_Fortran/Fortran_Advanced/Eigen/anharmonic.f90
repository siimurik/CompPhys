!========================================================
! Compile and execute with:
!   $ gfortran -O2 anharmonic.f90 -o an -llapack -lblas
!   $ ./an
!========================================================
program anharmonic_elevels
!========================================================
    implicit none
    integer, parameter                 :: P     = 1000
    integer, parameter                 :: LWORK = 3*P-1
    integer                            :: DIM, iostat
    double precision, dimension(P,P)   :: H, X, X4 ! Hamiltonian+Position Ops
    double precision, dimension(P)     :: E        ! energy eigenvalues
    double precision, dimension(LWORK) :: WORK
    double precision                   :: lambda
    integer                            :: i

    H = 0.0D0; X = 0.0D0; X4 = 0.0D0
    E = 0.0D0; WORK = 0.0D0

    do
        ! Prompt the user for input
        write(*, *) '# Enter Hilbert Space dimension:'

        ! Try to read the input as an integer
        read(*, *, iostat=iostat) DIM

        ! Check if the read was successful
        if (iostat .eq. 0) then
            write(*, *) '# You entered:', DIM
            exit   ! Exit the loop if input is successful
        else
            write(*, *) '# Warning: Input is not an integer. Please try again.'
        end if
    end do

    print *, '# Enter lambda:'
    read *, lambda
    print *, '# lambda = ', lambda
    
    ! Print Message:
    print *, '# ################################################'
    print *, '# Energy spectrum of anharmonic oscillator'
    print *, '# using matrix methods.'
    print *, '# Hilbert Space Dimension DIM = ', DIM
    print *, '# lambda coupling = ', lambda
    print *, '# ################################################'
    print *, '# Output: DIM   lambda   E_0   E_1   ....  E_{N-1}'
    print *, '# ------------------------------------------------'

    ! Calculate X^4 operator:
    call calculate_X4(X, X4, DIM)
    ! Calculate eigenvalues:
    call calculate_evs(H, X4, E, WORK, lambda, DIM)
    write(6,100) 'EV ', DIM, lambda, (E(i),i=1, DIM)
100 FORMAT(A3, I8, 20000G25.15) 

end program anharmonic_elevels

!========================================================
subroutine calculate_evs(H, X4, E, WORK, lambda, DIM)
!========================================================
    implicit none
    integer, parameter                  :: P     = 1000
    integer, parameter                  :: LWORK = 3*P-1
    double precision, dimension (P, P)  :: H, X4
    double precision, dimension (P)     :: E
    double precision, dimension (LWORK) :: WORK
    integer             :: DIM
    double precision    :: lambda
    character ( 1 )     :: JOBZ, UPLO
    integer             :: INFO, i, j
    call calculate_H(H, X4, lambda, DIM)
    JOBZ='V';   UPLO='U'
    call DSYEV( JOBZ, UPLO, DIM, H, P, E, WORK, LWORK, INFO)
    print *, '#**********************EVEC*******************'
    do j = 1, DIM
        write (6, 101) '# EVEC ', lambda, (H(i,j), i=1, DIM)
    enddo
    print *, '#**********************EVEC*******************'

101 FORMAT(A7, F8.3, 20000G20.6)

    ! If INFO is nonzero then we have an error
    if (INFO .ne. 0) then
        print *, ' DSYEV failed. INFO = ', INFO
        stop
    end if
end subroutine calculate_evs

!========================================================
subroutine calculate_H(H, X4, lambda, DIM)
!========================================================
    implicit none
    integer, parameter               :: P = 1000
    double precision, dimension(P,P) :: H, X4
    integer                          :: DIM
    double precision                 :: lambda
    integer                          :: i, j
    do j = 1, DIM
        do i = 1, DIM
            H(i,j)=lambda * X4(i,j)
        enddo
        H(j,j) = H(j,j) + DBLE(j) - 0.5D0 !E_n=n+1/ 2 ,n=j -1=>E_n=j -1/2
    enddo
    
    print *, '#**********************H**********************'
    do j = 1, DIM
        write (6 ,102) '# HH ', (H(i,j), i=1, DIM)
    enddo
    print *, '#**********************H**********************'

102 FORMAT(A5, 20000G16.6)
end subroutine calculate_H

!========================================================
subroutine calculate_X4(X, X4, DIM)
!========================================================
    implicit none
    integer, parameter               :: P = 1000
    double precision, dimension(P,P) :: X, X4
    double precision, allocatable, dimension(:,:) :: X2
    integer                          :: DIM
    integer                          :: i, j, m, n
    double precision, parameter      :: isqrt2 = 1.0D0/sqrt(2.0D0)
    
    ! Compute the position operator:
    X = 0.0D0
    
    ! Compute the nonzero elemen ts
    do i = 1, DIM
        n = i-1 ! indices0, ..., DIM-1
        ! The delta_{n ,m+1} term , i.e. m = n-1
        m = n-1 ! the energy level n -> i=n+1 , m-> j=m+1
        j = m+1
        if (j .ge. 1)   X(i,j) = isqrt2*sqrt(DBLE(m+1))
        ! The delta_{n, m-1} term, i.e. m = n+1
        m = n+1
        j = m+1
        if (j .le. DIM) X(i,j) = isqrt2*sqrt(DBLE(m))
    enddo

    ! Compute the Hamiltonian operator:
    ! Start with the X^4 operator:
    allocate(X2(P,P))
    X2 = MATMUL(X , X ) ! first X2, then X4:
    X4 = MATMUL(X2, X2)
end subroutine calculate_X4