program matrix_multiply

    ! include the LAPACK module
    !use lapack
    implicit none

    ! Define the two matrices that we want to multiply
    real :: A(1000,1000), B(1000,1000), C(1000,1000)
    integer :: i, j, n

    ! initialize the matrices with random values
    do i = 1, 1000
        do j = 1, 1000
            A(i,j) = (i-1) * 1000 + j
            B(i,j) = -((i-1) * 1000 + j)
        enddo
    enddo
    !call random_number(A)
    !call random_number(B)

    do i = 1, 1000
        do j = 1, 1000
            C(i,j) = 0.0
        enddo
    enddo

    ! Perform the matrix multiplication using the SGEMM function from LAPACK
    call sgemm('N', 'N', 1000, 1000, 1000, 1.0, A, 1000, B, 1000, 0.0, C, 1000)

    print *, "Top left corner of matrix A:"
    print 20, ((A(i,j), j = 1,min(1000,6)), i = 1,min(1000,6))
    print *, ""

    print *, "Top left corner of matrix B:"
    print 20, ((B(i,j),j = 1,min(1000,6)), i = 1,min(1000,6))
    print *, ""

20      format(6(F12.0,1x))

    print *, "Top left corner of matrix C:"
    print 30, ((C(i,j), j = 1,min(1000,6)), i = 1,min(1000,6))
    print *, ""

30      format(6(ES12.4,1x))

end program matrix_multiply
