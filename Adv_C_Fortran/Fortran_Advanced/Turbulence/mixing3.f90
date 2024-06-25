!===============================================================
! Compile and execute with
!   > gfortran mixing.f90 -o mix
!   > ./mix
!===============================================================
program main
    implicit none
    integer, parameter :: N = 512 ! Set your desired matrix size
    real(8), parameter :: pi = 3.14159265358979323846D0
    real(8) :: a = 200.0d0
    real(8), dimension(:, :), allocatable :: c1, c2
    integer :: i, j, k
    integer :: numCycles
    integer :: start_time, end_time, elapsed_time, rate
    real(8) :: elapsed_seconds

    ! Allocate the size of original matrix c1
    allocate(c1(N, N))
    allocate(c2(N, N))

    ! Initialize c2 with sine values
    do i = 1, N
        do j = 1, N
            c2(i, j) = sin(j*2.0d0*pi / dble(N))
        end do
    end do

    ! Get the number of mixing cycles
    write(*,*) "Input the number of mixing cycles:"
    read(*,*) numCycles

    ! Get start time
    call system_clock(count=start_time, count_rate=rate)

    ! Perform the iterations
    call iterate(numCycles, c2, c1, N)

    ! Get end time
    call system_clock(count=end_time)

    ! Calculate elapsed time
    elapsed_time = end_time - start_time
    elapsed_seconds = dble(elapsed_time) / dble(rate)

    ! Print matrix
    print *, "Computations completed."
    print *, ""
    print *, "Top left corner of matrix c1:"
    print 15, ((c1(i,j), j = 1, MIN(N,5)), i = 1, MIN(N,5))
    print *, ""
15  format (6(ES12.4,1x))

    ! Print formatted calculation time
    write (*,16) elapsed_seconds
16  format('Elapsed time:', E10.3, ' seconds.')

    ! Export values into a custom CSV file
    open(unit=1, file='output.csv', status='replace')
    do i = 1, N
        ! Assuming c1 is a 2D array with dimensions (N, M)
        do j = 1, N
            ! Add a comma after each value except the last one
            if (j /= N) then
                write(1, '(F20.10, A)', advance='no') c1(i, j), ','
            else
                write(1, '(F20.10)') c1(i, j)
            end if
        end do
    end do
    close(1)
    
    ! Free memory
    deallocate(c1)
    deallocate(c2)

contains
    subroutine iterate(no_of_times, matrix, result, M)
        integer, intent(in) :: no_of_times, M
        real(8), dimension(M, M), intent(inout) :: matrix
        real(8), dimension(M, M), intent(out) :: result
        !integer :: k, i, j
        !real(8) :: a = 200.0d0

        ! Initialize temporary matrix
        real(8), dimension(M, M) :: temp

        ! Iterate for the given number of times
        do k = 0, no_of_times - 1
            ! Perform mixing
            do i = 1, M
                do j = 1, M
                    if (mod(k, 2) == 0) then
                        temp(i, j) = matrix(i, mod(j + int(a * cos(i * pi * 2.0d0 / M)), M) + 1)
                    else
                        temp(i, j) = matrix(mod(i + int(a * cos(j * pi * 2.0d0 / M)), M) + 1, j)
                    end if
                end do
            end do
            ! Perform diffusion
            call p_diffuse(temp, matrix, mod(k, 2) == 0, M)
        end do

        ! Copy the final result
        result = matrix
    end subroutine iterate

    subroutine p_diffuse(input_matrix, output_matrix, in_x_direction, NN)
        real(8), dimension(NN, NN), intent(in) :: input_matrix
        real(8), dimension(NN, NN), intent(out) :: output_matrix
        logical, intent(in) :: in_x_direction
        integer, intent(in) :: NN
        !real(8) :: a = 200.0d0
        real(8) :: w1, w2
        !integer :: i, j

        ! Initialize the sine matrix
        real(8), dimension(NN, NN) :: c_init
        do i = 1, NN
            do j = 1, NN
                c_init(i, j) = sin(j * 2.0d0 * pi / dble(NN))
            end do
        end do

        ! Loop through each element in the matrix
        do i = 1, NN
            do j = 1, NN
                if (in_x_direction) then
                    ! If diffusion is in the x-direction
                    w1 = a * cos(dble(j-1.d0) * 2.0d0 * pi / dble(NN)) - floor(a * cos(dble(j-1.d0) * 2.0d0 * pi / dble(NN)))
                    w2 = 1.0d0 - w1

                    if (j == 1) then
                        output_matrix(i, j) = input_matrix(i, N) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                    else
                        output_matrix(i, j) = input_matrix(i, j-1) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                    end if
                else
                    ! If diffusion is in the y-direction
                    w1 = a * cos(dble(i-1.d0) * 2.0d0 * pi / dble(NN)) - floor(a * cos(dble(i-1.d0) * 2.0d0 * pi / dble(NN)))
                    w2 = 1.0d0 - w1

                    if (i == 1) then
                        output_matrix(i, j) = input_matrix(N, j) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                    else
                        output_matrix(i, j) = input_matrix(i-1, j) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                    end if
                end if
            end do
        end do
    end subroutine p_diffuse

end program main

!! Assuming c_in is a 2D array of size (N, N)
!integer, dimension(N) :: indices(N)
!real(kind=8), dimension(N, N) :: c_in
!real(kind=8), dimension(N*N) :: temp_array
!
!! Initialize an array with indices from 1 to N
!indices = [(i, i=1,N)]
!
!! Compute the sine values directly using vectorized operations
!temp_array = SIN((indices - 1) * 2.0d0 * pi / dble(N))
!
!! Reshape the temporary array to match the shape of c_in
!c_in = RESHAPE(temp_array, [N, N])
