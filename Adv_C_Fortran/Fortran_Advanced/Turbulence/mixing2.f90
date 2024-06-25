!===============================================================
! Compile and execute with
!   > gfortran mixing2.f90 -o mix2
!   > ./mix2
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
    call iterate(numCycles, c2, c1)

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
    open(unit=1, file='c1_values.csv', status='replace')
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
        subroutine iterate(no_of_times, matrix, result)
            integer, intent(in) :: no_of_times
            real(8), dimension(:, :), intent(inout) :: matrix
            real(8), dimension(:, :), intent(out)   :: result
            real(8), dimension(:, :), allocatable   :: temp
            !integer :: k, i, j
            !real(8) :: a = 200.0d0
            integer :: rows, cols

            ! Get the dimensions of the input matrix
            rows = size(matrix, dim=1)
            cols = size(matrix, dim=2)

            ! Allocate memory for temporary matrix
            allocate(temp(rows, cols))

            ! Iterate for the given number of times
            do k = 0, no_of_times - 1
                ! Perform mixing
                do i = 1, rows
                    do j = 1, cols
                        if (mod(k, 2) == 0) then
                            temp(i, j) = matrix(i, mod(j + int(a * cos(i * pi * 2.0d0 / cols)), rows) + 1)
                        else
                            temp(i, j) = matrix(mod(i + int(a * cos(j * pi * 2.0d0 / cols)), rows) + 1, j)
                        end if
                    end do
                end do
                ! Perform diffusion
                call p_diffuse(temp, matrix, mod(k, 2) == 0)
            end do

            ! Copy the final result
            result = matrix

            ! Deallocate c_init
            deallocate(temp)
        end subroutine iterate

        subroutine p_diffuse(input_matrix, output_matrix, in_x_direction)
            real(8), dimension(:, :), intent(in)  :: input_matrix
            real(8), dimension(:, :), intent(out) :: output_matrix
            logical, intent(in)                   :: in_x_direction
            real(8), dimension(:, :), allocatable :: c_init
            !real(8) :: a = 200.0d0
            real(8) :: w1, w2
            integer :: rows, cols

            ! Get the dimensions of the input matrix
            rows = size(input_matrix, dim=1)
            cols = size(input_matrix, dim=2)
        
            ! Allocate memory for c_init
            allocate(c_init(rows, cols))

            ! Initialize the sine matrix
            do i = 1, rows
                do j = 1, cols
                    c_init(i, j) = sin(j*2.0d0*pi/dble(rows))
                end do
            end do

            ! Loop through each element in the matrix
            do i = 1, rows
                do j = 1, cols
                    if (in_x_direction) then
                        ! If diffusion is in the x-direction
                        w1 = a * cos(j * 2.0d0 * pi / dble(cols)) - floor(a * cos(j * 2.0d0 * pi / dble(cols)))
                        w2 = 1.0d0 - w1

                        if (j == 1) then
                            output_matrix(i, j) = input_matrix(i, cols) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                        else
                            output_matrix(i, j) = input_matrix(i, j-1) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                        end if
                    else
                        ! If diffusion is in the y-direction
                        w1 = a * cos(i * 2.0d0 * pi / dble(cols)) - floor(a * cos(i * 2.0d0 * pi / dble(rows)))
                        w2 = 1.0d0 - w1

                        if (i == 1) then
                            output_matrix(i, j) = input_matrix(rows, j) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                        else
                            output_matrix(i, j) = input_matrix(i-1, j) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                        end if
                    end if
                end do
            end do

            ! Deallocate c_init
            deallocate(c_init)

        end subroutine p_diffuse

end program main
