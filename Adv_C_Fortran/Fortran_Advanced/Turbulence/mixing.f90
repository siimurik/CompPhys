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
    real(8), dimension(:, :), allocatable :: c_in, c_out
    integer :: i, j, k
    integer :: numCycles
    integer :: start_time, end_time, elapsed_time, rate
    real(8) :: elapsed_seconds

    ! Allocate the size of original matrix c1
    allocate(c_in(N, N))
    allocate(c_out(N, N))

    ! Initialize c2 with sine values
    do i = 1, N
        do j = 1, N
            c_in(i, j) = sin((j-1)*2.0d0*pi / dble(N))
        end do
    end do

    ! Get the number of mixing cycles
    write(*,*) "Input the number of mixing cycles:"
    read(*,*) numCycles

    ! Get start time
    call system_clock(count=start_time, count_rate=rate)

    ! Perform the iterations
    call iterate(numCycles, c_in, c_out)

    ! Get end time
    call system_clock(count=end_time)

    ! Calculate elapsed time
    elapsed_time = end_time - start_time
    elapsed_seconds = dble(elapsed_time) / dble(rate)

    ! Print matrix
    print *, "Computations completed."
    print *, "Top left corner of the matrix:"
    call print_matrix(c_out)
    print *, ""
!    print 15, ((c_out(i,j), j = 1, MIN(N,5)), i = 1, MIN(N,5))
!    print *, ""
!15  format (6(ES12.4,1x))

    ! Print formatted calculation time
    write (*,16) elapsed_seconds
16  format('Elapsed time:', E10.3, ' seconds.')

    ! Export values into a custom CSV file
    call export_matrix_to_csv(c_out, 'c1_values.csv')
    !open(unit=1, file='c1_values.csv', status='replace')
    !do i = 1, N
    !    ! Assuming c1 is a 2D array with dimensions (N, M)
    !    do j = 1, N
    !        ! Add a comma after each value except the last one
    !        if (j /= N) then
    !            write(1, '(F20.10, A)', advance='no') c_out(i, j), ','
    !        else
    !            write(1, '(F20.10)') c_out(i, j)
    !        end if
    !    end do
    !end do
    !close(1)

    ! Free memory
    deallocate(c_in)
    deallocate(c_out)

contains
    subroutine iterate(no_of_times, matrix, result)
        integer, intent(in) :: no_of_times
        real(8), dimension(:, :), intent(inout) :: matrix
        real(8), dimension(:, :), intent(out)   :: result
        real(8), dimension(:, :), allocatable   :: temp
        real(8), dimension(:), allocatable :: cos_values
        integer :: index! k, i, j
        !real(8) :: a = 200.0d0
        integer :: rows, cols

        ! Get the dimensions of the input matrix
        rows = size(matrix, dim=1)
        cols = size(matrix, dim=2)

        ! Allocate memory for temporary matrix
        allocate(temp(rows, cols))
        allocate(cos_values(rows))

        ! Precompute cosine values
        do i = 1, N
            cos_values(i) = cos((i - 1) * pi * 2.0 / N)
        end do

        ! Iterate for the given number of times
        do k = 0, no_of_times - 1
            ! Perform mixing
            do i = 1, N
                do j = 1, N
                    if (mod(k, 2) == 0) then
                        index = iand(j - 1 + int(a * cos_values(i)), N - 1)
                        temp(i, j) = matrix(i, index + 1)
                        !print *, index
                    else
                        index = iand(i - 1 + int(a * cos_values(j)), N - 1)
                        temp(i, j) = matrix(index + 1, j)
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
        deallocate(cos_values)
    end subroutine iterate

    subroutine p_diffuse(input_matrix, output_matrix, in_x_direction)
        real(8), dimension(:, :), intent(in)  :: input_matrix
        real(8), dimension(:, :), intent(out) :: output_matrix
        logical, intent(in)                   :: in_x_direction
        real(8), dimension(:, :), allocatable :: c_init
        real(8) :: w1, w2
        integer :: rows, cols
    
        rows = size(input_matrix, dim=1)
        cols = size(input_matrix, dim=2)
        a = 200.0d0
    
        allocate(c_init(rows, cols))
    
        ! Initialize the sine matrix
        do i = 1, rows
            do j = 1, cols
                c_init(i, j) = sin(dble(j-1) * 2.0d0 * pi / dble(N))
            end do
        end do
    
        ! Loop through each element in the matrix
        do i = 1, rows
            do j = 1, cols
                if (in_x_direction) then
                    w1 = a * cos(dble(j-1) * 2.0d0 * pi / dble(cols)) - floor(a * cos(dble(j-1) * 2.0d0 * pi / dble(cols)))
                    w2 = 1.0d0 - w1
    
                    if (j == 1) then
                        output_matrix(i, j) = input_matrix(i, cols) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                    else
                        output_matrix(i, j) = input_matrix(i, j-1) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                    end if
                else
                    w1 = a * cos(dble(i-1) * 2.0d0 * pi / dble(rows)) - floor(a * cos(dble(i-1) * 2.0d0 * pi / dble(rows)))
                    w2 = 1.0d0 - w1
    
                    if (i == 1) then
                        output_matrix(i, j) = input_matrix(rows, j) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                    else
                        output_matrix(i, j) = input_matrix(i-1, j) * w1 + input_matrix(i, j) * w2 + c_init(i, j)
                    end if
                end if
            end do
        end do
    
        deallocate(c_init)
    end subroutine p_diffuse
    
    subroutine print_matrix(matrix)
        implicit none
        real(8), dimension(:,:) :: matrix
        !integer :: i, j, N
        !N = size(matrix, 1)  ! Assuming the matrix is square (N x N)
        
        do i = 1, 5
            do j = 1, 5
                write(*, '(F10.6)', advance="no") matrix(i, j)
                if (j < 5) then
                    write(*, '(" ")', advance="no")
                end if
            end do
            print *
        end do
    end subroutine print_matrix

    subroutine export_matrix_to_csv(print_mat, filename)
        implicit none
        real(8) :: print_mat(:,:)
        !integer, parameter :: dim = size(print_mat, 1) ! Assuming matrix is N x M
        character(len=*), intent(in) :: filename
        integer :: dim
        character(20) :: fmt

        dim = size(print_mat, 1)
    
        ! Open the file for writing
        open(unit=1, file=filename, status='replace')
    
        ! Set the format for writing values
        fmt = '(E20.10, A)'
    
        do i = 1, dim
            do j = 1, dim
                ! Add a comma after each value except the last one
                if (j /= dim) then
                    write(1, fmt, advance='no') print_mat(i, j), ','
                else
                    write(1, '(E20.10)') print_mat(i, j)
                end if
            end do
        end do
    
        ! Close the file
        close(1)
    end subroutine export_matrix_to_csv
    

end program main
