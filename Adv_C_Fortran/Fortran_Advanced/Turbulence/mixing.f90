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
    real(8), dimension(:, :), allocatable :: c1
    integer :: i, j, k
    integer :: numCycles
    integer :: start_time, end_time, elapsed_time, rate
    real(8) :: elapsed_seconds

    ! Allocate the size of original matrix c1
    allocate(c1(N, N))

    ! Initialize c1 with sine values
    do i = 1, N
        do j = 1, N
            c1(i, j) = sin(j*2.0d0*pi / dble(N))
        end do
    end do

    ! Get the number of mixing cycles
    write(*,*) "Input the number of mixing cycles:"
    read(*,*) numCycles

    ! Get start time
    call system_clock(count=start_time, count_rate=RATE)

    ! Perform the iterations
    call iterate(numCycles, c1, N)

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

    ! Include necessary subroutines
    contains

        subroutine iterate(no_of_times, matrix, dim)
            implicit none
            integer, intent(in) :: no_of_times, dim
            !logical :: in_x_direction
            real(8), dimension(:,:) :: matrix
            real(8), dimension(size(matrix, dim=1), size(matrix, dim=2)) :: c
        
            a = 200.0d0
            !allocate(c(N, N))
            c = 0.0d0
        
            do k = 1, no_of_times
                do i = 1, dim
                    do j = 1, dim
                        if (MOD(k, 2) == 0) then
                            c(i, j) = matrix(i, MOD(j + INT(a * COS(i * 2.0d0 * pi / dim)), dim) + 1)
                        else
                            c(i, j) = matrix(MOD(i + INT(a * COS(j * 2.0d0 * pi / dim)), dim) + 1, j)
                        end if
                    end do
                end do
                ! Assuming p_diffuse is a subroutine that modifies 'c'
                call p_diffuse(c, MOD(k, 2) == 0, matrix)
                matrix = matrix + c
            end do
        
        end subroutine iterate

        subroutine p_diffuse(matrix, in_x_direction, new_matrix)
            real(8), dimension(:,:) :: matrix, new_matrix
            real(8), dimension(size(matrix, dim=1), size(matrix, dim=2) ) :: c_init
            logical, intent(in) :: in_x_direction
            !integer :: i, j
            !real(8), parameter :: pi = 3.14159265358979323846
            real(8) :: w1, w2
            integer :: rows, cols
        
            ! Get the dimensions of the input matrix
            rows = size(matrix, dim=1)
            cols = size(matrix, dim=2)
        
            ! Initialize new matrices with zeros
            new_matrix = 0.0D0
            c_init = 0.0D0
        
            ! Set the diffusion coefficient
            a = 200.0D0
        
            ! Loop through each element in the matrix
            do i = 1, rows
                do j = 1, cols
                    if (in_x_direction) then
                        ! If diffusion is in the x-direction
                        w1 = a * cos(j * 2.0D0 * pi / cols) - dble(floor(a * cos(j * 2.0D0 * pi / cols)))
                        w2 = 1.0D0 - w1
        
                        ! Handle the special case when j is 1
                        if (j == 1) then
                            ! Wrap around to the last column for the diffusion calculation
                            new_matrix(i,j) = matrix(i,cols) * w1 + matrix(i,j) * w2 + c_init(i,j)
                        else
                            ! Standard diffusion calculation
                            new_matrix(i,j) = matrix(i,j-1) * w1 + matrix(i,j) * w2 + c_init(i,j)
                        end if
                    else
                        ! If diffusion is in the y-direction
                        w1 = a * cos(i * 2.0D0 * pi / rows) - dble(floor(a * cos(i * 2.0D0 * pi / rows)))
                        w2 = 1.0D0 - w1
        
                        ! Handle the special case when i is 1
                        if (i == 1) then
                            ! Wrap around to the last row for the diffusion calculation
                            new_matrix(i,j) = matrix(rows,j) * w1 + matrix(i,j) * w2 + c_init(i,j)
                        else
                            ! Standard diffusion calculation
                            new_matrix(i,j) = matrix(i-1,j) * w1 + matrix(i,j) * w2 + c_init(i,j)
                        end if
                    end if
                end do
            end do
        
        end subroutine p_diffuse

end program main
