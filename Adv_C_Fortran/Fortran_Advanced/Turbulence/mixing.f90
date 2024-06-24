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
    !allocate(c2(N, N))

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
    !deallocate(c2)

    ! Include necessary subroutines
    contains

        subroutine iterate(no_of_times, matrix, dim)
            implicit none
            integer, intent(in) :: no_of_times, dim
            real(8), dimension(:, :), intent(inout) :: matrix
            real(8), dimension(:, :), allocatable :: c
            !integer :: k, i, j
            !real(8) :: a, pi
            !parameter (a = 200.0d0)
            !parameter (pi = 3.141592653589793d0)
        
            ! Allocate memory for temporary matrix
            allocate(c(dim, dim))
        
            ! Iterate over the number of times specified
            do k = 1, no_of_times
                ! Iterate over the matrix dimensions
                do i = 1, dim
                    do j = 1, dim
                        if (MOD(k, 2) == 0) then
                            c(i, j) = matrix(i, MOD(j - 1 + INT(a * COS(i * 2.0d0 * pi / dim)), dim) + 1)
                        else
                            c(i, j) = matrix(MOD(i - 1 + INT(a * COS(j * 2.0d0 * pi / dim)), dim) + 1, j)
                        end if
                    end do
                end do
                call p_diffuse(c, MOD(k, 2) == 0)
                do i = 1, dim
                    do j = 1, dim
                        matrix(i, j) = c(i, j)
                    end do
                end do
            end do
        
            ! Deallocate the temporary matrix
            deallocate(c)
        end subroutine iterate
        
        subroutine p_diffuse(matrix, in_x_direction)
            implicit none
            logical, intent(in) :: in_x_direction
            real(8), dimension(:, :), intent(inout) :: matrix
            integer :: rows, cols!, i, j 
            real(8), dimension(:, :), allocatable :: new_matrix, c_init
            real(8) :: w1, w2 !, a, pi
            !parameter (a = 200.0d0)
            !parameter (pi = 3.141592653589793d0)
        
            rows = size(matrix, 1)
            cols = size(matrix, 2)
        
            ! Allocate memory for temporary matrices
            allocate(new_matrix(rows, cols))
            allocate(c_init(rows, cols))
        
            ! Initialize the c_init matrix
            do i = 1, rows
                do j = 1, cols
                    c_init(i, j) = sin(j * pi * 2.0d0 / real(cols))
                end do
            end do
        
            ! Loop through each element in the matrix
            do i = 1, rows
                do j = 1, cols
                    if (in_x_direction) then
                        w1 = a * cos(j * 2.0d0 * pi / real(cols)) - floor(a * cos(j * 2.0d0 * pi / real(cols)))
                        w2 = 1.0d0 - w1
        
                        if (j == 1) then
                            new_matrix(i, j) = matrix(i, cols) * w1 + matrix(i, j) * w2 + c_init(i, j)
                        else
                            new_matrix(i, j) = matrix(i, j - 1) * w1 + matrix(i, j) * w2 + c_init(i, j)
                        end if
                    else
                        w1 = a * cos(i * 2.0d0 * pi / real(rows)) - floor(a * cos(i * 2.0d0 * pi / real(rows)))
                        w2 = 1.0d0 - w1
        
                        if (i == 1) then
                            new_matrix(i, j) = matrix(rows, j) * w1 + matrix(i, j) * w2 + c_init(i, j)
                        else
                            new_matrix(i, j) = matrix(i - 1, j) * w1 + matrix(i, j) * w2 + c_init(i, j)
                        end if
                    end if
                end do
            end do
        
            ! Update the matrix with the new values
            matrix = new_matrix
        
            ! Deallocate temporary matrices
            deallocate(new_matrix)
            deallocate(c_init)
        end subroutine p_diffuse
    

end program main
