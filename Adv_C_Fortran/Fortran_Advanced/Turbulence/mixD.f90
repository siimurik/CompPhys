program test_diffuse
    implicit none
    real(8), dimension(:,:), allocatable :: cin, cout
    real(8), dimension(:,:), allocatable :: temp
    real(8), dimension(:), allocatable :: cos_values
    real(8), parameter :: pi = 3.14159265358979323846d0
    real(8) :: a = 200.0d0
    integer :: i, j, k, N, index
    !logical :: direc

    N = 5
    allocate(cin(N, N))
    allocate(cout(N, N))

    ! Initialize input matrix
    do i = 1, N
        do j = 1, N
            cin(i, j) = dble((i-1) + (j-1))  ! Corrected array initialization
        end do
    end do

    ! Print the input matrix
    print *, 'Input Matrix:'
    call print_matrix(cin)

    !direc = .true.

    ! Initialize temporary matrix
    allocate(temp(N,N))
    allocate(cos_values(N))

    ! Precompute cosine values
    do i = 1, N
        cos_values(i) = cos((i - 1) * pi * 2.0 / N)
    end do

    ! Iterate for the given number of times
    do k = 0, 2 - 1
        ! Perform mixing
        do i = 1, N
            do j = 1, N
                if (mod(k, 2) == 0) then
                    index = iand(j - 1 + int(a * cos_values(i)), N - 1)
                    temp(i, j) = cin(i, index + 1)
                    print *, index
                else
                    index = iand(i - 1 + int(a * cos_values(j)), N - 1)
                    temp(i, j) = cin(index + 1, j)
                end if
            end do
        end do
        ! Perform diffusion
        call p_diffuse(temp, cin, mod(k, 2) == 0)
    end do

    ! Copy the final result
    cout = cin

    ! Print the output matrix
    print *, 'Output Matrix:'
    call print_matrix(cout)

    deallocate(cin)
    deallocate(cout)
    deallocate(temp)

contains

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
  
    N = size(matrix, 1)  ! Assuming the matrix is square (N x N)
  
    do i = 1, N
        do j = 1, N
            write(*, '(F10.5)', advance="no") matrix(i, j)
            if (j < N) then
                write(*, '(" ")', advance="no")
            end if
        end do
        print *
    end do
  end subroutine print_matrix

end program test_diffuse