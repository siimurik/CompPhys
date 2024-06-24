program mixD
    implicit none
    integer :: num_cycles, i, j
    real(8), dimension(:,:), allocatable :: c1
  
    ! Prompt the user for the number of mixing cycles
    print *, "Input the number of mixing cycles:"
    read *, num_cycles
  
    ! Initialize the matrix (adjust size as needed)
    allocate(c1(5,5))  ! Change dimensions to match your actual matrix size
    call initialize_matrix(c1)
  
    ! Print the initial state of the matrix
    print *, "Initial state of matrix c1:"
    call print_matrix(c1)
  
    ! Perform the mixing cycles
    do i = 1, num_cycles
      call p_diffuse(c1)
      print *, "State of matrix c1 after iteration ", i, ":"
      call print_matrix(c1)
    end do
  
    ! Deallocate the matrix
    deallocate(c1)
  
    print *, "Computations completed."
    
  contains
  
    subroutine initialize_matrix(matrix)
      real(8), dimension(:,:), intent(out) :: matrix
      !integer :: i, j
  
      ! Initialize matrix with some values (adjust as needed)
      do i = 1, size(matrix, 1)
        do j = 1, size(matrix, 2)
          matrix(i, j) = real(i + j, 8)
        end do
      end do
    end subroutine initialize_matrix
  
    subroutine print_matrix(matrix)
      real(8), dimension(:,:), intent(in) :: matrix
      !integer :: i, j
  
      do i = 1, size(matrix, 1)
        do j = 1, size(matrix, 2)
          write(*,'(F10.6)', advance='no') matrix(i, j)
          write(*,'(A)', advance='no') ' '
        end do
        print *
      end do
    end subroutine print_matrix
  
    subroutine p_diffuse(matrix)
      real(8), dimension(:,:), intent(inout) :: matrix
      real(8), dimension(size(matrix, 1), size(matrix, 2)) :: temp_matrix
      !integer :: i, j
  
      ! Example diffusion logic (replace with actual logic)
      temp_matrix = matrix
      do i = 2, size(matrix, 1) - 1
        do j = 2, size(matrix, 2) - 1
          matrix(i, j) = (temp_matrix(i-1, j) + temp_matrix(i+1, j) + &
                          temp_matrix(i, j-1) + temp_matrix(i, j+1)) / 4.0
        end do
      end do
  
      ! Handle boundaries (example logic, adjust as needed)
      do i = 1, size(matrix, 1)
        matrix(i, 1) = matrix(i, 2)
        matrix(i, size(matrix, 2)) = matrix(i, size(matrix, 2) - 1)
      end do
  
      do j = 1, size(matrix, 2)
        matrix(1, j) = matrix(2, j)
        matrix(size(matrix, 1), j) = matrix(size(matrix, 1) - 1, j)
      end do
    end subroutine p_diffuse
  
end program mixD
  