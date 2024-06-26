program benchmark
  implicit none
  integer, parameter :: n = 100000000
  real(8), dimension(n) :: array
  real(8) :: sum
  integer :: i

  ! Initialize the array with random numbers
  call random_seed()
  call random_number(array)

  sum = 0.0d0
  ! Sum the elements of the array
  do i = 1, n
    sum = sum + array(i)
  end do

  print *, 'Sum of array elements: ', sum
end program benchmark

