program benchmark
  use omp_lib
  implicit none
  integer, parameter :: n = 10000000
  real(8), allocatable :: a(:)
  real(8) :: sum
  integer :: i

  allocate(a(n))
  
  ! Initialize the array
  !$omp parallel do
  do i = 1, n
     a(i) = real(i, 8) * 0.5
  end do
  !$omp end parallel do

  ! Compute the sum of the array elements
  sum = 0.0
  !$omp parallel do reduction(+:sum)
  do i = 1, n
     sum = sum + a(i)
  end do
  !$omp end parallel do

  print *, 'Sum of array elements: ', sum

  deallocate(a)
end program benchmark

