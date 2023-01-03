program solve_system
  !use lapack
  implicit none

  integer :: N, NRHS, LDA, LDB, info
  double precision, allocatable :: A(:,:), B(:,:), X(:,:)

  ! Set the size of the system
  N = 3
  NRHS = 1
  LDA = N
  LDB = NRHS

  ! Allocate memory for the matrices
  allocate(A(LDA,N), B(LDB,NRHS), X(LDB,N))

  ! Set the values of the matrices
  A = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9], shape(A))
  B = reshape([1, 0, 1], shape(B))

  ! Solve the system of equations
  call dgesv(N, NRHS, A, LDA, B, LDB, X, LDB, info)

  ! Print the solution
  print *, X

  ! Deallocate memory
  deallocate(A, B, X)
end program solve_system
