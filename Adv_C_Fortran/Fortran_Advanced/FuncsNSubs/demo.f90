program demo
    implicit none
    integer, parameter :: m = 5
    real :: a(m), b(1:m), c(0:m-1)
    
    a = [1,2,3,4,5]       ! Indices 1-5
    b = [1,2,3,4,5]       ! Indices 1-5 
    c = [1,2,3,4,5]       ! Indices 0-4
    
    print *, a(1), a(m)    ! 1.0, 5.0
    print *, b(1), b(m)    ! 1.0, 5.0 
    print *, c(0), c(m-1)  ! 1.0, 5.0   ! #MindBlown
end program