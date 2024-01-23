module LinspaceModule

    contains

    function linspace(a, b, initArray) result(finArray)
        implicit none
        double precision, intent(in) :: a, b, initArray(:)
        double precision, allocatable :: finArray(:)
        integer :: num, i
        double precision :: step

        num = size(initArray)

        if  (num .le. 1) then
            print *, "Error: Size of array must be greater than 1."
            return
        end if

        step = (b-a)/dble(num-1.D0)

        allocate(finArray(num))
        do i = 1, num
            finArray(i) = a + (i-1)*step
            !print *, finArray
        end do

    end function linspace
end module LinspaceModule

program test_linspace
    implicit none
    use  LinspaceModule
    double precision :: a, b, c(10), x(10)
    !double precision, external :: linspace
    integer :: i

    call random_number(c)
    a = -2.D0
    b =  2.D0

    x = linspace(a, b, c)

    do i = 1, 10
        print *, x(i)
    end do

end program test_linspace