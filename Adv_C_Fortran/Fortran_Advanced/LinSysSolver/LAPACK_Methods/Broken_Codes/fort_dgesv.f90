!=========================================================
!   $ gfortran fort_dgesv.f90 -o fd -llapack -lblas
!---------------------------------------------------------
program main
    implicit none (type, external)
    external :: sgesv
    integer, parameter :: n = 200
    integer :: i, j
    double precision :: h, sum
    double precision, dimension(n+1,n+1) :: A
    double precision, dimension(n+1) :: y, x, u, b
    double precision :: pivot(n+1) ! Pivot indices (list of swap operations).
    integer :: rc ! Return code.
    double precision :: start_time, end_time, elapsed_time

    h = 2.0/n
    do i = 1, n+1
        x(i) = (i-1)*h
    end do
    A = 0.0
    y = 0.0

    A(1,1) = 1.0/h + h/2*(4.0-x(1))
    A(1,2) = -1.0/h
    y(1)   = h/2*(x(1)+5.0) - 1.0

    do i = 2, n
        A(i,i-1) = -1.0/h
        A(i,i)   =  2.0/h + h*(4.0-x(i))
        A(i,i+1) = -1.0/h
        y(i)     = h*(x(i)+5.0)
    end do
    A(n+1,n)   = -1/h
    A(n+1,n+1) = 1/h + h/2*(4-x(n+1))
    y(n+1)     = h/2*(x(n+1)+5) - 1

    do i = 1, n+1
        b(i) = y(i)
    enddo

    ! record the starting time
    call cpu_time(start_time)

    call sgesv(n+1, 1, A, n+1, pivot, b, n+1, rc)

    if (rc /= 0) then
        print '(a, i0)', 'Error: ', rc
        stop
    end if

    ! record the ending time
    call cpu_time(end_time)

    ! compute the elapsed time
    elapsed_time = end_time - start_time

    ! print the elapsed time
    write(*,'(A, F8.6)') "Elapsed time: ", elapsed_time, " seconds"

    ! plot the solution
    do i = 1, 5
        write(*,*) x(i), b(i)
    end do
end program main
