!  =============================================================================
!   $ gfortran dgtsv_solver.f90 -o tri -llapack
!   $ ./tri
! In GNU Plot:
!   $ plot 'out_2.csv' with lines
!  =============================================================================
program dgstv_solver
    implicit none
    integer          N, NRHS
    parameter        ( N = 5001, NRHS = 1 )
    integer          LDA, LDB, i
    parameter        ( LDA = N, LDB = N )
    integer          INFO
    double precision DL( N-1 ), DU( N-1 )
    double precision D( N )
    double precision A( LDA, N ), B( LDB, NRHS )
    double precision h
    double precision X(N)
    double precision start_time, end_time, elapsed_time

    ! Initializing 
    h = 2.D0/float(N-1)
    do i = 1, N
        X(i) = (i-1)*h
    end do
    DL = 0.D0; D = 0.D0; DU = 0.D0;
    B = 0.D0

    ! Boundary conditions
    D(1) =  1.D0/h + h/2.D0*(4.D0-X(1))
    B(1,1)   = h/2.D0*(X(1)+5.D0) - 1.D0

    ! Elements of the main diagonal D and 
    ! righthand side vector B
    do i = 2, N-1
        D(i) =  2.D0/h + h*(4.D0-X(i))
        B(i,1) = h*(X(i)+5.D0)
    end do

    ! Elements of off-diagonal elements
    do i = 1, N-1
        DL(i) = -1.D0/h
        DU(i) = -1.D0/h
    end do

    ! Final values
    D(N)    =  1.D0/h + h/2.D0*(4.D0-X(N))
    B(N,1)  =  h/2.D0*(X(N)+5.D0) - 1.D0

    !do i = 1, 3
    !    write(*,'(6f12.5)') D(i)
    !end do
    !write(*,*) '            ...'
    !do i = N-2, N
    !    write(*,'(6f12.5)') D(i)
    !end do

    ! record the starting time
    call cpu_time(start_time)

    ! Solve the equations A*X = B.
    call dgtsv( N, NRHS, DL, D, DU, B, LDB, INFO )

    ! record the ending time
    call cpu_time(end_time)

    ! checking for problems
    if ( INFO.GT.0 ) then
        write(*,*)'The diagonal element of the triangular factor of A,'
        write(*,*)'U(',INFO,',',INFO,') is zero, so that'
        write(*,*)'A is singular; the solution could not be computed.'
        stop
    end if

    ! printing out the first and last 3 elements
    write(*,*)'DGTSV Program Results'
    do i = 1, 3
        write(*,'(6f12.5)') X(i), B(i,1)
    end do
    write(*,*) '            ...'
    do i = N-2, N
        write(*,'(6f12.5)') X(i), B(i,1)
    end do

    ! Print the solution into a separate file
    open(unit = 11, file="out_2.csv")
    do i = 1, N
        write(11,*) X(i), B(i,1)
    end do
    close(11)

    ! compute the elapsed time
    elapsed_time = end_time - start_time

    ! print the elapsed time
    write(*,9) elapsed_time
9   format ( 'Elapsed time: ', F8.6, ' seconds.')

end program
