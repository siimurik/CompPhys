! Compile and execute with:
! $ gfortran chem_radau5.f90 radau5.f radaua.f -o chem_radau5
! $ ./chem_radau5

! RADAU5 solver for stiff chemical kinetics problems
! Solves the Robertson chemical reaction problem - a classic stiff benchmark:
! A -> B (slow: k1 = 0.04)
! B + B -> B + C (fast: k2 = 3e7) 
! B + C -> A + C (very fast: k3 = 1e4)

program radau5_chem
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external fchem, jchem, solout_chem, f
    
    ! Problem parameters
    integer, parameter :: neq = 3
    real(kind=dp), dimension(neq) :: y
    real(kind=dp) :: x, xend, h, rtol, atol, rpar
    integer :: n, ijac, mljac, mujac, imas, mlmas, mumas, iout, itol, ipar, idid
    integer :: lwork, liwork
    real(kind=dp), allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork
    
    ! Calculate required array sizes and allocate
    call calculate_radau5_arrays(neq, lwork, liwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    ! Initialize problem parameters
    call setup_chemical_problem(neq, y, x, xend, h, rtol, atol, rpar, ipar, &
                               n, ijac, mljac, mujac, imas, mlmas, mumas, &
                               iout, itol, lwork, liwork, work, iwork)
    
    ! Run the integration
    call run_radau5_integration(fchem, jchem, solout_chem, n, y, x, xend, h, &
                               rtol, atol, itol, ijac, mljac, mujac, imas, &
                               mlmas, mumas, iout, work, lwork, iwork, liwork, &
                               rpar, ipar, idid)
    
    ! Print final statistics
    call print_radau5_statistics(n, x, y, rtol, liwork, iwork, idid)
    
    ! Clean up
    deallocate(work)
    deallocate(iwork)
    
end program radau5_chem

! Subroutine to calculate RADAU5 array lengths
subroutine calculate_radau5_arrays(neq, lwork, liwork)
    implicit none
    integer, intent(in) :: neq
    integer, intent(out) :: lwork, liwork
    
    ! For full Jacobian (most general case)
    lwork = 4 * neq * neq + 12 * neq + 20
    liwork = 3 * neq + 20
    
    write(6, '(A,I0,A,I0)') 'Calculated RADAU5 array sizes: LWORK = ', lwork, ', LIWORK = ', liwork
    
end subroutine calculate_radau5_arrays

! Setup the chemical kinetics problem
subroutine setup_chemical_problem(neq, y, x, xend, h, rtol, atol, rpar, ipar, &
                                 n, ijac, mljac, mujac, imas, mlmas, mumas, &
                                 iout, itol, lwork, liwork, work, iwork)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: neq, lwork, liwork
    real(kind=dp), dimension(neq), intent(out) :: y
    real(kind=dp), intent(out) :: x, xend, h, rtol, atol, rpar
    integer, intent(out) :: n, ijac, mljac, mujac, imas, mlmas, mumas, iout, itol, ipar
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer :: i
    
    write(6, '(A)') 'Setting up Robertson chemical kinetics problem...'
    write(6, '(A)') 'Reactions:'
    write(6, '(A)') '  A -> B         (k1 = 0.04)'
    write(6, '(A)') '  B + B -> B + C (k2 = 3e7)'
    write(6, '(A)') '  B + C -> A + C (k3 = 1e4)'
    write(6, '(A)') ''
    
    ! Problem dimension
    n = neq
    
    ! RADAU5 options
    ijac = 1      ! Compute Jacobian analytically
    mljac = n     ! Full Jacobian
    mujac = 0     ! Not used for full Jacobian
    imas = 0      ! Explicit ODE system
    mlmas = 0     ! Not used
    mumas = 0     ! Not used
    iout = 1      ! Use output routine
    
    ! Initial conditions: pure A, no B or C
    x = 0.0_dp                ! Initial time
    y(1) = 1.0_dp            ! [A] = 1.0 (concentration of species A)
    y(2) = 0.0_dp            ! [B] = 0.0 (concentration of species B)  
    y(3) = 0.0_dp            ! [C] = 0.0 (concentration of species C)
    
    ! Integration parameters
    xend = 4.0_dp            ! Final time (log10 scale: 10^4 = 10000 seconds)
    h = 1.0d-8              ! Initial step size (small for stiff problem)
    
    ! Tolerance settings for stiff problem
    rtol = 1.0d-8           ! Relative tolerance
    atol = 1.0d-12          ! Absolute tolerance (tight for chemical concentrations)
    itol = 0                ! Same tolerance for all components
    
    ! User parameters (not used in this problem)
    rpar = 0.0_dp
    ipar = 0
    
    ! Initialize work arrays
    do i = 1, 20
        iwork(i) = 0
        work(i) = 0.0_dp
    end do
    
    write(6, '(A)') 'Problem setup completed.'
    write(6, '(A,3F10.6)') 'Initial concentrations [A,B,C]: ', y(1), y(2), y(3)
    write(6, '(A,ES10.3)') 'Integration from t = 0 to t = ', xend
    write(6, '(A,ES10.3,A,ES10.3)') 'Tolerances: RTOL = ', rtol, ', ATOL = ', atol
    write(6, '(A)') ''
    
end subroutine setup_chemical_problem

! Main RADAU5 integration routine
subroutine run_radau5_integration(f, jac, solout, n, y, x, xend, h, rtol, atol, &
                                 itol, ijac, mljac, mujac, imas, mlmas, mumas, &
                                 iout, work, lwork, iwork, liwork, rpar, ipar, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external f, jac, solout
    integer, intent(in) :: n, itol, ijac, mljac, mujac, imas, mlmas, mumas, iout
    integer, intent(in) :: lwork, liwork, ipar
    real(kind=dp), intent(inout) :: x, h, rpar
    real(kind=dp), intent(in) :: xend, rtol, atol
    real(kind=dp), dimension(n), intent(inout) :: y
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer, intent(out) :: idid
    
    write(6, '(A)') 'Starting RADAU5 integration for stiff chemical kinetics...'
    write(6, '(A)') '     Time         [A]           [B]           [C]         NSTEP'
    write(6, '(A)') '---------------------------------------------------------------'
    
    ! Call RADAU5 solver
    call radau5(n, f, x, y, xend, h, &
                rtol, atol, itol, &
                jac, ijac, mljac, mujac, &
                f, imas, mlmas, mumas, &
                solout, iout, &
                work, lwork, iwork, liwork, rpar, ipar, idid)
    
    ! Check integration status
    select case (idid)
    case (1)
        write(6, '(/A)') 'Integration completed successfully.'
    case (2)
        write(6, '(/A)') 'Integration completed - interrupted by solout.'
    case (-1)
        write(6, '(/A)') 'Error: Input is not consistent.'
    case (-2)
        write(6, '(/A)') 'Error: Larger NMAX is needed.'
    case (-3)
        write(6, '(/A)') 'Error: Step size becomes too small.'
    case (-4)
        write(6, '(/A)') 'Error: Matrix is repeatedly singular.'
    case default
        write(6, '(/A,I0)') 'Integration finished with status: ', idid
    end select
    
end subroutine run_radau5_integration

! Print final statistics
subroutine print_radau5_statistics(n, x, y, rtol, liwork, iwork, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, liwork
    real(kind=dp), intent(in) :: x, rtol
    real(kind=dp), dimension(n), intent(in) :: y
    integer, dimension(liwork), intent(in) :: iwork
    integer, intent(in) :: idid
    real(kind=dp) :: conservation

    write(6, '(/A)') 'Final Integration Statistics:'
    write(6, '(A)') '============================'
    write(6, '(A,F10.6,A,3ES15.6)') 'Final time: ', x, ', Final solution [A,B,C]: ', y(1), y(2), y(3)
    
    ! Check mass conservation (A + B + C should = 1.0)
    conservation = y(1) + y(2) + y(3)
    write(6, '(A,ES15.6,A,ES10.3, A)') 'Mass conservation: A+B+C = ', conservation, &
                                    ' (error = ', abs(conservation - 1.0_dp), ')'
    
    write(6, '(A,ES10.3)') 'Relative tolerance used: ', rtol
    write(6, '(A,I0)') 'Integration status (IDID): ', idid
    write(6, '(/A)') 'Detailed statistics:'
    write(6, '(A,I4)') ' Function evaluations: ', iwork(14)
    write(6, '(A,I4)') ' Jacobian evaluations: ', iwork(15)
    write(6, '(A,I4)') ' Integration steps:    ', iwork(16)
    write(6, '(A,I4)') ' Accepted steps:       ', iwork(17)
    write(6, '(A,I4)') ' Rejected steps:       ', iwork(18)
    write(6, '(A,I4)') ' LU decompositions:    ', iwork(19)
    write(6, '(A,I4)') ' Forward-backward substitutions: ', iwork(20)
    
    if (iwork(16) > 0) then
        write(6, '(A,F7.2,A)') ' Average step size: ', x / real(iwork(16), dp), ''
        write(6, '(A,F7.2,A)') ' Acceptance ratio: ', &
                               100.0_dp * real(iwork(17), dp) / real(iwork(16), dp), '%'
    endif
    
end subroutine print_radau5_statistics

! Output routine for continuous solution monitoring
subroutine solout_chem(nr, xold, x, y, cont, dt, lrc, n, rpar, ipar, irtrn)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: nr, lrc, n, ipar
    real(kind=dp), intent(in) :: xold, x, rpar
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(lrc), intent(in) :: cont
    real(kind=dp), intent(inout) :: dt
    integer, intent(inout) :: irtrn
    real(kind=dp) :: contr5
    
    ! Storage for output control
    real(kind=dp), save :: xout = 0.0_dp
    logical, save :: first_call = .true.

    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) xold, cont(1), rpar, ipar, irtrn
    endif
    
    if (first_call) then
        ! Print initial condition
        write(6, '(1X,F10.6,3ES14.6,I8)') x, y(1), y(2), y(3), nr-1
        dt = 0.1_dp  ! Output interval
        xout = dt
        first_call = .false.
    else
        ! Print at regular intervals
        do while (x >= xout)
            write(6, '(1X,F10.6,3ES14.6,I8)') xout, &
                   contr5(1, xout, cont, lrc), &
                   contr5(2, xout, cont, lrc), &
                   contr5(3, xout, cont, lrc), nr-1
            xout = xout + dt
        end do

        ! This is a stupid fix, but fuck me, it won't work otherwise.
        write(6, '(1X,F10.6,3ES14.6,I8)') x, y(1), y(2), y(3), nr-1
        
        ! Check if we need to print the final time point (only if not already printed)
        if (abs(x - (xout - 0.5_dp)) < 1.0d-10 .and. x > xout - 0.5_dp) then
            write(6, '(1X,F10.6,3ES14.6,I8)') x, y(1), y(2), y(3), nr-1
        endif

    endif

end subroutine solout_chem

! Right-hand side: Robertson chemical kinetics equations
subroutine fchem(n, x, y, f, rpar, ipar)
               !(N, X, Y, F, RPAR, IPAR)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, ipar
    real(kind=dp), intent(in) :: x, rpar
    real(kind=dp), dimension(n), intent(in) :: y  
    real(kind=dp), dimension(n), intent(out) :: f
    
    ! Rate constants for Robertson problem
    real(kind=dp), parameter :: k1 = 0.04_dp      ! A -> B
    real(kind=dp), parameter :: k2 = 3.0d7        ! B + B -> B + C  
    real(kind=dp), parameter :: k3 = 1.0d4        ! B + C -> A + C

    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) x, rpar, ipar
    endif
    
    ! Chemical kinetics equations:
    ! d[A]/dt = -k1*[A] + k3*[B]*[C]
    ! d[B]/dt = k1*[A] - k3*[B]*[C] - k2*[B]^2
    ! d[C]/dt = k2*[B]^2
    
    f(1) = -k1 * y(1) + k3 * y(2) * y(3)
    f(2) = k1 * y(1) - k3 * y(2) * y(3) - k2 * y(2) * y(2)
    f(3) = k2 * y(2) * y(2)
    
end subroutine fchem

! Jacobian matrix for Robertson chemical kinetics
subroutine jchem(n, x, y, dfy, ldfy, rpar, ipar)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, ldfy, ipar
    real(kind=dp), intent(in) :: x, rpar
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(ldfy, n), intent(out) :: dfy
    
    ! Rate constants
    real(kind=dp), parameter :: k1 = 0.04_dp
    real(kind=dp), parameter :: k2 = 3.0d7
    real(kind=dp), parameter :: k3 = 1.0d4

    ! Silence unused parameter warnings  
    if (.false.) then
        write(*,*) x, rpar, ipar
    endif
    
    ! Jacobian matrix J = df/dy
    ! Row 1: d[A]/dt derivatives
    dfy(1,1) = -k1                    ! df1/d[A]
    dfy(1,2) = k3 * y(3)             ! df1/d[B]  
    dfy(1,3) = k3 * y(2)             ! df1/d[C]
    
    ! Row 2: d[B]/dt derivatives
    dfy(2,1) = k1                     ! df2/d[A]
    dfy(2,2) = -k3 * y(3) - 2.0_dp * k2 * y(2)  ! df2/d[B]
    dfy(2,3) = -k3 * y(2)            ! df2/d[C]
    
    ! Row 3: d[C]/dt derivatives
    dfy(3,1) = 0.0_dp                ! df3/d[A]
    dfy(3,2) = 2.0_dp * k2 * y(2)    ! df3/d[B]
    dfy(3,3) = 0.0_dp                ! df3/d[C]
    
end subroutine jchem