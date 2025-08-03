! Compile and execute with:
! $ gfortran chem_rodas.f90 rodas.f90 decsol.f90 dc_decsol.f90 -o chem_rodas
! $ ./chem_rodas

! RODAS solver for stiff chemical oscillator problems
! Solves the Brusselator reaction system - a classic stiff oscillatory benchmark:
! A -> X (k1 = 1.0)
! 2X + Y -> 3X (k2 = 3.0) 
! B + X -> Y + D (k3 = 1.0)
! X -> E (k4 = 1.0)
! This creates the famous Brusselator oscillations in X and Y concentrations

program rodas_rlcd
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external fcircuit, jcircuit, solout_circuit, f
    
    ! Problem parameters
    integer, parameter :: neq = 2
    real(kind=dp), dimension(neq) :: y
    real(kind=dp) :: x, xend, h, rtol, atol, rpar
    integer :: n, ifcn, ijac, mljac, mujac, imas, mlmas, mumas, iout, itol, ipar, idid
    integer :: lwork, liwork, idfx
    real(kind=dp), allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork
    
    ! Calculate required array sizes and allocate
    call calculate_rodas_arrays(neq, lwork, liwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    ! Initialize problem parameters
    call setup_circuit_problem(neq, y, x, xend, h, rtol, atol, rpar, ipar, &
                                  n, ifcn, ijac, mljac, mujac, imas, mlmas, mumas, &
                                  iout, itol, idfx, lwork, liwork, work, iwork)
    
    ! Run the integration
    call run_rodas_integration(fcircuit, jcircuit, solout_circuit, n, y, x, xend, h, &
                              rtol, atol, itol, ifcn, ijac, mljac, mujac, imas, &
                              mlmas, mumas, iout, idfx, work, lwork, iwork, liwork, &
                              rpar, ipar, idid)
    
    ! Print final statistics
    call print_circuit_statistics(n, x, y, rtol, liwork, iwork, idid)
    
    ! Clean up
    deallocate(work)
    deallocate(iwork)
    
end program rodas_rlcd

! Subroutine to calculate RODAS array lengths
subroutine calculate_rodas_arrays(neq, lwork, liwork)
    implicit none
    integer, intent(in) :: neq
    integer, intent(out) :: lwork, liwork
    
    ! For full Jacobian RODAS requirements
    lwork = 2 * neq * neq + 14 * neq + 20
    liwork = neq + 20
    
    write(6, '(A,I0,A,I0)') 'Calculated RODAS array sizes: LWORK = ', lwork, ', LIWORK = ', liwork
    
end subroutine calculate_rodas_arrays

! Setup the Brusselator chemical oscillator problem
subroutine setup_circuit_problem(neq, y, x, xend, h, rtol, atol, rpar, ipar, &
                                     n, ifcn, ijac, mljac, mujac, imas, mlmas, mumas, &
                                     iout, itol, idfx, lwork, liwork, work, iwork)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: neq, lwork, liwork
    real(kind=dp), dimension(neq), intent(out) :: y
    real(kind=dp), intent(out) :: x, xend, h, rtol, atol, rpar
    integer, intent(out) :: n, ifcn, ijac, mljac, mujac, imas, mlmas, mumas, iout, itol, ipar, idfx
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer :: i
    
    write(6, '(A)') 'Setting up Brusselator chemical oscillator problem...'
    write(6, '(A)') 'Reactions:'
    write(6, '(A)') '  A -> X         (k1 = 1.0)'
    write(6, '(A)') '  2X + Y -> 3X   (k2 = 3.0)'
    write(6, '(A)') '  B + X -> Y + D (k3 = 1.0)'
    write(6, '(A)') '  X -> E         (k4 = 1.0)'
    write(6, '(A)') 'This creates oscillations in X and Y concentrations.'
    write(6, '(A)') ''
    
    ! Problem dimension
    n = neq
    
    ! RODAS options
    ifcn = 0      ! Problem is autonomous (no explicit time dependence)
    ijac = 1      ! Compute Jacobian analytically
    mljac = n     ! Full Jacobian
    mujac = 0     ! Not used for full Jacobian
    imas = 0      ! Explicit ODE system (no mass matrix)
    mlmas = 0     ! Not used
    mumas = 0     ! Not used
    iout = 1      ! Use output routine
    idfx = 0      ! df/dx not provided (autonomous system)
    
    ! Initial conditions for Brusselator oscillations
    x = 0.0_dp               ! Initial time
    y(1) = 1.5_dp           ! [X] initial concentration
    y(2) = 3.0_dp           ! [Y] initial concentration
    
    ! Integration parameters for observing oscillations
    xend = 20.0_dp          ! Final time (enough to see several oscillations)
    h = 1.0d-4              ! Initial step size
    
    ! Tolerance settings for stiff oscillatory problem
    rtol = 1.0d-8           ! Relative tolerance
    atol = 1.0d-10          ! Absolute tolerance
    itol = 0                ! Same tolerance for all components
    
    ! Brusselator parameters (A=1, B=3 gives nice oscillations)
    rpar = 3.0_dp           ! B parameter (A is hardcoded as 1.0)
    ipar = 0                ! Not used
    
    ! Initialize work arrays
    do i = 1, 20
        iwork(i) = 0
        work(i) = 0.0_dp
    end do
    
    write(6, '(A)') 'Problem setup completed.'
    write(6, '(A,2F10.6)') 'Initial concentrations [X,Y]: ', y(1), y(2)
    write(6, '(A,F8.3)') 'Integration from t = 0 to t = ', xend
    write(6, '(A,ES10.3,A,ES10.3)') 'Tolerances: RTOL = ', rtol, ', ATOL = ', atol
    write(6, '(A,F4.1)') 'Brusselator parameter B = ', rpar
    write(6, '(A)') ''
    
end subroutine setup_circuit_problem

! Main RODAS integration routine
subroutine run_rodas_integration(f, jac, solout, n, y, x, xend, h, rtol, atol, &
                                itol, ifcn, ijac, mljac, mujac, imas, mlmas, mumas, &
                                iout, idfx, work, lwork, iwork, liwork, rpar, ipar, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external f, jac, solout
    integer, intent(in) :: n, itol, ifcn, ijac, mljac, mujac, imas, mlmas, mumas, iout, idfx
    integer, intent(in) :: lwork, liwork, ipar
    real(kind=dp), intent(inout) :: x, h, rpar
    real(kind=dp), intent(in) :: xend, rtol, atol
    real(kind=dp), dimension(n), intent(inout) :: y
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer, intent(out) :: idid
    
    write(6, '(A)') 'Starting RODAS integration for Brusselator oscillator...'
    write(6, '(A)') '     Time         [X]           [Y]         NSTEP'
    write(6, '(A)') '----------------------------------------------------'
    
    ! Call RODAS solver
    call rodas(n, f, ifcn, x, y, xend, h, &
               rtol, atol, itol, &
               jac, ijac, mljac, mujac, f, idfx, &
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
    
end subroutine run_rodas_integration

! Print final statistics
subroutine print_circuit_statistics(n, x, y, rtol, liwork, iwork, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, liwork
    real(kind=dp), intent(in) :: x, rtol
    real(kind=dp), dimension(n), intent(in) :: y
    integer, dimension(liwork), intent(in) :: iwork
    integer, intent(in) :: idid
    
    write(6, '(/A)') 'Final Integration Statistics:'
    write(6, '(A)') '============================'
    write(6, '(A,F10.6,A,2ES15.6)') 'Final time: ', x, ', Final solution [X,Y]: ', y(1), y(2)
    
    write(6, '(A,ES10.3)') 'Relative tolerance used: ', rtol
    write(6, '(A,I0)') 'Integration status (IDID): ', idid
    write(6, '(/A)') 'Detailed statistics:'
    write(6, '(A,I6)') ' Function evaluations: ', iwork(14)
    write(6, '(A,I6)') ' Jacobian evaluations: ', iwork(15)
    write(6, '(A,I6)') ' Integration steps:    ', iwork(16)
    write(6, '(A,I6)') ' Accepted steps:       ', iwork(17)
    write(6, '(A,I6)') ' Rejected steps:       ', iwork(18)
    write(6, '(A,I6)') ' LU decompositions:    ', iwork(19)
    write(6, '(A,I6)') ' Linear solves:        ', iwork(20)
    
    if (iwork(16) > 0) then
        write(6, '(A,ES10.3)') ' Average step size: ', x / real(iwork(16), dp)
        write(6, '(A,F7.2,A)') ' Acceptance ratio: ', &
                               100.0_dp * real(iwork(17), dp) / real(iwork(16), dp), '%'
    endif
    
    write(6, '(/A)') 'Note: RODAS is particularly efficient for this oscillatory'
    write(6, '(A)') '      stiff problem due to its L-stable properties.'
    
end subroutine print_circuit_statistics

! Output routine for continuous solution monitoring
subroutine solout_circuit(nr, xold, x, y, cont, lrc, n, rpar, ipar, irtrn)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: nr, lrc, n, ipar
    real(kind=dp), intent(in) :: xold, x, rpar
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(lrc), intent(in) :: cont
    integer, intent(inout) :: irtrn
    real(kind=dp) :: contro
    
    ! Storage for output control
    real(kind=dp), save :: xout = 0.0_dp
    real(kind=dp), save :: dt = 0.5_dp  ! Output interval for oscillations
    logical, save :: first_call = .true.

    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) xold, rpar, ipar, irtrn
    endif
    
    if (first_call) then
        ! Print initial condition
        write(6, '(1X,F10.6,2ES14.6,I8)') x, y(1), y(2), nr-1
        xout = dt
        first_call = .false.
    else
        ! Print at regular intervals to capture oscillations
        do while (xout <= x)
            write(6, '(1X,F10.6,2ES14.6,I8)') xout, &
                   contro(1, xout, cont, lrc), &
                   contro(2, xout, cont, lrc), nr-1
            xout = xout + dt
        end do
    endif

end subroutine solout_circuit

! Right-hand side: Brusselator chemical oscillator equations
subroutine fcircuit(n, t, y, f, rpar, ipar)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, ipar
    real(kind=dp), intent(in) :: t, rpar
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(n), intent(out) :: f

    ! Circuit parameters
    real(kind=dp), parameter :: R = 100.0_dp, L = 0.1_dp, C = 1e-6_dp
    real(kind=dp) :: Vin, diode_on

    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) rpar, ipar
    endif

    ! Vin is a step input: 5V from t > 0
    Vin = merge(5.0_dp, 0.0_dp, t > 0.0_dp)

    ! y(1) = I(t), y(2) = V_C(t)
    ! Diode model: ideal switch, blocks when V_C < 0
    diode_on = merge(1.0_dp, 0.0_dp, y(2) > 0.0_dp)

    ! dI/dt and dVC/dt
    f(1) = (Vin - R * y(1) - y(2)) / L
    f(2) = diode_on * y(1) / C
end subroutine fcircuit


subroutine jcircuit(n, t, y, dfy, ldfy, rpar, ipar)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, ldfy, ipar
    real(kind=dp), intent(in) :: t, rpar
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(ldfy,n), intent(out) :: dfy

    real(kind=dp), parameter :: R = 100.0_dp, L = 0.1_dp, C = 1e-6_dp
    real(kind=dp) :: diode_on

    if (.false.) then
        write (*,*) t, rpar, ipar
    end if

    diode_on = merge(1.0_dp, 0.0_dp, y(2) > 0.0_dp)

    dfy(1,1) = -R / L      ! df1/dI
    dfy(1,2) = -1.0_dp / L ! df1/dVC
    dfy(2,1) = diode_on / C
    dfy(2,2) = 0.0_dp
end subroutine jcircuit
