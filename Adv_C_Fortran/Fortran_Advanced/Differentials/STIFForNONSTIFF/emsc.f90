program enhanced_stiffness_detector
    implicit none
    
    ! Maximum system size
    integer, parameter :: MAX_N = 10
    integer :: choice
    !double precision :: params(10)  ! Parameter array
    !character(len=50) :: system_name
    
    write(*,*) 'Enhanced Stiffness Detector'
    write(*,*) '=========================='
    write(*,*) '1. RLC Circuit'
    write(*,*) '2. Van der Pol Oscillator' 
    write(*,*) '3. Lorenz System'
    write(*,*) '4. Custom 2x2 System'
    write(*,*) 'Choose system (1-4): '
    read(*,*) choice
    
    select case(choice)
        case(1)
            call analyze_rc_circuit()
        case(2) 
            call analyze_van_der_pol()
        case(3)
            call analyze_lorenz()
        case(4)
            call analyze_custom_2x2()
        case default
            write(*,*) 'Invalid choice!'
    end select
    
end program enhanced_stiffness_detector

!==============================================================================
! RC Circuit Analysis
!==============================================================================
subroutine analyze_rc_circuit()
    implicit none
    integer, parameter :: n = 2
    double precision :: R, C, L
    double precision :: A(n,n)
    
    write(*,*) 'RC Circuit Analysis'
    write(*,*) 'Enter R (Resistance): '
    read(*,*) R
    write(*,*) 'Enter C (Capacitance): '
    read(*,*) C  
    write(*,*) 'Enter L (Inductance): '
    read(*,*) L
    
    call compute_rc_jacobian(R, C, L, A, n)
    call analyze_stiffness(A, n, 'RC Circuit')
    
end subroutine analyze_rc_circuit

subroutine compute_rc_jacobian(R, C, L, J, n)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: R, C, L
    double precision, intent(out) :: J(n,n)
    
    ! Jacobian for RC circuit system
    J(1,1) = -1.0d0/(R*C)
    J(1,2) = 1.0d0/C
    J(2,1) = -1.0d0/L  
    J(2,2) = 0.0d0
    
end subroutine compute_rc_jacobian

!==============================================================================
! Van der Pol Oscillator Analysis  
!==============================================================================
subroutine analyze_van_der_pol()
    implicit none
    integer, parameter :: n = 2
    double precision :: mu, x0, y0
    double precision :: A(n,n)
    
    write(*,*) 'Van der Pol Oscillator Analysis'
    write(*,*) 'Enter mu (damping parameter): '
    read(*,*) mu
    write(*,*) 'Enter x0 (linearization point): '
    read(*,*) x0
    write(*,*) 'Enter y0 (linearization point): '  
    read(*,*) y0
    
    call compute_vdp_jacobian(mu, x0, y0, A, n)
    call analyze_stiffness(A, n, 'Van der Pol')
    
end subroutine analyze_van_der_pol

subroutine compute_vdp_jacobian(mu, x, y, J, n)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: mu, x, y
    double precision, intent(out) :: J(n,n)
    
    ! Jacobian for Van der Pol oscillator
    J(1,1) = 0.0d0
    J(1,2) = 1.0d0
    J(2,1) = -2.0d0*mu*x*y - 1.0d0
    J(2,2) = mu*(1.0d0 - x*x)
    
end subroutine compute_vdp_jacobian

!==============================================================================
! Lorenz System Analysis
!==============================================================================
subroutine analyze_lorenz()
    implicit none
    integer, parameter :: n = 3  
    double precision :: sigma, rho, beta, x0, y0, z0
    double precision :: A(n,n)
    
    write(*,*) 'Lorenz System Analysis'
    write(*,*) 'Enter sigma (Prandtl number): '
    read(*,*) sigma
    write(*,*) 'Enter rho (Rayleigh number): '
    read(*,*) rho
    write(*,*) 'Enter beta (geometric factor): '
    read(*,*) beta
    write(*,*) 'Enter x0, y0, z0 (linearization point): '
    read(*,*) x0, y0, z0
    
    call compute_lorenz_jacobian(sigma, rho, beta, x0, y0, z0, A, n)
    call analyze_stiffness(A, n, 'Lorenz System')
    
end subroutine analyze_lorenz

subroutine compute_lorenz_jacobian(sigma, rho, beta, x, y, z, J, n)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: sigma, rho, beta, x, y, z  
    double precision, intent(out) :: J(n,n)
    
    ! Jacobian for Lorenz system
    J(1,1) = -sigma
    J(1,2) = sigma  
    J(1,3) = 0.0d0
    J(2,1) = rho - z
    J(2,2) = -1.0d0
    J(2,3) = -x
    J(3,1) = y
    J(3,2) = x
    J(3,3) = -beta
    
end subroutine compute_lorenz_jacobian

!==============================================================================
! Custom 2x2 System
!==============================================================================
subroutine analyze_custom_2x2()
    implicit none
    integer, parameter :: n = 2
    double precision :: A(n,n)
    integer :: i, j
    
    write(*,*) 'Custom 2x2 Jacobian Matrix'
    write(*,*) 'Enter matrix elements:'
    
    do i = 1, n
        do j = 1, n
            write(*,'(A,I1,A,I1,A)',advance='no') 'J(', i, ',', j, '): '
            read(*,*) A(i,j)
        end do
    end do
    
    call analyze_stiffness(A, n, 'Custom 2x2')
    
end subroutine analyze_custom_2x2

!==============================================================================
! General Stiffness Analysis Routine
!==============================================================================
subroutine analyze_stiffness(A, n, system_name)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: A(n,n)
    character(len=*), intent(in) :: system_name
    
    integer :: lda, ldvl, ldvr, lwork, info, i
    double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:), WORK(:)
    double precision, allocatable :: A_copy(:,:)
    double precision :: max_re, min_re, ratio, max_abs, condition_number
    
    lda = n
    ldvl = n  
    ldvr = n
    lwork = 4 * n
    
    ! Allocate arrays
    allocate(A_copy(n,n), WR(n), WI(n), VL(n,n), VR(n,n), WORK(lwork))
    
    ! Copy matrix (LAPACK destroys input)
    A_copy = A
    
    write(*,*) 
    write(*,*) '======================================='
    write(*,'(A,A)') 'Stiffness Analysis for: ', trim(system_name)
    write(*,*) '======================================='
    
    ! Display Jacobian matrix
    write(*,*) 'Jacobian Matrix:'
    do i = 1, n
        write(*,'(10F12.6)') A(i,:)
    end do
    
    ! Compute eigenvalues
    call dgeev('V', 'V', n, A_copy, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, lwork, info)
    
    if (info /= 0) then
        write(*,*) "Error computing eigenvalues: INFO =", info
        return
    end if
    
    ! Display eigenvalues
    write(*,*) 
    write(*,*) 'Eigenvalues:'
    do i = 1, n
        if (abs(WI(i)) < 1.0d-12) then
            write(*,'(A,I1,A,F15.8)') 'λ', i, ' = ', WR(i)
        else
            write(*,'(A,I1,A,F15.8,A,F15.8,A)') 'λ', i, ' = ', WR(i), ' + ', WI(i), 'i'
        end if
    end do
    
    ! Calculate stiffness metrics
    max_re = 0.0d0
    min_re = huge(1.0d0)
    max_abs = 0.0d0
    
    do i = 1, n
        max_abs = max(max_abs, sqrt(WR(i)**2 + WI(i)**2))
        if (abs(WR(i)) > 1.0d-12) then  ! Avoid division by very small numbers
            max_re = max(max_re, abs(WR(i)))
            min_re = min(min_re, abs(WR(i)))
        end if
    end do
    
    if (min_re < huge(1.0d0) .and. min_re > 1.0d-12) then
        ratio = max_re / min_re
        condition_number = max_abs / min_re
    else
        ratio = 0.0d0
        condition_number = 0.0d0
    end if
    
    write(*,*) 
    write(*,*) 'Stiffness Metrics:'
    write(*,'(A,ES15.6)') 'Maximum |Re(λ)|: ', max_re
    write(*,'(A,ES15.6)') 'Minimum |Re(λ)|: ', min_re  
    write(*,'(A,ES15.6)') 'Stiffness Ratio: ', ratio
    write(*,'(A,ES15.6)') 'Condition Number: ', condition_number
    
    write(*,*)
    if (ratio > 1000.0d0) then
        write(*,*) '>>> This system is likely STIFF.'
        write(*,*) '    Consider using implicit methods (BDF, Radau, etc.)'
    else if (ratio > 100.0d0) then
        write(*,*) '>>> This system has moderate stiffness.'  
        write(*,*) '    Semi-implicit methods may be beneficial.'
    else
        write(*,*) '>>> This system is likely NOT stiff.'
        write(*,*) '    Explicit methods (RK4, etc.) should work well.'
    end if
    
    ! Clean up
    deallocate(A_copy, WR, WI, VL, VR, WORK)
    
end subroutine analyze_stiffness