program gsl_integration
    use iso_c_binding
    implicit none
    
    ! GSL constants
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: GSL_SUCCESS = 0
    
    ! GSL function structure (matches GSL's gsl_function)
    type, bind(c) :: gsl_function
        type(c_funptr) :: function
        type(c_ptr) :: params
    end type
    
    ! Variables for integration
    real(dp) :: result, error
    real(dp) :: a, b  ! Integration limits
    real(dp) :: epsabs, epsrel  ! Absolute and relative error tolerances
    integer :: limit, key
    integer :: status
    type(c_ptr) :: workspace
    type(gsl_function) :: F
    
    ! External GSL functions
    interface
        function gsl_integration_workspace_alloc(n) bind(c, name='gsl_integration_workspace_alloc')
            use iso_c_binding
            integer(c_size_t), value :: n
            type(c_ptr) :: gsl_integration_workspace_alloc
        end function
        
        subroutine gsl_integration_workspace_free(w) bind(c, name='gsl_integration_workspace_free')
            use iso_c_binding
            type(c_ptr), value :: w
        end subroutine
        
        function gsl_integration_qag(f, a, b, epsabs, epsrel, limit, key, workspace, result, abserr) &
                bind(c, name='gsl_integration_qag')
            use iso_c_binding
            import :: gsl_function
            type(gsl_function) :: f
            real(c_double), value :: a, b, epsabs, epsrel
            integer(c_size_t), value :: limit
            integer(c_int), value :: key
            type(c_ptr), value :: workspace
            real(c_double) :: result, abserr
            integer(c_int) :: gsl_integration_qag
        end function
        
        function gsl_integration_qags(f, a, b, epsabs, epsrel, limit, workspace, result, abserr) &
                bind(c, name='gsl_integration_qags')
            use iso_c_binding
            import :: gsl_function
            type(gsl_function) :: f
            real(c_double), value :: a, b, epsabs, epsrel
            integer(c_size_t), value :: limit
            type(c_ptr), value :: workspace
            real(c_double) :: result, abserr
            integer(c_int) :: gsl_integration_qags
        end function
    end interface
    
    ! Set up integration parameters
    a = 0.0_dp        ! Lower limit
    b = 1.0_dp        ! Upper limit
    epsabs = 0.0_dp   ! Absolute error tolerance
    epsrel = 1.0e-7_dp ! Relative error tolerance
    limit = 1000      ! Maximum number of subintervals
    key = 1           ! Integration rule (15-point Gauss-Kronrod)
    
    ! Allocate workspace
    workspace = gsl_integration_workspace_alloc(int(limit, c_size_t))
    if (.not. c_associated(workspace)) then
        write(*,*) 'Failed to allocate GSL workspace'
        stop 1
    end if
    
    write(*,*) '================================================='
    write(*,*) '    GSL Integration Examples in Fortran'
    write(*,*) '================================================='
    write(*,*)
    
    ! Example 1: Simple polynomial function x^2
    write(*,*) 'Example 1: Integrating x^2 from 0 to 1'
    write(*,*) 'Analytical result: 1/3 = 0.333333...'
    write(*,*)
    
    F%function = c_funloc(polynomial_function)
    F%params = c_null_ptr
    
    status = gsl_integration_qag(F, real(a, c_double), real(b, c_double), &
                                real(epsabs, c_double), real(epsrel, c_double), &
                                int(limit, c_size_t), int(key, c_int), &
                                workspace, result, error)
    
    if (status == GSL_SUCCESS) then
        write(*,'(A,F12.8)') '  Numerical result: ', result
        write(*,'(A,E12.4)') '  Estimated error:  ', error
        write(*,'(A,F12.8)') '  Analytical diff:  ', abs(result - 1.0_dp/3.0_dp)
    else
        write(*,*) '  Integration failed with status:', status
    end if
    
    write(*,*)
    write(*,*) '-------------------------------------------------'
    write(*,*)
    
    ! Example 2: Exponential function e^(-x^2) - Gaussian
    write(*,*) 'Example 2: Integrating exp(-x^2) from 0 to 1'
    write(*,*) 'Using adaptive QAGS algorithm'
    write(*,*)
    
    F%function = c_funloc(gaussian_function)
    F%params = c_null_ptr
    
    status = gsl_integration_qags(F, real(a, c_double), real(b, c_double), &
                                 real(epsabs, c_double), real(epsrel, c_double), &
                                 int(limit, c_size_t), workspace, result, error)
    
    if (status == GSL_SUCCESS) then
        write(*,'(A,F12.8)') '  Numerical result: ', result
        write(*,'(A,E12.4)') '  Estimated error:  ', error
    else
        write(*,*) '  Integration failed with status:', status
    end if
    
    write(*,*)
    write(*,*) '-------------------------------------------------'
    write(*,*)
    
    ! Example 3: Oscillatory function sin(10*x)
    write(*,*) 'Example 3: Integrating sin(10*x) from 0 to 1'
    write(*,*) 'Analytical result: (1-cos(10))/10 = 0.155633...'
    write(*,*)
    
    F%function = c_funloc(oscillatory_function)
    F%params = c_null_ptr
    
    ! Use key=6 for oscillatory functions and higher precision
    epsrel = 1.0e-10_dp
    key = 6
    status = gsl_integration_qag(F, real(a, c_double), real(b, c_double), &
                                real(epsabs, c_double), real(epsrel, c_double), &
                                int(limit, c_size_t), int(key, c_int), &
                                workspace, result, error)
    
    if (status == GSL_SUCCESS) then
        write(*,'(A,F12.8)') '  Numerical result: ', result
        write(*,'(A,E12.4)') '  Estimated error:  ', error
        write(*,'(A,F12.8)') '  Analytical diff:  ', abs(result - (1.0_dp - cos(10.0_dp))/10.0_dp)
    else
        write(*,*) '  Integration failed with status:', status
    end if
    
    ! Clean up
    call gsl_integration_workspace_free(workspace)
    
    write(*,*)
    write(*,*) '================================================='
    write(*,*) 'Integration completed successfully!'
    write(*,*) '================================================='

contains

    ! Function 1: Simple polynomial x^2
    function polynomial_function(x, params) result(y) bind(c)
        use iso_c_binding
        real(c_double), value :: x
        type(c_ptr), value :: params
        real(c_double) :: y
        
        y = x * x
    end function polynomial_function
    
    ! Function 2: Gaussian exp(-x^2)
    function gaussian_function(x, params) result(y) bind(c)
        use iso_c_binding
        real(c_double), value :: x
        type(c_ptr), value :: params
        real(c_double) :: y
        
        y = exp(-x * x)
    end function gaussian_function
    
    ! Function 3: Oscillatory sin(10*x)
    function oscillatory_function(x, params) result(y) bind(c)
        use iso_c_binding
        real(c_double), value :: x
        type(c_ptr), value :: params
        real(c_double) :: y
        
        y = sin(10.0d0 * x)
    end function oscillatory_function

end program gsl_integration