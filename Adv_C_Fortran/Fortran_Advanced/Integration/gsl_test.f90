program main
    !use iso_c_binding
    !use gsl_qag

    implicit none

    ! Function to be integrated
    real(c_float) function f(x, params)
        real(c_float), intent(in) :: x
        type(c_ptr), intent(in) :: params
        f = x**2
    end function

    ! Integration workspace
    type(gsl_integration_workspace) :: workspace
    integer(c_int) :: status
    real(c_float) :: result, error

    ! Allocate workspace
    workspace = gsl_integration_workspace_alloc(10000)

    ! Set up parameters for integration
    real(c_float), parameter :: epsabs = 1.0e-6, epsrel = 1.0e-6
    integer(c_int), parameter :: key = 6

    ! Perform the integration
    call gsl_integration_qag(f, 0.0, 1.0, epsabs, epsrel, 10000, key, workspace, result, error, status)

    ! Print the result
    print *, "Result = ", result
    print *, "Error = ", error

    ! Deallocate workspace
    call gsl_integration_workspace_free(workspace)

end program
