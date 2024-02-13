!===============================================================
! Program to solve a 2 ODE system using Runge-Kutta Method
! User must supply derivatives
! dx1/dt = f1(t, x1, x2)    dx2/dt = f2(t, x1, x2)
! as real(kind=8) functions
! Output is written in file rk.dat
!===============================================================

program rk_solve
    implicit none
    !integer, parameter         :: P = 110000
    real(kind=8), allocatable  :: T(:), X1(:), X2(:), xhvec(:) ! Declare variables as allocatable
    real(kind=8)               :: Ti, Tf, X10, X20, dt
    integer                    :: Nt, i
    real(kind=8), EXTERNAL     :: xh, vh

    ! Input 
    print *, 'Runge-Kutta Method for 2-ODEs Integration'
    !print *, 'Enter Nt, Ti, Tf, X10, X20:'
    !read  *, Nt, Ti, Tf, X10, X20
    Nt = 50000
    Ti = 0.0D0
    Tf = 6.0D0
    X10 = 0.2D0
    X20 = 0.0D0
    print *, 'Nt =', Nt
    print *, 'Time: Initial Ti = ',  Ti, ' Final Tf = ',  Tf
    print *, '           X1(Ti)= ', X10, '    X2(Ti)= ', X20
    !if (Nt .gt. P) stop 'Nt > P'

    ! Allocate memory for arrays
    allocate(T(Nt), X1(Nt), X2(Nt), xhvec(Nt))

    ! The Calculation:
    call  RK(T, X1, X2, Ti, Tf, X10, X20, Nt)

    ! Output:
    open (unit = 11, file = 'rk.dat')
    do i = 1, Nt
        write (11, *) T(i), X1(i), X2(i)
    end do
    close(11)

    dt = (Tf-Ti)/(Nt-1)
    do i = 1, Nt
        xhvec(i) = xh(i*dt)
    end do

    open (unit = 12, file = 'err.dat')
    do i = 1, Nt
        write (12, *) T(i), abs(X1(i) - xhvec(i))
    end do
    close(12)

    ! Deallocate memory
    deallocate(T, X1, X2, xhvec)

end program rk_solve

!=====================================================
! The functions f1, f2(t, x1, x2) provided by the user
!=====================================================
function f1(t, x1, x2)
    implicit none
    real(kind=8) :: f1
    real(kind=8) :: t, x1, x2

    f1 = x2     ! dx1/dt = v = x2

end function f1
!-----------------------------------------------------
function f2(t, x1, x2)
    implicit none
    real(kind=8) :: f2
    real(kind=8) :: t, x1, x2

    f2 = -10.0D0*x1 ! harmonic oscillator 

end function f2

function xh(t)
    implicit none
    real(kind=8) :: xh, t
    real(kind=8) :: x0, v0, omega

    x0 = 0.2D0
    v0 = 0.0D0
    omega = sqrt(10.D0)
    xh = x0*cos(omega*t) + v0/omega*sin(omega*t)

end function xh

function vh(t)
    implicit none
    real(kind=8) :: vh, t
    real(kind=8) :: x0, v0, omega

    x0 = 0.2D0
    v0 = 0.0D0
    omega = sqrt(10.D0)
    vh = v0*cos(omega*t) - x0*omega*sin(omega*t)

end function vh

!==============================================================
! RK(T, X1, X2, Ti, Tf, X10, X20, Nt) is the driver
! for the Runge-Kutta integration routine RKSTEP
! Input: Initial and final times Ti, Tf 
!        Initial values at t = Tf   X10, X20
!        Number of steps of integration: Nt-1
!        Size of arrays T, X1, X2
! Output: real(kind=8) arrays T(Nt), X1(Nt), X2(Nt) where
! T(1) = Ti     X1(1) = X10             X2(1) = X20
!               X1(k) = X1(at t = T(k)) X2(k) = X2(at t = T(k))
! T(Nt) = Tf
!==============================================================
subroutine RK(T, X1, X2, Ti, Tf, X10, X20, Nt)
    implicit none
    integer, intent(in) :: Nt                       ! intent(in) for values being used but not modified
    real(kind=8), dimension(Nt), intent(inout) :: T, X1, X2 ! intent(inout) for arrays being modified
    real(kind=8),    intent(in) :: Ti, Tf, X10, X20         ! intent(in) for values being used but not modified
    real(kind=8)                :: dt
    real(kind=8)                :: Ts, X1S, X2S
    integer                     :: i
    
    ! Initialize variables
    dt = (Tf-Ti)/(Nt-1)
    T (1) = Ti
    X1(1) = X10
    X2(1) = X20
    Ts    = Ti
    X1S   = X10
    X2S   = X20

    ! Make RK steps: The arguments of RKSTEP
    ! are replaced with the new ones!
    do i = 2, Nt
        call RKSTEP(Ts, X1S, X2S, dt)
        T (i) = Ts
        X1(i) = X1S
        X2(i) = X2S
    end do
end subroutine RK

!=========================================================
! Subroutine RKSTEP(t, x1, x2, dt)
! Runge-Kutta Integration routine of ODE
! dx1/dt = f1(t, x1, x2) dx2/dt = f2(t, x1, x2)
! User must supply derivative functions:
!   real(kind=8) function f1(t, x1, x2)
!   real(kind=8) function f1(t, x1, x2)
! Given initial point (t, x1, x2) the routine advances it
! by time dt.
! Input:  Initial time t    and function values x1, x2
! Output: Final time   t+dt and function values x1, x2
! Careful!: values of t, x1, and x2 are overwritten
!=========================================================
subroutine RKSTEP(t, x1, x2, dt)
    implicit none
    real(kind=8) :: t, x1, x2, dt
    real(kind=8) :: f1, f2
    real(kind=8) :: k11, k12, k13, k14, k21, k22, k23, k24
    real(kind=8) :: h, h2, h6

    h = dt          ! h = dt, integartion step
    h2 = 0.5D0*h    ! h2 = h/2
    h6 = h/6.0D0      ! h6=h / 6
    k11 = f1(t   , x1         , x2     )
    k21 = f2(t   , x1         , x2     )
    k12 = f1(t+h2, x1+h2*k11, x2+h2*k21)
    k22 = f2(t+h2, x1+h2*k11, x2+h2*k21)
    k13 = f1(t+h2, x1+h2*k12, x2+h2*k22)
    k23 = f2(t+h2, x1+h2*k12, x2+h2*k22)
    k14 = f1(t+h , x1+h*k13 , x2+h*k23 )
    k24 = f2(t+h , x1+h*k13 , x2+h*k23 )
    t   = t + h
    x1  = x1 + h6*(k11 + 2.0D0*(k12+k13) + k14)
    x2  = x2 + h6*(k21 + 2.0D0*(k22+k23) + k24)

end subroutine RKSTEP