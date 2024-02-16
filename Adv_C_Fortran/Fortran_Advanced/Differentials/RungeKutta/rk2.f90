! Compile and execute with:
!   $ gfortran -O2 rk2.f90 rk2_g.f90 -o rk2
!=========================================================
! Program to solve a 4 ODE system using Runge-Kutta Method
! User mus supply derivatives
!
! dx1                          dx2
! --- = f1(t, x1, x2, x3, x4)  --- = f2(t, x1, x2, x3, x4)
! dt                           dt
!
! dx3                          dx4
! --- = f3(t, x1, x2, x3, x4)  --- = f4(t, x1, x2, x3, x4)
! dt                           dt
!
! as real(kind=8) functions.
! Output is written in rk2.dat
!=========================================================
program rk2_solve
    implicit none
    integer, parameter :: P = 1010000
    real(kind=8), allocatable :: T(:), X1(:), X2(:), V1(:), V2(:)
    real(kind=8)              :: Ti, Tf, X10, X20, V10, V20
    integer                   :: Nt, i
    real(kind=8)              :: k1, k2
    common / couplings /         k1, k2
    real(kind=8)              :: energy, E0, EF, DE
    
    ! Input
    print *, 'Runge-Kutta Method for 4-ODEs Integration'
    print *, 'Enter coupling constants:'
    print *, 'Try for example:'
    print *, '10.0, 0.0'
    read  *, k1, k2
    print *, 'k1 = ', k1, 'k2 = ', k2
    print *, 'Enter Nt, Ti, Tf, X10, X20, V10, V20:'
    print *, 'Try for example:'
    print *, '20000, 0.0, 0.2, 0.0, 0.0, 1.0, 1.0'
    read  *,        Nt, Ti, Tf, X10, X20, V10, V20
    print *, 'Nt = ', Nt
    print *, 'Time: Initial Ti = ', Ti, ' Final Tf = ', Tf
    print *, '           X1(Ti)= ',X10, '    X2(Ti)= ', X20
    print *, '           V1(Ti)= ',V10, '    V2(Ti)= ', V20
    
    allocate(T(Nt), X1(Nt), X2(Nt), V1(Nt), V2(Nt))

    ! The calculation:
    call RK4(T, X1, X2, V1, V2, Ti, Tf, X10, X20, V10, V20, Nt)
    ! Output:
    open(unit=11, file='rk2.dat')
    do i = 1, Nt
        write (11, *)  T(i), X1(i), X2(i), V1(i), V2(i), & 
                energy(T(i), X1(i), X2(i), V1(i), V2(i)) 
    end do
    close(11)

    ! Rutherford scattering angles:
    print *, 'v-angle: ', atan2(V2(Nt), V1(Nt))
    print *, 'b-angle: ', 2.0D0*atan(k1/(V10*V10*X20))
    E0 = energy(Ti   , X10   , X20   , V10   , V20   )     
    EF = energy(T(Nt), X1(Nt), X2(Nt), V1(Nt), V2(Nt)) 
    DE = abs(0.5D0*(EF-E0)/(EF+E0))
    print *, 'E0, EF, DE/E = ', E0, EF, DE
    
    deallocate(T, X1, X2, V1, V2)

end program rk2_solve

!========================================================
! The velocity functions f1, f2(t, x1, x2, x3, x4)
!========================================================
function f1(t, x1, x2, v1, v2)
    implicit none
    real(kind=8) :: f1
    real(kind=8) :: t, x1, x2, v1, v2

    f1 = v1     ! dx1/dt = v1

end function f1
!--------------------------------------------------------
function f2(t, x1, x2, v1, v2)
    implicit none
    real(kind=8) :: f2
    real(kind=8) :: t, x1, x2, v1, v2

    f2 = v2     ! dx2/dt = v2

end function f2
!=====================================================================
! RK4(T, X1, X2, V1, V2, Ti, Tf, X10, X20, V10, V20, Nt) is the driver
! for the Runge-Kutta integration routine RKSTEP
! Input : Initial and final times Ti, Tf
!         Initial values at t=Ti X10, X20, V10, V20
!         Number of integration steps: Nt-1
!         Size of arrays T, X1, X2, V1, V2
! Output : real arrays T(Nt), X1(Nt), X2(Nt),
!                             V1(Nt), V2(Nt) where
! T(1) = Ti X1(1) = X10           X2(1) = X20 
!           V1(1) = V10           V2(1) = V20
!           X1(k) = X1(at t=T(k)) X2(k) = X2(at t=T(k))
!           V1(k) = V1(at t=T(k)  V2(k) = V2(at t=T(k))
! T(Nt)= Tf
!=====================================================================
subroutine RK4(T, X1, X2, V1, V2, Ti, Tf, X10, X20, V10, V20, Nt)
    implicit none
    integer      :: Nt
    real(kind=8), dimension(Nt) :: T, X1, X2, V1, V2 
    real(kind=8) :: Ti, Tf 
    real(kind=8) :: X10, X20 
    real(kind=8) :: V10, V20
    real(kind=8) :: dt
    real(kind=8) :: TS, X1S, X2S ! Values of time and X1, X2, at given step 
    real(kind=8) ::     V1S, V2S
    integer      :: i

    ! Initialize variables:
    dt    = (Tf-Ti)/(Nt-1)
    T (1) = Ti
    X1(1) = X10; X2(1) = X20
    V1(1) = V10; V2(1) = V20
    TS    = Ti
    X1S   = X10; X2S   = X20 
    V1S   = V10; V2S   = V20 
    
    ! Make RK steps: The arguments of RKSTEP are
    ! replaced with new ones
    do i = 2, Nt
        call RKSTEP4(TS, X1S, X2S, V1S, V2S, dt)
        T (i) = TS
        X1(i) = X1S; X2(i) = X2S
        V1(i) = V1S; V2(i) = V2S
    end do
end subroutine RK4

!=============================================================
! Subroutine RKSTEP( t, x1, x2, dt)
! Runge-Kutta Integration routine of ODE
! dx1                          dx2
! --- = f1(t, x1, x2, x3, x4)  --- = f2(t, x1, x2, x3, x4)
! dt                           dt
!
! dx3                          dx4
! --- = f3(t, x1, x2, x3, x4)  --- = f4(t, x1, x2, x3, x4)
! dt                           dt  
! User must supply derivative functions:
! real function f1(t, x1, x2, x3, x4)
! real function f2(t, x1, x2, x3, x4)
! real function f3(t, x1, x2, x3, x4)
! real function f4(t, x1, x2, x3, x4)
! Given initial point (t, x1, x2) the routine advances it
! by time dt .
! Input  : Inital time t    and function values x1, x2, x3, x4
! Output : Final  time t+dt and function values x1, x2, x3, x4
! Careful!: values of t, x1, x2, x3, x4 are overwritten...
!=============================================================
subroutine RKSTEP4(t, x1, x2, x3, x4, dt)
    implicit none
    real(kind=8) :: t, x1, x2, x3, x4, dt
    real(kind=8) :: f1, f2, f3, f4
    real(kind=8) :: k11, k12, k13, k14, k21, k22, k23, k24
    real(kind=8) :: k31, k32, k33, k34, k41, k42, k43, k44
    real(kind=8) :: h, h2, h6
    h   = dt         ! h =dt , i n t e g r a t i o n s t e p
    h2  = 0.5D0 * h  ! h2=h / 2
    h6  = h/6.0D0   ! h6=h / 6
    k11 = f1( t   , x1         , x2         , x3         , x4          )
    k21 = f2( t   , x1         , x2         , x3         , x4          )
    k31 = f3( t   , x1         , x2         , x3         , x4          )
    k41 = f4( t   , x1         , x2         , x3         , x4          )
    k12 = f1( t+h2, x1+h2 * k11, x2+h2 * k21, x3+h2 * k31, x4+h2 * k41 )
    k22 = f2( t+h2, x1+h2 * k11, x2+h2 * k21, x3+h2 * k31, x4+h2 * k41 )
    k32 = f3( t+h2, x1+h2 * k11, x2+h2 * k21, x3+h2 * k31, x4+h2 * k41 )
    k42 = f4( t+h2, x1+h2 * k11, x2+h2 * k21, x3+h2 * k31, x4+h2 * k41 )
    k13 = f1( t+h2, x1+h2 * k12, x2+h2 * k22, x3+h2 * k32, x4+h2 * k42 )
    k23 = f2( t+h2, x1+h2 * k12, x2+h2 * k22, x3+h2 * k32, x4+h2 * k42 )
    k33 = f3( t+h2, x1+h2 * k12, x2+h2 * k22, x3+h2 * k32, x4+h2 * k42 )
    k43 = f4( t+h2, x1+h2 * k12, x2+h2 * k22, x3+h2 * k32, x4+h2 * k42 )
    k14 = f1( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43  )
    k24 = f2( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43  )
    k34 = f3( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43  )
    k44 = f4( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43  )
    t   = t+h
    x1  = x1+h6 * ( k11 + 2.0D0 * ( k12+k13 )+k14 )
    x2  = x2+h6 * ( k21 + 2.0D0 * ( k22+k23 )+k24 )
    x3  = x3+h6 * ( k31 + 2.0D0 * ( k32+k33 )+k34 )
    x4  = x4+h6 * ( k41 + 2.0D0 * ( k42+k43 )+k44 )
end subroutine RKSTEP4
