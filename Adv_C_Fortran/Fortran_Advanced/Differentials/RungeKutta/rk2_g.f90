!==============================================================
! The acceleration functions f3, f4(t, x1, x2, v1, v2) provided
! by the user
!==============================================================
! Free fall in constant gravitational field filled with
! g = -k2
function f3(t, x1, x2, v1, v2)
    implicit none
    real(kind=8) :: f3
    real(kind=8) :: t, x1, x2, v1, v2
    real(kind=8) ::            k1, k2
    common /couplings/         k1, k2
    
    f3 = 0.0D0  ! dx3/dt = dv1/dt = a1

end function f3
!--------------------------------------------------------------
function f4(t, x1, x2, v1, v2)
    implicit none
    real(kind=8) :: f4
    real(kind=8) :: t, x1, x2, v1, v2
    real(kind=8) ::            k1, k2
    common /couplings/         k1, k2
    
    f4 = -k1  ! dx3/dt = dv1/dt = a1

end function f4
!--------------------------------------------------------------
function energy(t, x1, x2, v1, v2)
    implicit none
    real(kind=8) :: energy
    real(kind=8) :: t, x1, x2, v1, v2
    real(kind=8) ::            k1, k2
    common /couplings/         k1, k2
    
    energy = 0.5D0*(v1*v1 + v2*v2) + k1*k2

end function energy