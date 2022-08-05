program main
    implicit none
    real(4) :: epsilon, x, xvana, fx, fxvana, xuus
    integer(4) :: count

    call loikaja(epsilon, x, xvana, fx, fxvana, xuus)
    print *, 'f (x) = ', fx, 'x = ',x

end program main

subroutine loikaja(epsilon, x, xvana, fx, fxvana, xuus)
    real(4) :: epsilon, x, xvana, fx, fxvana, xuus
    integer(4) :: count
    
    epsilon = 1e-6
    x       = 490
    xvana   = -2.0
    count   = 0
    do while (abs(x-xvana) >= epsilon)
        count   = count + 1
        fx      = 0.1*x     - 2.0E-4*x**2
        fxvana  = 0.1*xvana - 2.0E-4*xvana**2
        xuus    = x - fx*(x-xvana)/(fx-fxvana)      !Lõikajate meetod
        xvana   = x
        x       = xuus
    end do
    
end subroutine loikaja
