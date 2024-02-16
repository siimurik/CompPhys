program main
    ! Declaration of variables
    real    ::  y(4), work(27), mu_atm, m_keha, m_maa  
    integer :: iwork(5)
    common gm, r_maa, ro_atm_0, p1, p2
    external keplerf
    
    ! Opening output files
    open(20,file="fall.dat")
    open(21,file="fall_lv.dat")
    
    ! Constants
    pi     = 3.1415926536
    c_keha = 0.47
    r_keha = 2.
    m_keha = 2800.

    p_atm  = 101000.
    mu_atm = 29.
    r_gaas = 8314.
    t_atm  = 300.

    G     = 6.67408e-11
    m_maa = 5.9722e+24
    gm    = G*m_maa*1.e-9
    r_maa = 6378.

    s_keha   = pi*r_keha**2
    ro_atm_0 = p_atm*mu_atm/(r_gaas*t_atm)
    p2       = s_keha*c_keha*500./m_keha
    p1       = mu_atm*9814./(r_gaas*t_atm)

    ! Initial conditions
    neqn   = 4
    dt     = 1.
    nt     = 20000
    abserr = 1.e-8
    relerr = 1.e-8
    v0     = 7.7884079
    alpha  = 0.7
        
    y(1)  = 6578.0
    y(2)  = 0.
    y(3)  = -v0*sin(pi*alpha/180.)
    y(4)  =  v0*cos(pi*alpha/180.)
    t     = 0
    iflag = 1
        
    do i = 1, nt
        tout = t + dt
        call rkf45(keplerf, neqn, y, t, tout, relerr, abserr, iflag, work, iwork)

        if (iflag.ne.2) then
            iflag = 2
        endif

        v      = sqrt(y(3)**2 + y(4)**2)
        r      = sqrt(y(1)**2 + y(2)**2)
        h      = r - r_maa
        ro_atm = ro_atm_0 * exp(-p1 * h)
        ax     = -y(1) * gm / r**3 - p2 * ro_atm * v * y(3)
        ay     = -y(2) * gm / r**3 - p2 * ro_atm * v * y(4)
        a      = sqrt(ax**2 + ay**2)

        if (h.le.5.) then
            r_keha = 20.
            s_keha = pi * r_keha**2
            p2     = c_keha * s_keha * 500. / m_keha
        end if

        if (h.le.0.) stop

        write(20, '(10f15.5)') tout, y, h, v*1000, a*1000
        write(21, '(10(f15.5,"_"))') tout, y, h, v*1000, a*1000
    end do

    close(20)
    close(21)

end program main

subroutine keplerf(t, y, dy)
    real   :: y(4), dy(4)
    common gm, r_maa, ro_atm_0, p1, p2
    
    r      = sqrt(y(1)**2 + y(2)**2)
    vv     = sqrt(y(3)**2 + y(4)**2)
    h      = r - r_maa
    ro_atm = ro_atm_0 * exp(-p1 * h)
    dy(1)  =  y(3)
    dy(2)  =  y(4)
    dy(3)  = -y(1)*gm/r**3 - p2*ro_atm*vv*y(3)
    dy(4)  = -y(2)*gm/r**3 - p2*ro_atm*vv*y(4)
    return
end subroutine keplerf
