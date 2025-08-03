program main
    implicit none
    external f1, jac1
    integer :: i, iopar, iopt, iout, istate, itask, itol, iwork,    &
       leniw, lenrw, liw, lout, lrw, mband, meth, mf, miter,        &
       ml, mu, neq, nerr, nfe, nfea, nje, nout, nqu, nst             
    double precision atol, dtout, er, erm, ero, hu, rtol, rwork, t, &
       tout, tout1, y                                               
    dimension y(25), rwork(697), iwork(45)                          
    data lout/6/, tout1/1.39283880203d0/, dtout/2.214773875d0/      

    nerr = 0
    itol = 1
    rtol = 0.0d0
    atol = 1.0d-6
    lrw = 697
    liw = 45
    iopt = 0
    !----------------
    ! First problem
    !----------------
    neq = 2
    nout = 4
    write (lout,110) neq,itol,rtol,atol
    110  format(/' Demonstration program for DLSODE package'/// &
            ' Problem 1:  Van der Pol oscillator:'/             &
            '  xdotdot - 3*(1 - x**2)*xdot + x = 0, ',          &
            '   x(0) = 2, xdot(0) = 0'/                         &
            ' neq =',i2/                                        &
            ' itol =',i3,'   rtol =',d10.1,'   atol =',d10.1//)
    !
        do 195 meth = 1,2
        do 190 miter = 0,3
        mf = 10*meth + miter
        write (lout,120) mf
    120  format(///' Solution with mf =',i3//   &
             5x,'t               x               xdot       nq      h'//)
        t = 0.0d0
        y(1) = 2.0d0
        y(2) = 0.0d0
        itask = 1
        istate = 1
        tout = tout1
        ero = 0.0d0
        do iout = 1,nout
            call dlsode(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,    &
            iopt,rwork,lrw,iwork,liw,jac1,mf)
            hu = rwork(11)
            nqu = iwork(14)
            write (lout,140) t,y(1),y(2),nqu,hu
    140    format(d15.5,d16.5,d14.3,i5,d14.3)
            if (istate .lt. 0) go to 175
            iopar = iout - 2*(iout/2)
            if (iopar .ne. 0) go to 170
            er = abs(y(1))/atol
            ero = max(ero,er)
            if (er .gt. 1000.0d0) then
            write (lout,150)
    150      format(//' Warning: error exceeds 1000 * tolerance'//)
            nerr = nerr + 1
            endif
    170    tout = tout + dtout
        end do
    175  continue
        if (istate .lt. 0) nerr = nerr + 1
        nst = iwork(11)
        nfe = iwork(12)
        nje = iwork(13)
        lenrw = iwork(17)
        leniw = iwork(18)
        nfea = nfe
        if (miter .eq. 2) nfea = nfe - neq*nje
        if (miter .eq. 3) nfea = nfe - nje
        write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
    180  format(//' Final statistics for this run:'/    &
        ' rwork size =',i4,'   iwork size =',i4/     &
        ' number of steps =',i5/                     &
        ' number of f-s   =',i5/                     &
        ' (excluding J-s) =',i5/                     &
        ' number of J-s   =',i5/                     &
        ' error overrun =',d10.2)
    190  continue
    195  continue

end program main

subroutine f1 (neq, t, y, ydot)
    implicit none
    integer neq
    double precision t, y, ydot
    dimension y(neq), ydot(neq)
    ydot(1) = y(2)
    ydot(2) = 3.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
    return
end subroutine f1

subroutine jac1 (neq, t, y, ml, mu, pd, nrowpd)
    implicit none
    integer neq, ml, mu, nrowpd
    double precision t, y, pd
    dimension y(neq), pd(nrowpd,neq)
    pd(1,1) = 0.0d0
    pd(1,2) = 1.0d0
    pd(2,1) = -6.0d0*y(1)*y(2) - 1.0d0
    pd(2,2) = 3.0d0*(1.0d0 - y(1)*y(1))
    return
end subroutine jac1