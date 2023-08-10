program solbybisection
implicit none
        
    double precision, parameter :: xmin = 0.0, xmax = 2.0, max_err = 1E-08
    integer, parameter          :: iter = 10 ! , iter_max = 1000 
    double precision            :: sol, err, fcos
    integer                     :: iter_max
    external fcos

!    call bisect(xmin, xmax, fcos, sol, iter, err)
    call coolerBisect(xmin, xmax, fcos, sol, iter_max, max_err)   
     
    write (*, '(A, F10.8, A, ES12.6)') "sol = ", sol, " err = ", max_err
    write (*, '(A, I6)') "Iterations: ",  iter_max
end program solbybisection
    
function fcos(x)
    implicit none
    double precision :: fcos, x

    fcos = cos(x)

end function fcos
