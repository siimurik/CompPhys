subroutine coolerBisect(xmin, xmax, func, sol, iter_max, max_err)
    implicit none
    
    double precision :: xmin, xmax, sol, func
    integer          :: iter_max
    double precision :: x1, x2, res, err, max_err
    integer          :: i
    external func
    
    x1 = xmin
    x2 = xmax
    err = (xmax - xmin) / 2.0

    ! Check if there could be a unique solution in interval.
    res = func(x1) * func(x2)
    if (res.gt.0.0) then
        write(*, *) 'There is either no solution or more than one solution in this interval.'
        write(*, *) 'Try again with a different interval.'
        return
    endif

    i = 0
    do while (err > max_err) !.and. i < iter_max)
        sol = (x1 + x2) / 2.0
        res = func(sol) * func(x2)        
        if (res.gt.0.0) then
            x2 = sol ! The solution is between x1 and sol, shrink x2
        else
            x1 = sol ! The solution is between sol and x2, increase x1
        endif
        err = (x2 - x1) / 2.0
        i = i + 1
    end do
    iter_max = i
end subroutine coolerBisect

