subroutine bisect(xmin,xmax,func,sol,iter,err)
        implicit none
        
        ! On entry this routine must be provided with the interval [xmin,xmax]
        ! over which to search for a solution,
        ! the name of the function fund to search for
        ! and the number of iterations iter to apply.
        ! On exit this routine returns the solution in sol, and its error in err.
        
        real    :: xmin,xmax,sol,func
        integer :: iter
        real    :: x1,x2,res,err
        integer :: i
        external func
        
        x1 = xmin
        x2 = xmax
        ! Check if there could be a unique solution in interval.
        res = func(x1)*func(x2)
        if(res.gt.0.0) then
                write(*,*) 'There is either no solution or more than one solution in this interval.'
                write(*,*) 'Try again with a different interval.'
        endif
        do i = 1, iter
                sol = (x1+x2)/2.0 ! Calculate the mid-point of the interval
                res = func(sol)*func(x2)
                if (res.gt.0.0) then
                        x2 = sol ! The solution is between x1 and sol, shrink x2
                else
                        x1 = sol ! The solution is between sol and x2, increase x1
                endif
        enddo
        sol = (x2+x1)/2.0
        err = (x2-x1)/2.0
end subroutine bisect
