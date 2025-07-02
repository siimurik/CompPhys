FUNCTION fmin(x)
    INTEGER n,NP
    REAL fmin,x(*),fvec
    PARAMETER (NP=40)
    COMMON /newtv/ fvec(NP),n
    SAVE /newtv/
    ! USES funcv
    ! Returns f = 1/2 F·F at x. SUBROUTINE funcv(n,x,f) is a fixed-name, user-supplied
    ! routine that returns the vector of functions at x. The common block newtv communicates
    ! the function values back to newt.
    INTEGER i
    REAL sum
    call funcv(n,x,fvec)
    sum=0.
    do i=1,n
        sum=sum+fvec(i)**2
    enddo
    fmin=0.5*sum
    return
END FUNCTION