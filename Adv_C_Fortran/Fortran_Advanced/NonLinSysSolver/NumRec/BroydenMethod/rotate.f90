SUBROUTINE rotate(r,qt,n,np,i,a,b)
    INTEGER n,np,i
    REAL a,b,r(np,np),qt(np,np)
    ! Given n×n matrices r and qt of physical dimension np, carry out a Jacobi rotation on rows i
    ! and i + 1 of each matrix. a and b are the parameters of the rotation: cos θ = a/√a2 + b2,
    ! sin θ = b/√a2 + b2.
    INTEGER j
    REAL c,fact,s,w,y
    if(a.eq.0.)then !! Avoid unnecessary overflow or underflow.
        c = 0.
        s = sign(1.,b)
    else if (abs(a).gt.abs(b)) then
        fact = b/a
        c = sign(1./sqrt(1.+fact**2),a)
        s = fact*c
    else
        fact = a/b
        s    = sign(1./sqrt(1.+fact**2),b)
        c    = fact*s
    endif
    do j = i, n ! Premultiply r by Jacobi rotation.
        y = r(i,j)
        w = r(i+1,j)
        r(i,j)   = c*y - s*w
        r(i+1,j) = s*y + c*w
    enddo
    do j = 1, n !! Premultiply qt by Jacobi rotation.
        y=qt(i,j)
        w=qt(i+1,j)
        qt(i,j)=c*y-s*w
        qt(i+1,j)=s*y+c*w
    enddo
    return
END