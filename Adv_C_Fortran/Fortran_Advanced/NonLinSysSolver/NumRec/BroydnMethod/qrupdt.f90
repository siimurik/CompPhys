SUBROUTINE qrupdt(r,qt,n,np,u,v)
    INTEGER n,np
    REAL r(np,np),qt(np,np),u(np),v(np)
    ! USES rotate
    ! Given the QR decomposition of some n × n matrix, calculates the QR decomposition of
    ! the matrix Q·(R + u⊗v). The matrices r and qt have physical dimension np. Note that
    ! QT is input and returned in qt.
    INTEGER i,j,k
    do k = n, 1, -1 !! Find largest k such that u(k) .ne. 0.
        if (u(k).ne.0.) goto 1
    enddo
    k=1
1   do i=k-1,1,-1       !! Transform R + u⊗v to upper Hes-
        call rotate(r,qt,n,np,i,u(i),-u(i+1)) !! senberg.
        if (u(i).eq.0.) then
            u(i) = abs(u(i+1))
        else if (abs(u(i)).gt.abs(u(i+1))) then
            u(i)=abs(u(i))*sqrt(1.+(u(i+1)/u(i))**2)
        else
            u(i)=abs(u(i+1))*sqrt(1.+(u(i)/u(i+1))**2)
        endif
    enddo
    do j = 1, n
        r(1,j) = r(1,j) + u(1)*v(j)
    enddo
    do i = 1, k-1                   !! Transform upper Hessenberg matrix
        call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i)) !! to upper triangular.
    enddo
    return
END SUBROUTINE