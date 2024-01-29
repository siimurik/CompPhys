PROGRAM SVDExample
    INTEGER, PARAMETER :: np = 3
    REAL :: a(np, np), u(np, np), w(np), v(np, np), b(np), x(np)
    REAL :: wmax, wmin
    INTEGER :: i, j, n

    ! Set the size of the square matrix
    n = np

    ! Initialize matrix 'a' 
    a = reshape((/ 2.0,  3.0, -1.0, &
                   3.0,  5.0,  2.0, &
                   1.0, -1.0,  3.0 /),[np,np])
    ! Transpose the matrix
    a = transpose(a)

    ! Copy matrix 'a' into 'u'
    do i = 1, n
        do j = 1, n
            u(i, j) = a(i, j)
        end do
    end do

    ! Initialize vector b
    DATA b / 1.0, 1.0, 1.0 /  ! Adjust the values as needed

    ! Perform SVD on 'u'
    CALL svdcmp(u, n, n, np, np, w, v)

    ! Find the maximum singular value
    wmax = 0.0
    do j = 1, n
        if (w(j) .gt. wmax) wmax = w(j)
    end do

    ! Set the threshold for singular values
    wmin = wmax * 1.0e-6    ! 1.0D-12 for double

    ! Threshold singular values
    do j = 1, n
        if (w(j) < wmin) w(j) = 0.0
    end do

    ! Perform backsubstitution
    CALL svbksb(u, w, v, n, n, np, np, b, x)

    ! Print the results
    WRITE(*,*) 'Original Matrix A:'
    DO i = 1, n
        WRITE(*, '(6F10.6)') a(i, :)
    END DO

    WRITE(*,*) 'Singular Values:'
    WRITE(*, '(6F10.6)') w

    WRITE(*,*) 'Backsubstitution Result X:'
    WRITE(*, '(6F10.6)') x

END PROGRAM SVDExample


SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
    INTEGER, PARAMETER :: NMAX=500 !Maximum anticipated value of n.
    INTEGER m,mp,n,np
    REAL a(mp,np),v(np,np),w(np)
    ! USES pythag
    !Given a matrix a(1:m,1:n), with physical dimensions mp by np, this routine computes its
    !singular value decomposition, A = U · W · V T . The matrix U replaces a on output. The
    !diagonal matrix of singular values W is output as a vector w(1:n). The matrix V (not the
    !transpose V T ) is output as v(1:n,1:n).
    INTEGER i,its,j,jj,k,l,nm
    REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag
    g = 0.0! Householder reduction to bidiagonal form.
    scale = 0.0
    anorm = 0.0
    do i = 1, n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if(i.le.m)then
            do k=i,m
                scale=scale+abs(a(k,i))
            end do
            if(scale.ne.0.0)then
                do k=i,m
                    a(k,i)=a(k,i)/scale
                    s=s+a(k,i)*a(k,i)
                end do
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                do j=l,n
                    s=0.0
                    do k=i,m
                        s=s+a(k,i)*a(k,j)
                    end do
                    f=s/h
                    do k=i,m
                        a(k,j)=a(k,j)+f*a(k,i)
                    end do
                end do
                do k=i,m
                    a(k,i)=scale*a(k,i)
                end do
            end if
        end if
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if((i.le.m).and.(i.ne.n))then
            do k=l,n
                scale=scale+abs(a(i,k))
            end do
            if(scale.ne.0.0)then
                do k=l,n
                    a(i,k)=a(i,k)/scale
                    s=s+a(i,k)*a(i,k)
                end do
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                do k=l,n
                    rv1(k)=a(i,k)/h
                end do
                do j=l,m
                    s=0.0
                    do k=l,n
                        s=s+a(j,k)*a(i,k)
                    end do
                    do k=l,n
                        a(j,k)=a(j,k)+s*rv1(k)
                    end do
                end do
                do k=l,n
                    a(i,k)=scale*a(i,k)
                end do
            end if
        end if
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
    end do
    do i=n,1,-1 !Accumulation of right-hand transformations.
        if(i.lt.n)then
            if(g.ne.0.0)then
                do j=l,n !Double division to avoid possible underflow.
                    v(j,i)=(a(i,j)/a(i,l))/g
                end do
                do j=l,n
                    s=0.0
                    do k=l,n
                        s=s+a(i,k)*v(k,j)
                    end do
                    do k=l,n
                        v(k,j)=v(k,j)+s*v(k,i)
                    end do
                end do
            end if
            do j=l,n
                v(i,j)=0.0
                v(j,i)=0.0
            end do
        end if
        v(i,i)=1.0
        g=rv1(i)
        l=i
    end do
    do i=min(m,n),1,-1 !Accumulation of left-hand transformations.
        l=i+1
        g=w(i)
        do j=l,n
            a(i,j)=0.0
        end do
        if(g.ne.0.0)then
            g=1.0/g
            do j=l,n
                s=0.0
                do k=l,m
                    s=s+a(k,i)*a(k,j)
                end do
                f=(s/a(i,i))*g
                do k=i,m
                    a(k,j)=a(k,j)+f*a(k,i)
                end do
            end do
            do j=i,m
                a(j,i)=a(j,i)*g
            end do
        else
            do j= i,m
                a(j,i)=0.0
            end do
        endif
        a(i,i)=a(i,i)+1.0
    end do
    do k=n,1,-1 !Diagonalization of the bidiagonal form: Loop over
        !singular values, and over allowed iterations.
        do its=1,30
            do l=k,1,-1 !Test for splitting.
                nm=l-1 !Note that rv1(1) is always zero.
                if((abs(rv1(l))+anorm).eq.anorm) goto 2
                if((abs(w(nm))+anorm).eq.anorm) goto 1
            enddo
1           c=0.0 !Cancellation of rv1(l), if l > 1.
            s=1.0
            do i=l,k
                f=s*rv1(i)
                rv1(i)=c*rv1(i)
                if((abs(f)+anorm).eq.anorm) goto 2
                g=w(i)
                h=pythag(f,g)
                w(i)=h
                h=1.0/h
                c= (g*h)
                s=-(f*h)
                do  j=1,m
                    y=a(j,nm)
                    z=a(j,i)
                    a(j,nm)=(y*c)+(z*s)
                    a(j,i)=-(y*s)+(z*c)
                end do
            end do
2           z=w(k)
            if(l.eq.k)then !Convergence.
                if(z.lt.0.0)then !Singular value is made nonnegative.
                    w(k)=-z
                    do j=1,n
                        v(j,k)=-v(j,k)
                    end do
                end if
                goto 3
            end if
            if(its.eq.30) stop 'no convergence in svdcmp'
            x=w(l) !Shift from bottom 2-by-2 minor.
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            g=pythag(f,1.0)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c=1.0 !Next QR transformation:
            s=1.0
            do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=pythag(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f= (x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                do jj=1,n
                    x=v(jj,j)
                    z=v(jj,i)
                    v(jj,j)= (x*c)+(z*s)
                    v(jj,i)=-(x*s)+(z*c)
                end do
                z=pythag(f,h)
                w(j)=z !Rotation can be arbitrary if z = 0.
                if(z.ne.0.0)then
                    z=1.0/z
                    c=f*z
                    s=h*z
                endif
                f= (c*g)+(s*y)
                x=-(s*g)+(c*y)
                do jj=1,m
                    y=a(jj,j)
                    z=a(jj,i)
                    a(jj,j)= (y*c)+(z*s)
                    a(jj,i)=-(y*s)+(z*c)
                end do
            end do
            rv1(l)=0.0
            rv1(k)=f
            w(k)=x
        end do 
3       continue
    end do
    return
END SUBROUTINE

FUNCTION pythag(a,b)
    REAL :: a, b, pythag
    ! Computes (a2 + b2)1/2 without destructive underflow or overflow.
    REAL :: absa, absb
    absa = abs(a)
    absb = abs(b)
    if (absa .gt. absb) then
        pythag = absa*sqrt(1.0 + (absb/absa)**2)
    else
        if (absb .eq. 0.0) then 
        pythag = 0.0
    else
        pythag = absb*sqrt(1.0 + (absa/absb)**2)
        endif
    endif
    return
END FUNCTION pythag

SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
    INTEGER, PARAMETER :: NMAX=500 ! Maximum anticipated value of n.
    INTEGER :: m,mp,n,np
    REAL    :: b(mp),u(mp,np),v(np,np),w(np),x(np)
    !Solves A · X = B for a vector X, where A is specified by the arrays u, w, v as returned by
    !'svdcmp'. m and n are the logical dimensions of a, and will be equal for square matrices. 
    !mp and np are the physical dimensions of a. b(1:m) is the input right-hand side. x(1:n) is
    !the output solution vector. No input quantities are destroyed, so the routine may be called
    !sequentially with different b’s.
    INTEGER :: i,j,jj
    REAL    :: s,tmp(NMAX)
    do j=1,n    ! Calculate U^T B.
        s = 0.0
        if(w(j) .ne. 0.0) then ! Nonzero result only if wj is nonzero.
            do i = 1, m
                s = s + u(i,j)*b(i)
            end do
            s = s/w(j) ! This is the divide by wj .
        end if
        tmp(j) = s
    end do
    do j = 1, n ! Matrix multiply by V to get answer.
        s = 0.0
        do jj=1,n
            s = s + v(j,jj)*tmp(jj)
        end do
        x(j) = s
    end do
    return
END SUBROUTINE