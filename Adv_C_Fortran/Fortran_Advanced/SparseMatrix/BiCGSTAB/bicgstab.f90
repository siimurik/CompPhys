PROGRAM BiconjugateGradientExample
    implicit none
    INTEGER, PARAMETER :: np = 4
    DOUBLE PRECISION :: tol, threshold, err
    INTEGER :: itol, itmax, iter
    DOUBLE PRECISION, DIMENSION(np, np) :: a
    DOUBLE PRECISION, DIMENSION(np) :: b, x
    integer :: n, nmax, i
    
    ! Set up the linear system: A * x = b
    DOUBLE PRECISION, DIMENSION(np*np) :: sa
    INTEGER, DIMENSION(np*np) :: ija

    ! Declare the COMMON block variables with fixed-size arrays
    !COMMON /mat/ sa, ija
    
    n = np
    
    ! Define the matrix A in sparse format
    a = reshape((/  4.D0, 2.D0, 0.D0, 1.D0, &
                    3.D0, 0.D0, 0.D0, 2.D0, &
                    0.D0, 1.D0, 1.D0, 1.D0, &
                    0.D0, 2.D0, 1.D0, 0.D0 /),[np,np])
  
    ! Set the right-hand side vector b
    DATA b /-1.D0, -0.5D0, -1.D0, 2.D0/
  
    ! Set the convergence parameters
    x = 1.D0
    x = matmul(A, x)
    itol = 1   ! Residual norm
    tol = 1.0D-6
    threshold = tol
    itmax = 100
    nmax = itmax

    ! Call the dsprsin subroutine to convert a dense matrix into sparse storage
    !CALL dsprsin(a, n, N, threshold, nmax, sa, ija)

    ! Call the linbcg subroutine
    CALL linbcg(n, b, x, itol, tol, itmax, iter, err)
    ! Define the right-hand side vector b
    !b = [-1, -0.5, -1, 2]

    ! Call the bicgstab subroutine
    !call bicgstab(zeros(4), b, matvec, [1.0e-5, 1000], x, exit_code)

    ! Print the result
    print *, "Solution x:"
    do i = 1, size(x)
        print *, x(i)
    end do
  
  END PROGRAM BiconjugateGradientExample
  

SUBROUTINE linbcg(n,b,x,itol,tol,itmax,iter,err)
    INTEGER            :: iter,itmax,itol,n
    DOUBLE PRECISION   :: err,tol,b(*),x(*) !Double precision is a good idea in this rou-
    INTEGER, PARAMETER :: NMAX=1024         !tine.
    DOUBLE PRECISION   :: EPS=1.d-14        
    !USES atimes,asolve,snrm
    !Solves A · x = b for x(1:n), given b(1:n), by the iterative biconjugate gradient method.
    !On input x(1:n) should be set to an initial guess of the solution (or all zeros); itol is
    !1,2,3, or 4, specifying which convergence test is applied (see text); itmax is the maximum
    !number of allowed iterations; and tol is the desired convergence tolerance. On output,
    !x(1:n) is reset to the improved solution, iter is the number of iterations actually taken,
    !and err is the estimated error. The matrix A is referenced only through the user-supplied
    !routines atimes, which computes the product of either A or its transpose on a vector; and
    !e · x = b or A
    !e T · x = b for some preconditioner matrix A
    !e (possibly
    !asolve, which solves A
    !the trivial diagonal part of A).
    INTEGER j
    DOUBLE PRECISION :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm
    DOUBLE PRECISION :: p(NMAX),pp(NMAX),r(NMAX),rr(NMAX),z(NMAX),zz(NMAX),snrm
    iter = 0 !Calculate initial residual.
    call atimes(n,x,r,0) !Input to atimes is x(1:n), output is r(1:n);
    do j=1,n!the final 0 indicates that the matrix (not
        r(j)=b(j)-r(j) !its transpose) is to be used.        
        rr(j)=r(j)
    end do
    ! call atimes(n,r,rr,0)
    !Uncomment this line to get the “minimum
    if(itol.eq.1) then  !residual” variant of the algorithm.
        bnrm = snrm(n,b,itol)
        call asolve(n,r,z,0) ! Input to asolve is r(1:n), output is z(1:n);
    else if (itol.eq.2) then !the final 0 indicates that the matrix A
        call asolve(n,b,z,0) !(not its transpose) is to be used.
        bnrm=snrm(n,z,itol)
        call asolve(n,r,z,0)
    else if (itol.eq.3.or.itol.eq.4) then
        call asolve(n,b,z,0)
        bnrm=snrm(n,z,itol)
        call asolve(n,r,z,0)
        znrm=snrm(n,z,itol)
    else
        stop 'illegal itol in linbcg'
    end if
    100 if (iter.le.itmax) then !Main loop.
        iter=iter+1         !eT .
        call asolve(n,rr,zz,1)!Final 1 indicates use of transpose matrix A
        bknum=0.d0
        do j=1,n!Calculate coefficient bk and direction vectors
            bknum=bknum+z(j)*rr(j) !p and pp.
        end do
        if(iter.eq.1) then
            do j=1,n
                p(j)=z(j)
                pp(j)=zz(j)
            end do
        else
            bk=bknum/bkden
            do j=1,n
                p(j)=bk*p(j)+z(j)
                pp(j)=bk*pp(j)+zz(j)
            end do
        end if
        bkden=bknum!Calculate coefficient ak, new iterate x, and
        call atimes(n,p,z,0)!new residuals r and rr.
        akden = 0.d0
        do j=1,n
            akden=akden+z(j)*pp(j)
        end do
        ak=bknum/akden
        call atimes(n,pp,zz,1)
        do j=1,n
            x(j)=x(j)+ak*p(j)
            r(j)=r(j)-ak*z(j)
            rr(j)=rr(j)-ak*zz(j)
        end do  !e ·z = r and check stopping criterion.
        call asolve(n,r,z,0) !Solve A
        if(itol.eq.1)then
            err=snrm(n,r,itol)/bnrm
        else if(itol.eq.2)then
            err=snrm(n,z,itol)/bnrm
        else if(itol.eq.3.or.itol.eq.4)then
            zm1nrm=znrm
            znrm=snrm(n,z,itol)
            if(abs(zm1nrm-znrm).gt.EPS*znrm) then
                dxnrm=abs(ak)*snrm(n,p,itol)
                err=znrm/abs(zm1nrm-znrm)*dxnrm
            else
                err=znrm/bnrm!Error may not be accurate, so loop again.
                goto 100
            end if
            xnrm=snrm(n,x,itol)
            if(err.le.0.5d0*xnrm) then
                err=err/xnrm
            else
                err=znrm/bnrm!Error may not be accurate, so loop again.
                goto 100
            end if
        end if
        write (*,*) ' iter=',iter,' err=',err
    if(err.gt.tol) goto 100
    end if
    return
END SUBROUTINE


SUBROUTINE dsprsax(sa,ija,x,b,n)
    INTEGER n,ija(*)
    DOUBLE PRECISION b(n),sa(*),x(n)
    !Multiply a matrix in row-index sparse storage arrays sa and ija by a vector x(1:n), giving
    !a vector b(1:n).
    INTEGER i,k
    if (ija(1).ne.n+2) stop 'mismatched vector and matrix in sprsax'
    do i=1,n
        b(i)=sa(i)*x(i) !Start with diagonal term.
        do k=ija(i),ija(i+1)-1 ! Loop over off-diagonal terms.
            b(i)=b(i)+sa(k)*x(ija(k))
        end do
    end do
    return
END SUBROUTINE dsprsax

SUBROUTINE dsprstx(sa,ija,x,b,n)
    INTEGER n,ija(*)
    DOUBLE PRECISION b(n),sa(*),x(n)
    !Multiply the transpose of a matrix in row-index sparse storage arrays sa and ija by a
    !vector x(1:n), giving a vector b(1:n).
    INTEGER i,j,k
    if (ija(1).ne.n+2) stop 'mismatched vector and matrix in sprstx'
    do i=1,n!Start with diagonal terms.
        b(i)=sa(i)*x(i)
    end do
    do i=1,n!Loop over off-diagonal terms.
        do k=ija(i),ija(i+1)-1
            j=ija(k)
            b(j)=b(j)+sa(k)*x(i)
        end do
    end do
    return
END SUBROUTINE dsprstx

SUBROUTINE atimes(n, x, r, itrnsp)
    INTEGER :: n,itrnsp,ija,NMAX
    DOUBLE PRECISION :: x(n),r(n), sa
    PARAMETER (NMAX=1000)
    COMMON /mat/ sa(NMAX),ija(NMAX) ! The matrix is stored somewhere.
    !USES dsprsax,dsprstx
    if (itrnsp.eq.0) then
        call dsprsax(sa,ija,x,r,n)
    else
        call dsprstx(sa,ija,x,r,n)
    end if
    return
END SUBROUTINE atimes

SUBROUTINE asolve(n, b, x, itrnsp)
    INTEGER :: n,itrnsp,ija,NMAX,i
    DOUBLE PRECISION :: x(n), b(n), sa
    PARAMETER (NMAX=1000)
    COMMON /mat/ sa(NMAX),ija(NMAX) ! The matrix is stored somewhere.
    do i=1,n
        x(i)=b(i)/sa(i)
    end do
    return
END SUBROUTINE asolve

FUNCTION snrm(n, sx, itol)
    INTEGER n,itol,i,isamax
    DOUBLE PRECISION sx(n),snrm
    ! Compute one of two norms for a vector sx(1:n), as signaled by itol. Used by linbcg.
    if (itol.le.3)then
        snrm=0.
        do i=1,n!Vector magnitude norm.
            snrm=snrm+sx(i)**2
        end do
        snrm=sqrt(snrm)
    else
        isamax=1
        do i=1,n !Largest component norm.
            if(abs(sx(i)) .gt. abs(sx(isamax))) isamax=i
        end do
        snrm=abs(sx(isamax))
    end if
    return
END FUNCTION snrm

SUBROUTINE dsprsin(a,n,np,thresh,nmax,sa,ija)
    INTEGER n,nmax,np,ija(nmax)
    DOUBLE PRECISION thresh,a(np,np),sa(nmax)
    !Converts a square matrix a(1:n,1:n) with physical dimension np into row-indexed sparse
    !storage mode. Only elements of a with magnitude ≥thresh are retained. Output is in
    !two linear arrays with physical dimension nmax (an input parameter): sa(1:) contains
    !array values, indexed by ija(1:). The logical sizes of sa and ija on output are both
    !ija(ija(1)-1)-1 (see text).
    INTEGER i,j,k
    do j=1,n !Store diagonal elements.
        sa(j)=a(j,j)
    end do
    ija(1)=n+2  !Index to 1st row off-diagonal element, if any.
    k=n+1
    do i=1,n !Loop over rows.
        do j=1,n !Loop over columns.
            if(abs(a(i,j)).ge.thresh)then
                if(i.ne.j)then
                !Store off-diagonal elements and their columns.
                k=k+1
                if(k.gt.nmax) stop 'nmax too small in sprsin'
                sa(k)=a(i,j)
                ija(k)=j
                end if
            end if
        end do
        ija(i+1)=k+1
        !As each row is completed, store index to next.
    end do
    return
END SUBROUTINE dsprsin