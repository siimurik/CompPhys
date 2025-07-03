FUNCTION snrm(n, sx, itol)
    INTEGER n,itol,i,isamax
    DOUBLE PRECISION sx(n),snrm
    ! Compute one of two norms for a vector sx(1:n), as signaled by itol. Used by linbcg.
    if (itol.le.3)then
        snrm = 0.
        do i = 1,n !! Vector magnitude norm.
            snrm = snrm + sx(i)**2
        end do
        snrm = sqrt(snrm)
    else
        isamax = 1
        do i = 1, n !! Largest component norm.
            if(abs(sx(i)) .gt. abs(sx(isamax))) isamax=i
        end do
        snrm = abs(sx(isamax))
    end if
    return
END FUNCTION snrm