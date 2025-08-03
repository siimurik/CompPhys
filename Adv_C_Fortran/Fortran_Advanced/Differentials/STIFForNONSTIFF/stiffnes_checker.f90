program stiffness_detector
    implicit none
    integer :: n, lda, ldvl, ldvr, lwork, info, i
    double precision, allocatable :: A(:,:), WR(:), WI(:), VL(:,:), VR(:,:), WORK(:)
    double precision :: max_re, min_re, ratio
    double precision :: R, L, C

    ! Set matrix size (you can change this)
    n = 2
    lda = n
    ldvl = n
    ldvr = n
    lwork = 4 * n

    ! Allocate arrays
    allocate(A(n,n), WR(n), WI(n), VL(n,n), VR(n,n), WORK(lwork))

    ! Define Jacobian matrix (replace with real system)
    !A = reshape([ -1000.0d0,   0.0d0,    0.0d0, &
    !                0.0d0,  -1.0d0,    0.0d0, &
    !                0.0d0,   0.0d0,  -0.1d0 ], shape(A))
    R = 100.D0
    C = 1.0D-3
    L = 1.D0

    A = reshape([ -1.d0/(R*C), 1.d0/C, &
                  -1.d0/L,   0.d0], shape(A))

    ! Compute eigenvalues and (optional) eigenvectors using LAPACK
    call dgeev('V', 'V', n, A, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, lwork, info)

    if (info /= 0) then
        print *, "Error computing eigenvalues: INFO =", info
        stop
    end if

    ! Calculate stiffness ratio (max(Re)/min(Re))
    max_re = abs(WR(1))
    min_re = abs(WR(1))
    do i = 2, n
        if (abs(WR(i)) > max_re) max_re = abs(WR(i))
        if (abs(WR(i)) < min_re .and. WR(i) /= 0.0d0) min_re = abs(WR(i))
    end do

    ratio = max_re / min_re

    print *, "Eigenvalue real parts: ", WR
    print *, "Stiffness ratio R = ", ratio

    if (ratio > 1000.0d0) then
        print *, ">>> This system is likely STIFF."
    else
        print *, ">>> This system is likely NOT stiff."
    end if

    ! Deallocate arrays (clean-up)
    deallocate(A, WR, WI, VL, VR, WORK)
end program stiffness_detector
