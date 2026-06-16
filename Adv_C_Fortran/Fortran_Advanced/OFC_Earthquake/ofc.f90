! Compile: gfortran -fopenmp -O2 ofc.f90 -o ofc
! Run:     ./ofc_omp
! Outputs:
!   ofc_stress.dat      - final L x L stress field (space-separated, readable by np.loadtxt)
!   ofc_avalanches.dat  - one integer per earthquake: total sites that fired
!   ofc_footprint.dat   - L x L binary mask (0/1) of the largest avalanche recorded
!
! Author: Siim Erik Pugal, June 2026

program ofc_model_openmp
    use omp_lib
    implicit none

    ! ------------------------------------------------------------------ !
    ! Parameters - tune these                                            !
    ! ------------------------------------------------------------------ !
    integer,  parameter :: L        = 256       ! grid side 
    integer,  parameter :: n_quakes = 25000000  ! need >> L^2 for SOC
    real,     parameter :: alpha     = 0.20   ! stress transf frac (keep < 0.25)
    real,     parameter :: threshold = 1.0    ! firing threshold

    ! ------------------------------------------------------------------ !
    ! Variables                                                          !
    ! ------------------------------------------------------------------ !
    real,    allocatable :: grid(:,:)         ! stress field
    real,    allocatable :: shed(:,:)         ! stress-to-shed buffer
    integer, allocatable :: footprint(:,:)    ! binary mask of largest avalanche
    integer, allocatable :: best_fp(:,:)      ! running best footprint

    real    :: max_stress
    integer :: i, j, step
    logical :: earthquake_active

    ! Avalanche diagnostics
    integer :: avalanche_size           ! sites fired in current earthquake
    integer :: wave_count               ! sites fired in current relaxation wave
    integer :: best_size                ! largest avalanche seen so far

    real    :: elapsed_seconds
    integer :: start_time, end_time, elapsed_time, rate

    ! ------------------------------------------------------------------ !
    ! Allocation                                                         !
    ! ------------------------------------------------------------------ !
    allocate(grid(L,L), shed(L,L), footprint(L,L), best_fp(L,L))

    ! ------------------------------------------------------------------ !
    ! Initialisation                                                     !
    ! ------------------------------------------------------------------ !
    call random_seed()
    call random_number(grid)   ! uniform random stress in [0,1)

    footprint = 0
    best_fp   = 0
    best_size = 0

    ! Open avalanche log (written incrementally - one line per quake)
    open(unit=20, file='ofc_avalanches.dat', status='replace')

    print "(A,I4,A,I4,A,I9,A,I3,A)", "OFC model -", L, " x", L, &
            " grid,", n_quakes, " quakes,", &
            omp_get_max_threads(), " OpenMP threads."
    print "(A, F4.2, A, F4.2)", "alpha = ", alpha, &
                        "  threshold = ", threshold

    ! ------------------------------------------------------------------ !
    ! Main simulation loop                                               !
    ! ------------------------------------------------------------------ !
    call system_clock(count=start_time, count_rate=rate)
    do step = 1, n_quakes

        ! ---- 1. Global driving ---------------------------------------- !
        ! Raise every site uniformly so the maximum just reaches threshold.
        max_stress = maxval(grid)

        !$omp parallel do private(i,j) schedule(static)
        do j = 1, L
            do i = 1, L
                grid(i,j) = grid(i,j) + (threshold - max_stress)
            end do
        end do
        !$omp end parallel do

        ! ---- 2. Avalanche relaxation ----------------------------------- !
        avalanche_size   = 0
        footprint        = 0
        earthquake_active = .true.

        do while (earthquake_active)
            earthquake_active = .false.
            wave_count        = 0
            shed              = 0.0

            ! Phase 1: identify and reset all sites at or above threshold.
            ! wave_count accumulates across threads via reduction.
            !$omp parallel do private(i,j) &
            !$omp& reduction(.or.:earthquake_active) &
            !$omp& reduction(+:wave_count) schedule(static)
            do j = 1, L
                do i = 1, L
                    if (grid(i,j) >= threshold) then
                        shed(i,j)     = grid(i,j)
                        grid(i,j)     = 0.0
                        footprint(i,j) = 1
                        earthquake_active = .true.
                        wave_count        = wave_count + 1
                    end if
                end do
            end do
            !$omp end parallel do

            avalanche_size = avalanche_size + wave_count

            ! Phase 2: redistribute stress to the four neighbours.
            ! Atomic updates prevent race conditions when multiple sites
            ! shed stress into the same neighbour simultaneously.
            if (earthquake_active) then
                !$omp parallel do private(i,j) schedule(static)
                do j = 1, L
                    do i = 1, L
                        if (shed(i,j) > 0.0) then
                            if (i > 1) then
                                !$omp atomic
                                grid(i-1,j) = grid(i-1,j) + alpha * shed(i,j)
                            end if
                            if (i < L) then
                                !$omp atomic
                                grid(i+1,j) = grid(i+1,j) + alpha * shed(i,j)
                            end if
                            if (j > 1) then
                                !$omp atomic
                                grid(i,j-1) = grid(i,j-1) + alpha * shed(i,j)
                            end if
                            if (j < L) then
                                !$omp atomic
                                grid(i,j+1) = grid(i,j+1) + alpha * shed(i,j)
                            end if
                        end if
                    end do
                end do
                !$omp end parallel do
            end if

        end do  ! end relaxation while loop

        ! ---- 3. Record diagnostics ------------------------------------ !
        write(20,*) avalanche_size

        ! Keep the footprint of the largest avalanche seen so far
        if (avalanche_size > best_size) then
            best_size = avalanche_size
            best_fp   = footprint
        end if

        ! Progress report every 10 % of n_quakes
        if (mod(step, n_quakes/10) == 0) then
            print '(A,I8,A,I8,A,I6)', &
                "  step ", step, " / ", n_quakes, &
                "   largest so far:", best_size
        end if

    end do  ! end main driving loop

    call system_clock(count=end_time)
    elapsed_time = end_time - start_time
    elapsed_seconds = real(elapsed_time) / real(rate)

    close(20)
    print *, "Avalanche log written to ofc_avalanches.dat"

    ! ------------------------------------------------------------------ !
    ! Output 1: final stress field                                       !
    ! Row i = i-th row of the Fortran array.                             !
    ! np.loadtxt('ofc_stress.dat') will give a (L,L) array directly.     !
    ! ------------------------------------------------------------------ !
    print *, "Writing stress field to ofc_stress.dat..."
    open(unit=10, file='ofc_stress.dat', status='replace')
    do i = 1, L
        write(10,*) (grid(i,j), j = 1, L)   ! implied-do: space-separated row
    end do
    close(10)

    ! ------------------------------------------------------------------ !
    ! Output 2: largest-avalanche footprint (binary 0/1 matrix)          !
    ! ------------------------------------------------------------------ !
    print *, "Writing avalanche footprint to ofc_footprint.dat..."
    open(unit=30, file='ofc_footprint.dat', status='replace')
    do i = 1, L
        write(30,*) (best_fp(i,j), j = 1, L)
    end do
    close(30)

    print '(A,I0,A,F5.2,A)', &
        "Done. Largest avalanche: ", best_size, " sites  (", &
        100.0*best_size/(L*L), "% of grid)"

    print *, ""
    print "(A, F10.3, A)", "Elapsed time:", elapsed_seconds, " s."
    ! ------------------------------------------------------------------ !
    ! Clean up                                                           !
    ! ------------------------------------------------------------------ !
    deallocate(grid, shed, footprint, best_fp)

end program ofc_model_openmp