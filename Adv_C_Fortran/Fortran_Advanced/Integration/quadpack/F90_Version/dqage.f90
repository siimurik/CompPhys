    subroutine dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr, &
    neval,ier,alist,blist,rlist,elist,iord,last)
!***begin prologue  dqage
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
!***keywords  automatic integrator, general-purpose,
!             integrand examinator, globally adaptive,
!             gauss-kronrod
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral   i = integral of f over (a,b),
!            hopefully satisfying following claim for accuracy
!            abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
!***description

!        computation of a definite integral
!        standard fortran subroutine
!        double precision version

!        parameters
!         on entry
!            f      - double precision
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.

!            a      - double precision
!                     lower limit of integration

!            b      - double precision
!                     upper limit of integration

!            epsabs - double precision
!                     absolute accuracy requested
!            epsrel - double precision
!                     relative accuracy requested
!                     if  epsabs.le.0
!                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.

!            key    - integer
!                     key for choice of local integration rule
!                     a gauss-kronrod pair is used with
!                          7 - 15 points if key.lt.2,
!                         10 - 21 points if key = 2,
!                         15 - 31 points if key = 3,
!                         20 - 41 points if key = 4,
!                         25 - 51 points if key = 5,
!                         30 - 61 points if key.gt.5.

!            limit  - integer
!                     gives an upperbound on the number of subintervals
!                     in the partition of (a,b), limit.ge.1.

!         on return
!            result - double precision
!                     approximation to the integral

!            abserr - double precision
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)

!            neval  - integer
!                     number of integrand evaluations

!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier.gt.0 abnormal termination of the routine
!                             the estimates for result and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value
!                             of limit.
!                             however, if this yields no improvement it
!                             is rather advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined(e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the
!                             subranges. if possible, an appropriate
!                             special-purpose integrator should be used
!                             which is designed for handling the type of
!                             difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             (epsabs.le.0 and
!                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                             result, abserr, neval, last, rlist(1) ,
!                             elist(1) and iord(1) are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.

!            alist   - double precision
!                      vector of dimension at least limit, the first
!                       last  elements of which are the left
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)

!            blist   - double precision
!                      vector of dimension at least limit, the first
!                       last  elements of which are the right
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)

!            rlist   - double precision
!                      vector of dimension at least limit, the first
!                       last  elements of which are the
!                      integral approximations on the subintervals

!            elist   - double precision
!                      vector of dimension at least limit, the first
!                       last  elements of which are the moduli of the
!                      absolute error estimates on the subintervals

!            iord    - integer
!                      vector of dimension at least limit, the first k
!                      elements of which are pointers to the
!                      error estimates over the subintervals,
!                      such that elist(iord(1)), ...,
!                      elist(iord(k)) form a decreasing sequence,
!                      with k = last if last.le.(limit/2+2), and
!                      k = limit+1-last otherwise

!            last    - integer
!                      number of subintervals actually produced in the
!                      subdivision process

!***references  (none)
!***routines called  d1mach,dqk15,dqk21,dqk31,
!                    dqk41,dqk51,dqk61,dqpsrt
!***end prologue  dqage

    double precision :: a,abserr,alist,area,area1,area12,area2,a1,a2,b, &
    blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach, &
    epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f, &
    resabs,result,rlist,uflow
    integer :: ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval, &
    nrmax

    dimension alist(limit),blist(limit),elist(limit),iord(limit), &
    rlist(limit)

    external f

!            list of major variables
!            -----------------------

!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                      (alist(i),blist(i))
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest
!                       error estimate
!           errmax    - elist(maxerr)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision


!           machine dependent constants
!           ---------------------------

!           epmach  is the largest relative spacing.
!           uflow  is the smallest positive magnitude.

!***first executable statement  dqage
    epmach = d1mach(4)
    uflow = d1mach(1)

!           test on validity of parameters
!           ------------------------------

    ier = 0
    neval = 0
    last = 0
    result = 0.0d+00
    abserr = 0.0d+00
    alist(1) = a
    blist(1) = b
    rlist(1) = 0.0d+00
    elist(1) = 0.0d+00
    iord(1) = 0
    if(epsabs <= 0.0d+00 .AND. &
    epsrel < dmax1(0.5d+02*epmach,0.5d-28)) ier = 6
    if(ier == 6) go to 999

!           first approximation to the integral
!           -----------------------------------

    keyf = key
    if(key <= 0) keyf = 1
    if(key >= 7) keyf = 6
    neval = 0
    if(keyf == 1) call dqk15(f,a,b,result,abserr,defabs,resabs)
    if(keyf == 2) call dqk21(f,a,b,result,abserr,defabs,resabs)
    if(keyf == 3) call dqk31(f,a,b,result,abserr,defabs,resabs)
    if(keyf == 4) call dqk41(f,a,b,result,abserr,defabs,resabs)
    if(keyf == 5) call dqk51(f,a,b,result,abserr,defabs,resabs)
    if(keyf == 6) call dqk61(f,a,b,result,abserr,defabs,resabs)
    last = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1) = 1

!           test on accuracy.

    errbnd = dmax1(epsabs,epsrel*dabs(result))
    if(abserr <= 0.5d+02*epmach*defabs .AND. abserr > errbnd) ier = 2
    if(limit == 1) ier = 1
    if(ier /= 0 .OR. (abserr <= errbnd .AND. abserr /= resabs) &
     .OR. abserr == 0.0d+00) go to 60

!           initialization
!           --------------


    errmax = abserr
    maxerr = 1
    area = result
    errsum = abserr
    nrmax = 1
    iroff1 = 0
    iroff2 = 0

!           main do-loop
!           ------------

    do 30 last = 2,limit
    
    !           bisect the subinterval with the largest error estimate.
    
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        if(keyf == 1) call dqk15(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 2) call dqk21(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 3) call dqk31(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 4) call dqk41(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 5) call dqk51(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 6) call dqk61(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 1) call dqk15(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf == 2) call dqk21(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf == 3) call dqk31(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf == 4) call dqk41(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf == 5) call dqk51(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf == 6) call dqk61(f,a2,b2,area2,error2,resabs,defab2)
    
    !           improve previous approximations to integral
    !           and error and test for accuracy.
    
        neval = neval+1
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1 == error1 .OR. defab2 == error2) go to 5
        if(dabs(rlist(maxerr)-area12) <= 0.1d-04*dabs(area12) &
         .AND. erro12 >= 0.99d+00*errmax) iroff1 = iroff1+1
        if(last > 10 .AND. erro12 > errmax) iroff2 = iroff2+1
        5 rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
        if(errsum <= errbnd) go to 8
    
    !           test for roundoff error and eventually set error flag.
    
        if(iroff1 >= 6 .OR. iroff2 >= 20) ier = 2
    
    !           set error flag in the case that the number of subintervals
    !           equals limit.
    
        if(last == limit) ier = 1
    
    !           set error flag in the case of bad integrand behaviour
    !           at a point of the integration range.
    
        if(dmax1(dabs(a1),dabs(b2)) <= (0.1d+01+0.1d+03* &
        epmach)*(dabs(a2)+0.1d+04*uflow)) ier = 3
    
    !           append the newly-created intervals to the list.
    
        8 if(error2 > error1) go to 10
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 20
        10 alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
    
    !           call subroutine dqpsrt to maintain the descending ordering
    !           in the list of error estimates and select the subinterval
    !           with the largest error estimate (to be bisected next).
    
        20 call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
    ! ***jump out of do-loop
        if(ier /= 0 .OR. errsum <= errbnd) go to 40
    30 END DO

!           compute final result.
!           ---------------------

    40 result = 0.0d+00
    do 50 k=1,last
        result = result+rlist(k)
    50 END DO
    abserr = errsum
    60 if(keyf /= 1) neval = (10*keyf+1)*(2*neval+1)
    if(keyf == 1) neval = 30*neval+15
    999 return
    end subroutine dqage
