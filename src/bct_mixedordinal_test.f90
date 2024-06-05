


subroutine estimate_bct_ordinal_test(postZmean, postZcov, P, numcorr, K, numG, BHat, sdHat, CHat, XtXi, samsize0, &
    burnin, Ntot, Njs, Xgroups, Ygroups, C_quantiles, sigma_quantiles, B_quantiles, BDrawsStore, &
    sigmaDrawsStore, CDrawsStore, sdMH, ordinal_in, Cat_in, maxCat, gLiuSab, seed, nuggetscale)
!
    implicit none

    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    integer(i6), intent(in) :: P, numcorr, K, numG, samsize0, Ntot, seed, maxCat, Njs(numG), burnin
    real(r15), intent(in)   :: BHat(numG,K,P), sdHat(numG,P), CHat(numG,P,P), XtXi(numG,K,K), Ygroups(numG,Ntot,P), &
                               Xgroups(numG,Ntot,K), sdMH(numG,P), ordinal_in(numG,P), nuggetscale, &
                               Cat_in(numG,P)
    real(r15), intent(out)  :: postZmean(numcorr,1), postZcov(numcorr,numcorr), B_quantiles(numG,K,P,3), &
                               C_quantiles(numG,P,P,3), sigma_quantiles(numG,P,3), BDrawsStore(samsize0,numG,K,P), &
                               sigmaDrawsStore(samsize0,numG,P), CDrawsStore(samsize0,numG,P,P), gLiuSab(samsize0,numG,P)
    real(r15)               :: BDraws(numG,K,P), sigmaDraws(numG,P), CDraws(numG,P,P), Ccan(P,P), Ds(P,P), dummyPP2(P,P), &
                            CcanInv(P,P), SS1(P,P), rnunif(1), errorMatj(P,P), & !Ccurr(P,P), CcurrInv(P,P),
                            sigma_can(P), aa, bb, SigmaMat(P,P), R_MH, epsteps(P,P), diffmat(Ntot,P), & !logR_MH,
                            varz1, varz2, varz1z2Plus, varz1z2Min, Cinv(P,P), Zcorr_sample(samsize0,numcorr), &
                            acceptSigma(numG,P), covBeta(P*K,P*K), betaDrawj(1,P*K), dummyPP(P,P), &
                            dummy3(samsize0), dummy2(samsize0), meanO(P*K), para(((P*K)*((P*K)+3)/2 + 1)), SS2(P,P), &
                            gLiuSab_curr(numG,P), Cnugget(P,P), acceptLS(numG,P), sdMHg(numG,P), Wgroups(numG,Ntot,P), &
                            Wdummy(numG,P,Ntot,maxCat), alphaMat(numG,maxCat+1,P)
    integer(i6)             :: s1, g1, i1, corrteller, c1, c2, p1, p2, k1, errorflag, lower_int, median_int, upper_int, &
                               iseed, ordinal(numG,P), Cat(numG,P)
!
    !set seed
    iseed = seed
!
    !initial posterior draws start at MLEs
    BDraws = BHat
    sigmaDraws = sdHat(1:numG,:)
    CDraws = CHat
    meanO = 0
    gLiuSab = 0
    sigmaDrawsStore = 0
    gLiuSab_curr = 1.0
    !

!
!
    !start Gibbs sampler
    do s1 = 1,samsize0

        corrteller = 0

        do g1 = 1,numG

            !draw B
            SigmaMat = matmul(matmul(diag(sigmaDraws(g1,:),P),CDraws(g1,:,:)),diag(sigmaDraws(g1,:),P))
            call kronecker(K,P,XtXi(g1,:,:),SigmaMat,covBeta)

            call setgmn(meanO,covBeta,P*K,para)
            call GENMN(para,betaDrawj(1,1:(P*K)),P*K,iseed)
            do p1 = 1,P
                BDraws(g1,:,p1) = betaDrawj(1,((p1-1)*K+1):(p1*K)) + BHat(g1,:,p1)
            end do
            !BDraws = BHat
            !write(*,*)'burnin BHat',s1,g1,2
            !write(*,*)BHat(1,1,1)
            !write(*,*)'burnin BDraws',s1,g1,3
            !write(*,*)BDraws(g1,1,:)

            !draw candidate draw for the correlation matrix
            diffmat(1:Njs(g1),1:P) = Ygroups(g1,1:Njs(g1),1:P) - matmul(Xgroups(g1,1:Njs(g1),1:K),BDraws(g1,1:K,1:P))
            !write(*,*)'diffmat'
            !write(*,*)diffmat(1,:)
            errorMatj = matmul(transpose(diffmat(1:Njs(g1),1:P)),diffmat(1:Njs(g1),1:P))
            Ds = diag(1/sqrt(diagonals(errorMatj,P)),P)
            diffmat(1:Njs(g1),1:P) = matmul(diffmat(1:Njs(g1),1:P),Ds) !diffmat is now epsilon in LD
            epsteps = matmul(transpose(diffmat(1:Njs(g1),1:P)),diffmat(1:Njs(g1),1:P))
            SS1 = matmul(matmul(diag(1/sigmaDraws(g1,:),P),epsteps),diag(1/sigmaDraws(g1,:),P))
            call FINDInv(SS1,SS2,P,errorflag)
            call gen_wish(SS2,Njs(g1)-P-1,dummyPP2,P,iseed)
            call FINDInv(dummyPP2,dummyPP,P,errorflag)
            Ccan = matmul(matmul(diag(1/sqrt(diagonals(dummyPP,P)),P),dummyPP),diag(1/sqrt(diagonals(dummyPP,P)),P))
            !write(*,*)'Ccan'
            !write(*,*)Ccan
            call FINDInv(Ccan,CcanInv,P,errorflag)
            Cinv = CcanInv
            CDraws(g1,:,:) = Ccan(:,:)
            do i1 = 1,P-1 !keep Fisher z transformed posterior draws of rho's
                Zcorr_sample(s1,(corrteller+1):(corrteller+P-i1)) = .5*log((1+CDraws(g1,(1+i1):P,i1))/ &
                    (1-CDraws(g1,(1+i1):P,i1)))
                corrteller = corrteller + (P - i1)
            end do
!
            do p1 = 1,P
                bb = sum(errorMatj(p1,:)*Cinv(p1,:)/sigmaDraws(g1,:)) - errorMatj(p1,p1)*Cinv(p1,p1)/sigmaDraws(g1,p1)
                aa = Cinv(p1,p1)*errorMatj(p1,p1)
                sigma_can(p1) = rnormal(iseed)
                sigma_can(p1) = sigma_can(p1)*sdMH(g1,p1) + sigmaDraws(g1,p1) !random walk
                R_MH = exp((-real(Njs(g1))+1.0)*(log(sigma_can(p1))-log(sigmaDraws(g1,p1)) ) &
                       -.5*aa*(1.0/sigma_can(p1)**2 - 1.0/sigmaDraws(g1,p1)**2) &
                       -bb*(1.0/sigma_can(p1) - 1.0/sigmaDraws(g1,p1)))
                !call random_number(rnunif)
                rnunif = runiform ( iseed )
                if(rnunif(1) < R_MH .and. sigma_can(p1)>0) then
                    sigmaDraws(g1,p1) = sigma_can(p1)
                    acceptSigma(g1,p1) = acceptSigma(g1,p1) + 1
                end if
            end do
        end do
!
        !store posterior draws
        BDrawsStore(s1,1:numG,1:K,1:P) = BDraws(1:numG,1:K,1:P)
        sigmaDrawsStore(s1,1:numG,1:P) = sigmaDraws(1:numG,1:P)
        CDrawsStore(s1,1:numG,1:P,1:P) = CDraws(1:numG,1:P,1:P)
    end do
!
    ! compute posterior mean
    do c1=1,numcorr
        dummy2(:) = Zcorr_sample(:,c1)
        dummy3 = dummy2
        call piksrt(samsize0,dummy3)
        postZmean(c1,1) = dummy3(int(samsize0*.5))
    end do
!
    ! compute posterior covariance matrix
    do c1=1,numcorr
        do c2=c1,numcorr
            call robust_covest(samsize0, Zcorr_sample(1:samsize0,c1), Zcorr_sample(1:samsize0,c2), postZmean(c1,1), &
                postZmean(c2,1), varz1, varz2, varz1z2Plus, varz1z2Min)
            postZcov(c1,c2) = (varz1*varz2)**.5 * (varz1z2Plus - varz1z2Min)/(varz1z2Plus + varz1z2Min)
            postZcov(c2,c1) = postZcov(c1,c2)
        end do
    end do
!
    ! compute posterior quantiles
    lower_int = int(samsize0*.025)
    median_int = int(samsize0*.5)
    upper_int = int(samsize0*.975)
    C_quantiles = 0
    do g1=1,numG
        do p1=1,P
            !for the sigma's
            dummy2(:) = sigmaDrawsStore(:,g1,p1)
            dummy3=dummy2
            call piksrt(samsize0,dummy3)
            sigma_quantiles(g1,p1,1) = dummy3(lower_int)
            sigma_quantiles(g1,p1,2) = dummy3(median_int)
            sigma_quantiles(g1,p1,3) = dummy3(upper_int)
            !for the beta coefficients
            do k1=1,K
                dummy2(:) = BDrawsStore(:,g1,k1,p1)
                dummy3 = dummy2
                call piksrt(samsize0,dummy3)
                B_quantiles(g1,k1,p1,1) = dummy3(lower_int)
                B_quantiles(g1,k1,p1,2) = dummy3(median_int)
                B_quantiles(g1,k1,p1,3) = dummy3(upper_int)
            end do
            if(p1>1)then
                do p2=1,p1-1
                    dummy2(:) = CDrawsStore(:,g1,p1,p2)
                    dummy3 = dummy2
                    call piksrt(samsize0,dummy3)
                    C_quantiles(g1,p1,p2,1) = dummy3(lower_int)
                    C_quantiles(g1,p1,p2,2) = dummy3(median_int)
                    C_quantiles(g1,p1,p2,3) = dummy3(upper_int)
                end do
            end if
        end do
    end do


contains


function runiform ( iseed )

!*****************************************************************************80
!
!! RUNIFORM returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity iseed is an integer variable.
!
!    This routine implements the recursion
!
!      iseed = ( 16807 * iseed ) mod ( 2^31 - 1 )
!      runiform = iseed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial iseed is 12345, then the first three computations are
!
!      Input     Output      RUNIFORM
!      iseed      iseed
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 8 ) iseed, the "iseed" value, which should
!    NOT be 0. On output, iseed has been updated.
!
!    Output, real ( kind = 8 ) RUNIFORM, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer, parameter :: r15 = selected_real_kind(15)
  integer, parameter :: i6 = selected_int_kind(6)

  integer ( kind = i6 ), parameter :: i4_huge = 2147483647
  integer ( kind = i6 ) k
  real ( kind = r15 ) runiform
  integer ( kind = i6 ) iseed

  k = iseed / 127773

  iseed = 16807 * ( iseed - k * 127773 ) - k * 2836

  if ( iseed < 0 ) then
    iseed = iseed + i4_huge
  end if

  runiform = real ( iseed, kind = r15 ) * 4.656612875D-10

return
end function


function stvaln ( p )

!*****************************************************************************80
!
!! STVALN provides starting values for the inverse of the normal distribution.
!
!  Discussion:
!
!    The routine returns an X for which it is approximately true that
!      P = CUMNOR(X),
!    that is,
!      P = Integral ( -infinity < U <= X ) exp(-U*U/2)/sqrt(2*PI) dU.
!
!  Reference:
!
!    William Kennedy, James Gentle,
!    Statistical Computing,
!    Marcel Dekker, NY, 1980, page 95,
!    QA276.4 K46
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the probability whose normal deviate
!    is sought.
!
!    Output, real ( kind = 8 ) STVALN, the normal deviate whose probability
!    is approximately P.
!
  implicit none

  integer, parameter :: r15 = selected_real_kind(15)
  integer, parameter :: i6 = selected_int_kind(6)

  !real ( kind = 8 ) eval_pol
  real ( kind = r15 ) p
  real ( kind = r15 ) sgn
  real ( kind = r15 ) stvaln
  real ( kind = r15 ), parameter, dimension(0:4) :: xden = (/ &
    0.993484626060D-01, &
    0.588581570495D+00, &
    0.531103462366D+00, &
    0.103537752850D+00, &
    0.38560700634D-02 /)
  real ( kind = r15 ), parameter, dimension(0:4) :: xnum = (/ &
    -0.322232431088D+00, &
    -1.000000000000D+00, &
    -0.342242088547D+00, &
    -0.204231210245D-01, &
    -0.453642210148D-04 /)
  real ( kind = r15 ) y
  real ( kind = r15 ) z

  if ( p <= 0.5D+00 ) then

    sgn = -1.0D+00
    z = p

  else

    sgn = 1.0D+00
    z = 1.0D+00 - p

  end if

  y = sqrt ( -2.0D+00 * log ( z ) )
  stvaln = y + eval_pol ( xnum, 4, y ) / eval_pol ( xden, 4, y )
  stvaln = sgn * stvaln

  return
end function


function eval_pol ( a, n, x )

!*****************************************************************************80
!
!! EVAL_POL evaluates a polynomial at X.
!
!  Discussion:
!
!    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
!
!  Modified:
!
!    15 December 1999
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(0:N), coefficients of the polynomial.
!
!    Input, integer ( kind = 4 ) N, length of A.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) EVAL_POL, the value of the polynomial at X.
!
  implicit none

  integer, parameter :: r15 = selected_real_kind(15)
  integer, parameter :: i6 = selected_int_kind(6)

  integer ( kind = i6 ) n

  real ( kind = r15 ) a(0:n)
  real ( kind = r15 ) eval_pol
  integer ( kind = i6 ) i
  real ( kind = r15 ) term
  real ( kind = r15 ) x

  term = a(n)
  do i = n - 1, 0, -1
    term = term * x + a(i)
  end do

  eval_pol = term

  return
end function eval_pol



function cumnor ( arg )

!*****************************************************************************
! The original code was modified
!
!
!! CUMNOR computes the cumulative normal distribution.
!
!  Discussion:
!
!    This function evaluates the normal distribution function:
!
!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo
!
!    This transportable program uses rational functions that
!    theoretically approximate the normal distribution function to
!    at least 18 significant decimal digits.  The accuracy achieved
!    depends on the arithmetic system, the compiler, the intrinsic
!    functions, and proper selection of the machine dependent
!    constants.
!
!  Author:
!
!    William Cody
!    Mathematics and Computer Science Division
!    Argonne National Laboratory
!    Argonne, IL 60439
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev approximations for the error function,
!    Mathematics of Computation,
!    1969, pages 631-637.
!
!    William Cody,
!    Algorithm 715:
!    SPECFUN - A Portable FORTRAN Package of Special Function Routines
!    and Test Drivers,
!    ACM Transactions on Mathematical Software,
!    Volume 19, Number 1, 1993, pages 22-32.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the upper limit of integration.
!
!    Output, real ( kind = 8 ) the Normal density CDF.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) EPS, the argument below which anorm(x)
!    may be represented by 0.5 and above which  x*x  will not underflow.
!    A conservative value is the largest machine number X
!    such that   1.0D+00 + X = 1.0D+00   to machine precision.
!
  implicit none

  integer, parameter :: r15 = selected_real_kind(15)
  integer, parameter :: i6 = selected_int_kind(6)

  real ( kind = r15 ), parameter, dimension ( 5 ) :: a = (/ &
    2.2352520354606839287D+00, &
    1.6102823106855587881D+02, &
    1.0676894854603709582D+03, &
    1.8154981253343561249D+04, &
    6.5682337918207449113D-02 /)
  real ( kind = r15 ) arg
  real ( kind = r15 ), parameter, dimension ( 4 ) :: b = (/ &
    4.7202581904688241870D+01, &
    9.7609855173777669322D+02, &
    1.0260932208618978205D+04, &
    4.5507789335026729956D+04 /)
  real ( kind = r15 ), parameter, dimension ( 9 ) :: c = (/ &
    3.9894151208813466764D-01, &
    8.8831497943883759412D+00, &
    9.3506656132177855979D+01, &
    5.9727027639480026226D+02, &
    2.4945375852903726711D+03, &
    6.8481904505362823326D+03, &
    1.1602651437647350124D+04, &
    9.8427148383839780218D+03, &
    1.0765576773720192317D-08 /)
  real ( kind = r15 ) cumnor
  real ( kind = r15 ), parameter, dimension ( 8 ) :: d = (/ &
    2.2266688044328115691D+01, &
    2.3538790178262499861D+02, &
    1.5193775994075548050D+03, &
    6.4855582982667607550D+03, &
    1.8615571640885098091D+04, &
    3.4900952721145977266D+04, &
    3.8912003286093271411D+04, &
    1.9685429676859990727D+04 /)
  real ( kind = r15 ) del
  real ( kind = r15 ) eps
  integer ( kind = i6 ) i
  real ( kind = r15 ), parameter, dimension ( 6 ) :: p = (/ &
    2.1589853405795699D-01, &
    1.274011611602473639D-01, &
    2.2235277870649807D-02, &
    1.421619193227893466D-03, &
    2.9112874951168792D-05, &
    2.307344176494017303D-02 /)
  real ( kind = r15 ), parameter, dimension ( 5 ) :: q = (/ &
    1.28426009614491121D+00, &
    4.68238212480865118D-01, &
    6.59881378689285515D-02, &
    3.78239633202758244D-03, &
    7.29751555083966205D-05 /)
  real ( kind = r15 ), parameter :: root32 = 5.656854248D+00
  real ( kind = r15 ), parameter :: sixten = 16.0D+00
  real ( kind = r15 ), parameter :: sqrpi = 3.9894228040143267794D-01
  real ( kind = r15 ), parameter :: thrsh = 0.66291D+00
  real ( kind = r15 ) x
  real ( kind = r15 ) xden
  real ( kind = r15 ) xnum
  real ( kind = r15 ) y
  real ( kind = r15 ) xsq
!
!  Machine dependent constants
!
  eps = epsilon ( 1.0D+00 ) * 0.5D+00

  x = arg
  y = abs ( x )

  if ( y <= thrsh ) then
!
!  Evaluate  anorm  for  |X| <= 0.66291
!
    if ( eps < y ) then
      xsq = x * x
    else
      xsq = 0.0D+00
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do
    cumnor = x * ( xnum + a(4) ) / ( xden + b(4) )
    cumnor = 0.5D+00 + cumnor
!
!  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
!
  else if ( y <= root32 ) then

    xnum = c(9) * y
    xden = y
    do i = 1, 7
      xnum = ( xnum + c(i) ) * y
      xden = ( xden + d(i) ) * y
    end do
    cumnor = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = aint ( y * sixten ) / sixten
    del = ( y - xsq ) * ( y + xsq )
    cumnor = exp ( - xsq * xsq * 0.5D+00 ) * exp ( -del * 0.5D+00 ) * cumnor

    if ( 0.0D+00 < x ) then
       cumnor = 1D+00 - cumnor
    end if
!
!  Evaluate ANORM for sqrt(32) < |X|.
!
  else

    cumnor = 0.0D+00
    xsq = 1.0D+00 / ( x * x )
    xnum = p(6) * xsq
    xden = xsq
    do i = 1, 4
      xnum = ( xnum + p(i) ) * xsq
      xden = ( xden + q(i) ) * xsq
    end do

    cumnor = xsq * ( xnum + p(5) ) / ( xden + q(5) )
    cumnor = ( sqrpi - cumnor ) / y
    xsq = aint ( x * sixten ) / sixten
    del = ( x - xsq ) * ( x + xsq )
    cumnor = exp ( - xsq * xsq * 0.5D+00 ) &
      * exp ( - del * 0.5D+00 ) * cumnor

    if ( 0.0D+00 < x ) then
        cumnor = 1D+00 - cumnor
    end if

  end if

  if ( cumnor < tiny ( cumnor ) ) then
    cumnor = 0.0D+00
  end if

  return
 end function



 function dinvnr ( p )

!*****************************************************************************80
!
!! DINVNR computes the inverse of the normal distribution.
!
!  Discussion:
!
!    This routine returns X such that
!
!      CUMNOR(X) = P,
!
!    that is, so that
!
!      P = integral ( -oo <= T <= X ) exp(-U*U/2)/sqrt(2*PI) dU
!
!    The rational function is used as a
!    starting value for the Newton method of finding roots.
!
!  Reference:
!
!    William Kennedy, James Gentle,
!    Statistical Computing,
!    Marcel Dekker, NY, 1980,
!    QA276.4 K46
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P.
!
!    Output, real ( kind = 8 ) DINVNR, the argument X for which the
!    Normal CDF has the value P.
!
  implicit none

  integer, parameter :: r15 = selected_real_kind(15)
  integer, parameter :: i6 = selected_int_kind(6)

  real ( kind = r15 ) cum
  real ( kind = r15 ) dinvnr
  real ( kind = r15 ) dx
  real ( kind = r15 ), parameter :: eps = 1.0D-13
  integer ( kind = i6 ) i
  integer ( kind = i6 ), parameter :: maxit = 100
  real ( kind = r15 ) p
  real ( kind = r15 ) pp
  real ( kind = r15 ), parameter :: r2pi = 0.3989422804014326D+00
  real ( kind = r15 ) strtx
  !real ( kind = r15 ) stvaln
  real ( kind = r15 ) xcur

  pp = min ( p, 1-p )
  strtx = stvaln(pp)
  xcur = strtx
!
!  Newton iterations.
!
  do i = 1, maxit

    cum = cumnor(xcur)
    dx = ( cum - pp ) / ( r2pi * exp ( -0.5D+00 * xcur * xcur ) )
    xcur = xcur - dx

    if ( abs ( dx / xcur ) < eps ) then
      if ( p <= 1-p ) then
        dinvnr = xcur
      else
        dinvnr = -xcur
      end if
      return
    end if

  end do

  if ( p <= 1-p ) then
    dinvnr = strtx
  else
    dinvnr = -strtx
  end if

  return
end function dinvnr



subroutine kronecker(dimA,dimB,A,B,AB)
!
    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)
!
    integer(i6), intent(in)  :: dimA, dimB
    real(r15), intent(in)    :: A(dimA,dimA), B(dimB,dimB) !dummy arguments
    real(r15), intent(out)   :: AB(dimA*dimB,dimA*dimB) !output matrix of the kronecker product
    integer(i6)              :: i,j !loop counters
!
    do i=1,dimA
        do j=1,dimA
            AB((1+dimB*(i-1)):(dimB+dimB*(i-1)),(1+dimB*(j-1)):(dimB+dimB*(j-1))) = A(i,j)*B(:,:)
        end do
    end do
!
end subroutine kronecker



FUNCTION sdot(N,SX,INCX,SY,INCY)
 !
 !  -- Reference BLAS level1 routine (version 3.8.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2017
 !
    implicit none
 !
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)
 !
 !     .. Scalar Arguments ..
    INTEGER(i6) INCX,INCY,N
 !     ..
 !     .. Array Arguments ..
    REAL(r15) SX(*),SY(*)
 !    ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
    REAL(r15) STEMP, sdot
    INTEGER(i6) I,IX,IY,M,MP1
 !     ..
 !     .. Intrinsic Functions ..
    INTRINSIC mod
 !     ..
    stemp = 0.0e0
    sdot = 0.0e0
    IF (n.LE.0) RETURN
        IF (incx.EQ.1 .AND. incy.EQ.1) THEN
 !
 !        code for both increments equal to 1
 !
 !
 !        clean-up loop
 !
        m = mod(n,5)
        IF (m.NE.0) THEN
            DO i = 1,m
                stemp = stemp + sx(i)*sy(i)
            END DO
            IF (n.LT.5) THEN
                sdot = stemp
                RETURN
            END IF
        END IF
        mp1 = m + 1
        DO i = mp1,n,5
            stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) + sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
        END DO
    ELSE
 !
 !        code for unequal increments or equal increments
 !          not equal to 1
 !
        ix = 1
        iy = 1
        IF (incx.LT.0) ix = (-n+1)*incx + 1
            IF (incy.LT.0) iy = (-n+1)*incy + 1
                DO i = 1,n
                    stemp = stemp + sx(ix)*sy(iy)
                    ix = ix + incx
                    iy = iy + incy
                END DO
            END IF
            sdot = stemp
        RETURN
END FUNCTION sdot


function diag(A, n)

    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)
!
    integer(i6) n,i
    real(r15) A(n), check(n,n)
    real(r15) diag(n,n)

    check = 0
    do i=1,n
        check(i,i)=A(i)
    end do
    diag(:,:)=check(:,:)

    return

end function diag


function diagonals(A, n)

    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    integer(i6) n,i
    real(r15) A(n,n), diagonals(n), check(n)

    do i=1,n
        check(i)= A(i,i)
    enddo
    diagonals(:)=check(:)
     return
end function diagonals



!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
    IMPLICIT NONE

    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    !Declarations
    INTEGER(i6), INTENT(IN)  :: n
    REAL(r15), INTENT(IN)    :: matrix(n,n)  !Input matrix
    INTEGER(i6), INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
    REAL(r15), INTENT(OUT)   :: inverse(n,n) !Inverted matrix

    LOGICAL :: FLAG = .TRUE.
    INTEGER(i6) :: i, j, k
    REAL(r15) :: m
    REAL(r15), DIMENSION(n,2*n) :: augmatrix !augmented matrix

    inverse = eye(n)
    !Augment input matrix with an identity matrix
    DO i = 1, n
        DO j = 1, 2*n
            IF (j <= n ) THEN
                augmatrix(i,j) = matrix(i,j)
            ELSE IF ((i+n) == j) THEN
                augmatrix(i,j) = 1
            Else
                augmatrix(i,j) = 0
            ENDIF
        END DO
    END DO

    !Reduce augmented matrix to upper traingular form
    DO k =1, n-1
        IF (augmatrix(k,k) == 0) THEN
            FLAG = .FALSE.
            DO i = k+1, n
                IF (augmatrix(i,k) /= 0) THEN
                    DO j = 1,2*n
                        augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                    END DO
                    FLAG = .TRUE.
                    EXIT
                ENDIF
                IF (FLAG .EQV. .FALSE.) THEN
!                    PRINT*, "Matrix is non - invertible"
                    inverse = 0
                    errorflag = -1
                    return
                ENDIF
            END DO
        ENDIF
        DO j = k+1, n
            m = augmatrix(j,k)/augmatrix(k,k)
            DO i = k, 2*n
                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
            END DO
        END DO
    END DO

    !Test for invertibility
    DO i = 1, n
        IF (augmatrix(i,i) == 0) THEN
!            PRINT*, "Matrix is non - invertible"
            inverse = 0
            errorflag = -1
            return
        ENDIF
    END DO

    !Make diagonal elements as 1
    DO i = 1 , n
        m = augmatrix(i,i)
        DO j = i , (2 * n)
            augmatrix(i,j) = (augmatrix(i,j) / m)
        END DO
    END DO

    !Reduced right side half of augmented matrix to identity matrix
    DO k = n-1, 1, -1
        DO i = 1, k
            m = augmatrix(i,k+1)
            DO j = k, (2*n)
                augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
            END DO
        END DO
    END DO

    !store answer
    DO i =1, n
        DO j = 1, n
            inverse(i,j) = augmatrix(i,j+n)
        END DO
    END DO
    errorflag = 0
END SUBROUTINE FINDinv



recursive function det(a,n,permanent) result(accumulation)
    ! setting permanent to 1 computes the permanent.
    ! setting permanent to -1 computes the determinant.

    IMPLICIT NONE

    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    integer(i6), intent(in) :: n, permanent
    real(r15), dimension(n,n), intent(in) :: a
    real(r15), dimension(n-1, n-1) :: b
    real(r15) :: accumulation
    integer(i6) :: i, sgn

    if (n .eq. 1) then
      accumulation = a(1,1)
    else
      accumulation = 0
      sgn = 1
      do i=1, n
        b(:, :(i-1)) = a(2:, :i-1)
        b(:, i:) = a(2:, i+1:)
        accumulation = accumulation + sgn * a(1, i) * det(b, n-1, permanent)
        sgn = sgn * permanent
      enddo
    endif

end function det



function rnormal (iseed)

!*****************************************************************************80
!
!! RNORMAL returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!
!    Output, real ( kind = 8 ) RNORMAL, a normally distributed
!    random value.
!
  implicit none

  integer, parameter :: r15 = selected_real_kind(15)
  integer, parameter :: i6 = selected_int_kind(6)

  real ( kind = r15 ) r1
  real ( kind = r15 ) r2
  real ( kind = r15 ) r3
  real ( kind = r15 ) rnormal
  real ( kind = r15 ), parameter :: pi = 3.141592653589793D+00
  !real ( kind = r15 ) GG
  real ( kind = r15 ) x
  integer ( kind = i6 ) iseed

  !nseed = 1344
  r1 = runiform(iseed)
  r2 = runiform(iseed)
  r3 = runiform(iseed)
  !PRINT*, iseed, r1, r2, r3
  !call random_number(GG)
  !r1 = GG
  !call random_number(GG)
  !r2 = GG
  x = sqrt ( - 2.0D+00 * log ( r3 ) ) * cos ( 2.0D+00 * pi * r2 )

  rnormal = x

  return
end function rnormal



function eye(n)

    implicit none

    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    integer(i6) i,n
    real(r15) eye(n,n)
    real(r15) check(n,n)

    check=0
    do i=1,n
        check(i,i)= 1
    enddo

    eye(:,:)=check(:,:)
    return

end function eye



subroutine gen_wish(A,nu,B,P,iseed)
!
    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)
!
    !Declare local variables

    integer(i6), intent (in)  :: nu,P,iseed
    real(r15), intent (in)    :: A(P,P)
    real(r15), intent (out)   :: B(P,P)
    real(r15)                 :: RNmat(nu,P),para((P*(P+3)/2) + 1),m0(P)
    integer(i6)               :: i
!
    !sample from Wishart distribution as in Press (2005, p. 109)
    m0=0

    call setgmn(m0,A,P,para)

    do i=1,nu

        call GENMN(para,RNmat(i,:),P,iseed)

    end do

    B = matmul(transpose(RNmat),RNmat)
!
end subroutine gen_wish


!
subroutine robust_covest(m, betas1, betas2, mn1, mn2, varb1, varb2, varb1b2Plus, varb1b2Min)

    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)
!
    !Declare local variables
    integer(i6), intent(in)  :: m
    real(r15), intent(in)    :: betas1(m), betas2(m), mn1, mn2
    real(r15), intent(out)   :: varb1, varb2, varb1b2Plus, varb1b2Min

    real(r15)                :: dummy1(m), dummy2(m), Phi075, xxx
    integer(i6)              :: mmin, i
!
    xxx=0.75
    Phi075 = dinvnr(xxx)
    mmin = 0
!
    !robust variance estimators of beta1 and beta2

    dummy1=abs(betas1 - mn1)
    call piksrt(m,dummy1)
    do i=1,m
        if(dummy1(i)>0) then
            mmin = i
            exit
        end if
    end do
    varb1 = ((dummy1(int(mmin+(m-mmin)*.5)) + dummy1(int(mmin+(m-mmin)*.5+1)))*.5/Phi075)**2.0
    dummy1=abs(betas2 - mn2)
    call piksrt(m,dummy1)
    do i=1,m
        if(dummy1(i)>0) then
            mmin = i
            exit
        end if
    end do
    varb2 = ((dummy1(int(mmin+(m-mmin)*.5)) + dummy1(int(mmin+(m-mmin)*.5+1)))*.5/Phi075)**2.0
!
    !robust variance estimators of beta1 + beta2
    dummy2 = betas1 + betas2
    dummy1=abs(dummy2 - mn1 - mn2)
    call piksrt (m,dummy1)
    do i=1,m
        if(dummy1(i)>0) then
            mmin = i
            exit
        end if
    end do
    varb1b2Plus = ((dummy1(int(mmin+(m-mmin)*.5)) + dummy1(int(mmin+(m-mmin)*.5+1)))*.5/Phi075)**2.0
!
    !robust variance estimators of beta1 - beta2
    dummy2 = betas1 - betas2
    dummy1= abs(dummy2 - mn1 + mn2)
    call piksrt(m,dummy1)
    do i=1,m
        if(dummy1(i)>0) then
            mmin = i
            exit
        end if
    end do
    varb1b2Min = ((dummy1(int(mmin+(m-mmin)*.5)) + dummy1(int(mmin+(m-mmin)*.5+1)))*.5/Phi075)**2.0
!
end subroutine



SUBROUTINE piksrt(n,arr)

  implicit none
!
  integer, parameter :: r15 = selected_real_kind(15)
  integer, parameter :: i6 = selected_int_kind(6)

  integer(i6):: n, i,j
  real(r15):: arr(n), a

  do j=2, n
    a=arr(j)
    do i=j-1,1,-1
      if (arr(i)<=a) goto 10
      arr(i+1)=arr(i)
    end do
  i=0
10  arr(i+1)=a
  end do

  return

END SUBROUTINE



SUBROUTINE setgmn(meanv,covm,p,parm)
!**********************************************************************
!
!     SUBROUTINE SETGMN( MEANV, COVM, P, PARM)
!            SET Generate Multivariate Normal random deviate
!
!
!                              Function
!
!
!      Places P, MEANV, and the Cholesky factoriztion of COVM
!      in GENMN.
!
!
!                              Arguments
!
!
!     MEANV --> Mean vector of multivariate normal distribution.
!                                        REAL MEANV(P)
!
!     COVM   <--> (Input) Covariance   matrix    of  the  multivariate
!                 normal distribution
!                 (Output) Destroyed on output
!                                        REAL !OVM(P,P)
!
!     P     --> Dimension of the normal, or length of MEANV.
!                                        INTEGER P
!
!     PARM <-- Array of parameters needed to generate multivariate norma
!                deviates (P, MEANV and !holesky decomposition of
!                !OVM).
!                1 : 1                - P
!                2 : P + 1            - MEANV
!                P+2 : P*(P+3)/2 + 1  - !holesky decomposition of !OVM
!                                             REAL PARM(P*(P+3)/2 + 1)
!
!**********************************************************************

      implicit none
!
      integer, parameter :: r15 = selected_real_kind(15)
      integer, parameter :: i6 = selected_int_kind(6)

!     .. Scalar Arguments ..
      INTEGER(i6) p
!     ..
!     .. Array Arguments ..
      REAL(r15) covm(p,p),meanv(p),parm(p*(p+3)/2+1)
!     ..
!     .. Local Scalars ..
      INTEGER(i6) i,icount,info,j
!     ..
!     .. External Subroutines ..

!     ..
!     .. Executable Statements ..
!
!
!     TEST THE INPUT
!
      IF (.NOT. (p.LE.0)) GO TO 10
!      WRITE (*,*) 'P nonpositive in SETGMN'
!      WRITE (*,*) 'Value of P: ',p
!      STOP 'P nonpositive in SETGMN'

   10 parm(1) = p
!
!     PUT P AND MEANV INTO PARM
!
      DO 20,i = 2,p + 1
          parm(i) = meanv(i-1)
   20 CONTINUE
!
!      Cholesky decomposition to find A s.t. trans(A)*(A) = COVM
!
      CALL spofa(covm,p,p,info)
      IF (.NOT. (info.NE.0)) GO TO 30
!      WRITE (*,*) ' !OVM not positive definite in SETGMN'
!      STOP ' COVM not positive definite in SETGMN'

   30 icount = p + 1
!
!     PUT UPPER HALF OF A, WHICH IS NOW THE !HOLESKY FA!TOR, INTO PARM
!          !OVM(1,1) = PARM(P+2)
!          !OVM(1,2) = PARM(P+3)
!                    :
!          !OVM(1,P) = PARM(2P+1)
!          !OVM(2,2) = PARM(2P+2)  ...
!
      DO 50,i = 1,p
          DO 40,j = i,p
              icount = icount + 1
              parm(icount) = covm(i,j)
   40     CONTINUE
   50 CONTINUE
      RETURN
!
END SUBROUTINE



SUBROUTINE genmn(parm,x,p,iseed)
  !**********************************************************************
  !
  !     SUBROUTINE GENMN(PARM,X,WORK)
  !              GENerate Multivariate Normal random deviate
   !
   !
   !                              Arguments
   !
  !
  !     PARM --> Parameters needed to generate multivariate normal
  !               deviates (MEANV and Cholesky decomposition of
  !               COVM). Set by a previous call to SETGMN.
  !               1 : 1                - size of deviate, P
  !               2 : P + 1            - mean vector
  !               P+2 : P*(P+3)/2 + 1  - upper half of cholesky
  !                                       decomposition of cov matrix
  !                                             DOUBLE PRECISION PARM(*)
  !
  !     X    <-- Vector deviate generated.
  !                                             DOUBLE PRECISION X(P)
  !
  !     WORK <--> Scratch array
  !                                             DOUBLE PRECISION WORK(P)
  !
  !
  !                              Method
  !
  !
  !     1) Generate P independent standard normal deviates - Ei ~ N(0,1)
  !
  !     2) Using Cholesky decomposition find A s.t. trans(A)*A = COVM
  !
  !     3) trans(A)E + MEANV ~ N(MEANV,!OVM)
  !
  !**********************************************************************
  !
        implicit none
  !
        integer, parameter :: r15 = selected_real_kind(15)
        integer, parameter :: i6 = selected_int_kind(6)
  !
  !     .. Array Arguments ..
        integer(i6), intent(in) :: p, iseed
        real(r15), intent(in)   :: parm(p*(p+3)/2 + 1)
        real(r15)               :: work(p)
        real(r15), intent(out)  :: x(p)
  !     ..
  !     .. Local Scalars ..
        real(r15) ae
        INTEGER(i6) i,icount,j
  !    ..
  !     .. External Functions ..
  !     ..
  !     .. Executable Statements ..
  !
  !     Generate P independent normal deviates - WORK ~ N(0,1)
  !
        DO 10,i = 1,p
            work(i) = rnormal(iseed)
     10 CONTINUE
        DO 30,i = 1,p
  !
  !     PARM (P+2 : P*(P+3)/2 + 1) contains A, the Cholesky
  !      decomposition of the desired covariance matrix.
  !          trans(A)(1,1) = PARM(P+2)
  !          trans(A)(2,1) = PARM(P+3)
  !          trans(A)(2,2) = PARM(P+2+P)
  !          trans(A)(3,1) = PARM(P+4)
  !          trans(A)(3,2) = PARM(P+3+P)
  !          trans(A)(3,3) = PARM(P+2-1+2P)  ...
  !
  !     trans(A)*WORK + MEANV ~ N(MEANV,COVM)
  !
            icount = 0
            ae = 0.0
            DO 20,j = 1,i
                icount = icount + j - 1
                ae = ae + parm(i+ (j-1)*p-icount+p+1)*work(j)
     20     CONTINUE
            x(i) = ae + parm(i+1)
     30 CONTINUE
        RETURN
  !
END SUBROUTINE




subroutine spofa(a,lda,n,info)

      implicit none
!
      integer, parameter :: r15 = selected_real_kind(15)
      integer, parameter :: i6 = selected_int_kind(6)
!
      integer(i6) lda,n,info
      real(r15) a(lda,n)
!
!     spofa factors a real symmetric positive definite matrix.
!
!     spofa is usually called by spoco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for spo!o) = (1 + 18/n)*(time for spofa) .
!
!     on entry
!
!        a       real(lda, n)
!                the symmetric matrix to be factored.  only the
!                diagonal and upper triangle are used.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix  r  so that  a = trans(r)*r
!                where  trans(r)  is the transpose.
!                the stri!t lower triangle is unaltered.
!                if  info .ne. 0 , the fa!torization is not !omplete.
!
!        info    integer
!                = 0  for normal return.
!                = k  signals an error !ondition.  the leading minor
!                     of order  k  is not positive definite.
!
!     linpa!k.  this version dated 08/14/78 .
!     !leve moler, university of new mexi!o, argonne national lab.
!
!     subroutines and functions
!
!     blas sdot
!     fortran sqrt
!
!     internal variables
!
      real(r15) t
      real(r15) s
      integer(i6) j,jm1,k
!     begin block with ...exits to 40
!
         do 30 j = 1, n
            info = j
            s = 0.0e0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               t = a(k,j) - sdot(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
   10       continue
   20       continue
            s = a(j,j) - s
!     ......exit
            if (s .le. 0.0e0) go to 40
            a(j,j) = sqrt(s)
   30    continue
         info = 0
   40 continue
      return
end SUBROUTINE spofa


  end subroutine estimate_bct_ordinal_test





