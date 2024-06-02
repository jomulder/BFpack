
module rkinds
   use, intrinsic :: iso_c_binding
   private
   integer, parameter, public :: rint = c_int
   integer, parameter, public :: rdp = c_double
end module


subroutine estimate_bct_ordinal(postZmean, postZcov, P, numcorr, K, numG, BHat, sdHat, CHat, XtXi, samsize0, &
    burnin, Ntot, Njs, Xgroups, Ygroups, C_quantiles, sigma_quantiles, B_quantiles, BDrawsStore, &
    sigmaDrawsStore, CDrawsStore, sdMH, ordinal_in, Cat_in, maxCat, gLiuSab, seed, nuggetscale)
!
    use rkinds, only: rint, rdp
!
    implicit none
!
!    integer, parameter :: r15 = selected_real_kind(15)
!    integer, parameter :: i6 = selected_int_kind(6)
!
    integer(rint), intent(in) :: P, numcorr, K, numG, samsize0, burnin, Ntot, maxCat, seed, Njs(numG)
    real(rdp), intent(in) ::  BHat(numG,K,P), sdHat(numG,P), CHat(numG,P,P), XtXi(numG,K,K), Cat_in(numG,P), &
                              sdMH(numG,P), Xgroups(numG,Ntot,K), Ygroups(numG,Ntot,P), ordinal_in(numG,P), &
                              nuggetscale
    real(rdp), intent(inout)::  postZmean(numcorr,1), postZcov(numcorr,numcorr), B_quantiles(numG,K,P,3), &
                              C_quantiles(numG,P,P,3), sigma_quantiles(numG,P,3), BDrawsStore(samsize0,numG,K,P), &
                              sigmaDrawsStore(samsize0,numG,P), CDrawsStore(samsize0,numG,P,P), &
                              gLiuSab(samsize0,numG,P)
    real(rdp) ::  BDraws(numG,K,P), CDraws(numG,P,P), sigmaDraws(numG,P), meanMat(Ntot,P), SigmaMatDraw(P,P), &
                  R_MH, covBeta(K*P,K*P), Ds(P,P), Ccan(P,P), CcanInv(P,P), Ccurr(P,P), epsteps(P,P), &
                  SS1(P,P), SS1inv(P,P), rnunif(1), errorMatj(P,P), sigma_can(P), aa, bb, &
                  betaDrawj(1,P*K), acceptSigma(numG,P), dummyPP(P,P), dummyPPinv(P,P), &
                  varz1, varz2, varz1z2Plus, varz1z2Min, Cnugget(P,P), SigmaInv(P,P), sdMHg(numG,P), gLiuSab_can, &
                  Wgroups(numG,Ntot,P), alphaMin, alphaMax, Cinv(P,P), Bmean(K,P), acceptLS(numG,P), &
                  alphaMat(numG,maxCat+1,P), Wdummy(numG,P,Ntot,maxCat), condMean, condVar, &
                  Zcorr_sample(samsize0,numcorr), dummy3(samsize0), dummy2(samsize0), &
                  diffmat(Ntot,P), meanO(P*K), para((P*K)*((P*K)+3)/2 + 1), randraw, gLiuSab_curr(numG,P)
    integer(rint) ::s1, g1, i1, corrteller, Cat(numG,P), ordinal(numG,P), &
                  c1, c2, p1, Yi1Categorie, tellers(numG,maxCat,P), k1, p2, iseed, errorflag, &
                  lower_int, median_int, upper_int
!
!   set seed
    iseed = seed
!
    !initial posterior draws
    BDraws = BHat
    sigmaDraws = sdHat
    CDraws = CHat
    meanO = 0.0
    gLiuSab_curr = 1.0
!
    do g1=1,numG
        do p1=1,P
            ordinal(g1,p1) = int(ordinal_in(g1,p1))
            Cat(g1,p1) = int(Cat_in(g1,p1))
        end do
    end do
    do p1=1,P
        do g1=1,numG
            !initial values
            if(ordinal(g1,p1)==1) then
                sigmaDraws(g1,p1) = 1.0
                sigma_quantiles(g1,p1,1) = 1.0
                sigma_quantiles(g1,p1,2) = 1.0
                sigma_quantiles(g1,p1,3) = 1.0
            end if
        end do
    end do
!
    !define nugget matrix to avoid approximate nonpositive definite correlation matrices for candidates
    Cnugget = nuggetscale
    do p1=1,P
        Cnugget(p1,p1) = 1.0
    end do
!
    !count number of accepted draws for R (over all groups)
    acceptSigma = 0.0
    acceptLS = 0.0
    sdMHg = .1 !for gLiuBanhatti parameter
!
    !initial values for latent W's corresponding to ordinal DVs
    Wgroups = Ygroups
    Wdummy = 0.0
!
    !initial values of boundary values alpha to link between ordinal Y and continuous latent W
    alphaMat = 0.0
    alphaMat(:,1,:) = -1e10  !alpha0
    alphaMat(:,2,:) = 0.0      !alpha1
    do p1=1,P
        do g1=1,numG
            if(ordinal(g1,p1)>0) then
                do c1=3,Cat(g1,p1)
                    alphaMat(g1,c1,p1) = .3*(real(c1)-2.0)
                end do
                alphaMat(g1,Cat(g1,p1)+1,p1) = 1e10
            end if
        end do
    end do

    !test write
    Ccan = 0
    CDrawsStore(1,1,:,:) = 1
    CDrawsStore = 1
    CDrawsStore(1,1,:,:) - Ccan(:,:)

contains



subroutine robust_covest(m, betas1, betas2, mn1, mn2, varb1, varb2, varb1b2Plus, varb1b2Min)

    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

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

  integer, parameter :: r15 = selected_real_kind(15)
  integer, parameter :: i6 = selected_int_kind(6)

  integer(i6) :: n, i,j
  real(r15)   :: arr(n), a

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


function eye(n)

    implicit none

    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    integer(i6):: i,n
    real(r15):: eye(n,n)
    real(r15):: check(n,n)

    check=0
    do i=1,n
        check(i,i)= 1
    enddo

    eye(:,:)=check(:,:)
    return

end function eye



subroutine kronecker(dimA,dimB,A,B,AB)
!
    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)
!
    integer(i6), intent(in) :: dimA, dimB
    real(r15), intent(in)   :: A(dimA,dimA), B(dimB,dimB) !dummy arguments
    real(r15), intent(out)  :: AB(dimA*dimB,dimA*dimB) !output matrix of the kronecker product
    integer(i6)             :: i,j !loop counters
!
    do i=1,dimA
        do j=1,dimA
            AB((1+dimB*(i-1)):(dimB+dimB*(i-1)),(1+dimB*(j-1)):(dimB+dimB*(j-1))) = A(i,j)*B(:,:)
        end do
    end do
!
end subroutine kronecker


function diag(A, n)

    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    integer(i6) :: n,i
    real(r15)   :: A(n), check(n,n)
    real(r15)   :: diag(n,n)

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

    integer(i6) :: n,i
    real(r15)   :: A(n,n), diagonals(n), check(n)

    do i=1,n
        check(i)= A(i,i)
    enddo
    diagonals(:)=check(:)

    return

end function diagonals



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
!
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
!
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
end function runiform


recursive function det(a,n,permanent) result(accumulation)
    ! setting permanent to 1 computes the permanent.
    ! setting permanent to -1 computes the determinant.

    implicit none
!
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


subroutine gen_wish(A,nu,B,P,iseed)
!
    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)
!
    !Declare local variables

    integer(i6), intent (in)    :: nu,P,iseed
    real(r15), intent (in)    :: A(P,P)
    real(r15), intent (out)   :: B(P,P)
    real(r15)                 :: RNmat(nu,P),para((P*(P+3)/2) + 1),m0(P)
    integer(i6)                 :: i
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
!     .. Scalar Arguments ..
      implicit none
!
      integer, parameter :: r15 = selected_real_kind(15)
      integer, parameter :: i6 = selected_int_kind(6)

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
END SUBROUTINE setgmn



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
  !     .. Array Arguments ..

        implicit none
!
        integer, parameter :: r15 = selected_real_kind(15)
        integer, parameter :: i6 = selected_int_kind(6)

        integer(i6), intent(in) :: p, iseed
        real(r15), intent(in) :: parm(p*(p+3)/2 + 1)
        real(r15)             :: work(p)
        real(r15), intent(out):: x(p)
  !     ..
  !     .. Local Scalars ..
        real(r15) ae
        integer(i6) :: i,icount,j
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
END SUBROUTINE genmn


subroutine spofa(a,lda,n,info)

      implicit none
!
      integer, parameter :: r15 = selected_real_kind(15)
      integer, parameter :: i6 = selected_int_kind(6)

      integer(i6) ::lda,n,info
      real(r15) :: a(lda,n)

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
      integer(i6) :: j,jm1,k
!     begin block with ...exits to 40
!
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



FUNCTION sdot(N,SX,INCX,SY,INCY)
 !
 !  -- Reference BLAS level1 routine (version 3.8.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2017
 !
 !     .. Scalar Arguments ..

    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    INTEGER(i6) :: INCX,INCY,N
 !     ..
  !    .. Array Arguments ..
    REAL(r15) :: SX(*),SY(*)
    real ( kind = r15 ) sdot
 !    ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
    REAL(r15) :: STEMP
    INTEGER(i6) :: I,IX,IY,M,MP1
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



subroutine compute_condMeanVar(welke,dimIn,meanIn,covmIn,obsIn,condMean,condVar)

    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    integer(i6), intent(in)         :: welke, dimIn
    real (kind = r15), intent(in)   :: meanIn(dimIn), covmIn(dimIn,dimIn), obsIn(dimIn)
    real (kind = r15), intent(out)  :: condMean, condVar
    real (kind = r15)               :: dummy3(1,1), dummy2(dimIn-1,1), S12(1,dimIn-1), S22(dimIn-1,dimIn-1), &
                                       S22inv(dimIn-1,dimIn-1), meanLocal(dimIn,1)
    integer                         :: errorflag
!
    meanLocal(1:dimIn,1) = meanIn(1:dimIn)
!
    S12(1,1:(welke-1)) = covmIn(welke,1:(welke-1))
    S12(1,welke:(dimIn-1)) = covmIn(welke,(welke+1):dimIn)
    S22(1:(welke-1),1:(welke-1)) = covmIn(1:(welke-1),1:(welke-1))
    S22(1:(welke-1),welke:(dimIn-1)) = covmIn(1:(welke-1),(welke+1):dimIn)
    S22(welke:(dimIn-1),1:(welke-1)) = covmIn((welke+1):dimIn,1:(welke-1))
    S22(welke:(dimIn-1),welke:(dimIn-1)) = covmIn((welke+1):dimIn,(welke+1):dimIn)
    call FINDInv(S22,S22inv,dimIn-1,errorflag)
    dummy2(1:(welke-1),1) = obsIn(1:(welke-1)) - meanLocal(1:(welke-1),1)
    dummy2(welke:(dimIn-1),1) = obsIn((welke+1):dimIn) - meanLocal((welke+1):dimIn,1)
    dummy3 = matmul(matmul(S12,S22inv),dummy2)
    condMean = meanLocal(welke,1) + dummy3(1,1) !conditional mean
    dummy3 = matmul(matmul(S12,S22inv),transpose(S12))
    condVar = covmIn(welke,welke) - dummy3(1,1) !conditional variance
!
end subroutine compute_condMeanVar


subroutine inverse_prob_sampling(condMean,condVar,LBtrue,UBtrue,LB,UB,condDraw,iseed)

    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    integer(i6), intent(in)           :: LBtrue, UBtrue, iseed
    real ( kind = r15 ), intent(in)   :: condMean, condVar, LB, UB
    real ( kind = r15 ), intent(out)  :: condDraw
    real ( kind = r15 )               :: xdraw
    real ( kind = r15 )               :: LBstand, UBstand, yUB, yLB, rnIPS, &
                                         pi, machPres, rrand
    integer(i6)                          teller
    logical                           :: uppie

    parameter(pi=3.141592653)
    uppie = .false.
!
    machPres = 1e-6
    !normalize bounds
    UBstand = (UB - condMean)/sqrt(condVar)
    LBstand = (LB - condMean)/sqrt(condVar)
    !yUB = anordf (UBstand)
    yUB = alnorm ( UBstand, uppie )
    !yLB = anordf (LBstand)
    yLB = alnorm ( LBstand, uppie )
    teller = 0
601    rrand = runiform(iseed)
    rnIPS = rrand*(yUB - yLB) + yLB
    if(LBtrue+UBtrue==0) then !unconstrained sampling
        rrand = rnormal(iseed)
        condDraw = rrand*sqrt(condVar) + condMean
    else if(abs(rnIPS) > machPres .and. abs(rnIPS-1) > machPres) then
        !inverse probability sampling
        xdraw = dinvnr ( rnIPS )
        condDraw = xdraw * sqrt(condVar) + condMean
    else if(UBstand>-4.0 .and. LBstand<4.0) then !IPS must be redone
        teller = teller + 1
        go to 601
    else
        if(condMean > UB) then
            condDraw = UB
        else if(condMean < LB) then
            condDraw = LB
        end if
    end if
!
end subroutine inverse_prob_sampling



function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none
!
  integer, parameter :: r15 = selected_real_kind(15)
  integer, parameter :: i6 = selected_int_kind(6)

  real ( kind = r15 ), parameter :: a1 = 5.75885480458D+00
  real ( kind = r15 ), parameter :: a2 = 2.62433121679D+00
  real ( kind = r15 ), parameter :: a3 = 5.92885724438D+00
  real ( kind = r15 ) alnorm
  real ( kind = r15 ), parameter :: b1 = -29.8213557807D+00
  real ( kind = r15 ), parameter :: b2 = 48.6959930692D+00
  real ( kind = r15 ), parameter :: c1 = -0.000000038052D+00
  real ( kind = r15 ), parameter :: c2 = 0.000398064794D+00
  real ( kind = r15 ), parameter :: c3 = -0.151679116635D+00
  real ( kind = r15 ), parameter :: c4 = 4.8385912808D+00
  real ( kind = r15 ), parameter :: c5 = 0.742380924027D+00
  real ( kind = r15 ), parameter :: c6 = 3.99019417011D+00
  real ( kind = r15 ), parameter :: con = 1.28D+00
  real ( kind = r15 ), parameter :: d1 = 1.00000615302D+00
  real ( kind = r15 ), parameter :: d2 = 1.98615381364D+00
  real ( kind = r15 ), parameter :: d3 = 5.29330324926D+00
  real ( kind = r15 ), parameter :: d4 = -15.1508972451D+00
  real ( kind = r15 ), parameter :: d5 = 30.789933034D+00
  real ( kind = r15 ), parameter :: ltone = 7.0D+00
  real ( kind = r15 ), parameter :: p = 0.398942280444D+00
  real ( kind = r15 ), parameter :: q = 0.39990348504D+00
  real ( kind = r15 ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( kind = r15 ), parameter :: utzero = 18.66D+00
  real ( kind = r15 ) x
  real ( kind = r15 ) y
  real ( kind = r15 ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 &
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
end function alnorm



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
!
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
!
  integer, parameter :: r15 = selected_real_kind(15)
  integer, parameter :: i6 = selected_int_kind(6)

  !real ( kind = r15 ) eval_pol
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
end function stvaln


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
!
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
!
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
 end function cumnor


!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
SUBROUTINE FINDinv(matrix, inverse, n, errorflag)
    implicit none
!
    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    !Declarations
    INTEGER(i6), INTENT(IN) :: n
    REAL(r15), INTENT(IN) :: matrix(n,n)  !Input matrix
    INTEGER(i6), INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
    REAL(r15), INTENT(OUT) :: inverse(n,n) !Inverted matrix

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



end subroutine estimate_bct_ordinal

