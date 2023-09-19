


subroutine estimate_bct_ordinal(postZmean, postZcov, P, numcorr, K, numG, BHat, sdHat, CHat, XtXi, samsize0, &
    burnin, Ntot, Xgroups, Ygroups, C_quantiles, sigma_quantiles, B_quantiles, BDrawsStore, sigmaDrawsStore, &
    CDrawsStore, sdMH, ordinal, Cat, maxCat, seed)
!
    implicit none
!
    integer, intent(in) :: P, numcorr, K, numG, samsize0, burnin, Ntot, ordinal(P), Cat(numG,P), &
                           maxCat, seed
    real(8), intent(in) :: BHat(numG,K,P), sdHat(numG,P), CHat(numG,P,P), XtXi(numG,K,K), &
                           sdMH(numG,P), Xgroups(numG,Ntot,K), Ygroups(numG,Ntot,P)
    real(8), intent(out):: postZmean(numcorr,1), postZcov(numcorr,numcorr), B_quantiles(numG,K,P,3), &
                           C_quantiles(numG,P,P,3), sigma_quantiles(numG,P,3), BDrawsStore(samsize0,numG,K,P), &
                           sigmaDrawsStore(samsize0,numG,P), CDrawsStore(samsize0,numG,P,P)
    real(8) :: BDraws(numG,K,P), CDraws(numG,P,P), sigmaDraws(numG,P), SSj(P,P), meanMat(Ntot,P), SigmaMatDraw(P,P), &
               R_MH, cholSigmaPK(K*P,K*P), covBeta(K*P,K*P), Ds(P,P), Ccan(P,P), CcanInv(P,P), Ccurr(P,P), &
               CcurrInv(P,P), SS1(P,P), SS1inv(P,P), rnunif(1), errorMatj(P,P), sigma_can(P), aa, bb, cc, dummy1(1),  &
               telft, telct, f1, f0, f_1, UB, LB, alpha1, beta1, rr_can(1), betaDrawj(1,P*K), acceptSigma(numG,P), &
               C_can(P,P), C_canInv(P,P), bounds(2), MHpart1, MHpart2, MHpart3, dummyPP(P,P), dummyPPinv(P,P), &
               Wp(P,P), epsteps(P,P), logR_MH_part1, logR_MH_part2, logR_MH_part3, logR_MH_part4, &
               varz1, varz2, varz1z2Plus, varz1z2Min, Cnugget(P,P), SigmaInv(P,P), sdMHg(numG,P), gLiuSab_can, &
               Wgroups(numG,Ntot,P), dummyVar, alphaMin, alphaMax, Cinv(P,P), Bmean(K,P), acceptLS(numG,P), &
               alphaMat(numG,maxCat+1,P), Wdummy(numG,P,Ntot,maxCat), condMean, condVar, gLiuSab(numG,P), &
               ones(samsize0,1), Zcorr_sample(samsize0,numcorr), dummy3(samsize0), dummy2(samsize0), diffmat(Ntot,P), &
               meanO(P*K), para(((P*K)*((P*K)+3)/2 + 1))
    integer :: s1, g1, rankP, acceptC(numG), i1, dummyI, testPr, rr1, rr2, nutarget, a0, iter1, corrteller, &
               c1, c2, p1, Yi1Categorie, tellers(numG,maxCat,P), k1, k2, CatUpper, p2, iperm1(samsize0), iseed, &
               errorflag, Njs(numG)
!
!   set seed
    iseed = seed
!
    !initial posterior draws
    BDraws = BHat
    sigmaDraws = sdHat
    CDraws = CHat
    meanO = 0
    gLiuSab = 1
    iseed = seed
    Njs = Ntot
!
    !initial prior C
    nutarget = P + 1 !in the case of a marginally uniform prior
    do p1=1,P
        !initial values
        if(ordinal(p1)==1) then
            sigmaDraws(:,p1) = 1
            sigma_quantiles(:,p1,1) = 1
            sigma_quantiles(:,p1,2) = 1
            sigma_quantiles(:,p1,3) = 1
            BDraws(:,:,p1) = 0
        end if
    end do
    !specify candidate prior for R
    Wp = eye(P)
    a0 = nutarget !note that this implies an inv-W(Wp,a0) candidate prior.
                  !the parameterization in Remark 9 of Liu&Daniels seems a bit odd.
!
    !define nugget matrix to avoid approximate nonpositive definite correlation matrices for candidates
    Cnugget = Wp*1e-4
!
    telft = 0
    telct = 0
!
    !count number of accepted draws for R (over all groups)
    acceptC = 0
    acceptSigma = 0
    acceptLS = 0
    sdMHg = .1 !for gLiuBanhatti parameter
    dummyVar = 1.0
!
    !initial values for latent W's corresponding to ordinal DVs
    Wgroups = Ygroups
    Wdummy = 0
    ones = 1
!
    !initial values of boundary values alpha to link between ordinal Y and continuous latent W
    alphaMat = 0
    alphaMat(:,1,:) = -1e10  !alpha0
    alphaMat(:,2,:) = 0      !alpha1
    do p1=1,P
        if(ordinal(p1)>0) then
            do g1=1,numG
                do c1=3,Cat(g1,p1)
                    alphaMat(g1,c1,p1) = .3*(c1-2.0)
                end do
                alphaMat(g1,Cat(g1,p1)+1,p1) = 1e10
            end do
        end if
    end do
!

!    write(*,*)'sampling for burn-in period'
    !start Gibbs sampler
    do s1 = 1,burnin
        corrteller = 0
        tellers = 0
        do g1 = 1,numG
            !compute means of latent W's for all observations
            meanMat(1:Njs(g1),1:P) = matmul(Xgroups(g1,1:Njs(g1),1:K),BDraws(g1,1:K,1:P))
            Ccurr = CDraws(g1,:,:)
            SigmaMatDraw = matmul(matmul(diag(sigmaDraws(g1,:),P),Ccurr),diag(sigmaDraws(g1,:),P))
!
            !draw latent W's for the ordinal Y's
            !compute mean vector for

            do p1=1,P
                if(ordinal(p1)>0) then
                    do i1=1,Njs(g1)
                        Yi1Categorie = Ygroups(g1,i1,p1)
                        call compute_condMeanVar(p1,P,meanMat(i1,1:P),SigmaMatDraw, &
                            Wgroups(g1,i1,1:P),condMean,condVar)
                        select case (Yi1Categorie)
                            case(1)
                                call inverse_prob_sampling(condMean,condVar,0,1,alphaMat(g1,1,p1), &
                                    alphaMat(g1,2,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,1,p1) = tellers(g1,1,p1) + 1
                                Wdummy(g1,p1,tellers(g1,1,p1),1) = Wgroups(g1,i1,p1)
                            case(2)
                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,2,p1), &
                                    alphaMat(g1,3,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,2,p1) = tellers(g1,2,p1) + 1
                                Wdummy(g1,p1,tellers(g1,2,p1),2) = Wgroups(g1,i1,p1)
                            case(3)
                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,3,p1), &
                                    alphaMat(g1,4,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,3,p1) = tellers(g1,3,p1) + 1
                                Wdummy(g1,p1,tellers(g1,3,p1),3) = Wgroups(g1,i1,p1)
                            case(4)
                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,4,p1), &
                                    alphaMat(g1,5,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,4,p1) = tellers(g1,4,p1) + 1
                                Wdummy(g1,p1,tellers(g1,4,p1),4) = Wgroups(g1,i1,p1)
                            case(5)
                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,5,p1), &
                                    alphaMat(g1,6,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,5,p1) = tellers(g1,5,p1) + 1
                                Wdummy(g1,p1,tellers(g1,5,p1),5) = Wgroups(g1,i1,p1)
                            case(6)
                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,6,p1), &
                                    alphaMat(g1,7,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,6,p1) = tellers(g1,6,p1) + 1
                                Wdummy(g1,p1,tellers(g1,6,p1),6) = Wgroups(g1,i1,p1)
                            case(7)
                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,7,p1), &
                                    alphaMat(g1,8,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,7,p1) = tellers(g1,7,p1) + 1
                                Wdummy(g1,p1,tellers(g1,7,p1),7) = Wgroups(g1,i1,p1)
                            case(8)
                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,8,p1), &
                                    alphaMat(g1,9,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,8,p1) = tellers(g1,8,p1) + 1
                                Wdummy(g1,p1,tellers(g1,8,p1),8) = Wgroups(g1,i1,p1)
                            case(9)
                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,9,p1), &
                                    alphaMat(g1,10,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,9,p1) = tellers(g1,9,p1) + 1
                                Wdummy(g1,p1,tellers(g1,9,p1),9) = Wgroups(g1,i1,p1)
                            case(10)
                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,10,p1), &
                                    alphaMat(g1,11,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,10,p1) = tellers(g1,10,p1) + 1
                                Wdummy(g1,p1,tellers(g1,10,p1),10) = Wgroups(g1,i1,p1)
                            case(11)
                                call inverse_prob_sampling(condMean,condVar,1,0,alphaMat(g1,11,p1), &
                                    alphaMat(g1,12,p1),Wgroups(g1,i1,p1),iseed)
                                tellers(g1,11,p1) = tellers(g1,11,p1) + 1
                                Wdummy(g1,p1,tellers(g1,11,p1),11) = Wgroups(g1,i1,p1)
                         end select
                    end do
!
                    !draw boundary's in alphaMat
                    if(Cat(g1,p1)>2) then
                        do c1=3,Cat(g1,p1)
                            alphaMin = maxval(Wdummy(g1,p1,1:tellers(g1,c1-1,p1),c1-1))
                            alphaMax = minval(Wdummy(g1,p1,1:tellers(g1,c1,p1),c1))
                            alphaMat(g1,c1,p1) = runiform(iseed)*(alphaMax-alphaMin)*.9999999+alphaMin+.0000001
                        end do
                    end if
                end if
!
            end do


            ! HOLE

            Bmean(1:K,1:P) = matmul(matmul(XtXi(g1,:,:),transpose(Xgroups(g1,1:Njs(g1),1:K))), &
                Wgroups(g1,1:Njs(g1),1:P))
            call kronecker(K,P,XtXi(g1,:,:),SigmaMatDraw,covBeta)

            call setgmn(meanO,covBeta,P*K,para)
            call GENMN(para,betaDrawj(1,1:(P*K)),P*K,iseed)
            do p1 = 1,P
                BDraws(g1,:,p1) = betaDrawj(1,((p1-1)*K+1):(p1*K)) + Bmean(1:K,p1)
            end do

            !draw R using method of Liu and Daniels (LD, 2006)
            !draw candidate R
            diffmat(1:Njs(g1),1:P) = Wgroups(g1,1:Njs(g1),1:P) - matmul(Xgroups(g1,1:Njs(g1),1:K), &
                BDraws(g1,1:K,1:P))
            errorMatj = matmul(transpose(diffmat(1:Njs(g1),1:P)),diffmat(1:Njs(g1),1:P))
            Ds = diag(1/sqrt(diagonals(errorMatj,P)),P)
            diffmat(1:Njs(g1),1:P) = matmul(diffmat(1:Njs(g1),1:P),Ds) !diffmat is now epsilon in LD
            epsteps = matmul(transpose(diffmat(1:Njs(g1),1:P)),diffmat(1:Njs(g1),1:P))
            SS1 = matmul(matmul(diag(1/sigmaDraws(g1,:),P),epsteps),diag(1/sigmaDraws(g1,:),P))
            call FINDInv(SS1,SS1inv,P,errorflag)
            call gen_wish(SS1inv,Njs(g1),dummyPP,P,iseed)
            call FINDInv(dummyPP,dummyPPinv,P,errorflag)
	    	    Ccan = matmul(matmul(diag(1/sqrt(diagonals(dummyPPinv,P)),P),dummyPPinv), &
	    	        diag(1/sqrt(diagonals(dummyPPinv,P)),P))
!	    	Ccan = Ccan + Cnugget
            call FINDInv(Ccan,CcanInv,P,errorflag)
            call FINDInv(Ccurr,CcurrInv,P,errorflag)
            !target prior of Barnard et al. with nu = p + kappa0
            !proposal prior
            logR_MH_part3 = (-.5*real(P+1))*(log(det(Ccurr,P,-1))-log(det(Ccan,P,-1)))
            R_MH = exp(logR_MH_part3)

            !note that if Wp is not the identity matrix we have to
            !compute sum(diagonals(matmul(Wp,RcurrInv)))
            rnunif = runiform(iseed)
            if(rnunif(1) < R_MH) then
                CDraws(g1,:,:) = Ccan(:,:)
                acceptC(g1) = acceptC(g1) + 1
                Cinv = CcanInv
            else
                Cinv = CcurrInv
            end if


        end do
!
!        write(60,*)Wgroups(1,6,3),Wgroups(1,1,3),Wgroups(1,3,3),Wgroups(1,5,3)
!        write(61,*)alphaMat(1,2,3),alphaMat(1,3,3),alphaMat(1,4,3)
!
!        if(modulo(s1,burnin/10)==0) then
!            write(*,*)100*s1/burnin,' percent done'
!        end if
    end do


!!!!!!!!!!!!
!!!!
!!!! ERROR WHEN 'SELECT CASE' PART IS ADDED
!!!!
!!!!!!!!!!!!

!
!
!
contains
!
!
!


function eye(n)

    integer i,n
    real(8) eye(n,n)
    real(8) check(n,n)

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
    integer, intent(in) :: dimA, dimB
    real(8), intent(in)    :: A(dimA,dimA), B(dimB,dimB) !dummy arguments
    real(8), intent(out)   :: AB(dimA*dimB,dimA*dimB) !output matrix of the kronecker product
    integer             :: i,j !loop counters
!
    do i=1,dimA
        do j=1,dimA
            AB((1+dimB*(i-1)):(dimB+dimB*(i-1)),(1+dimB*(j-1)):(dimB+dimB*(j-1))) = A(i,j)*B(:,:)
        end do
    end do
!
end subroutine kronecker



function diag(A, n)

    integer n,i
    real(8) A(n), check(n,n)
    real(8) diag(n,n)


    check = 0
    do i=1,n
        check(i,i)=A(i)
    end do
    diag(:,:)=check(:,:)
    return
end function diag


function diagonals(A, n)

    integer n,i
    real(8) A(n,n), diagonals(n), check(n)

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

  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) rnormal
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  !real ( kind = 8 ) GG
  real ( kind = 8 ) x
  integer ( kind = 4 ) iseed

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

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) runiform
  integer ( kind = 4 ) iseed

  k = iseed / 127773

  iseed = 16807 * ( iseed - k * 127773 ) - k * 2836

  if ( iseed < 0 ) then
    iseed = iseed + i4_huge
  end if

  runiform = real ( iseed, kind = 8 ) * 4.656612875D-10

return
end function


recursive function det(a,n,permanent) result(accumulation)
    ! setting permanent to 1 computes the permanent.
    ! setting permanent to -1 computes the determinant.

    integer, intent(in) :: n, permanent
    real(8), dimension(n,n), intent(in) :: a
    real(8), dimension(n-1, n-1) :: b
    real(8) :: accumulation
    integer :: i, sgn
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
    !Declare local variables

    integer, intent (in)    :: nu,P,iseed
    real(8), intent (in)    :: A(P,P)
    real(8), intent (out)   :: B(P,P)
    real(8)                 :: RNmat(nu,P),para((P*(P+3)/2) + 1),m0(P)
    integer                 :: i
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
      INTEGER p
!     ..
!     .. Array Arguments ..
      REAL(8) covm(p,p),meanv(p),parm(p*(p+3)/2+1)
!     ..
!     .. Local Scalars ..
      INTEGER i,icount,info,j
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
  !     .. Array Arguments ..
        integer, intent(in) :: p, iseed
        real(8), intent(in) :: parm(p*(p+3)/2 + 1)
        real(8)             :: work(p)
        real(8), intent(out):: x(p)
  !     ..
  !     .. Local Scalars ..
        real(8) ae
        INTEGER i,icount,j
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

      integer lda,n,info
      real(8) a(lda,n)

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
      real(8) t
      real(8) s
      integer j,jm1,k
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



REAL(8) FUNCTION sdot(N,SX,INCX,SY,INCY)
 !
 !  -- Reference BLAS level1 routine (version 3.8.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2017
 !
 !     .. Scalar Arguments ..
    INTEGER INCX,INCY,N
 !     ..
  !    .. Array Arguments ..
    REAL(8) SX(*),SY(*)
 !    ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
    REAL(8) STEMP
    INTEGER I,IX,IY,M,MP1
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

    integer, intent(in)            :: welke, dimIn
    real (kind = 8), intent(in)    :: meanIn(dimIn), covmIn(dimIn,dimIn), obsIn(dimIn)
    real (kind = 8), intent(out)   :: condMean, condVar
    real (kind = 8)                :: dummy3(1,1), dummy2(dimIn-1,1), S12(1,dimIn-1), S22(dimIn-1,dimIn-1), &
                                      S22inv(dimIn-1,dimIn-1), ZS(dimIn,1), meanLocal(dimIn,1)
    integer                        :: errorflag
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

    integer, intent(in)            :: LBtrue, UBtrue, iseed
    real ( kind = 8 ), intent(in)  :: condMean, condVar, LB, UB
    real ( kind = 8 ), intent(out) :: condDraw
    real ( kind = 8 )              :: xdraw
    real ( kind = 8 )              :: LBstand, UBstand, yUB, yLB, rnIPS, diffUBLB, &
                                      bb, cc, Discr, xi, Zstar, &
                                      lambda, pi, machPres, rrand
    logical                        :: uppie

    parameter(pi=3.141592653)
    uppie = .true.
!
    machPres = 1e-6
    !normalize bounds
    UBstand = (UB - condMean)/sqrt(condVar)
    LBstand = (LB - condMean)/sqrt(condVar)
    !yUB = anordf (UBstand)
    yUB = alnorm ( UBstand, uppie )
    !yLB = anordf (LBstand)
    yLB = alnorm ( LBstand, uppie )
    rrand = runiform(iseed)
601     rnIPS = rrand*(yUB - yLB) + yLB
    if(LBtrue+UBtrue==0) then !unconstrained sampling
        rrand = rnormal(iseed)
        condDraw = rrand*sqrt(condVar) + condMean
    else if(abs(rnIPS) > machPres .and. abs(rnIPS-1) > machPres) then
        !inverse probability sampling
        !call normal_01_cdf_inv ( rnIPS, xdraw )
        xdraw = dinvnr ( rnIPS )
        condDraw = xdraw * sqrt(condVar) + condMean
        write(*,*),rnIPS,xdraw
    else if(UBstand>-5.0 .and. LBstand<5.0) then !IPS must be redone
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


subroutine ipv(condMean,condVar,LBtrue,UBtrue,LB,UB,condDraw,iseed)

    implicit none

    integer, intent(in)            :: LBtrue, UBtrue, iseed
    real ( kind = 8 ), intent(in)  :: condMean, condVar, LB, UB
    real ( kind = 8 ), intent(out) :: condDraw
    real ( kind = 8 )              :: xdraw
    real ( kind = 8 )              :: LBstand, UBstand, yUB, yLB, rnIPS, diffUBLB, &
                                      bb, cc, Discr, xi, Zstar, &
                                      lambda, pi, machPres, rrand
    logical                        :: uppie

    parameter(pi=3.141592653)
    uppie = .true.
!
    machPres = 1e-6
    !normalize bounds
    UBstand = (UB - condMean)/sqrt(condVar)
    LBstand = (LB - condMean)/sqrt(condVar)

    yUB = alnorm ( UBstand, uppie )
    yLB = alnorm ( LBstand, uppie )
    rrand = runiform(iseed)
    rnIPS = rrand * (yUB - yLB) + yLB
    xdraw = dinvnr ( rnIPS )
!    call normal_01_cdf_inv ( rnIPS, xdraw )
!    condDraw = xdraw * sqrt(condVar) + condMean
!
end subroutine ipv



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

  real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
  real ( kind = 8 ) alnorm
  real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
  real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
  real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
  real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
  real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
  real ( kind = 8 ), parameter :: con = 1.28D+00
  real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
  real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
  real ( kind = 8 ), parameter :: ltone = 7.0D+00
  real ( kind = 8 ), parameter :: p = 0.398942280444D+00
  real ( kind = 8 ), parameter :: q = 0.39990348504D+00
  real ( kind = 8 ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( kind = 8 ), parameter :: utzero = 18.66D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

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
end



subroutine normal_01_cdf_inv ( p, x )

!*****************************************************************************80
!
!! NORMAL_01_CDF_INV inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10^16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2015
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.  If P is outside this range, an
!    "infinite" value will be returned.
!
!    Output, real ( kind = 8 ) X, the normal deviate value
!    with the property that the probability of a standard normal deviate being
!    less than or equal to the value is P.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = 8 ), parameter :: const1 = 0.180625D+00
  real ( kind = 8 ), parameter :: const2 = 1.6D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) r8poly_value_horner
  real ( kind = 8 ), parameter :: split1 = 0.425D+00
  real ( kind = 8 ), parameter :: split2 = 5.0D+00
  real ( kind = 8 ) x

  if ( p <= 0.0D+00 ) then
    x = - huge ( x )
    return
  end if

  if ( 1.0D+00 <= p ) then
    x = huge ( x )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    x = q * r8poly_value_horner ( 7, a, r ) &
          / r8poly_value_horner ( 7, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then

      x = huge ( x )

    else

      r = sqrt ( - log ( r ) )

      if ( r <= split2 ) then

        r = r - const2
        x = r8poly_value_horner ( 7, c, r ) &
          / r8poly_value_horner ( 7, d, r )

      else

        r = r - split2
        x = r8poly_value_horner ( 7, e, r ) &
          / r8poly_value_horner ( 7, f, r )

      end if

    end if

    if ( q < 0.0D+00 ) then
      x = -x
    end if

  end if

  return
end



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

  real ( kind = 8 ) cum
  real ( kind = 8 ) dinvnr
  real ( kind = 8 ) dx
  real ( kind = 8 ), parameter :: eps = 1.0D-13
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: maxit = 100
  real ( kind = 8 ) p
  real ( kind = 8 ) pp
  real ( kind = 8 ), parameter :: r2pi = 0.3989422804014326D+00
  real ( kind = 8 ) strtx
  !real ( kind = 8 ) stvaln
  real ( kind = 8 ) xcur

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

  !real ( kind = 8 ) eval_pol
  real ( kind = 8 ) p
  real ( kind = 8 ) sgn
  real ( kind = 8 ) stvaln
  real ( kind = 8 ), parameter, dimension(0:4) :: xden = (/ &
    0.993484626060D-01, &
    0.588581570495D+00, &
    0.531103462366D+00, &
    0.103537752850D+00, &
    0.38560700634D-02 /)
  real ( kind = 8 ), parameter, dimension(0:4) :: xnum = (/ &
    -0.322232431088D+00, &
    -1.000000000000D+00, &
    -0.342242088547D+00, &
    -0.204231210245D-01, &
    -0.453642210148D-04 /)
  real ( kind = 8 ) y
  real ( kind = 8 ) z

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

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(0:n)
  real ( kind = 8 ) eval_pol
  integer ( kind = 4 ) i
  real ( kind = 8 ) term
  real ( kind = 8 ) x

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

  real ( kind = 8 ), parameter, dimension ( 5 ) :: a = (/ &
    2.2352520354606839287D+00, &
    1.6102823106855587881D+02, &
    1.0676894854603709582D+03, &
    1.8154981253343561249D+04, &
    6.5682337918207449113D-02 /)
  real ( kind = 8 ) arg
  real ( kind = 8 ), parameter, dimension ( 4 ) :: b = (/ &
    4.7202581904688241870D+01, &
    9.7609855173777669322D+02, &
    1.0260932208618978205D+04, &
    4.5507789335026729956D+04 /)
  real ( kind = 8 ), parameter, dimension ( 9 ) :: c = (/ &
    3.9894151208813466764D-01, &
    8.8831497943883759412D+00, &
    9.3506656132177855979D+01, &
    5.9727027639480026226D+02, &
    2.4945375852903726711D+03, &
    6.8481904505362823326D+03, &
    1.1602651437647350124D+04, &
    9.8427148383839780218D+03, &
    1.0765576773720192317D-08 /)
  real ( kind = 8 ) cumnor
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    2.2266688044328115691D+01, &
    2.3538790178262499861D+02, &
    1.5193775994075548050D+03, &
    6.4855582982667607550D+03, &
    1.8615571640885098091D+04, &
    3.4900952721145977266D+04, &
    3.8912003286093271411D+04, &
    1.9685429676859990727D+04 /)
  real ( kind = 8 ) del
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter, dimension ( 6 ) :: p = (/ &
    2.1589853405795699D-01, &
    1.274011611602473639D-01, &
    2.2235277870649807D-02, &
    1.421619193227893466D-03, &
    2.9112874951168792D-05, &
    2.307344176494017303D-02 /)
  real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
    1.28426009614491121D+00, &
    4.68238212480865118D-01, &
    6.59881378689285515D-02, &
    3.78239633202758244D-03, &
    7.29751555083966205D-05 /)
  real ( kind = 8 ), parameter :: root32 = 5.656854248D+00
  real ( kind = 8 ), parameter :: sixten = 16.0D+00
  real ( kind = 8 ), parameter :: sqrpi = 3.9894228040143267794D-01
  real ( kind = 8 ), parameter :: thrsh = 0.66291D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xden
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) xsq
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


!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
    IMPLICIT NONE
    !Declarations
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: matrix(n,n)  !Input matrix
    INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
    REAL(8), INTENT(OUT) :: inverse(n,n) !Inverted matrix

    LOGICAL :: FLAG = .TRUE.
    INTEGER :: i, j, k
    REAL(8) :: m
    REAL(8), DIMENSION(n,2*n) :: augmatrix !augmented matrix

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

