

! rngfuncs.f90
module rngfuncs
implicit none
public
! See the R header `R_ext/Random.h` for details
interface
   ! double unif_rand(void);
   function unif_rand() bind(c,name="unif_rand")
      use, intrinsic :: iso_c_binding, only: c_double
      real(c_double) :: unif_rand
   end function

   ! void GetRNGstate(void);
   subroutine getrngstate() bind(c,name="GetRNGstate")
   end subroutine

   ! void PutRNGstate(void);
   subroutine putrngstate() bind(c,name="PutRNGstate")
   end subroutine
end interface
end module


module rkinds0
   use, intrinsic :: iso_c_binding
   use, intrinsic :: iso_fortran_env
   private
   integer, parameter, public :: rint = int32   ! Using int32 from iso_fortran_env
   integer, parameter, public :: rdp = real64   ! Using real64 from iso_fortran_env
   ! Using real64 from iso_fortran_env
end module


subroutine estimate_bct_ordinal_test(postZmean, postZcov, P, numcorr, K, numG, BHat, sdHat, CHat, XtXi, samsize0, &
    burnin, Ntot, Njs_in, Xgroups, Ygroups, C_quantiles, sigma_quantiles, B_quantiles, BDrawsStore, &
    sigmaDrawsStore, CDrawsStore, sdMH, ordinal_in, Cat_in, maxCat, gLiuSab, seed, nuggetscale, WgroupsStore, &
    meanMatMeanStore, SigmaMatDrawStore, CheckStore)
!
    use rkinds0, only: rint, rdp
    use rngfuncs
!
    implicit none
!
    integer(rint), intent(in) ::P, numcorr, K, numG, samsize0, burnin, Ntot, maxCat, seed
    real(rdp), intent(in) ::  BHat(numG,K,P), sdHat(numG,P), CHat(numG,P,P), XtXi(numG,K,K), Cat_in(numG,P), &
                              sdMH(numG,P), Xgroups(numG,Ntot,K), Ygroups(numG,Ntot,P), ordinal_in(numG,P), &
                              nuggetscale, Njs_in(numG,1)
    real(rdp), intent(out)::  postZmean(numcorr,1), postZcov(numcorr,numcorr), B_quantiles(numG,K,P,3), &
                              C_quantiles(numG,P,P,3), sigma_quantiles(numG,P,3), BDrawsStore(samsize0,numG,K,P), &
                              sigmaDrawsStore(samsize0,numG,P), CDrawsStore(samsize0,numG,P,P), &
                              gLiuSab(samsize0,numG,P), WgroupsStore(samsize0,numG,Ntot,P), &
                              meanMatMeanStore(samsize0,Ntot,P), SigmaMatDrawStore(samsize0,P,P), &
                              CheckStore(samsize0,numG,Ntot,P,3*P+2+3)
    real(rdp) ::  BDraws(numG,K,P), CDraws(numG,P,P), sigmaDraws(numG,P), meanMat(Ntot,P), SigmaMatDraw(P,P), &
                  R_MH, covBeta(K*P,K*P), Ds(P,P), Ccan(P,P), CcanInv(P,P), Ccurr(P,P), epsteps(P,P), &
                  SS1(P,P), SS1inv(P,P), rnunif(1), errorMatj(P,P), sigma_can(P), aa, bb, &
                  betaDrawj(1,P*K), acceptSigma(numG,P), dummyPP(P,P), dummyPPinv(P,P), &
                  varz1, varz2, varz1z2Plus, varz1z2Min, Cnugget(P,P), SigmaInv(P,P), sdMHg(numG,P), gLiuSab_can, &
                  Wgroups(numG,Ntot,P), alphaMin, alphaMax, Cinv(P,P), Bmean(K,P), acceptLS(numG,P), &
                  alphaMat(numG,maxCat+1,P), Wdummy(numG,P,Ntot,maxCat), condMean, condVar, &
                  Zcorr_sample(samsize0,numcorr), dummy3(samsize0), dummy2(samsize0), &
                  diffmat(Ntot,P), meanO(P*K), para((P*K)*((P*K)+3)/2 + 1), randraw, gLiuSab_curr(numG,P)
    integer(rint) ::s1, g1, i1, corrteller, Cat(numG,P), ordinal(numG,P), Njs(numG), &
                  c1, c2, p1, Yi1Categorie, tellers(numG,maxCat,P), k1, p2, iseed, errorflag, &
                  lower_int, median_int, upper_int
!
!    call getrngstate

    do i1 = 1, p
        CDrawsStore(1,1,1,i1) = unif_rand()
    end do

!    call putrngstate


end subroutine estimate_bct_ordinal_test





