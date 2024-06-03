
!module rkinds
!   use, intrinsic :: iso_c_binding
!   private
!   integer, parameter, public :: rint = c_int
!   integer, parameter, public :: rdp = c_double
!end module


subroutine estimate_bct_ordinal(postZmean, postZcov, P, numcorr, K, numG, BHat, sdHat, CHat, XtXi, samsize0, &
    burnin, Ntot, Njs_in, Xgroups, Ygroups, C_quantiles, sigma_quantiles, B_quantiles, BDrawsStore, &
    sigmaDrawsStore, CDrawsStore, sdMH, ordinal_in, Cat_in, maxCat, gLiuSab, seed, nuggetscale)
!
!    use rkinds, only: rint, rdp
!
    implicit none
!
    integer, intent(in) :: P, numcorr, K, numG, samsize0, burnin, Ntot, maxCat, seed
    double precision, intent(in) ::  BHat(numG,K,P), sdHat(numG,P), CHat(numG,P,P), XtXi(numG,K,K), Cat_in(numG,P), &
                              sdMH(numG,P), Xgroups(numG,Ntot,K), Ygroups(numG,Ntot,P), ordinal_in(numG,P), &
                              nuggetscale, Njs_in(numG,1)
    double precision, intent(inout)::  postZmean(numcorr,1), postZcov(numcorr,numcorr), B_quantiles(numG,K,P,3), &
                              C_quantiles(numG,P,P,3), sigma_quantiles(numG,P,3), BDrawsStore(samsize0,numG,K,P), &
                              sigmaDrawsStore(samsize0,numG,P), CDrawsStore(samsize0,numG,P,P), &
                              gLiuSab(samsize0,numG,P)
    double precision ::  BDraws(numG,K,P), CDraws(numG,P,P), sigmaDraws(numG,P), meanMat(Ntot,P), SigmaMatDraw(P,P), &
                  R_MH, covBeta(K*P,K*P), Ds(P,P), Ccan(P,P), CcanInv(P,P), Ccurr(P,P), epsteps(P,P), &
                  SS1(P,P), SS1inv(P,P), rnunif(1), errorMatj(P,P), sigma_can(P), aa, bb, &
                  betaDrawj(1,P*K), acceptSigma(numG,P), dummyPP(P,P), dummyPPinv(P,P), &
                  varz1, varz2, varz1z2Plus, varz1z2Min, Cnugget(P,P), SigmaInv(P,P), sdMHg(numG,P), gLiuSab_can, &
                  Wgroups(numG,Ntot,P), alphaMin, alphaMax, Cinv(P,P), Bmean(K,P), acceptLS(numG,P), &
                  alphaMat(numG,maxCat+1,P), Wdummy(numG,P,Ntot,maxCat), condMean, condVar, &
                  Zcorr_sample(samsize0,numcorr), dummy3(samsize0), dummy2(samsize0), &
                  diffmat(Ntot,P), meanO(P*K), para((P*K)*((P*K)+3)/2 + 1), randraw, gLiuSab_curr(numG,P)
    integer     ::s1, g1, i1, corrteller, Cat(numG,P), ordinal(numG,P), Njs(numG), &
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
        Njs(g1) = int(Njs_in(g1,1))
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
    alphaMat(:,1,:) = -1e10    !alpha0
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

    !
    !
    ! removed part
    !
    !

    !test write
    Ccan = 1
    CDrawsStore(1,1,:,:) = Ccan(:,:)

end subroutine estimate_bct_ordinal

