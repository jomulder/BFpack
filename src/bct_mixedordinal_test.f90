
!module rkinds
 !  use, intrinsic :: iso_c_binding
  ! private
  ! integer, parameter, public :: rint = c_int
 !  integer, parameter, public :: rdp = c_double
!end module

subroutine estimate_bct_ordinal_test(P, numcorr, numG, Ntot, Njs, samsize0, Ygroups, seed, ordinal, Cat, &
    CheckStore, postZmean, postZcov, sigma_quantiles, nuggetscale, Cnugget, maxCat)
!
!    use rkinds, only: rint, rdp
!
    implicit none

    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    integer(i6), intent(in) :: P, numcorr, numG, Ntot, samsize0, seed, maxCat
    integer(i6), intent(inout) :: Njs(numG), ordinal(numG,P), Cat(numG,P)
    real(r15), intent(inout)   :: Ygroups(numG,Ntot,P), sigma_quantiles(numG,P,3), nuggetscale
    real(r15), intent(inout)  :: postZmean(numcorr,1), postZcov(numcorr,numcorr), CheckStore(samsize0,numG,P), &
                            Cnugget(P,P)
    real(r15)               :: alphaMat(numG,maxCat+1,P), sigmaDraws(numG,P), Wgroups(numG,Ntot,P), &
                            Wdummy(numG,P,Ntot,maxCat), &
                            CDraws(numG,P,P), Ccan(P,P), Ds(P,P), dummyPP2(P,P), &
                            CcanInv(P,P), SS1(P,P), rnunif(1), errorMatj(P,P), &
                            sigma_can(P), aa, bb, SigmaMat(P,P), R_MH, epsteps(P,P), diffmat(Ntot,P), &
                            varz1, varz2, varz1z2Plus, varz1z2Min, Cinv(P,P), Zcorr_sample(samsize0,numcorr), &
                            acceptSigma(numG,P), dummyPP(P,P), &
                            dummy3(samsize0), dummy2(samsize0), SS2(P,P), &
                            gLiuSab_curr(numG,P), acceptLS(numG,P), sdMHg(numG,P)
    integer(i6)           :: c1, p1, s1, g1, i1, corrteller, iseed
!
!   set seed
    iseed = seed
!
    do p1=1,P
        do g1=1,numG
            !initial values
            if(ordinal(g1,p1)==1.0) then
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
    if(maxCat>1) then
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
    end if
!
  postZmean(1:numcorr,1) = (/1.3,.2,5.9/)

  end subroutine estimate_bct_ordinal_test





