
module rkinds2
   use, intrinsic :: iso_c_binding
   private
   integer, parameter, public :: rint = c_int
   integer, parameter, public :: rdp = c_double
end module
!
!subroutine estimate_bct_ordinal_test(P, numcorr, numG, Ntot, Njs, samsize0, Ygroups, seed, ordinal, Cat, &
!    CheckStore, postZmean, postZcov, sigma_quantiles, nuggetscale, Cnugget, maxCat)
!
subroutine estimate_bct_ordinal_test(P, numG, Ntot, Ygroups, Cnugget)
!
    use rkinds2, only: rint, rdp
!    use, intrinsic ::  iso_c_binding
!
    implicit none

!    integer, parameter :: r15 = selected_real_kind(15)
!    integer, parameter :: i6 = selected_int_kind(6)

    integer(rint), intent(in)    :: P, numG, Ntot
    real(rdp), intent(in)    :: Ygroups(numG,Ntot,P)
    real(rdp), intent(out)   :: Cnugget(P,P)
    real(rdp)                :: nuggetscale
    integer(rint)                :: p1
!
!
    !define nugget matrix to avoid approximate nonpositive definite correlation matrices for candidates
    nuggetscale = 0.995_rdp
    Cnugget = nuggetscale
    do p1=1,P
        Cnugget(p1,p1) = 1.0_rdp
    end do

end subroutine estimate_bct_ordinal_test





