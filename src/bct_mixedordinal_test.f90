
!module rkinds
 !  use, intrinsic :: iso_c_binding
  ! private
  ! integer, parameter, public :: rint = c_int
 !  integer, parameter, public :: rdp = c_double
!end module
!
!subroutine estimate_bct_ordinal_test(P, numcorr, numG, Ntot, Njs, samsize0, Ygroups, seed, ordinal, Cat, &
!    CheckStore, postZmean, postZcov, sigma_quantiles, nuggetscale, Cnugget, maxCat)
!
subroutine estimate_bct_ordinal_test(P, numG, Ntot, Ygroups, Cnugget)
!
!    use rkinds, only: rint, rdp
    use, intrinsic ::  iso_c_binding
!
    implicit none

!    integer, parameter :: r15 = selected_real_kind(15)
!    integer, parameter :: i6 = selected_int_kind(6)

    integer(c_int), intent(in) :: P, numG, Ntot
    real(c_double), intent(in)   :: Ygroups(numG,Ntot,P)
    real(c_double), intent(out)   :: Cnugget(P,P)
    real(c_double)               :: nuggetscale
    integer(c_int)           :: p1
!
!
    !define nugget matrix to avoid approximate nonpositive definite correlation matrices for candidates
    nuggetscale = 0.995
    Cnugget = nuggetscale
    do p1=1,P
        Cnugget(p1,p1) = 1.0
    end do

end subroutine estimate_bct_ordinal_test





