

module rkinds0
   use, intrinsic :: iso_c_binding
   use, intrinsic :: iso_fortran_env
   private
   integer, parameter, public :: rint = int32   ! Using int32 from iso_fortran_env
   integer, parameter, public :: rdp = real64   ! Using real64 from iso_fortran_env
   ! Using real64 from iso_fortran_env
end module


subroutine estimate_bct_ordinal_test(P, numG, Ntot, Ygroups, Cnugget)
!
    use rkinds0, only: rint, rdp
!
    implicit none

    integer(rint), intent(in)    :: P, numG, Ntot
    real(rdp), intent(in)    :: Ygroups(numG,Ntot,P)
    real(rdp), intent(out)   :: Cnugget(P,P)
    real(rdp)                :: nuggetscale
    integer(rint)                :: p1
!
    !define nugget matrix to avoid approximate nonpositive definite correlation matrices for candidates
    !nuggetscale = 0.995_rdp
    !Cnugget = nuggetscale
    !do p1=1,P
    !    Cnugget(p1,p1) = 1.0_rdp
    !end do

end subroutine estimate_bct_ordinal_test





