
!
module rkindss
   use, intrinsic :: iso_c_binding !c_int c_double
   use, intrinsic :: iso_fortran_env !int32 real64
   private
   integer, parameter, public :: rint = c_int   ! Using int32 from iso_fortran_env
   integer, parameter, public :: rdp = c_double   ! Using real64 from iso_fortran_env
   ! Using real64 from iso_fortran_env
end module


subroutine test_kind(X, Y)

    use rkindss, only: rint, rdp

    implicit none

    real(rdp), intent(in) ::  X
    real(rdp), intent(out) ::  Y

    Y = X + 1.0_rdp

end subroutine test_kind


