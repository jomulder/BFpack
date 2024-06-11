
! rngfuncs.f90
module rngfuncs1

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

module rkinds1
   use, intrinsic :: iso_c_binding
   use, intrinsic :: iso_fortran_env
   private
   integer, parameter, public :: rint = int32   ! Using int32 from iso_fortran_env
   integer, parameter, public :: rdp = real64   ! Using real64 from iso_fortran_env
   ! Using real64 from iso_fortran_env
end module


subroutine draw_ju(P,drawscorr,samsize,numcorrgroup,Fisher)
    ! Fortran implementation of the algorithm proposed by Joe (2006)

    use rkinds1, only: rint, rdp
    use rngfuncs1

    implicit none

    integer(rint), intent(in) :: P, samsize, numcorrgroup, Fisher
    integer(rint)             :: s1,r1, r2, i1, i2, k1, corrIndex(P,P), teldummy,&
                               t1, t2, error1
    real (rdp)              :: corrMat(P,P),draw1(1),&
                               R2inv(P,P), vec1(P,1), vec2(P,1),&
                               dummy11(1,1), dummy12(1,1), dummy22(1,1),&
                               Di1i2, preinv(P,P)
    real(rdp), intent (out) :: drawscorr(samsize,numcorrgroup)
    real(rdp)               :: alpha

!========================================================================================!

    ! create corrIndex matrix
    teldummy = 1

    do r1=2,P
        do r2=1,r1-1
            corrIndex(r1,r2) = teldummy
            corrIndex(r2,r1) = teldummy
            teldummy = teldummy + 1
        end do
    end do

    do s1=1,samsize

        ! create identity matrix
        do t1=1,P
            do t2=1,P
                if (t1==t2) then
                    corrMat(t1,t2)=1.0
                else
                    corrMat(t1,t2)=0.0
                ENDIF
            end do
        end do
        do r1 = 1,P-1
            alpha=P/2.0
            draw1 = random_beta(alpha, alpha, .true.)
            draw1 = draw1*2.0-1.0
            corrMat(r1,r1+1) = draw1(1)
            corrMat(r1+1,r1) = corrMat(r1,r1+1)
            drawscorr(s1,corrIndex(r1+1,r1)) = corrMat(r1,r1+1)
        end do
        R2inv(:,:) = 0
        preinv(:,:)= 0
        do r1 = 3,P
            do r2 = 1,P-r1+1
                i1 = r2
                i2 = r2+r1-1
                k1 = i2 - i1
                !draw partial correlations
                alpha = .5*(P+1-k1)
                draw1 = random_beta(alpha, alpha, .true.)
                draw1=draw1*2-1.0
                !rbeta(1,.5*(dim+1-k),.5*(dim+1-k))*2-1
                vec1(1:(k1-1),1) = corrMat(i1,(i1+1):(i1+k1-1))
                vec2(1:(k1-1),1) = corrMat(i2,(i1+1):(i1+k1-1))
                preinv(1:(i2-i1-1),1:(i2-i1-1)) = ((corrMat((i1+1):(i2-1),(i1+1):(i2-1))))
        !        R2inv(1:(i2-i1-1),1:(i2-i1-1)) = inverse(preinv(1:(i2-i1-1),1:(i2-i1-1)),(i2-i1-1))
                call FINDinv(preinv(1:(i2-i1-1),1:(i2-i1-1)),R2inv(1:(i2-i1-1),1:(i2-i1-1)),(i2-i1-1),error1)

                dummy11 = matmul(matmul(transpose(vec1(1:(k1-1),:)),R2inv(1:(i2-i1-1),1:(i2-i1-1))),vec1(1:(k1-1),:))
                dummy22 = matmul(matmul(transpose(vec2(1:(k1-1),:)),R2inv(1:(i2-i1-1),1:(i2-i1-1))),vec2(1:(k1-1),:))
                dummy12 = matmul(matmul(transpose(vec1(1:(k1-1),:)),R2inv(1:(i2-i1-1),1:(i2-i1-1))),vec2(1:(k1-1),:))
                Di1i2 = sqrt((1-dummy11(1,1))*(1-dummy22(1,1)))

                corrMat(i1,i2) = dummy12(1,1) + Di1i2*draw1(1)
                corrMat(i2,i1) = corrMat(i1,i2)

                drawscorr(s1,corrIndex(i1,i2)) = corrMat(i1,i2)
                end do
            end do
        end do
    if(Fisher==1) then
        drawscorr(1:samsize,1:numcorrgroup) = .5*log((1.0+drawscorr(1:samsize,1:numcorrgroup)) &
                                        /(1.0-drawscorr(1:samsize,1:numcorrgroup)))
    end if


contains


! Subroutine to find the inverse of a square matrix
SUBROUTINE FINDinv(matrix, inverse, n, errorflag)

    implicit none

    !Declarations
    INTEGER(rint), INTENT(IN) :: n
    REAL(rdp), INTENT(IN) :: matrix(n,n)  !Input matrix
    INTEGER(rint), INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
    REAL(rdp), INTENT(OUT) :: inverse(n,n) !Inverted matrix

    integer :: ipiv(n), info, lwork
    real(rdp) :: work(n)

    external :: dgetrf, dgetri

    errorflag = 0

    inverse = matrix
    call dgetrf(n,n,inverse,n,ipiv,info)
    if (info > 0) then
        inverse = 0
        errorflag = -1
        return
    end if

    lwork = n
    call dgetri(n,inverse,n,ipiv,work,lwork,info)
    if (info > 0) then
        inverse = 0
        errorflag = -1
        return
    end if

END SUBROUTINE FINDinv



FUNCTION random_beta(aa, bb, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     Author: Alan Miller
!             CSIRO Division of Mathematical & Information Sciences
!             Private Bag 10, Clayton South MDC
!             Clayton 3169, Victoria, Australia
!     Phone: (+61) 3 9545-8016      Fax: (+61) 3 9545-8080
!     e-mail: amiller @ bigpond.net.au

! FUNCTION GENERATES A RANDOM VARIATE IN [0,1]
! FROM A BETA DISTRIBUTION WITH DENSITY
! PROPORTIONAL TO BETA**(AA-1) * (1-BETA)**(BB-1).
! USING CHENG'S LOG LOGISTIC METHOD.

!     AA = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)
!     BB = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)

    implicit none

    REAL(rdp), INTENT(IN)    :: aa, bb
    LOGICAL, INTENT(IN) :: first
    !INTEGER(rint), INTENT(IN) :: iseed
    REAL ( kind = rdp )   :: fn_val

    !     Local variables
    REAL(rdp), PARAMETER  :: aln4 = 1.3862944, one=1.0, two=2.0, &
                             vlarge = HUGE(1.0), vsmall = TINY(1.0), zero = 0.0
    REAL ( kind = rdp ) :: a, b, g, r, s, x, y, z
    REAL ( kind = rdp ), SAVE        :: d, f, h, t, c
    LOGICAL, SAVE     :: swap



    IF (first) THEN                        ! Initialization, if necessary
    a = aa
    b = bb
    swap = b > a
    IF (swap) THEN
        g = b
        b = a
        a = g
    END IF
    d = a/b
    f = a+b
    IF (b > one) THEN
        h = SQRT((two*a*b - f)/(f - two))
        t = one
    ELSE
        h = b
        t = one/(one + (a/(vlarge*b))**b)
    END IF
    c = a+h
    END IF

    DO
    r = unif_rand()
    x = unif_rand()
    !print*, r,x

    s = r*r*x
    IF (r < vsmall .OR. s <= zero) CYCLE
    IF (r < t) THEN
        x = LOG(r/(one - r))/h
        y = d*EXP(x)
        z = c*x + f*LOG((one + d)/(one + y)) - aln4
        IF (s - one > z) THEN
        IF (s - s*z > one) CYCLE
        IF (LOG(s) > z) CYCLE
        END IF
        fn_val = y/(one + y)
    ELSE
        IF (4.0*s > (one + one/d)**f) CYCLE
        fn_val = one
    END IF
    EXIT
    END DO

    IF (swap) fn_val = one - fn_val
    RETURN
END FUNCTION random_beta



function eye(n)

    implicit none

    integer(rint):: i,n
    real(rdp):: eye(n,n)
    real(rdp):: check(n,n)

    check=0
    do i=1,n
        check(i,i)= 1
    enddo

    eye(:,:) = check(:,:)
    return

end function eye


end subroutine draw_ju






