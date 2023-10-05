
subroutine draw_ju(P,drawscorr,samsize,numcorrgroup,Fisher,seed)
    ! Fortran implementation of the algorithm proposed by Joe (2006)

    implicit none

    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    integer(i6), intent(in) :: P, samsize, numcorrgroup, Fisher, seed
    integer(i6)             :: s1,r1, r2, i1, i2, k1, corrIndex(P,P), teldummy,&
                               t1, t2, iseed, error1
    real (r15)              :: corrMat(P,P),draw1(1),&
                               R2inv(P,P), vec1(P,1), vec2(P,1),&
                               dummy11(1,1), dummy12(1,1), dummy22(1,1),&
                               Di1i2, preinv(P,P)
    real(r15), intent (out) :: drawscorr(samsize,numcorrgroup)
    real(r15)               :: alpha

!========================================================================================!

    !set seed
    !call RANDOM_SEED(size=nn)
    !allocate(iseed(nn))
    !iseed(:) = seed
    !call RANDOM_SEED(put=iseed)
    iseed = seed

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
            draw1 = random_beta(alpha, alpha, .true., iseed)
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
                draw1 = random_beta(alpha, alpha, .true., iseed)
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



FUNCTION random_beta(aa, bb, first, iseed) RESULT(fn_val)

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

    integer, parameter :: r15 = selected_real_kind(15)
    integer, parameter :: i6 = selected_int_kind(6)

    REAL(r15), INTENT(IN)    :: aa, bb
    LOGICAL, INTENT(IN) :: first
    INTEGER(i6), INTENT(IN) :: iseed
    REAL ( kind = r15 )   :: fn_val

    !     Local variables
    REAL(r15), PARAMETER  :: aln4 = 1.3862944, one=1.0, two=2.0, &
                             vlarge = HUGE(1.0), vsmall = TINY(1.0), zero = 0.0
    REAL ( kind = r15 ) :: a, b, g, r, s, x, y, z
    REAL ( kind = r15 ), SAVE        :: d, f, h, t, c
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
    r = runiform(iseed)
    x = runiform(iseed)
    !print*, r,x
    !CALL RANDOM_NUMBER(r)
    !CALL RANDOM_NUMBER(x)
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
end function


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



end subroutine draw_ju






