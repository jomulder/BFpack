
subroutine draw_ju(P,drawscorr,samsize,numcorrgroup,Fisher,seed)
    ! Fortran implementation of algorithm from on Joe (2006)
    implicit none
    integer                             :: s1,r1, r2, i1, i2, k1, corrIndex(P,P), teldummy,&
                                           t1, t2, seed, nn
    integer, intent(in)                 :: P, samsize, numcorrgroup, Fisher
    real (8)                            :: corrMat(P,P),draw1(1),&
                                           R2inv(P,P), vec1(P,1), vec2(P,1),&
                                           dummy11(1,1), dummy12(1,1), dummy22(1,1),&
                                           Di1i2, &
                                           preinv(P,P)
    integer, allocatable, dimension(:)  :: iseed
    real(8), intent (out)               :: drawsCorr(samsize,numcorrgroup)
    real(4)                             :: alpha

!========================================================================================!
    
    !set seed
    call RANDOM_SEED(size=nn)
    allocate(iseed(nn))
    iseed(:)=seed
    call RANDOM_SEED(put=iseed)
    
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
            draw1=draw1*2.0-1.0
            corrMat(r1,r1+1) = draw1(1)
            corrMat(r1+1,r1) = corrMat(r1,r1+1)
            drawsCorr(s1,corrIndex(r1+1,r1)) = corrMat(r1,r1+1)
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
                R2inv(1:(i2-i1-1),1:(i2-i1-1)) = inverse(preinv(1:(i2-i1-1),1:(i2-i1-1)),(i2-i1-1))
                dummy11 = matmul(matmul(transpose(vec1(1:(k1-1),:)),R2inv(1:(i2-i1-1),1:(i2-i1-1))),vec1(1:(k1-1),:))
                dummy22 = matmul(matmul(transpose(vec2(1:(k1-1),:)),R2inv(1:(i2-i1-1),1:(i2-i1-1))),vec2(1:(k1-1),:))
                dummy12 = matmul(matmul(transpose(vec1(1:(k1-1),:)),R2inv(1:(i2-i1-1),1:(i2-i1-1))),vec2(1:(k1-1),:))
                Di1i2 = sqrt((1-dummy11(1,1))*(1-dummy22(1,1)))
                
                corrMat(i1,i2) = dummy12(1,1) + Di1i2*draw1(1)
                corrMat(i2,i1) = corrMat(i1,i2)
                
                drawsCorr(s1,corrIndex(i1,i2)) = corrMat(i1,i2)
                end do
            end do
        end do
    if(Fisher==1) then
        drawsCorr(1:samsize,1:numcorrgroup) = .5*log((1.0+drawsCorr(1:samsize,1:numcorrgroup)) &
                                        /(1.0-drawsCorr(1:samsize,1:numcorrgroup)))
    end if
contains



function inverse(a,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================

integer n
double precision a(n,n), inverse(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    inverse(i,k) = x(i)
  end do
  b(k)=0.0
end do
end 




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

    REAL, INTENT(IN)    :: aa, bb
    LOGICAL, INTENT(IN) :: first
    REAL                :: fn_val

    !     Local variables
    REAL, PARAMETER  :: aln4 = 1.3862944, one=1.0, two=2.0, vlarge = HUGE(1.0), vsmall = TINY(1.0), zero = 0.0
    REAL             :: a, b, g, r, s, x, y, z
    REAL, SAVE       :: d, f, h, t, c
    LOGICAL, SAVE    :: swap



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
    CALL RANDOM_NUMBER(r)
    CALL RANDOM_NUMBER(x)
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

end subroutine draw_ju




subroutine compute_rcet(numE,drawsIn,wIn,delta,rcEt,samsize)
!estimates the density at 0 via a histogram estimate, e.g., mean(abs(draws)<delta)/(2*delta), with default delta=.1
    
    implicit none
    ! numE number of eqully constraint for numEQ(h1)
    ! modified drawscorr matrix

    integer, intent(in)             :: numE, samsize
    real(8), intent(in)             :: drawsIn(samsize,numE), wIn(numE), delta
    real(8), intent(out)            :: rcEt
    integer                         :: c1, i1
    real                            :: dummyvec1(samsize), checkvec1(samsize)

    dummyvec1 = 1
    do c1=1,numE
        checkvec1 = 1
        do i1=1,samsize
            if(abs(drawsIn(i1,c1)-wIn(c1))>delta) then
                checkvec1(i1) = 0
            end if
        end do
        dummyvec1 = dummyvec1 * checkvec1
    end do
    rcEt = 1.0/(2*delta)**real(numE)*sum(dummyvec1)/real(samsize)

end subroutine compute_rcet



subroutine compute_rcet2(numE,drawsIn,wIn,delta,rcEt,meanOut,covmOut,samsize,numcorr)
!estimates the density at 0 via a histogram estimate, e.g., mean(abs(draws)<delta)/(2*delta), with default delta=.1.
!and compute conditional mean and covariance matrix of unconstrained parameters.    
    implicit none
    
    integer, intent(in)                 :: numE, samsize, numcorr
    real(8), intent(in)                 :: drawsIn(samsize,numcorr), wIn(numE), delta
    real(8), intent(out)                :: rcEt, meanOut(numcorr-numE), covmOut(numcorr-numE,numcorr-numE)
    integer                             :: c1, i1, check1, tel1,m1
    real(8)                             :: dummyvec1(samsize), drawsIE(samsize,numcorr), &
                                           meanDummy(numcorr,1), covmDummy(numcorr,numcorr)
    real(8), allocatable                :: diffs(:,:), ones(:,:)
    
    tel1 = 0
    dummyvec1 = 1
    check1 = 0
    do i1=1,samsize
        check1 = 1
        do c1=1,numE
            if(abs(drawsIn(i1,c1)-wIn(c1))>delta) then
                check1 = 0
                exit
            end if
        end do
        if(check1==1) then
            tel1 = tel1 + 1
            drawsIE(tel1,:) = drawsIn(i1,:)
        end if
    end do
    allocate(diffs(tel1,numcorr),ones(tel1,1))
    rcEt = 1.0/(2*delta)**real(numE)*real(tel1)/real(samsize)
    do m1=1, numcorr
        meanDummy(m1,1)=sum(drawsIE(1:tel1,m1))/tel1
    end do 
    meanOut(:) = meanDummy((numE+1):numcorr,1)
    ones = 1
    diffs = drawsIE(1:tel1,1:numcorr) - matmul(ones,transpose(meanDummy))
    covmDummy(1:numcorr,1:numcorr) = matmul(transpose (diffs), diffs)/real(tel1)
    !covmDummy = covmDummy/real(tel1)
    covmOut(1:(numcorr-numE),1:(numcorr-numE)) = covmDummy((numE+1):numcorr,(numE+1):numcorr)
!    
end subroutine compute_rcet2





