!! Caspar: For now, I have commented out all code, because it is buggy
!!         and prevents the R-package from loading. After debugging, see
!!         the fortran_functions.R file for how to write wrappers for Fortran
!!         functions.
!!
!!   FORTRAN 95 syntax for the BCT method (Mulder & Gelissen, 2019) for testing equality and order constrained
!!   Hypotheses on measures of association.
!!
!!
!module global_variables
!!
!    implicit none
!!
!    !Declare global variables
!    real, allocatable ::    Rcon(:,:), CONin(:,:), YXin(:,:), Xgroups(:,:,:), BHat(:,:,:), sdHat(:,:),&
!                            CHat(:,:,:), SS(:,:,:), SigmaHat(:,:,:), Ygroups(:,:,:), XtXi(:,:,:), &
!                            Rcomplement(:,:)
!    integer ::              P, numG, Ntot, numcorr, jointunif, margunif, samsize, samsize0, conTot, &
!                            K, numH, intercept, numcorrgroup
!    integer, allocatable :: Njs(:), numIneq(:), numEq(:)
!!
!end module global_variables
!
!
!program BCT
!!
!    use IMSL_LIBRARIES
!    use global_variables
!!
!    implicit none
!!
!    ! Variables
!    integer ::              iseed, i1, i2, index1, index2, c1, telcorr, priorDraws, g1, same1, s1, &
!                            p1, j1, categorie, maxCat, simulation, m1, iterPostEst, testValue, checkDependent, &
!                            priorApprox, beginning0, numIneqTot, numEqTot, h1, teldummy, r1, r2, K1, &
!                            h, i, numER, numIR, j, rowrank, whichH, intdummy, header, tellerComp, perpteller
!    integer, allocatable :: tellerj(:), whichOrdinal(:), numCat(:,:), numf(:)
!    real ::                 ct, ft, ctApprox, alpha1, beta1, sdBeta, ctprop, ftprop, delta, ffe, pi, dummy6(6), &
!                            teller_c, dumE
!    real, allocatable ::    priorMean(:), priorCov(:,:), relComp(:), drawsJU(:,:), &
!                            rcE(:,:), rcI(:,:), Rh1(:,:), drawsTrans(:,:), wh1(:,:), statX(:,:), means(:,:), &
!                            covm(:,:), diffs(:,:), covmAll(:,:), meansAll(:,:), transMatrix(:,:), &
!                            meanI(:), covmI(:,:), wh1IE(:,:), Rh1IE(:,:), postmeanJU(:,:), postmeanMU(:,:), &
!                            postcovJU(:,:), postcovMU(:,:), diff1(:,:), DummyPP(:,:), transR(:,:), &
!                            rfE(:,:), rfI(:,:), BFtu(:,:), theta(:), thetacov(:,:), ff(:), &
!                            bet(:), betacov(:,:), tempmatrix(:,:),invbetacov(:,:), invDiag(:,:), IE(:,:), &
!                            B(:,:), constant(:), transcon(:), Rr(:,:), ERr1(:,:), IRr(:,:), ones(:,:), &
!                            RcomplementDummy(:,:), RconDummy(:,:), TransPerp(:,:), eyeC(:,:), dummyRow(:), &
!                            transMatrix1(:,:), BFmatrix(:,:), evaldraws(:,:), hits(:), meanX(:,:), sdX(:)
!    character ::            charH(1)
!!
!800 format(20F8.5)
!900 format(20F8.4)
!!
!    !open input files
!    open(10,file='BCT_input.txt')
!    open(20,file='data.txt')
!    open(30,file='BCT_output_relComp.txt',status='replace')
!    open(40,file='BCT_output_relFit.txt',status='replace')
!    open(50,file='BCT_output.txt',status='replace')
!!
!    read(10,*)
!    read(10,*)
!    read(10,*)P, K1, intercept, numG, Ntot, header
!    K = K1 + intercept
!    allocate(whichOrdinal(P))
!    numcorr = .5*P*(P-1)*numG
!    numcorrgroup = .5*P*(P-1)
!    allocate(transMatrix(numcorr,numcorr),postcovMU(numcorr,numcorr),postcovJU(numcorr,numcorr),&
!             postmeanJU(numcorr,1),postmeanMU(numcorr,1))
!    read(10,*)
!    read(10,*)
!    read(10,*)whichOrdinal
!    read(10,*)
!    read(10,*)
!    read(10,*)
!    read(10,*)numH
!    read(10,*)
!    read(10,*)
!    allocate(numEq(numH+1),numIneq(numH+1),relComp(numH),rcE(numH+1,1),rcI(numH+1,1),rfE(numH+1,1),&
!            rfI(numH+1,1),BFtu(numH+1,1))
!    numEq = 0
!    numIneq = 0
!    do i1 = 1,numH
!        read(10,*)numEq(i1),numIneq(i1)
!    end do
!    numEqTot = sum(numEq)
!    numIneqTot = sum(numIneq)
!    conTot = numEqTot + numIneqTot
!    allocate(CONin(conTot,6), Rcon(conTot,numcorr+1), RcomplementDummy(numIneqTot,numcorr+1))
!    read(10,*)
!    read(10,*)
!    read(10,*)
!    Rcon = 0
!    RcomplementDummy = 0
!    tellerComp = 0
!    !whichCorrTest = 0
!    h1 = 1
!    whichH = numEq(h1)+numIneq(h1)
!    do i1 = 1,conTot
!        read(10,*)CONin(i1,:)
!        if(CONin(i1,2) < CONin(i1,3)) then
!            intdummy = CONin(i1,2)
!            CONin(i1,2) = CONin(i1,3)
!            CONin(i1,3) = intdummy
!        end if
!        if(CONin(i1,4)==0) then
!            index1 = .5*(CONin(i1,2)-1.0)*(CONin(i1,2)-2.0) + CONin(i1,3) + numcorrgroup*(CONin(i1,1)-1.0)
!            Rcon(i1,index1) = CONin(i1,5)
!            Rcon(i1,numcorr+1) = CONin(i1,6)*CONin(i1,5)
!!            whichCorrTest(index1) = 1
!        else
!	    if(CONin(i1,5) < CONin(i1,6)) then
!            	intdummy = CONin(i1,5)
!            	CONin(i1,5) = CONin(i1,6)
!            	CONin(i1,6) = intdummy
!            end if
!            index1 = .5*(CONin(i1,2)-1.0)*(CONin(i1,2)-2.0) + CONin(i1,3) + numcorrgroup*(CONin(i1,1)-1.0)
!            index2 = .5*(CONin(i1,5)-1.0)*(CONin(i1,5)-2.0) + CONin(i1,6) + numcorrgroup*(CONin(i1,4)-1.0)
!            Rcon(i1,index1) = 1
!            Rcon(i1,index2) = -1
!        end if
!        !collect inequality constraints for complement hypothesis
!        if(numEq(h1)==0) then
!            tellerComp = tellerComp + 1
!            RcomplementDummy(tellerComp,:) = Rcon(i1,:)
!        end if
!        if(i1==whichH) then
!            read(10,*)
!            h1 = h1 + 1
!            whichH = whichH + numEq(h1) + numIneq(h1)
!        end if
!    end do
!    !construct matrix with inequality constraints for complement hypothesis
!    numIneq(numH+1) = tellerComp
!    allocate(Rcomplement(tellerComp,numcorr+1))
!    Rcomplement(:,:) = RcomplementDummy(1:tellerComp,:)
!
!    !transform the constant vector to Fisher-transformed scale
!    Rcon(:,numcorr+1) = .5*log((1.0+Rcon(:,numcorr+1))/(1.0-Rcon(:,numcorr+1)))
!    Rcomplement(:,numcorr+1) = .5*log((1.0+Rcomplement(:,numcorr+1))/(1.0-Rcomplement(:,numcorr+1)))
!!
!    jointunif = 1
!    margunif = 0
!    read(10,*)
!    read(10,*)
!    read(10,*)iseed, samsize, samsize0
!    delta = .1
!    allocate(drawsJU(samsize,numcorr),ones(samsize,1), statX(15,numcorrgroup), means(numcorrgroup,1), &
!             covm(numcorrgroup,numcorrgroup), diffs(samsize,numcorrgroup),covmAll(numcorr,numcorr), &
!             meansAll(numcorr,1))
!    ones = 1
!!
!    !set iseed
!    call rnset(iseed)
!!
!    rcE = 1
!    rcI = 1
!!
!    !use joint uniform prior for the correlations in the correlation matrices
!    drawsJU = 0
!    call draw_JU(drawsJU(1:samsize,1:numcorrgroup))
!    drawsJU(1:samsize,1:numcorrgroup) = .5*log((1.0+drawsJU(1:samsize,1:numcorrgroup)) &
!                                        /(1.0-drawsJU(1:samsize,1:numcorrgroup)))
!    teldummy = numcorrgroup
!    if(numG>1) then
!        do g1=2,numG
!            drawsJU(1:(samsize-g1+1),(teldummy+1):(teldummy+numcorrgroup)) = drawsJU(g1:samsize,1:numcorrgroup)
!            drawsJU((samsize-g1+2):samsize,(teldummy+1):(teldummy+numcorrgroup)) = drawsJU(1:(g1-1),1:numcorrgroup)
!            teldummy = teldummy + numcorrgroup
!        end do
!    end if
!    call uvsta(drawsJU(1:samsize,1:numcorrgroup),statX)
!    means(:,1) = statX(1,:)
!    diffs = drawsJU(1:samsize,1:numcorrgroup) - matmul(ones,transpose(means))
!    covm = 0
!    covm(1:numcorrgroup,1:numcorrgroup) = (diffs.tx.diffs)/samsize
!    covmAll = 0
!    do g1=1,numG
!        covmAll(((g1-1)*numcorrgroup+1):(g1*numcorrgroup),(numcorrgroup*(g1-1)+1):(g1*numcorrgroup)) = &
!            covm(1:numcorrgroup,1:numcorrgroup)
!        meansAll(((g1-1)*numcorrgroup+1):(g1*numcorrgroup),1) = means(1:numcorrgroup,1)
!    end do
!!
!    !Allocate estimates of parameters (theta) and their covariance matrix (thetacov)
!    allocate(theta(numcorr+1))
!    allocate(thetacov(numcorr,numcorr))
!
!    theta(1:numcorr) = meansAll(1:numcorr,1)
!    theta(numcorr+1) = -1  !We add an extra -1 in theta so that Rtheta>0 --> [R|r]theta>0, where [R|r] is an augmented matrix.
!!
!    thetacov = covmAll
!!
!    teldummy = 0
!    if(numcorr==1) then
!        do h1=1,numH
!            if(numIneq(h1)>0) then
!                allocate(Rh1(numIneq(h1),numcorr),drawsTrans(samsize,numIneq(h1)),wh1(samsize,numIneq(h1)),hits(samsize))
!                wh1 = 0
!                Rh1(1:numIneq(h1),1:numcorr) = Rcon((teldummy+1):(teldummy+numIneq(h1)),1:numcorr)
!                drawsTrans = matmul(drawsJU,transpose(Rh1))
!                hits = 1
!                do s1=1,samsize
!                    do c1=1,numIneq(h1)
!                        if(drawsTrans(s1,c1) < Rcon(c1,numcorr+1)) then
!                            hits(s1) = 0
!                            exit
!                        end if
!                    end do
!                end do
!                rcI(h1,1) = sum(hits)/samsize
!!
!                deallocate(Rh1,drawsTrans,wh1,hits)
!                teldummy = teldummy + numIneq(h1)
!            end if
!            if(numEq(h1)>0) then
!                rcE(h1,1) = 2*exp(2*Rcon(teldummy+1,numcorr+1)) / (exp(2*Rcon(teldummy+1,numcorr+1))+1)**2
!                teldummy = teldummy + numEq(h1)
!            end if
!        end do
!    end if
!    if(numcorr>1) then
!        do h1=1,numH
!            if(numEq(h1)>0 .and. numIneq(h1)==0) then
!                !compute rel.comp. of E
!                allocate(Rh1(numEq(h1),numcorr),drawsTrans(samsize,numEq(h1)),wh1(numEq(h1),1))
!                Rh1(1:numEq(h1),1:numcorr) = Rcon((teldummy+1):(teldummy+numEq(h1)),1:numcorr)
!                wh1(:,1) = Rcon((teldummy+1):(teldummy+numEq(h1)),numcorr+1)
!                drawsTrans = matmul(drawsJU,transpose(Rh1))
!                call compute_rcEt(numEq(h1),drawsTrans,wh1(:,1),delta,rcE(h1,1))
!                deallocate(Rh1,drawsTrans,wh1)
!                teldummy = teldummy + numEq(h1)
!            end if
!            if(numIneq(h1)>0 .and. numEq(h1)==0) then
!                !compute rel.comp. of I
!                allocate(Rh1(numIneq(h1),numcorr+1),evaldraws(samsize,numIneq(h1)),hits(samsize))
!                Rh1(1:numIneq(h1),1:(numcorr+1)) = Rcon((teldummy+1):(teldummy+numIneq(h1)),1:(numcorr+1))
!                evaldraws = matmul(drawsJU,transpose(Rh1(:,1:numcorr))) ! - Rh1(:,numcorr+1)
!                hits = 1
!                do s1=1,samsize
!                    do c1=1,numIneq(h1)
!                        if(evaldraws(s1,c1) < Rh1(c1,numcorr+1)) then
!                            hits(s1) = 0
!                            exit
!                        end if
!                    end do
!                end do
!                rcI(h1,1) = sum(hits)/samsize
!                teldummy = teldummy + numIneq(h1)
!                deallocate(Rh1,evaldraws,hits)
!            end if
!            if(numIneq(h1)>0 .and. numEq(h1)>0) then
!                !first compute rel.comp. of E
!                allocate(drawsTrans(samsize,numcorr),wh1(numEq(h1),1), TransPerp(numcorr,numcorr), &
!                    meanI(numcorr-numEq(h1)),covmI(numcorr-numEq(h1),numcorr-numEq(h1)), &
!                    eyeC(numcorr,numcorr), dummyRow(numcorr), transMatrix1(numEq(h1),numcorr))
!                transMatrix = 0
!                transMatrix(1:numEq(h1),1:numcorr) = Rcon((teldummy+1):(teldummy+numEq(h1)),1:numcorr)
!                transMatrix1(1:numEq(h1),1:numcorr) = Rcon((teldummy+1):(teldummy+numEq(h1)),1:numcorr)
!                wh1(:,1) = Rcon((teldummy+1):(teldummy+numEq(h1)),numcorr+1)
!    !
!                !add independent vectors to make the transformation one-to-one
!                eyeC = 0
!                do r1=1,numcorr
!                    eyeC(r1,r1) = 1
!                end do
!                TransPerp = eyeC - matmul(transpose(transMatrix1), matmul(.i.(matmul(transMatrix1, &
!                    transpose(transMatrix1))), transMatrix1))
!                transMatrix(numEq(h1)+1,1:numcorr) = TransPerp(1,1:numcorr)
!                perpteller = 2
!                do r1=2,numcorr
!                    dummyRow(1:numcorr) = TransPerp(r1,1:numcorr)
!                    same1 = 0
!                    do r2=1,r1-1
!                        if(sum(abs(TransPerp(r2,1:numcorr) - dummyRow(1:numcorr))) < .00001) then
!                            same1 = 1
!                        end if
!                    end do
!                    if(same1 == 0) then
!                        transMatrix(numEq(h1)+perpteller,1:numcorr) = dummyRow(1:numcorr)
!                        perpteller = perpteller + 1
!                    end if
!                end do
!
!                drawsTrans = matmul(drawsJU,transpose(transMatrix))
!                call compute_rcEt2(numEq(h1),drawsTrans,wh1(:,1),delta,rcE(h1,1),meanI,covmI)
!                teldummy = teldummy + numEq(h1)
!    !
!                !now compute rel.comp. of I|E
!                allocate(Rh1(numIneq(h1),numcorr),wh1IE(numIneq(h1),1),Rh1IE(numIneq(h1),numcorr-numEq(h1)+1))
!                Rh1(1:numIneq(h1),1:numcorr) = Rcon((teldummy+1):(teldummy+numIneq(h1)),1:numcorr)
!                !compute inequality restriction matrix based on transformed parameters.
!                Rh1 = Rh1.xi.transMatrix
!                wh1IE = Rh1(:,1:numEq(h1)).x.wh1(:,:)
!                Rh1IE(1:numIneq(h1),1:(numcorr-numEq(h1))) = Rh1(1:numIneq(h1),(numEq(h1)+1):numcorr)
!                Rh1IE(1:numIneq(h1),numcorr-numEq(h1)+1) = Rcon((teldummy+1):(teldummy+numIneq(h1)),numcorr+1) - &
!                    wh1IE(:,1)
!		! The Bain function (Gu, Hoijtink, Mulder, Rosseel, 2019) is used for computing the probability
!                call compute_prob_bain(0,numIneq(h1),numcorr-numEq(h1),Rh1IE(1:0,:),Rh1IE,meanI,covmI,dumE,rcI(h1,1))
!
!                teldummy = teldummy + numIneq(h1)
!                deallocate(Rh1,wh1IE,Rh1IE,drawsTrans,wh1,meanI,covmI)
!            end if
!        end do
!    end if
!
!    !compute relative complexity for the complement hypothesis
!    if(tellerComp > 0) then
!        do h1=1,numH
!            if(numIneq(h1)>0 .and. numEq(h1)==0) then
!                rcI(numH+1,1) = rcI(numH+1,1) - rcI(h1,1)
!            end if
!        end do
!    end if
!
!    write(*,*)'JU prior done'
!    write(*,*)'relative complexities for the hypotheses'
!    write(*,800)rcE(1:(numH+1),1)*rcI(1:(numH+1),1)
!    write(*,*)
!!
!    write(30,*)
!    write(30,*)
!    write(30,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(30,*)'!                                                                                   !'
!    write(30,*)'!                                   OUTPUT BCT                                      !'
!    write(30,*)'!                                                                                   !'
!    write(30,*)'!     When using this software please refer to Mulder & Gelissen. Bayes factor      !'
!    write(30,*)'!     testing of equality and order constraints on measures of association in       !'
!    write(30,*)'!     social research (under review).                                               !'
!    write(30,*)'!                                                                                   !'
!    write(30,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(30,*)
!    write(30,*)
!    write(30,*)'relative complexity using the joint uniform prior'
!    write(30,*)
!    write(30,*)'rc      rcE     rcI'
!    write(30,*)
!    do h1=1,numH
!        write(30,"(A)",advance='no')' Hypothesis H'
!        write(30,"(I)")h1
!        write(30,800)rcE(h1,1)*rcI(h1,1),rcE(h1,1),rcI(h1,1)
!        write(30,*)
!    end do
!    if(rcI(numH+1,1)>0) then
!        write(30,"(A)",advance='no')' Complement hypothesis (H'
!        write(30,"(I)",advance='no')numH+1
!        write(30,"(A)")')*'
!        write(30,800)rcE(numH+1,1)*rcI(numH+1,1),rcE(numH+1,1),rcI(numH+1,1)
!        write(30,*)
!        write(30,*)'* To include the complement hypothesis it is important that'
!        write(30,*)'the hypotheses with only inequality constraints are nonnested.'
!    else
!        write(30,*)
!        write(30,*)'  Note that the complement hypothesis is not included because the'
!        write(30,*)'  order-constrained hypotheses cover the complete parameter space.'
!    end if
!!
!    !read data
!    allocate(YXin(Ntot,K1+P+1),Njs(numG),Ygroups(numG,Ntot,P),Xgroups(numG,Ntot,K),tellerj(numG),numCat(numG,P),&
!             BHat(numG,K,P),sdHat(numG,P),CHat(numG,P,P),XtXi(numG,K,K),SS(numG,P,P),SigmaHat(numG,P,P),&
!             DummyPP(P,P))
!    Njs = 0
!    tellerj = 1
!    Ygroups = 1
!    Xgroups = 1
!    YXin = 999
!    if(header .ne. 0) then
!        read(20,*)
!    end if
!    if(numG>1) then
!        do i1 = 1,Ntot
!            read(20,*)YXin(i1,:)
!            do i2 = 1,numG
!                if(YXin(i1,K1+P+1)==i2) then
!                    Njs(i2) = Njs(i2) + 1
!                    Ygroups(i2,tellerj(i2),1:P) = YXin(i1,1:P)
!                    Xgroups(i2,tellerj(i2),(intercept+1):(intercept+K1)) = YXin(i1,(P+1):(P+K1))
!                    tellerj(i2) = tellerj(i2)+1
!                end if
!            end do
!        end do
!    else
!        Njs = Ntot
!        do i1 = 1,Ntot
!            read(20,*)YXin(i1,1:(K1+P))
!            Ygroups(1,i1,1:P) = YXin(i1,1:P)
!            Xgroups(1,i1,(intercept+1):(intercept+K1)) = YXin(i1,(P+1):(P+K1))
!        end do
!    end if
!
!    if(K1 > 0) then !standardize covariates
!        deallocate(statX,ones)
!        allocate(statX(15,K1),meanX(1,K1),sdX(K1))
!        do g1 = 1,numG
!            allocate(ones(Njs(g1),1))
!            ones = 1
!            call uvsta(Xgroups(g1,1:Njs(g1),(intercept+1):(intercept+K1)),statX)
!            meanX(1,:) = statX(1,:)
!            sdX(:) = statX(3,:)
!            Xgroups(g1,1:Njs(g1),(intercept+1):(intercept+K1)) = Xgroups(g1,1:Njs(g1),(intercept+1):(intercept+K1)) &
!                                                                    - matmul(ones,meanX)
!            Xgroups(g1,1:Njs(g1),(intercept+1):(intercept+K1)) = matmul(Xgroups(g1,1:Njs(g1),&
!                                                                    (intercept+1):(intercept+K1)),diag(1/sdX(:)))
!            deallocate(ones)
!        end do
!    end if
!!
!    !compute ML estimates
!    do g1=1,numG
!        allocate(diff1(Njs(g1),P))
!        XtXi(g1,:,:) = .i.(Xgroups(g1,1:Njs(g1),1:K).tx.Xgroups(g1,1:Njs(g1),1:K))
!        BHat(g1,:,:) = XtXi(g1,:,:).xt.Xgroups(g1,1:Njs(g1),1:K).x.Ygroups(g1,1:Njs(g1),1:P)
!        diff1 = Ygroups(g1,1:Njs(g1),1:P) - (Xgroups(g1,1:Njs(g1),1:K).x.BHat(g1,1:K,1:P))
!        SS(g1,:,:) = diff1.tx.diff1
!        SigmaHat(g1,:,:) = SS(g1,:,:)/(Njs(g1)-K)
!        sdHat(g1,:) = sqrt(diagonals(SigmaHat(g1,:,:)))
!        DummyPP(:,:) = diag(1/sdHat(g1,:))
!        CHat(g1,:,:) = DummyPP(:,:).x.SigmaHat(g1,:,:).x.DummyPP(:,:)
!        deallocate(diff1)
!    end do
!!!
!    !estimate mean and covariance matrix of Fisher transformed correlations
!    if(sum(whichOrdinal)==0) then !all DVs are continuous
!        call estimate_postMeanCov_FisherZ(postmeanJU, postcovJU, 1)
!        write(*,*)'posterior has been estimated based on the JU-prior'
!        write(*,*)
!    else
!        numCat = 0
!        do p1=1,P
!            if(whichOrdinal(p1)>0) then
!                do g1=1,numG
!                    numCat(g1,p1) = maxval(Ygroups(g1,1:Njs(g1),p1))
!                end do
!            end if
!        end do
!        maxCat = maxval(numCat)
!!
!        write(*,*)
!        write(*,*)'Posterior sampling...'
!        write(*,*)
!        call estimate_postMeanCov_FisherZordinal(whichOrdinal, numCat, maxCat, postmeanJU, postcovJU, 1)
!        write(*,*)
!        write(*,*)'...is finished'
!        write(*,*)
!!
!    end if
!!
!!   the computation of the relative fit is an adaptation of a function
!!   in BaIn (Gu, Hoijtink, Mulder, Rosseel, 2019).
!!
!    !compute relative fit based on posterior under the JU-prior
!    rfE = 1
!    rfI = 1
!!
!    theta(1:numcorr) = postmeanJU(1:numcorr,1)
!    theta(numcorr+1) = -1  !We add an extra -1 in theta so that Rtheta>0 --> [R|r]theta>0, where [R|r] is an augmented matrix.
!!
!    thetacov = postcovJU
!!
!    teldummy = 0
!    do h = 1,numH ! numH means the number of hypotheses under consideration.
!        numER = numEq(h)
!        numIR = numIneq(h)
!        call compute_prob_bain(numEq(h),numIneq(h),numcorr,Rcon((teldummy+1):(teldummy+numER),1:(numcorr+1)),&
!            Rcon((teldummy+numER+1):(teldummy+numER+numIR),1:(numcorr+1)),theta,thetacov,rfE(h,1),rfI(h,1))
!        teldummy = teldummy + numER + numIR
!    end do
!!
!    !compute relative fit for the complement hypothesis
!    if(tellerComp > 0) then
!        do h1=1,numH
!            if(numIneq(h1)>0 .and. numEq(h1)==0) then
!                rfI(numH+1,1) = rfI(numH+1,1) - rfI(h1,1)
!            end if
!        end do
!    end if
!!
!    write(40,*)
!    write(40,*)
!    write(40,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(40,*)'!                                                                                   !'
!    write(40,*)'!                                   OUTPUT BCT                                      !'
!    write(40,*)'!                                                                                   !'
!    write(40,*)'!     When using this software please refer to Mulder & Gelissen. Bayes factor      !'
!    write(40,*)'!     testing of equality and order constraints on measures of association in       !'
!    write(40,*)'!     social research (under review).                                               !'
!    write(40,*)'!                                                                                   !'
!    write(40,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(40,*)
!    write(40,*)
!    write(40,*)'relative fit using the joint uniform prior and/or the marginally uniform prior'
!    write(40,*)
!    write(40,*)'rf      rfE     rfI'
!    write(40,*)
!    do h=1,numH
!        write(40,"(A)",advance='no')' Hypothesis H'
!        write(40,"(I)")h
!        write(40,800)rfE(h,1)*rfI(h,1),rfE(h,1),rfI(h,1)
!        write(40,*)
!    end do
!    if(rcI(numH+1,1)>0) then
!        write(40,"(A)",advance='no')' Complement hypothesis (H'
!        write(40,"(I)",advance='no')numH+1
!        write(40,"(A)")')*'
!        write(40,800)rfE(numH+1,1)*rfI(numH+1,1),rfE(numH+1,1),rfI(numH+1,1)
!        write(40,*)
!        write(40,*)'* To include the complement hypothesis it is important that'
!        write(40,*)'the hypotheses with only inequality constraints are nonnested.'
!    else
!        write(40,*)
!        write(40,*)'Note that the complement hypothesis is not included because the'
!        write(40,*)'order-constrained hypotheses cover the complete parameter space.'
!    end if
!!!
!    close(40)
!    close(20)
!    close(10)
!!!
!!    !compute BFs and PMPs based on compute relative complexities and relative fits
!    do h=1,numH+1
!        BFtu(h,1) = rfE(h,1)*rfI(h,1)/(rcE(h,1)*rcI(h,1))
!    end do
!    allocate(BFmatrix(numH+1,numH+1))
!    do h=1,numH+1
!        do h1=1,numH+1
!            BFmatrix(h,h1) = BFtu(h,1)/BFtu(h1,1)
!        end do
!    end do
!    charH = 'H'
!
!    write(50,*)
!    write(50,*)
!    write(50,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(50,*)'!                                                                                   !'
!    write(50,*)'!                                   OUTPUT BCT                                      !'
!    write(50,*)'!                                                                                   !'
!    write(50,*)'!     When using this software please refer to Mulder & Gelissen. Bayes factor      !'
!    write(50,*)'!     testing of equality and order constraints on measures of association in       !'
!    write(50,*)'!     social research (under review).                                               !'
!    write(50,*)'!                                                                                   !'
!    write(50,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(50,*)
!    write(50,*)
!    write(50,*)'Posterior probabilities for the hypotheses (assuming equal prior probabilities)'
!    write(50,*)
!    if(rcI(numH+1,1)>0) then
!        do h=1,numH
!            write(50,"(A)",advance='no')' Hypothesis H'
!            write(50,"(I)")h
!            write(50,900)BFtu(h,1)/sum(BFtu(:,1))!,BFtu(h,2)/sum(BFtu(:,2))
!            write(50,*)
!        end do
!!
!        write(50,"(A)",advance='no')' Complement hypothesis (H'
!        write(50,"(I)",advance='no')numH+1
!        write(50,"(A)")')*'
!        write(50,900)BFtu(numH+1,1)/sum(BFtu(:,1))!,BFtu(numH+1,2)/sum(BFtu(:,2))
!        write(50,*)
!        write(50,*)'* To include the complement hypothesis it is important that'
!        write(50,*)'the hypotheses with only inequality constraints are nonnested.'
!        write(50,*)
!        write(50,*)
!        write(50,*)
!        write(50,*)'Bayes factors between each pair of hypotheses'
!        write(50,*)
!        write(50,"(A)",advance='no')'    '
!        do h=1,numH
!            write(50,"(A)",advance='no')'H'
!            write(50,"(I)",advance='no')h
!            write(50,"(A)",advance='no')'      '
!        end do
!        write(50,"(A)",advance='no')'H'
!        write(50,"(I)")numH+1
!        do h=1,numH+1
!             write(50,"(A)",advance='no')'H'
!             write(50,"(I)",advance='no')h
!             write(50,900)BFmatrix(h,1:(numH+1))
!        end do
!    else
!        do h=1,numH
!            write(50,"(A)",advance='no')' Hypothesis H'
!            write(50,"(I)")h
!            write(50,900)BFtu(h,1)/sum(BFtu(1:numH,1))!,BFtu(h,2)/sum(BFtu(:,2))
!            write(50,*)
!        end do
!        write(50,*)
!        write(50,*)
!        write(50,*)'Bayes factors between each pair of hypotheses'
!        write(50,*)
!        write(50,"(A)",advance='no')'    '
!        do h=1,numH-1
!            write(50,"(A)",advance='no')'H'
!            write(50,"(I)",advance='no')h
!            write(50,"(A)",advance='no')'      '
!        end do
!        write(50,"(A)",advance='no')'H'
!        write(50,"(I)")numH
!        do h=1,numH
!             write(50,"(A)",advance='no')'H'
!             write(50,"(I)",advance='no')h
!             write(50,900)BFmatrix(h,1:numH)
!        end do
!    end if
!    close(50)
!
!end program BCT
!!
!!
!!
!!draw rho's from joint uniform distribution
!!use algorithm of Joe (2006)
!subroutine draw_JU(drawsCorr)
!!
!    use IMSL_LIBRARIES
!    use global_variables
!    implicit none
!!
!    !Declare local variables
!!    integer, intent (in) ::
!    real, intent (out) :: drawsCorr(samsize,numcorrgroup)
!    integer ::            r1, r2, s1, i1, i2, k1, corrIndex(P,P), teldummy
!    real ::               dummy1, corrMat(P,P), draw1(1), &
!                          alpha, R2inv(P,P), vec1(P,1), vec2(P,1), &
!                          dummy11(1,1), dummy12(1,1), dummy22(1,1), Di1i2
!!
!!    open(60,file='priordraws.txt',status='replace')
!!
!    teldummy = 1
!    do r1=2,P
!        do r2=1,r1-1
!            corrIndex(r1,r2) = teldummy
!            corrIndex(r2,r1) = teldummy
!            teldummy = teldummy + 1
!        end do
!    end do
!
!    do s1=1,samsize
!        corrMat = eye(P)
!        do r1 = 1,P-1
!            alpha = P/2.0
!            call rnbet(alpha,alpha,draw1)
!            corrMat(r1,r1+1) = draw1(1)*2.0-1.0
!            corrMat(r1+1,r1) = corrMat(r1,r1+1)
!            drawsCorr(s1,corrIndex(r1+1,r1)) = corrMat(r1,r1+1)
!        end do
!
!        R2inv(:,:) = 0
!        do r1 = 3,P
!            do r2 = 1,P-r1+1
!                i1 = r2
!                i2 = r2+r1-1
!                k1 = i2 - i1
!                !draw partial correlations
!                alpha = .5*(P+1-k1)
!                call rnbet(alpha,alpha,draw1)
!                draw1 = draw1*2-1.0
!                !rbeta(1,.5*(dim+1-k),.5*(dim+1-k))*2-1
!                vec1(1:(k1-1),1) = corrMat(i1,(i1+1):(i1+k1-1))
!                vec2(1:(k1-1),1) = corrMat(i2,(i1+1):(i1+k1-1))
!                R2inv(1:(i2-i1-1),1:(i2-i1-1)) = .i.(corrMat((i1+1):(i2-1),(i1+1):(i2-1)))
!                dummy11 = BLINF(R2inv(1:(i2-i1-1),1:(i2-i1-1)),vec1(1:(k1-1),1),vec1(1:(k1-1),1))
!                dummy22 = BLINF(R2inv(1:(i2-i1-1),1:(i2-i1-1)),vec2(1:(k1-1),1),vec2(1:(k1-1),1))
!                dummy12 = BLINF(R2inv(1:(i2-i1-1),1:(i2-i1-1)),vec1(1:(k1-1),1),vec2(1:(k1-1),1))
!                Di1i2 = sqrt((1-dummy11(1,1))*(1-dummy22(1,1)))
!                corrMat(i1,i2) = dummy12(1,1) + Di1i2*draw1(1)
!                corrMat(i2,i1) = corrMat(i1,i2)
!                drawsCorr(s1,corrIndex(i1,i2)) = corrMat(i1,i2)
!            end do
!        end do
!    end do
!!
!end subroutine draw_JU
!
!
!
!subroutine compute_rcEt(numE,drawsIn,wIn,delta,rcEt)
!!estimates the density at 0 via a histogram estimate, e.g., mean(abs(draws)<delta)/(2*delta), with default delta=.2
!
!    use IMSL_LIBRARIES
!    use global_variables
!    implicit none
!
!    integer, intent(in) :: numE
!    real, intent(in) ::    drawsIn(samsize,numE), wIn(numE), delta
!    real, intent(out) ::   rcEt
!    integer ::             c1, i1
!    real ::                dummyvec1(samsize), checkvec1(samsize)
!
!    dummyvec1 = 1
!    do c1=1,numE
!        checkvec1 = 1
!        do i1=1,samsize
!            if(abs(drawsIn(i1,c1)-wIn(c1))>delta) then
!                checkvec1(i1) = 0
!            end if
!        end do
!        dummyvec1 = dummyvec1 * checkvec1
!    end do
!    rcEt = 1.0/(2*delta)**real(numE)*sum(dummyvec1)/real(samsize)
!
!end subroutine compute_rcEt
!
!
!subroutine compute_rcEt2(numE,drawsIn,wIn,delta,rcEt,meanOut,covmOut)
!!estimates the density at 0 via a histogram estimate, e.g., mean(abs(draws)<delta)/(2*delta), with default delta=.2.
!!and compute mean and covariance matrix of unconstrained parameters.
!
!    use IMSL_LIBRARIES
!    use global_variables
!    implicit none
!
!    integer, intent(in) :: numE
!    real, intent(in) ::    drawsIn(samsize,numcorr), wIn(numE), delta
!    real, intent(out) ::   rcEt, meanOut(numcorr-numE), covmOut(numcorr-numE,numcorr-numE)
!    integer ::             c1, i1, check1, tel1
!    real ::                dummyvec1(samsize), checkvec1(samsize), drawsIE(samsize,numcorr), &
!                           statC(15,numcorr), meanDummy(numcorr,1), covmDummy(numcorr,numcorr)
!    real, allocatable ::   diffs(:,:), ones(:,:)
!
!    tel1 = 0
!    dummyvec1 = 1
!    check1 = 0
!    do i1=1,samsize
!        check1 = 1
!        do c1=1,numE
!            if(abs(drawsIn(i1,c1)-wIn(c1))>delta) then
!                check1 = 0
!                exit
!            end if
!        end do
!        if(check1==1) then
!            tel1 = tel1 + 1
!            drawsIE(tel1,:) = drawsIn(i1,:)
!        end if
!    end do
!    allocate(diffs(tel1,numcorr),ones(tel1,1))
!    rcEt = 1.0/(2*delta)**real(numE)*real(tel1)/real(samsize)
!
!    call uvsta(drawsIE(1:tel1,1:numcorr),statC)
!    meanDummy(:,1) = statC(1,:)
!    meanOut(:) = meanDummy((numE+1):numcorr,1)
!    ones = 1
!    diffs = drawsIE(1:tel1,1:numcorr) - matmul(ones,transpose(meanDummy))
!    covmDummy(1:numcorr,1:numcorr) = (diffs.tx.diffs)/real(tel1)
!    covmOut(1:(numcorr-numE),1:(numcorr-numE)) = covmDummy((numE+1):numcorr,(numE+1):numcorr)
!!
!end subroutine compute_rcEt2
!
!!
!!
!
!subroutine inverse_prob_sampling(condMean,condVar,LBtrue,UBtrue,LB,UB,condDraw)
!
!    use IMSL_LIBRARIES
!    use global_variables
!    implicit none
!
!    integer, intent(in) :: LBtrue, UBtrue
!    real, intent(in) :: condMean, condVar, LB, UB
!    real, intent(out) :: condDraw
!
!    real :: LBstand, UBstand, yUB, yLB, rnIPS, diffUBLB, bb, cc, Discr, xi, Zstar, lambda, pi, &
!            machPres
!    parameter(pi=3.141592653)
!!
!    machPres = 1e-6
!    !normalize bounds
!    UBstand = (UB - condMean)/sqrt(condVar)
!    LBstand = (LB - condMean)/sqrt(condVar)
!    yUB = anordf(UBstand)
!    yLB = anordf(LBstand)
!601     rnIPS = d_rnunf()*(anordf(UBstand) - anordf(LBstand)) + anordf(LBstand)
!    if(LBtrue+UBtrue==0) then !unconstrained sampling
!        condDraw = d_rnnof()*sqrt(condVar) + condMean
!    else if(abs(rnIPS) > machPres .and. abs(rnIPS-1) > machPres) then !inverse probability sampling
!        condDraw = anorin(rnIPS)*sqrt(condVar) + condMean
!    else if(UBstand>-2.33 .and. LBstand<2.33) then !IPS must be redone
!        go to 601
!    else
!        if(condMean > UB) then
!            condDraw = UB
!        else if(condMean < LB) then
!            condDraw = LB
!        end if
!    end if
!!
!end subroutine inverse_prob_sampling
!!
!
!
!subroutine estimate_postMeanCov_FisherZ(postZmean, postZcov, JUprior)
!!
!    use IMSL_LIBRARIES
!    use global_variables
!    implicit none
!
!    real, intent(out)   :: postZmean(numcorr,1), postZcov(numcorr,numcorr)
!    integer, intent(in) :: JUprior
!    real                :: BDraws(numG,K,P), sigmaDraws(numG,P), CDraws(numG,P,P), SSj(P,P), Wp(P,P),&
!                           cholSigmaPK(P*K,P*K), Ds(P,P), Ccan(P,P), CcanInv(P,P), Ccurr(P,P), CcurrInv(P,P),&
!                           SS1(P,P), rnunif(1), sdMH(numG,P), errorMatj(P,P), sigma_can(P), aa, bb, cc, dummy1(1), &
!                           f1, f0, f_1, UB, LB, alpha1, beta1, rr_can(1), RDrawsPrior(numG,P,P), &
!                           R_can(P,P), R_canInv(P,P), MHpart1, MHpart2, MHpart3, SigmaMat(P,P), R_MH, &
!                           epsteps(P,P), logR_MH_part1, logR_MH_part2, logR_MH_part3, logR_MH_part4, &
!                           varz1, varz2, varz1z2Plus, varz1z2Min, Cinv(P,P), acceptSigma(numG,P), dummyDraw1(samsize0), &
!                           detR, covBeta(P*K,P*K), betaDrawj(1,P*K), DummyPP(P,P), dummyDraw2(samsize0), &
!                           B_quantiles(numG,K,P,3), C_quantiles(numG,P,P,3), sigma_quantiles(numG,P,3), &
!                           Cnugget(P,P), c21
!    integer             :: s1, g1, rankP, acceptC(numG), i1, dummyI, testPr, rr1, rr2, nutarget, a0, iter1, corrteller,&
!                           c1, c2, p1, p2, k1
!    integer, allocatable:: iperm1(:)
!    real, allocatable   :: diffmat(:,:), ones(:,:),Zcorr_sample(:,:), dummy3(:), dummy2(:), &
!                           BDrawsStore(:,:,:,:), sigmaDrawsStore(:,:,:), CDrawsStore(:,:,:,:)
!!
!    allocate(ones(samsize0,1), Zcorr_sample(samsize0,numcorr), dummy3(samsize0), dummy2(samsize0), iperm1(samsize0),&
!             BDrawsStore(samsize0,numG,K,P), sigmaDrawsStore(samsize0,numG,P), CDrawsStore(samsize0,numG,P,P))
!!
!    !initial posterior draws
!    BDraws = BHat
!    sigmaDraws = sdHat
!    CDraws = CHat
!!
!800 format(20F8.3)
!!    open(60,file='posteriordraws1.txt',status='replace')
!    open(70,file='BCT_estimates.txt',status='replace')
!!
!    !initial prior C
!    nutarget = P + 1 !in the case of a marginally uniform prior
!    !specify candidate prior for R
!    Wp = eye(P)
!    a0 = nutarget !note that this implies an inv-W(Wp,a0) candidate prior.
!                  !the parameterization in Remark 9 of Liu&Daniels seems a bit odd.
!    Cnugget = Wp*1e-4
!!
!    !count number of accepted draws for R (over all groups)
!    acceptC = 0
!    acceptSigma = 0
!    sdMH = .1
!!
!    allocate(diffmat(Ntot,P))
!    !start Gibbs sampler
!    do s1 = 1,samsize0
!!
!        corrteller = 0
!        do g1 = 1,numG
!
!            !draw B
!            SigmaMat = matmul(matmul(diag(sigmaDraws(g1,:)),CDraws(g1,:,:)),diag(sigmaDraws(g1,:)))
!            call kronecker(K,P,XtXi(g1,:,:),SigmaMat,covBeta)
!            call chfac(covBeta,rankP,cholSigmaPK)
!            call rnmvn(cholSigmaPK,betaDrawj)
!            do p1=1,P
!                BDraws(g1,:,p1) = betaDrawj(1,((p1-1)*K+1):(p1*K)) + BHat(g1,:,p1)
!            end do
!
!            !draw C using method of Liu and Daniels (LD, 2006)
!            !draw candidate C from flat prior p(C)=1
!            diffmat(1:Njs(g1),1:P) = Ygroups(g1,1:Njs(g1),1:P) - matmul(Xgroups(g1,1:Njs(g1),1:K),BDraws(g1,1:K,1:P))
!            errorMatj = matmul(transpose(diffmat(1:Njs(g1),1:P)),diffmat(1:Njs(g1),1:P))
!            Ds = diag(1/sqrt(diagonals(errorMatj)))
!    	    diffmat(1:Njs(g1),1:P) = matmul(diffmat(1:Njs(g1),1:P),Ds) !diffmat is now epsilon in LD
!    	    epsteps = matmul(transpose(diffmat(1:Njs(g1),1:P)),diffmat(1:Njs(g1),1:P))
!    	    SS1 = matmul(matmul(diag(1/sigmaDraws(g1,:)),epsteps),diag(1/sigmaDraws(g1,:)))
!    	    SS1 = .i.SS1
!    	    call gen_wish(SS1,Njs(g1)-P-1,dummyPP)
!    	    dummyPP = .i.dummyPP
!    	    dummyPP = dummyPP !+ Cnugget
!    	    Ccan = matmul(matmul(diag(1/sqrt(diagonals(dummyPP))),dummyPP),diag(1/sqrt(diagonals(dummyPP))))
!    	    Ccan = Ccan !+ Cnugget
!    	    CcanInv = .i.Ccan
!    	    CDraws(g1,:,:) = Ccan(:,:) !always accepted because candidate was drawn from joint uniform prior
!    	    Cinv = CcanInv
!	        !End option 2
!!
!            do i1 = 2,P !keep Fisher z transformed posterior draws of rho's
!                Zcorr_sample(s1,(corrteller+1):(corrteller+i1-1)) = .5*log((1+CDraws(g1,i1,1:(i1-1)))/(1-CDraws(g1,i1,1:(i1-1))))
!                corrteller = corrteller + i1 - 1
!            end do
!!
!!            !draw sigma's
!            do p1 = 1,P
!                bb = sum(errorMatj(p1,:)*Cinv(p1,:)/sigmaDraws(g1,:))-errorMatj(p1,p1)*Cinv(p1,p1)/sigmaDraws(g1,p1)
!                aa = Cinv(p1,p1)*errorMatj(p1,p1)
!                sigma_can(:) = sigmaDraws(g1,:)
!                sigma_can(p1) = rnnof()
!                sigma_can(p1) = sigma_can(p1)*sdMH(g1,p1) + sigmaDraws(g1,p1) !random walk
!                R_MH = exp((-real(Njs(g1))-1.0)*(log(sigma_can(p1)-log(sigmaDraws(g1,p1)))) &
!                       -.5*aa*(sigma_can(p1)**(-2) - sigmaDraws(g1,p1)**(-2)) &
!                       -bb*(sigma_can(p1)**(-1) - sigmaDraws(g1,p1)**(-1)))
!                if(rnunif(1) < R_MH .and. sigma_can(p1)>0) then
!                    sigmaDraws(g1,p1) = sigma_can(p1)
!                    acceptSigma(g1,p1) = acceptSigma(g1,p1) + 1
!                end if
!            end do
!
!        end do
!!        write(60,800)CDraws(1,:,:)
!!
!        !store posterior draws
!        BDrawsStore(s1,1:numG,1:K,1:P) = BDraws(1:numG,1:K,1:P)
!        sigmaDrawsStore(s1,1:numG,1:P) = sigmaDraws(1:numG,1:P)
!        CDrawsStore(s1,1:numG,1:P,1:P) = CDraws(1:numG,1:P,1:P)
!
!        write(60,*)CDraws(1,2,1),CDraws(2,2,1),CDraws(3,2,1),CDraws(4,2,1)
!!
!    end do
!!
!    !compute medians Fisher transformed draws
!    do c1=1,numcorr
!        dummy2(:) = Zcorr_sample(:,c1)
!        call svrgp(dummy2,dummy3,iperm1)
!        postZmean(c1,1) = dummy3(int(samsize0*.5))
!    end do
!!
!    !robust covariance matrix estimation Fisher transformed draws
!    postZcov = 0
!    do c1=1,numcorr
!        do c2=c1,numcorr
!            call robust_covest(samsize0, Zcorr_sample(1:samsize0,c1), Zcorr_sample(1:samsize0,c2), &
!                postZmean(c1,1), postZmean(c2,1), varz1, varz2, varz1z2Plus, varz1z2Min)
!            postZcov(c1,c2) = (varz1*varz2)**.5 * (varz1z2Plus - varz1z2Min)/(varz1z2Plus + varz1z2Min)
!            postZcov(c2,c1) = postZcov(c1,c2)
!        end do
!    end do
!!
!    !compute medians, lower bounds CI, upperbounds CI
!    do g1=1,numG
!        C_quantiles(g1,:,:,1) = eye(P)
!        C_quantiles(g1,:,:,2) = eye(P)
!        C_quantiles(g1,:,:,3) = eye(P)
!        do p1=1,P
!            !for the sigma's
!            dummy2(:) = sigmaDrawsStore(:,g1,p1)
!            call svrgp(dummy2,dummy3,iperm1)
!            sigma_quantiles(g1,p1,1) = dummy3(int(samsize0*.025))
!            sigma_quantiles(g1,p1,2) = dummy3(int(samsize0*.5))
!            sigma_quantiles(g1,p1,3) = dummy3(int(samsize0*.975))
!            !for the beta coefficients
!            do k1=1,K
!                dummy2(:) = BDrawsStore(:,g1,k1,p1)
!                call svrgp(dummy2,dummy3,iperm1)
!                B_quantiles(g1,k1,p1,1) = dummy3(int(samsize0*.025))
!                B_quantiles(g1,k1,p1,2) = dummy3(int(samsize0*.5))
!                B_quantiles(g1,k1,p1,3) = dummy3(int(samsize0*.975))
!            end do
!            if(p1>1)then
!                do p2=1,p1-1
!                    dummy2(:) = CDrawsStore(:,g1,p1,p2)
!                    call svrgp(dummy2,dummy3,iperm1)
!                    C_quantiles(g1,p1,p2,1) = dummy3(int(samsize0*.025))
!                    C_quantiles(g1,p1,p2,2) = dummy3(int(samsize0*.5))
!                    C_quantiles(g1,p1,p2,3) = dummy3(int(samsize0*.975))
!                end do
!            end if
!        end do
!    end do
!!
!    !plot results
!        write(70,*)
!    write(70,*)
!    write(70,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(70,*)'!                                                                                   !'
!    write(70,*)'!                                   OUTPUT BCT                                      !'
!    write(70,*)'!                                                                                   !'
!    write(70,*)'!     When using this software please refer to Mulder & Gelissen. Bayes factor      !'
!    write(70,*)'!     testing of equality and order constraints on measures of association in       !'
!    write(70,*)'!     social research (under review).                                               !'
!    write(70,*)'!                                                                                   !'
!    write(70,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(70,*)
!    write(70,*)
!    write(70,*)'Estimates were obtained under the unconstrained model'
!    write(70,*)
!    write(70,*)'Correlation matrix'
!    write(70,*)
!    do g1 = 1,numG
!        write(70,*)'Population',g1
!        write(70,*)
!        write(70,*)'lower bound of 95%-CI'
!        do p1=1,P
!            write(70,800)C_quantiles(g1,p1,1:p1,1)
!        end do
!        write(70,*)
!        write(70,*)'median'
!        do p1=1,P
!            write(70,800)C_quantiles(g1,p1,1:p1,2)
!        end do
!        write(70,*)
!        write(70,*)'upper bound of 95%-CI'
!        do p1=1,P
!            write(70,800)C_quantiles(g1,p1,1:p1,3)
!        end do
!        write(70,*)
!        write(70,*)
!    end do
!    write(70,*)'B-matrix with intercepts and regression coefficients'
!    write(70,*)
!    do g1 = 1,numG
!        write(70,*)'Population',g1
!        write(70,*)
!        write(70,*)'lower bound of 95%-CI'
!        do k1=1,K
!            write(70,800)B_quantiles(g1,k1,1:P,1)
!        end do
!        write(70,*)
!        write(70,*)'median'
!        do k1=1,K
!            write(70,800)B_quantiles(g1,k1,1:P,2)
!        end do
!        write(70,*)
!        write(70,*)'upper bound of 95%-CI'
!        do k1=1,K
!            write(70,800)B_quantiles(g1,k1,1:P,3)
!        end do
!        write(70,*)
!        write(70,*)
!    end do
!    write(70,*)'standard deviations'
!    write(70,*)
!    do g1 = 1,numG
!        write(70,*)'Population',g1
!        write(70,*)
!        write(70,*)'lower bound of 95%-CI'
!        write(70,800)sigma_quantiles(g1,1:P,1)
!        write(70,*)
!        write(70,*)'median'
!        write(70,800)sigma_quantiles(g1,1:P,2)
!        write(70,*)
!        write(70,*)'upper bound of 95%-CI'
!        write(70,800)sigma_quantiles(g1,1:P,3)
!        write(70,*)
!        write(70,*)
!    end do
!!
!!
!end subroutine estimate_postMeanCov_FisherZ
!!
!!
!!
!!   Gibbs sampler to estimate posterior mean of Z scores when some DVs are ordinal.
!!
!subroutine estimate_postMeanCov_FisherZordinal(ordinal, Cat, maxCat, postZmean, postZcov, JUprior)
!!
!    use IMSL_LIBRARIES
!    use global_variables
!    implicit none
!
!    integer, intent(in) :: ordinal(P), Cat(numG,P), maxCat, JUprior
!    real, intent(out)   :: postZmean(numcorr,1), postZcov(numcorr,numcorr)
!    real :: BDraws(numG,K,P), CDraws(numG,P,P), sigmaDraws(numG,P), SSj(P,P), meanMat(Ntot,P), SigmaMatDraw(P,P), R_MH, &
!            cholSigmaPK(K*P,K*P), covBeta(K*P,K*P), Ds(P,P), Ccan(P,P), CcanInv(P,P), Ccurr(P,P), CcurrInv(P,P),&
!            SS1(P,P), rnunif(1), sdMH(numG,P), errorMatj(P,P), sigma_can(P), aa, bb, cc, dummy1(1), acceptSigma(numG,P), &
!            telft, telct, f1, f0, f_1, UB, LB, alpha1, beta1, rr_can(1), betaDrawj(1,P*K), &
!            C_can(P,P), C_canInv(P,P), bounds(2), MHpart1, MHpart2, MHpart3, dummyPP(P,P), &
!            Wp(P,P), epsteps(P,P), logR_MH_part1, logR_MH_part2, logR_MH_part3, logR_MH_part4, &
!            varz1, varz2, varz1z2Plus, varz1z2Min, Cnugget(P,P), SigmaInv(P,P), sdMHg(numG,P), gLiuSab_can, &
!            Wgroups(numG,Ntot,P), dummyVar, alphaMin, alphaMax, Cinv(P,P), Bmean(K,P), acceptLS(numG,P), &
!            alphaMat(numG,maxCat+1,P), Wdummy(numG,P,Ntot,maxCat), condMean, condVar, &
!            B_quantiles(numG,K,P,3), C_quantiles(numG,P,P,3), sigma_quantiles(numG,P,3), gLiuSab(numG,P)
!    real, allocatable :: diffmat(:,:), ones(:,:),Zcorr_sample(:,:), dummy3(:), dummy2(:), &
!                         BDrawsStore(:,:,:,:), sigmaDrawsStore(:,:,:), CDrawsStore(:,:,:,:)
!    integer           :: s1, g1, rankP, acceptC(numG), i1, dummyI, testPr, rr1, rr2, nutarget, a0, iter1, corrteller, &
!                         c1, c2, p1, Yi1Categorie, tellers(numG,maxCat,P), burnin, k1, k2, CatUpper, p2
!    integer, allocatable :: iperm1(:)
!!
!800 format(30F8.3)
!    open(60,file='posteriordraws1.txt',status='replace')
!    open(61,file='posteriordraws2.txt',status='replace')
!    open(62,file='posteriordraws3.txt',status='replace')
!!    open(63,file='posteriordraws4.txt',status='replace')
!    open(70,file='BCT_estimates.txt',status='replace')
!!
!    allocate(ones(samsize0,1),Zcorr_sample(samsize0,numcorr), dummy3(samsize0), dummy2(samsize0), iperm1(samsize0), &
!             BDrawsStore(samsize0,numG,K,P), sigmaDrawsStore(samsize0,numG,P), CDrawsStore(samsize0,numG,P,P))
!!
!    !initial posterior draws
!    BDraws = BHat
!    sigmaDraws = sdHat
!    CDraws = CHat
!    gLiuSab = 1
!!
!    !initial prior C
!    nutarget = P + 1 !in the case of a marginally uniform prior
!    do p1=1,P
!        !initial values
!        if(ordinal(p1)==1) then
!            sigmaDraws(:,p1) = 1
!            sigma_quantiles(:,p1,1) = 1
!            sigma_quantiles(:,p1,2) = 1
!            sigma_quantiles(:,p1,3) = 1
!            BDraws(:,:,p1) = 0
!        end if
!    end do
!    !specify candidate prior for R
!    Wp = eye(P)
!    a0 = nutarget !note that this implies an inv-W(Wp,a0) candidate prior.
!                  !the parameterization in Remark 9 of Liu&Daniels seems a bit odd.
!!
!    !define nugget matrix to avoid approximate nonpositive definite correlation matrices for candidates
!    Cnugget = Wp*1e-4
!!
!    telft = 0
!    telct = 0
!!
!    !count number of accepted draws for R (over all groups)
!    acceptC = 0
!    acceptSigma = 0
!    acceptLS = 0
!    sdMH = .2 !for error standard deviations
!    sdMHg = .1 !for gLiuBanhatti parameter
!    dummyVar = 1.0
!!
!    !initial values for latent W's corresponding to ordinal DVs
!    Wgroups = Ygroups
!    Wdummy = 0
!    ones = 1
!!
!    !initial values of boundary values alpha to link between ordinal Y and continuous latent W
!    alphaMat = 0
!    alphaMat(:,1,:) = -1e10  !alpha0
!    alphaMat(:,2,:) = 0      !alpha1
!    do p1=1,P
!        if(ordinal(p1)>0) then
!            do g1=1,numG
!                do c1=3,Cat(g1,p1)
!                    alphaMat(g1,c1,p1) = .3*(c1-2.0)
!                end do
!                alphaMat(g1,Cat(g1,p1)+1,p1) = 1e10
!            end do
!        end if
!    end do
!!
!    allocate(diffmat(Ntot,P))
!    burnin = 500
!!
!    write(*,*)'sampling for burn-in period'
!    !start Gibbs sampler
!    do s1 = 1,burnin
!        corrteller = 0
!        tellers = 0
!        do g1 = 1,numG
!            !compute means of latent W's for all observations
!            meanMat(1:Njs(g1),1:P) = matmul(Xgroups(g1,1:Njs(g1),1:K),BDraws(g1,1:K,1:P))
!            Ccurr = CDraws(g1,:,:)
!            SigmaMatDraw = matmul(matmul(diag(sigmaDraws(g1,:)),Ccurr),diag(sigmaDraws(g1,:)))
!!
!            !draw latent W's for the ordinal Y's
!            !compute mean vector for
!            do p1=1,P
!                if(ordinal(p1)>0) then
!                    do i1=1,Njs(g1)
!                        Yi1Categorie = Ygroups(g1,i1,p1)
!                        call compute_condMeanVar(p1,P,meanMat(i1,1:P),SigmaMatDraw,Wgroups(g1,i1,1:P),condMean,condVar)
!                        select case (Yi1Categorie)
!                            case(1)
!                                call inverse_prob_sampling(condMean,condVar,0,1,alphaMat(g1,1,p1), &
!                                    alphaMat(g1,2,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,1,p1) = tellers(g1,1,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,1,p1),1) = Wgroups(g1,i1,p1)
!                            case(2)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,2,p1), &
!                                    alphaMat(g1,3,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,2,p1) = tellers(g1,2,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,2,p1),2) = Wgroups(g1,i1,p1)
!                            case(3)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,3,p1), &
!                                    alphaMat(g1,4,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,3,p1) = tellers(g1,3,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,3,p1),3) = Wgroups(g1,i1,p1)
!                            case(4)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,4,p1), &
!                                    alphaMat(g1,5,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,4,p1) = tellers(g1,4,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,4,p1),4) = Wgroups(g1,i1,p1)
!                            case(5)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,5,p1), &
!                                    alphaMat(g1,6,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,5,p1) = tellers(g1,5,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,5,p1),5) = Wgroups(g1,i1,p1)
!                            case(6)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,6,p1), &
!                                    alphaMat(g1,7,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,6,p1) = tellers(g1,6,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,6,p1),6) = Wgroups(g1,i1,p1)
!                            case(7)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,7,p1), &
!                                    alphaMat(g1,8,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,7,p1) = tellers(g1,7,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,7,p1),7) = Wgroups(g1,i1,p1)
!                            case(8)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,8,p1), &
!                                    alphaMat(g1,9,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,8,p1) = tellers(g1,8,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,8,p1),8) = Wgroups(g1,i1,p1)
!                            case(9)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,9,p1), &
!                                    alphaMat(g1,10,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,9,p1) = tellers(g1,9,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,9,p1),9) = Wgroups(g1,i1,p1)
!                            case(10)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,10,p1), &
!                                    alphaMat(g1,11,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,10,p1) = tellers(g1,10,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,10,p1),10) = Wgroups(g1,i1,p1)
!                            case(11)
!                                call inverse_prob_sampling(condMean,condVar,1,0,alphaMat(g1,11,p1), &
!                                    alphaMat(g1,12,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,11,p1) = tellers(g1,11,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,11,p1),11) = Wgroups(g1,i1,p1)
!                         end select
!                    end do
!!
!                    !draw boundary's in alphaMat
!                    if(Cat(g1,p1)>2) then
!                        do c1=3,Cat(g1,p1)
!                            alphaMin = maxval(Wdummy(g1,p1,1:tellers(g1,c1-1,p1),c1-1))
!                            alphaMax = minval(Wdummy(g1,p1,1:tellers(g1,c1,p1),c1))
!                            alphaMat(g1,c1,p1) = d_rnunf()*(alphaMax-alphaMin)*.9999999+alphaMin+.0000001
!                        end do
!                    end if
!                end if
!!
!            end do
!!
!            !draw B
!            Bmean(1:K,1:P) = matmul(matmul(XtXi(g1,:,:),transpose(Xgroups(g1,1:Njs(g1),1:K))),Wgroups(g1,1:Njs(g1),1:P))
!            call kronecker(K,P,XtXi(g1,:,:),SigmaMatDraw,covBeta)
!            call chfac(covBeta,rankP,cholSigmaPK)
!            call rnmvn(cholSigmaPK,betaDrawj)
!            do p1=1,P
!                BDraws(g1,:,p1) = betaDrawj(1,((p1-1)*K+1):(p1*K)) + Bmean(1:K,p1)
!            end do
!!
!            !draw R using method of Liu and Daniels (LD, 2006)
!            !draw candidate R
!            diffmat(1:Njs(g1),1:P) = Wgroups(g1,1:Njs(g1),1:P) - matmul(Xgroups(g1,1:Njs(g1),1:K),BDraws(g1,1:K,1:P))
!            errorMatj = matmul(transpose(diffmat(1:Njs(g1),1:P)),diffmat(1:Njs(g1),1:P))
!            Ds = diag(1/sqrt(diagonals(errorMatj)))
!            diffmat(1:Njs(g1),1:P) = matmul(diffmat(1:Njs(g1),1:P),Ds) !diffmat is now epsilon in LD
!            epsteps = matmul(transpose(diffmat(1:Njs(g1),1:P)),diffmat(1:Njs(g1),1:P))
!            SS1 = matmul(matmul(diag(1/sigmaDraws(g1,:)),epsteps),diag(1/sigmaDraws(g1,:)))
!            SS1 = .i.SS1
!            !proposal send using a flat prior for the covariance matrix
!            call gen_wish(SS1,Njs(g1)-P-1,dummyPP)
!            dummyPP = .i.dummyPP
!	    	Ccan = matmul(matmul(diag(1/sqrt(diagonals(dummyPP))),dummyPP),diag(1/sqrt(diagonals(dummyPP))))
!            CcanInv = .i.Ccan
!            !always accept because the proposal prior is the same as the target prior
!            CDraws(g1,:,:) = Ccan(:,:)
!            Cinv = CcanInv
!!
!            !draw sigma's
!            do p1 = 1,P
!                if(ordinal(p1)==0) then
!                    bb = sum(errorMatj(p1,:)*Cinv(p1,:)/sigmaDraws(g1,:))-errorMatj(p1,p1)*Cinv(p1,p1)/sigmaDraws(g1,p1)
!                    aa = Cinv(p1,p1)*errorMatj(p1,p1)
!                    sigma_can(:) = sigmaDraws(g1,:)
!                    sigma_can(p1) = rnnof()
!                    sigma_can(p1) = sigma_can(p1)*sdMH(g1,p1) + sigmaDraws(g1,p1) !random walk
!                    R_MH = exp((-real(Njs(g1))-1.0)*(log(sigma_can(p1)-log(sigmaDraws(g1,p1)))) &
!                           -.5*aa*(sigma_can(p1)**(-2) - sigmaDraws(g1,p1)**(-2)) &
!                           -bb*(sigma_can(p1)**(-1) - sigmaDraws(g1,p1)**(-1)))
!                    if(rnunif(1) < R_MH .and. sigma_can(p1)>0) then
!                        sigmaDraws(g1,p1) = sigma_can(p1)
!                        acceptSigma(g1,p1) = acceptSigma(g1,p1) + 1
!                    end if
!                end if
!            end do
!
!            !Draw parameter extended parameter by Liu and Sabatti (2001) via random walk
!            SigmaInv = matmul(matmul(diag(1/sigmaDraws(g1,:)),Cinv),diag(1/sigmaDraws(g1,:)))
!            do p1 = 1,P
!                if(ordinal(p1)>0) then !draw gLiuSab(g1,p1)
!                    aa = errorMatj(p1,p1)*SigmaInv(p1,p1)/2
!                    bb = sum(errorMatj(p1,:)*SigmaInv(p1,:)*gLiuSab(g1,:)) - errorMatj(p1,p1)*SigmaInv(p1,p1)*gLiuSab(g1,p1)
!                    gLiuSab_can = rnnof()
!                    gLiuSab_can = gLiuSab_can*sdMHg(g1,p1) + gLiuSab(g1,p1) ! random (moon) walk
!                    R_MH = exp((K + Cat(g1,p1) - 2 + Njs(g1) - 1)*(log(gLiuSab_can) - log(gLiuSab(g1,p1))) &
!                                -aa*(gLiuSab_can**2 - gLiuSab(g1,p1)**2) - bb*(gLiuSab_can - gLiuSab(g1,p1)))
!                    if(rnunif(1) < R_MH .and. gLiuSab_can>0) then
!                        gLiuSab(g1,p1) = gLiuSab_can
!                        acceptLS(g1,p1) = acceptLS(g1,p1) + 1
!                        !update the other parameter through the parameter transformation g(x) = g * x
!                        BDraws(g1,1:K,p1) = BDraws(g1,1:K,p1)*gLiuSab(g1,p1)
!                        alphaMat(g1,3:Cat(g1,p1),p1) = alphaMat(g1,3:Cat(g1,p1),p1)*gLiuSab(g1,p1)
!                        Wgroups(g1,1:Njs(g1),p1) = Wgroups(g1,1:Njs(g1),p1)*gLiuSab(g1,p1)
!                    end if
!                end if
!            end do
!        end do
!!
!        write(60,*)Wgroups(1,6,3),Wgroups(1,1,3),Wgroups(1,3,3),Wgroups(1,5,3)
!        write(61,*)alphaMat(1,2,3),alphaMat(1,3,3),alphaMat(1,4,3)
!
!        if(modulo(s1,burnin/10)==0) then
!            write(*,*)100*s1/burnin,' percent done'
!        end if
!    end do
!!
!    write(*,*)
!    write(*,*)'sampling from stationary distribution'
!!
!    !start sampling under (stationary) posterior distribution
!    do s1 = 1,samsize0
!!
!        corrteller = 0
!        tellers = 0
!        do g1 = 1,numG
!            !compute means of latent W's for all observations
!            meanMat(1:Njs(g1),1:P) = matmul(Xgroups(g1,1:Njs(g1),1:K),BDraws(g1,1:K,1:P))
!            SigmaMatDraw = matmul(matmul(diag(sigmaDraws(g1,:)),CDraws(g1,:,:)),diag(sigmaDraws(g1,:)))
!!
!            !draw latent W's for the ordinal Y's
!            !compute mean vector for
!            do p1=1,P
!                if(ordinal(p1)>0) then
!                    do i1=1,Njs(g1)
!                        Yi1Categorie = Ygroups(g1,i1,p1)
!                        call compute_condMeanVar(p1,P,meanMat(i1,1:P),SigmaMatDraw,Wgroups(g1,i1,1:P),condMean,condVar)
!                        select case (Yi1Categorie)
!                            case(1)
!                                call inverse_prob_sampling(condMean,condVar,0,1,alphaMat(g1,1,p1), &
!                                    alphaMat(g1,2,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,1,p1) = tellers(g1,1,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,1,p1),1) = Wgroups(g1,i1,p1)
!                            case(2)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,2,p1), &
!                                    alphaMat(g1,3,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,2,p1) = tellers(g1,2,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,2,p1),2) = Wgroups(g1,i1,p1)
!                            case(3)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,3,p1), &
!                                    alphaMat(g1,4,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,3,p1) = tellers(g1,3,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,3,p1),3) = Wgroups(g1,i1,p1)
!                            case(4)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,4,p1), &
!                                    alphaMat(g1,5,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,4,p1) = tellers(g1,4,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,4,p1),4) = Wgroups(g1,i1,p1)
!                            case(5)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,5,p1), &
!                                    alphaMat(g1,6,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,5,p1) = tellers(g1,5,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,5,p1),5) = Wgroups(g1,i1,p1)
!                            case(6)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,6,p1), &
!                                    alphaMat(g1,7,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,6,p1) = tellers(g1,6,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,6,p1),6) = Wgroups(g1,i1,p1)
!                            case(7)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,7,p1), &
!                                    alphaMat(g1,8,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,7,p1) = tellers(g1,7,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,7,p1),7) = Wgroups(g1,i1,p1)
!                            case(8)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,8,p1), &
!                                    alphaMat(g1,9,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,8,p1) = tellers(g1,8,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,8,p1),8) = Wgroups(g1,i1,p1)
!                            case(9)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,9,p1), &
!                                    alphaMat(g1,10,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,9,p1) = tellers(g1,9,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,9,p1),9) = Wgroups(g1,i1,p1)
!                            case(10)
!                                call inverse_prob_sampling(condMean,condVar,1,1,alphaMat(g1,10,p1), &
!                                    alphaMat(g1,11,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,10,p1) = tellers(g1,10,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,10,p1),10) = Wgroups(g1,i1,p1)
!                            case(11)
!                                call inverse_prob_sampling(condMean,condVar,1,0,alphaMat(g1,11,p1), &
!                                    alphaMat(g1,12,p1),Wgroups(g1,i1,p1))
!                                tellers(g1,11,p1) = tellers(g1,11,p1) + 1
!                                Wdummy(g1,p1,tellers(g1,11,p1),11) = Wgroups(g1,i1,p1)
!                         end select
!                    end do
!!
!                    !draw boundary's in alphaMat
!                    if(Cat(g1,p1)>2) then
!                        do c1=3,Cat(g1,p1)
!                            alphaMin = maxval(Wdummy(g1,p1,1:tellers(g1,c1-1,p1),c1-1))
!                            alphaMax = minval(Wdummy(g1,p1,1:tellers(g1,c1,p1),c1))
!                            alphaMat(g1,c1,p1) = d_rnunf()*(alphaMax-alphaMin)*.9999998+alphaMin+.0000001
!                        end do
!                    end if
!                end if
!!
!            end do
!!
!            !draw B
!            Bmean(1:K,1:P) = XtXi(g1,:,:).xt.Xgroups(g1,1:Njs(g1),1:K).x.Wgroups(g1,1:Njs(g1),1:P)
!            call kronecker(K,P,XtXi(g1,:,:),SigmaMatDraw,covBeta)
!            call chfac(covBeta,rankP,cholSigmaPK)
!            call rnmvn(cholSigmaPK,betaDrawj)
!            do p1=1,P
!                BDraws(g1,:,p1) = betaDrawj(1,((p1-1)*K+1):(p1*K)) + Bmean(1:K,p1)
!            end do
!!
!            !draw R using method of Liu and Daniels (LD, 2006)
!            !draw candidate R
!            diffmat(1:Njs(g1),1:P) = Wgroups(g1,1:Njs(g1),1:P) - matmul(Xgroups(g1,1:Njs(g1),1:K),BDraws(g1,1:K,1:P))
!            errorMatj = matmul(transpose(diffmat(1:Njs(g1),1:P)),diffmat(1:Njs(g1),1:P))
!            Ds = diag(1/sqrt(diagonals(errorMatj)))
!            diffmat(1:Njs(g1),1:P) = matmul(diffmat(1:Njs(g1),1:P),Ds) !diffmat is now epsilon in LD
!            epsteps = matmul(transpose(diffmat(1:Njs(g1),1:P)),diffmat(1:Njs(g1),1:P))
!            SS1 = matmul(matmul(diag(1/sigmaDraws(g1,:)),epsteps),diag(1/sigmaDraws(g1,:)))
!            SS1 = .i.SS1
!            !proposal send using a flat prior for the covariance matrix
!            call gen_wish(SS1,Njs(g1)-P-1,dummyPP)
!            dummyPP = .i.dummyPP
!	    	Ccan = matmul(matmul(diag(1/sqrt(diagonals(dummyPP))),dummyPP),diag(1/sqrt(diagonals(dummyPP))))
!            CcanInv = .i.Ccan
!            !Metropolis-Hastings not needed because proposal prior is same as target prior
!            !logR_MH_part3 = (-.5*real(P+1))*(log(det(Ccurr))-log(det(Ccan))) !proposal prior
!            !R_MH = exp(logR_MH_part3)
!            !Thus, always accept because the proposal prior is the same as the target prior
!            CDraws(g1,:,:) = Ccan(:,:)
!            Cinv = CcanInv
!    !
!            do i1 = 2,P !store Fisher z transformed posterior draws of rho's
!                Zcorr_sample(s1,(corrteller+1):(corrteller+i1-1)) = .5*log((1+CDraws(g1,i1,1:(i1-1)))/(1-CDraws(g1,i1,1:(i1-1))))
!                corrteller = corrteller + i1 - 1
!            end do
!    !
!            !draw sigma's
!            do p1 = 1,P
!                if(ordinal(p1)==0) then
!                    bb = sum(errorMatj(p1,:)*Cinv(p1,:)/sigmaDraws(g1,:))-errorMatj(p1,p1)*Cinv(p1,p1)/sigmaDraws(g1,p1)
!                    aa = Cinv(p1,p1)*errorMatj(p1,p1)
!                    sigma_can(:) = sigmaDraws(g1,:)
!                    sigma_can(p1) = rnnof()
!                    sigma_can(p1) = sigma_can(p1)*sdMH(g1,p1) + sigmaDraws(g1,p1) !random walk
!                    R_MH = exp((-real(Njs(g1))-1.0)*(log(sigma_can(p1)-log(sigmaDraws(g1,p1)))) &
!                           -.5*aa*(sigma_can(p1)**(-2) - sigmaDraws(g1,p1)**(-2)) &
!                           -bb*(sigma_can(p1)**(-1) - sigmaDraws(g1,p1)**(-1)))
!                    if(rnunif(1) < R_MH .and. sigma_can(p1)>0) then
!                        sigmaDraws(g1,p1) = sigma_can(p1)
!                        acceptSigma(g1,p1) = acceptSigma(g1,p1) + 1
!                    end if
!                end if
!            end do
!
!            !Draw parameter extended parameter by Liu and Sabatti (2001) via random walk
!            SigmaInv = matmul(matmul(diag(1/sigmaDraws(g1,:)),Cinv),diag(1/sigmaDraws(g1,:)))
!            do p1 = 1,P
!                if(ordinal(p1)>0) then !draw gLiuSab(g1,p1)
!                    aa = errorMatj(p1,p1)*SigmaInv(p1,p1)/2
!                    bb = sum(errorMatj(p1,:)*SigmaInv(p1,:)*gLiuSab(g1,:)) - errorMatj(p1,p1)*SigmaInv(p1,p1)*gLiuSab(g1,p1)
!                    gLiuSab_can = rnnof()
!                    gLiuSab_can = gLiuSab_can*sdMHg(g1,p1) + gLiuSab(g1,p1) ! random (moon) walk
!                    R_MH = exp((K + Cat(g1,p1) - 2 + Njs(g1) - 1)*(log(gLiuSab_can) - log(gLiuSab(g1,p1))) &
!                                -aa*(gLiuSab_can**2 - gLiuSab(g1,p1)**2) - bb*(gLiuSab_can - gLiuSab(g1,p1)))
!                    if(rnunif(1) < R_MH .and. gLiuSab_can>0) then
!                        gLiuSab(g1,p1) = gLiuSab_can
!                        acceptLS(g1,p1) = acceptLS(g1,p1) + 1
!                        !update the other parameter through the parameter transformation g(x) = g * x
!                        BDraws(g1,1:K,p1) = BDraws(g1,1:K,p1)*gLiuSab(g1,p1)
!                        alphaMat(g1,3:Cat(g1,p1),p1) = alphaMat(g1,3:Cat(g1,p1),p1)*gLiuSab(g1,p1)
!                        Wgroups(g1,1:Njs(g1),p1) = Wgroups(g1,1:Njs(g1),p1)*gLiuSab(g1,p1)
!                    end if
!                end if
!            end do
!!
!        end do
!        BDrawsStore(s1,1:numG,1:K,1:P) = BDraws(1:numG,1:K,1:P)
!        sigmaDrawsStore(s1,1:numG,1:P) = sigmaDraws(1:numG,1:P)
!        CDrawsStore(s1,1:numG,1:P,1:P) = CDraws(1:numG,1:P,1:P)
!!
!        write(60,*)Wgroups(1,6,3),Wgroups(1,1,3),Wgroups(1,3,3),Wgroups(1,5,3)
!        write(61,*)alphaMat(1,2,3),alphaMat(1,3,3),alphaMat(1,4,3)
!        write(62,*)Zcorr_sample(s1,1:3)
!!
!        if(modulo(s1,samsize0/10)==0) then
!            write(*,*)100*s1/samsize0,' percent done'
!        end if
!    end do
!!
!    !compute medians Fisher transformed draws
!    do c1=1,numcorr
!        dummy2(:) = Zcorr_sample(:,c1)
!        call svrgp(dummy2,dummy3,iperm1)
!        postZmean(c1,1) = dummy3(int(samsize0*.5))
!    end do
!!
!    !robust covariance matrix estimation Fisher transformed draws
!    postZcov = 0
!    do c1=1,numcorr
!        do c2=c1,numcorr
!            call robust_covest(samsize0, Zcorr_sample(1:samsize0,c1), Zcorr_sample(1:samsize0,c2), &
!                postZmean(c1,1), postZmean(c2,1), varz1, varz2, varz1z2Plus, varz1z2Min)
!            postZcov(c1,c2) = (varz1*varz2)**.5 * (varz1z2Plus - varz1z2Min)/(varz1z2Plus + varz1z2Min)
!            postZcov(c2,c1) = postZcov(c1,c2)
!        end do
!    end do
!!
!    !compute medians, lower bounds CI, upperbounds CI
!    do g1=1,numG
!        C_quantiles(g1,:,:,1) = eye(P)
!        C_quantiles(g1,:,:,2) = eye(P)
!        C_quantiles(g1,:,:,3) = eye(P)
!        do p1=1,P
!            !for the sigma's
!            dummy2(:) = sigmaDrawsStore(:,g1,p1)
!            call svrgp(dummy2,dummy3,iperm1)
!            sigma_quantiles(g1,p1,1) = dummy3(int(samsize0*.025))
!            sigma_quantiles(g1,p1,2) = dummy3(int(samsize0*.5))
!            sigma_quantiles(g1,p1,3) = dummy3(int(samsize0*.975))
!            !for the beta coefficients
!            do k1=1,K
!                dummy2(:) = BDrawsStore(:,g1,k1,p1)
!                call svrgp(dummy2,dummy3,iperm1)
!                B_quantiles(g1,k1,p1,1) = dummy3(int(samsize0*.025))
!                B_quantiles(g1,k1,p1,2) = dummy3(int(samsize0*.5))
!                B_quantiles(g1,k1,p1,3) = dummy3(int(samsize0*.975))
!            end do
!            if(p1>1)then
!                do p2=1,p1-1
!                    dummy2(:) = CDrawsStore(:,g1,p1,p2)
!                    call svrgp(dummy2,dummy3,iperm1)
!                    C_quantiles(g1,p1,p2,1) = dummy3(int(samsize0*.025))
!                    C_quantiles(g1,p1,p2,2) = dummy3(int(samsize0*.5))
!                    C_quantiles(g1,p1,p2,3) = dummy3(int(samsize0*.975))
!                end do
!            end if
!        end do
!    end do
!!
!    !plot results
!        write(70,*)
!    write(70,*)
!    write(70,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(70,*)'!                                                                                   !'
!    write(70,*)'!                            OUTPUT BCT_estimates                                   !'
!    write(70,*)'!                                                                                   !'
!    write(70,*)'!       When using this software please refer to Mulder & Gelissen. Bayesian        !'
!    write(70,*)'!       testing of equality and order constraints on measures of association        !'
!    write(70,*)'!       in social research (under review).                                          !'
!    write(70,*)'!                                                                                   !'
!    write(70,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(70,*)
!    write(70,*)
!    write(70,*)'Estimates were obtained under the unconstrained model'
!    write(70,*)
!    write(70,*)'Correlation matrix'
!    write(70,*)
!    do g1 = 1,numG
!        write(70,*)'Group',g1
!        write(70,*)
!        write(70,*)'lower bound of 95%-CI'
!        do p1=1,P
!            write(70,800)C_quantiles(g1,p1,1:p1,1)
!        end do
!        write(70,*)
!        write(70,*)'median'
!        do p1=1,P
!            write(70,800)C_quantiles(g1,p1,1:p1,2)
!        end do
!        write(70,*)
!        write(70,*)'upper bound of 95%-CI'
!        do p1=1,P
!            write(70,800)C_quantiles(g1,p1,1:p1,3)
!        end do
!        write(70,*)
!        write(70,*)
!    end do
!    write(70,*)'B-matrix with intercepts and regression coefficients'
!    write(70,*)
!    do g1 = 1,numG
!        write(70,*)'Group',g1
!        write(70,*)
!        write(70,*)'lower bound of 95%-CI'
!        do k1=1,K
!            write(70,800)B_quantiles(g1,k1,1:P,1)
!        end do
!        write(70,*)
!        write(70,*)'median'
!        do k1=1,K
!            write(70,800)B_quantiles(g1,k1,1:P,2)
!        end do
!        write(70,*)
!        write(70,*)'upper bound of 95%-CI'
!        do k1=1,K
!            write(70,800)B_quantiles(g1,k1,1:P,3)
!        end do
!        write(70,*)
!        write(70,*)
!    end do
!    write(70,*)'standard deviations'
!    write(70,*)
!    do g1 = 1,numG
!        write(70,*)'Group',g1
!        write(70,*)
!        write(70,*)'lower bound of 95%-CI'
!        write(70,800)sigma_quantiles(g1,1:P,1)
!        write(70,*)
!        write(70,*)'median'
!        write(70,800)sigma_quantiles(g1,1:P,2)
!        write(70,*)
!        write(70,*)'upper bound of 95%-CI'
!        write(70,800)sigma_quantiles(g1,1:P,3)
!        write(70,*)
!        write(70,*)
!    end do
!!
!!
!end subroutine estimate_postMeanCov_FisherZordinal
!
!
!!
!!   Calculate robust estimates of variances of beta1, beta2, beta1+beta2, and beta1-beta2
!!
!    subroutine robust_covest(m, betas1, betas2, mn1, mn2, varb1, varb2, varb1b2Plus, varb1b2Min)
!!
!    use IMSL_LIBRARIES
!    implicit none
!!
!    !Declare local variables
!    integer, intent(in) :: m
!    real, intent(in)    :: betas1(m), betas2(m), mn1, mn2
!    real, intent(out)   :: varb1, varb2, varb1b2Plus, varb1b2Min
!
!    real                :: dummy1(m), dummy2(m), Phi075
!    integer             :: iperm(m), mmin, i
!!
!    Phi075 = anorin(.75)
!!
!    !robust variance estimators of beta1 and beta2
!    call svrgp(abs(betas1 - mn1), dummy1, iperm)
!    do i=1,m
!        if(dummy1(i)>0) then
!            mmin = i
!            exit
!        end if
!    end do
!    varb1 = ((dummy1(int(mmin+(m-mmin)*.5)) + dummy1(int(mmin+(m-mmin)*.5+1)))*.5/Phi075)**2.0
!    call svrgp(abs(betas2 - mn2), dummy1, iperm)
!    do i=1,m
!        if(dummy1(i)>0) then
!            mmin = i
!            exit
!        end if
!    end do
!    varb2 = ((dummy1(int(mmin+(m-mmin)*.5)) + dummy1(int(mmin+(m-mmin)*.5+1)))*.5/Phi075)**2.0
!!
!    !robust variance estimators of beta1 + beta2
!    dummy2 = betas1 + betas2
!    call svrgp(abs(dummy2 - mn1 - mn2), dummy1, iperm)
!    do i=1,m
!        if(dummy1(i)>0) then
!            mmin = i
!            exit
!        end if
!    end do
!    varb1b2Plus = ((dummy1(int(mmin+(m-mmin)*.5)) + dummy1(int(mmin+(m-mmin)*.5+1)))*.5/Phi075)**2.0
!!
!    !robust variance estimators of beta1 - beta2
!    dummy2 = betas1 - betas2
!    call svrgp(abs(dummy2 - mn1 + mn2), dummy1, iperm)
!    do i=1,m
!        if(dummy1(i)>0) then
!            mmin = i
!            exit
!        end if
!    end do
!    varb1b2Min = ((dummy1(int(mmin+(m-mmin)*.5)) + dummy1(int(mmin+(m-mmin)*.5+1)))*.5/Phi075)**2.0
!!
!    end subroutine robust_covest
!!
!
!
!subroutine kronecker(dimA,dimB,A,B,AB)
!!
!    use global_variables
!    implicit none
!!
!    integer, intent(in) :: dimA, dimB
!    real, intent(in)    :: A(dimA,dimA), B(dimB,dimB) !dummy arguments
!    real, intent(out)   :: AB(dimA*dimB,dimA*dimB) !output matrix of the kronecker product
!    integer             :: i,j !loop counters
!!
!    do i=1,dimA
!        do j=1,dimA
!            AB((1+dimB*(i-1)):(dimB+dimB*(i-1)),(1+dimB*(j-1)):(dimB+dimB*(j-1))) = A(i,j)*B(:,:)
!        end do
!    end do
!!
!end subroutine kronecker
!
!
!
!!draw randomly from Wishart distribution with scale matrix A and nu degrees of freedom
!!
!subroutine gen_wish(A,nu,B)
!!
!    use IMSL_LIBRARIES
!    use global_variables
!    implicit none
!!
!    !Declare local variables
!    real, intent (in)    :: A(P,P)
!    integer, intent (in) :: nu
!    real, intent (out)   :: B(P,P)
!    integer              :: irank
!    real                 :: chA(P,P), RNmat(nu,P)
!!
!    !sample from Wishart distribution as in Press (2005, p. 109)
!    call chfac(A, irank, chA)
!    call rnmvn(chA, RNmat)
!    B = RNmat.tx.RNmat
!!
!end subroutine gen_wish
!
!
!
!subroutine Mrank(numIR, numSP, rowrank, IRr, transR, constant, transcon)
!!
!   use IMSL_LIBRARIES
!   implicit none
!!
!   integer, intent(in)                                     :: numIR, numSP
!   integer, intent(out)                                    :: rowrank
!   Integer                                                 :: i, j, k
!   real                                                    :: temp
!   integer, allocatable, dimension(:)                      :: rownumber
!   real, intent(inout)                                     :: constant(numIR), transcon(numIR)
!   real, allocatable, dimension(:,:)                       :: A, Ainv
!   real, intent(inout)                                     :: IRr(numIR,numSP+1), transR(numIR,numIR)
!!
!   if (numIR>numSP) then
!      allocate (A(numIR,numIR), Ainv(numIR,numIR), rownumber(numIR))
!   else
!      allocate (A(numIR,numSP), Ainv(numIR,numSP), rownumber(numIR))
!   end if
!!
!   A=0
!   Ainv=0
!   transR=0
!   rownumber=1
!   do i = 1,numIR
!      A(i,1:numSP)=IRr(i,1:numSP)
!      Ainv(i,i)=1
!      if (i>1) then
!         rownumber(i)=rownumber(i-1)+1
!      end if
!   end do
!   rowrank=numIR
!   do i=1,numIR
!      do k=i+1,numIR
!         if (A(i,i).EQ.0 .and. A(k,i).NE.0) then
!            do j=1,max(numIR,numSP)
!               temp=A(i,j)
!               A(i,j)=A(k,j)
!               A(k,j)=temp
!               temp=Ainv(i,j)
!               Ainv(i,j)=Ainv(k,j)
!               Ainv(k,j)=temp
!            end do
!            do j=1,(numSP+1)
!               temp=IRr(i,j)
!               IRr(i,j)=IRr(k,j)
!               IRr(k,j)=temp
!            end do
!               temp=rownumber(i)
!               rownumber(i)=rownumber(k)
!               rownumber(k)=temp
!         end if
!      end do
!      if (A(i,i).NE.0) then
!         temp=A(i,i)
!         do j=1,max(numIR,numSP)
!            A(i,j)=A(i,j)/temp
!            Ainv(i,j)=Ainv(i,j)/temp
!         end do
!         do k=1,numIR
!            if (k.NE.i) then
!               temp=A(k,i)
!               do j=1,max(numIR,numSP)
!                  A(k,j)=A(k,j)-A(i,j)*temp
!                  Ainv(k,j)=Ainv(k,j)-Ainv(i,j)*temp
!               end do
!            end if
!         end do
!      end if
!   end do
!   do i=1,numIR
!      if (sum(abs(A(i,1:numSP))).EQ.0) then
!         rowrank=rowrank-1
!      end if
!   end do
!   do i=1,numIR
!      do k=i+1,numIR
!         if (sum(abs(A(i,1:numSP))).EQ.0 .and. sum(abs(A(k,1:numSP))).NE.0) then
!            do j=1,max(numIR,numSP)
!               temp=A(i,j)
!               A(i,j)=A(k,j)
!               A(k,j)=temp
!               temp=Ainv(i,j)
!               Ainv(i,j)=Ainv(k,j)
!               Ainv(k,j)=temp
!            end do
!            do j=1,(numSP+1)
!               temp=IRr(i,j)
!               IRr(i,j)=IRr(k,j)
!               IRr(k,j)=temp
!            end do
!               temp=rownumber(i)
!               rownumber(i)=rownumber(k)
!               rownumber(k)=temp
!         end if
!      end do
!   end do
!
!   do i=1,numIR
!      do j=1,numIR
!         transR(i,j)=Ainv(i,rownumber(j))
!      end do
!
!      transcon(i)=sum(Ainv(i,1:numIR)*constant(1:numIR))
!   end do
!!
!deallocate(A,Ainv,rownumber)
!end subroutine Mrank
!
!
!
!subroutine compute_HPDinterval(samsize,sample,perc,LB,UB)
!
!    use IMSL_LIBRARIES
!    implicit none
!
!    integer, intent(in) :: samsize
!    real, intent(in)    :: sample(samsize), perc
!    real, intent(out)   :: LB, UB
!!
!    real                :: sampleSort(samsize), length(samsize)
!    integer             :: iperm(samsize), welkeLB, welkeUB, i1
!!
!    iperm = 1
!    call svrgp(sample,sampleSort,iperm)
!!
!    length = 99999
!    welkeLB = 1
!    welkeUB = samsize*perc/100.0
!!
!    do i1=1,samsize - welkeUB
!        length(i1) = sampleSort(welkeUB) - sampleSort(welkeLB)
!        welkeLB = welkeLB + 1
!        welkeUB = welkeUB + 1
!    end do
!
!end subroutine compute_HPDinterval
!
!
!
!subroutine compute_condMeanVar(welke,dimIn,meanIn,covmIn,obsIn,condMean,condVar)
!
!    use IMSL_LIBRARIES
!    implicit none
!
!    integer, intent(in) :: welke, dimIn
!    real, intent(in)    :: meanIn(dimIn), covmIn(dimIn,dimIn), obsIn(dimIn)
!    real, intent(out)   :: condMean, condVar
!    real                :: dummy3(1,1), dummy2(dimIn-1,1), S12(1,dimIn-1), S22(dimIn-1,dimIn-1), &
!                           S22inv(dimIn-1,dimIn-1), ZS(dimIn,1), meanLocal(dimIn,1)
!!
!    meanLocal(1:dimIn,1) = meanIn(1:dimIn)
!!
!    S12(1,1:(welke-1)) = covmIn(welke,1:(welke-1))
!    S12(1,welke:(dimIn-1)) = covmIn(welke,(welke+1):dimIn)
!    S22(1:(welke-1),1:(welke-1)) = covmIn(1:(welke-1),1:(welke-1))
!    S22(1:(welke-1),welke:(dimIn-1)) = covmIn(1:(welke-1),(welke+1):dimIn)
!    S22(welke:(dimIn-1),1:(welke-1)) = covmIn((welke+1):dimIn,1:(welke-1))
!    S22(welke:(dimIn-1),welke:(dimIn-1)) = covmIn((welke+1):dimIn,(welke+1):dimIn)
!    S22inv = .i.S22
!    dummy2(1:(welke-1),1) = obsIn(1:(welke-1)) - meanLocal(1:(welke-1),1)
!    dummy2(welke:(dimIn-1),1) = obsIn((welke+1):dimIn) - meanLocal((welke+1):dimIn,1)
!    dummy3 = matmul(matmul(S12,S22inv),dummy2)
!    condMean = meanLocal(welke,1) + dummy3(1,1) !conditional mean
!    dummy3 = matmul(matmul(S12,S22inv),transpose(S12))
!    condVar = covmIn(welke,welke) - dummy3(1,1) !conditional variance
!!
!end subroutine compute_condMeanVar
!
!
!
!
!
!
!
!
!
!
!
!
