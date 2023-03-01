MODULE VFI
  USE datetime_module
  USE MPIMod
  USE Param
  USE NL
  USE pdefMod
  IMPLICIT NONE
  
CONTAINS

  !
  !
  ! Stationary
  !
  !
  SUBROUTINE compute()
    REAL(wp) :: errPols, errFuture
    INTEGER :: iter

    iter = 1
    errPols = 1.0_wp
    errFuture = 1.0_wp
    DO WHILE (iter <= maxIter .AND. MAX(errFuture, errPols) > epsTol)
      CALL iteratePolicies(errPols)
      CALL iterateFuture(errFuture)

      IF (workerId == 0 .AND. MOD(iter, 2000) == 0) THEN
        WRITE (OUTPUT_UNIT, "(A,I5,2(A,E10.3))") "Iteration ", iter, &
          " errPols = ", errPols, " errFuture = ", errFuture
        FLUSH (OUTPUT_UNIT)
      END IF
      iter = iter + 1
    END DO
  END SUBROUTINE compute

  SUBROUTINE iteratePolicies(errVal)
    REAL(wp), INTENT(OUT) :: errVal
    INTEGER :: bIx, bPrIx
    REAL(wp) :: theMax, theSum, dStarHere, CdZero, CdInterior
    REAL(wp), DIMENSION(bSz) :: theExps
    REAL(wp), DIMENSION(bSz, bSz) :: uVal

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(bIx,bPrIx,dStarHere,CdZero,CdInterior)
    DO bIx = 1,bSz
    DO bPrIx = 1,bSz
      dStarHere = ((1.0_wp - kappa * q0(bPrIx)) * pay * bGrid(bIx) ) / &
        (meanZ * (1.0_wp - gamm2) * gamm0 * gamm1)

      IF (dStarHere >= 1.0_wp) THEN
        dStar(bIx, bPrIx) = dStarHere
        dStarHere = 1.0_wp
      ELSE IF (dStarHere <= 0.0_wp) THEN
        dStar(bIx, bPrIx) = dStarHere
        dStarHere = 0.0_wp
      ELSE
        dStarHere = dStarHere**(1.0_wp / (gamm1 - 1.0_wp))
        dStar(bIx, bPrIx) = dStarHere
      END IF

      CdZero = meanZ - pay * bGrid(bIx) &
        + q0(bPrIx) * (bGrid(bPrIx) - bGrid(bIx) * (1.0_wp - delta))
      CdInterior = meanZ * phi(dStarHere) - pay * (1.0_wp - dStarHere) * bGrid(bIx) &
        + q0(bPrIx) * (bGrid(bPrIx) - bGrid(bIx) * (1.0_wp - delta + kappa * pay * dStarHere))

      IF (withService) THEN
        CdZero = CdZero - exoLoanPay 
        CdInterior = CdInterior - exoLoanPay 
      END IF

      IF (CdInterior > CdZero) THEN
        d(bIx, bPrIx) = dStarHere
        C(bIx, bPrIx) = CdInterior
      ELSE
        d(bIx, bPrIx) = 0.0_wp
        C(bIx, bPrIx) = CdZero
      END IF
    END DO
    END DO

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(bIx, bPrIx)
    DO bIx = 1,bSz
    DO bPrIx = 1,bSz
      IF (C(bIx, bPrIx) <= 0_wp) THEN
        uVal(bIx, bPrIx) = veryNeg
      ELSE
        uVal(bIx, bPrIx) = uu(C(bIx, bPrIx)) + beta * V0(bPrIx)
      END IF
    END DO
    END DO

    !$OMP PARALLEL DO PRIVATE(bIx,theExps,theMax,theSum,bPrIx)
    DO bIx = 1,bSz
      theMax = MAXVAL( uVal(bIx, :) ) ! , uVal(bIx, :) > veryNeg )
      theExps = EXP( (uVal(bIx, :) - theMax) / tasteB )
      theSum = SUM( theExps )

      IF (theSum > 0_wp) THEN
        polPr(bIx, :) = theExps / theSum
        V1(bIx) = theMax + tasteB * (LOG(theSum) - eulerMascheroni)
      ELSE
        polPr(bIx, :) = 0_wp
        polPr(bIx, bIx) = 1.0_wp
        V1(bIx) = veryNeg
      END IF
    END DO

    errVal = MAXVAL(ABS(V1-V0))
    V0 = V1
  END SUBROUTINE iteratePolicies

  SUBROUTINE iterateFuture(errVal)
    REAL(wp), INTENT(OUT) :: errVal
    INTEGER :: bIx, bPrIx
    REAL(wp) :: qTmp, dPrVal 

    !$OMP PARALLEL DO PRIVATE(bIx,bPrIx,qTmp,dPrVal) 
    DO bIx = 1,bSz
      qTmp = 0_wp
      DO bPrIx = 1,bSz
        dPrVal = d(bIx, bPrIx)
        qTmp = qTmp + polPr(bIx, bPrIx) * ( pay * (1.0_wp - dPrVal) &
          + (1.0_wp - delta + kappa * pay * dPrVal) * q0(bPrIx) )
      END DO
      q1(bIx) = qTmp / rf
    END DO

    errVal = MAXVAL(ABS(q1 - q0))
    q0 = q1
  END SUBROUTINE iterateFuture

  SUBROUTINE expectedPolicies()
    INTEGER :: bIx

    !$OMP PARALLEL DO PRIVATE(bIx)
    DO bIx = 1,bSz
      EC(bIx) = DOT_PRODUCT( polPr(bIx, :), C(bIx, :) )
      Ed(bIx) = DOT_PRODUCT( polPr(bIx, :), d(bIx, :) )
      EbPr(bIx) = DOT_PRODUCT( polPr(bIx, :), bGrid ) 

      ECsrv(bIx) = DOT_PRODUCT( polPrSrv(bIx, :), Csrv(bIx, :) )
      EdSrv(bIx) = DOT_PRODUCT( polPrSrv(bIx, :), dSrv(bIx, :) )
      EbPrSrv(bIx) = DOT_PRODUCT( polPrSrv(bIx, :), bGrid ) 
    END DO
  END SUBROUTINE expectedPolicies

  !
  !
  ! Backward
  !
  !
  SUBROUTINE setupTerminal()
    INTEGER :: sIx, iIx

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(iIx,sIx)
    DO iIx = 1,muIsz
    DO sIx = 1,muSsz
      sirV1(sIx, iIx, :) = V0srv
      sirq1(sIx, iIx, :) = q0srv
      sirEbPr(sIx, iIx, :) = EbPrSrv
      sird(sIx, iIx, :, :) = dSrv
      sirPolPr(sIx, iIx, :, :) = polPrSrv
    END DO
    END DO

    sirEL = 0_wp
    sirLix = 1_sip

    sirV0 = sirV1
    sirq0 = sirq1
  END SUBROUTINE setupTerminal

  SUBROUTINE goBackwards(hIx)
    INTEGER, INTENT(IN) :: hIx
    INTEGER :: tIx
    INTEGER :: ixx, sIx, iIx, bIx, bPrIx
    INTEGER(sip) :: Lix
    REAL(wp) :: cSS, cII, deaths, nSS, nII
    REAL(wp) :: qTmp, dPrVal, dStarHere, CdZero, CdInterior
    REAL(wp) :: theSum, theMax
    REAL(wp), DIMENSION(Lsz) :: Cvals, dvals, vvals
    REAL(wp), DIMENSION(bSz) :: theExps, ww, cctmp

    tIx = hIx + before

    IF (workerId == 0)  THEN 
      WRITE (OUTPUT_UNIT, "(A,I6,A,F8.4)") &
        TAB // "Backwards ", hIx, " with wedge ", histWedge(tIx)
      FLUSH (OUTPUT_UNIT)
    END IF

    sirV0 = sirV1
    sirq0 = sirq1

    ! Bond Prices
    !$OMP PARALLEL DO PRIVATE(ixx,sIx,iIx,&
    !$OMP& cSS,cII,bIx,qTmp,bPrIx,Lix,dPrVal,deaths,nSS,nII)
    DO ixx = 1,chunk
      sIx = longToWide(workerId * chunk + ixx, 1) 
      iIx = longToWide(workerId * chunk + ixx, 2) 
      bIx = longToWide(workerId * chunk + ixx, 3) 

      cSS = muSgrid(sIx)
      cII = muIgrid(iIx)

      qTmp = 0.0_wp
      DO bPrIx = 1,bSz
        IF (sirPolPr(sIx, iIx, bIx, bPrIx) > 1.0D-12) THEN
          Lix = sirLix(sIx, iIx, bIx, bPrIx)
          dPrVal = sird(sIx, iIx, bIx, bPrIx)
          CALL sirLOM(histWedge(tIx + 1) * thetaLsq(Lix), cSS, cII, deaths, nSS, nII) 
          qTmp = qTmp + sirPolPr(sIx, iIx, bIx, bPrIx) * ( &
            pay * (1.0_wp - dPrVal) + (1.0_wp - delta + kappa * pay * dPrVal) &
            * iQi(nSS, nII, bPrIx) )
        END IF
      END DO
      buffq(ixx) = qTmp / rf
    END DO

    CALL MPI_ALLGATHER(buffq, chunk, MPI_DOUBLE_PRECISION, buffBIG, chunk, &
      MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
    CALL mpiBarrier()

    !$OMP PARALLEL DO PRIVATE(sIx,iIx,bIx)
    DO ixx = 1,spaceSz 
      sIx = longToWide(ixx, 1) 
      iIx = longToWide(ixx, 2) 
      bIx = longToWide(ixx, 3) 

      sirq1(sIx, iIx, bIx) = buffBIG(ixx)
    END DO
 
    ! Values and policies
    !$OMP PARALLEL DO PRIVATE(ixx, &
    !$OMP& sIx,iIx,bIx,cSS,cII,bPrIx,Lix,nSS,nII,qTmp,&
    !$OMP& dStarHere,CdZero,CdInterior,Cvals,dvals,vvals,ww,&
    !$OMP& cctmp,theMax,theSum,theExps,deaths)
    DO ixx = 1,chunk
      sIx = longToWide(workerId * chunk + ixx, 1) 
      iIx = longToWide(workerId * chunk + ixx, 2) 
      bIx = longToWide(workerId * chunk + ixx, 3) 

      cSS = muSgrid(sIx)
      cII = muIgrid(iIx)

      DO bPrIx = 1,bSz
      DO Lix = 1_sip,Lsz
        CALL sirLOM(histWedge(tIx) * thetaLsq(Lix), cSS, cII, deaths, nSS, nII) 
        qTmp = iQi(nSS, nII, bPrIx, .TRUE.)

        dStarHere = ((1.0_wp - kappa * qTmp) * pay * bGrid(bIx) ) / &
          (histZ(hIx) * (1.0_wp - thetaY * L(Lix))**drs * (1.0_wp - gamm2) * gamm0 * gamm1)

        IF (dStarHere >= 1.0_wp) THEN
          dStarHere = 1.0_wp
        ELSE IF (dStarHere <= 0.0_wp) THEN
          dStarHere = 0.0_wp
        ELSE
          dStarHere = dStarHere**(1.0_wp / (gamm1 - 1.0_wp))
        END IF

        CdZero = histZ(hIx) * (1.0_wp - thetaY * L(Lix))**drs + exoCash(hIx) &
          - pay * bGrid(bIx) + qTmp * (bGrid(bPrIx) - bGrid(bIx) * (1.0_wp - delta))
          
        CdInterior = histZ(hIx) * (1.0_wp - thetaY * L(Lix))**drs * phi(dStarHere) + exoCash(hIx) & 
          - pay * (1.0_wp - dStarHere) * bGrid(bIx) + qTmp &
          * (bGrid(bPrIx) - bGrid(bIx) * (1.0_wp - delta + kappa * pay * dStarHere))

        IF (CdInterior > CdZero) THEN
          dvals(Lix) = dStarHere
          Cvals(Lix) = CdInterior
        ELSE
          dvals(Lix) = 0.0_wp
          Cvals(Lix) = CdZero
        END IF

        IF (Cvals(Lix) > 0.0_wp) THEN
          vvals(Lix) = uu(Cvals(Lix)) - chi * deaths + beta * iVi(nSS, nII, bPrIx)
        ELSE
          vvals(Lix) = veryNeg
        END IF
      END DO ! Lix

      ! Best L for b'
      ww(bPrIx) = MAXVAL(vvals)
      Lix = MAXLOC(vvals, 1)
      cctmp(bPrIx) = Cvals(Lix) 
      buffL((ixx-1) * bSz + bPrIx) = Lix
      buffd((ixx-1) * bSz + bPrIx) = dvals(Lix)
      END DO ! bPrIx

      ! Best b'
      theMax = MAXVAL(ww)
      theExps = EXP((ww - theMax) / tasteB)
      theSum = SUM(theExps)
      IF (theSum > 0.0_wp) THEN
        theExps = theExps / theSum
        buffb((ixx-1)*bSz+1:ixx*bSz) = theExps
        buffv(ixx) = theMax + tasteB * (LOG(theSum) - eulerMascheroni)
        buffc(ixx) = DOT_PRODUCT( theExps, cctmp )
      ELSE
        buffb((ixx-1)*bSz+1:ixx*bSz) = 0.0_wp
        buffv(ixx) = veryNeg
        buffc(ixx) = veryNeg
      END IF
    END DO ! ixx

    CALL MPI_ALLGATHER(buffv, chunk, MPI_DOUBLE_PRECISION, buffBIG, chunk, &
      MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
    CALL mpiBarrier()
    !$OMP PARALLEL DO PRIVATE(sIx,iIx,bIx)
    DO ixx = 1,spaceSz 
      sIx = longToWide(ixx, 1) 
      iIx = longToWide(ixx, 2) 
      bIx = longToWide(ixx, 3) 

      sirV1(sIx, iIx, bIx) = buffBIG(ixx)
    END DO

    CALL MPI_ALLGATHER(buffc, chunk, MPI_DOUBLE_PRECISION, buffBIG, chunk, &
      MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
    CALL mpiBarrier()
    !$OMP PARALLEL DO PRIVATE(sIx,iIx,bIx)
    DO ixx = 1,spaceSz 
      sIx = longToWide(ixx, 1) 
      iIx = longToWide(ixx, 2) 
      bIx = longToWide(ixx, 3) 

      sirEC(sIx, iIx, bIx) = buffBIG(ixx)
    END DO

    CALL MPI_ALLGATHER(buffb, chunk*bSz, MPI_DOUBLE_PRECISION, buffHUGE, chunk*bSz, &
      MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
    CALL mpiBarrier()
    !$OMP PARALLEL DO PRIVATE(sIx,iIx,bIx)
    DO ixx = 1,spaceSz 
      sIx = longToWide(ixx, 1) 
      iIx = longToWide(ixx, 2) 
      bIx = longToWide(ixx, 3) 

      sirPolPr(sIx, iIx, bIx, :) = buffHUGE((ixx-1)*bSz+1:ixx*bSz)
    END DO

    CALL MPI_ALLGATHER(buffd, chunk*bSz, MPI_DOUBLE_PRECISION, buffHUGE, chunk*bSz, &
      MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
    CALL mpiBarrier()
    !$OMP PARALLEL DO PRIVATE(sIx,iIx,bIx)
    DO ixx = 1,spaceSz 
      sIx = longToWide(ixx, 1) 
      iIx = longToWide(ixx, 2) 
      bIx = longToWide(ixx, 3) 

      sird(sIx, iIx, bIx, :) = buffHUGE((ixx-1)*bSz+1:ixx*bSz)
    END DO

    CALL MPI_ALLGATHER(buffl, chunk*bSz, MPI_INTEGER1, buffHUGEi, chunk*bSz, &
      MPI_INTEGER1, MPI_COMM_WORLD, mpiErr)
    CALL mpiBarrier()
    !$OMP PARALLEL DO PRIVATE(sIx,iIx,bIx)
    DO ixx = 1,spaceSz 
      sIx = longToWide(ixx, 1) 
      iIx = longToWide(ixx, 2) 
      bIx = longToWide(ixx, 3) 

      sirLix(sIx, iIx, bIx, :) = buffHUGEi((ixx-1)*bSz+1:ixx*bSz)
    END DO

    ! Expected Policies
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(sIx,iIx,bIx,theExps)
    DO bIx = 1,bSz
    DO sIx = 1,muSsz
    DO iIx = 1,muIsz
      theExps = sirPolPr(sIx, iIx, bIx, :)
      sirEbPr(sIx, iIx, bIx) = DOT_PRODUCT( theExps, bGrid )
      sirEL(sIx, iIx, bIx) = DOT_PRODUCT( &
        theExps, L(sirLix(sIx, iIx, bIx, :)) )
      sirEd(sIx, iIx, bIx) = DOT_PRODUCT( &
        theExps, sird(sIx, iIx, bIx, :) )
    END DO ! iIx
    END DO ! sIx
    END DO ! bIx

    IF (workerId == 0) THEN
      CALL saveSirPeriod(hIx)
      CALL tstamp()
    END IF
  END SUBROUTINE goBackwards

  !
  !
  ! Forward
  !
  !
  SUBROUTINE goForwardsStat(tIx)
    INTEGER, INTENT(IN) :: tIx
    REAL(wp) :: bNow, bPrime, theQ, dStarHere, CdZero, CdInterior

    bNow = histB(tIx)

    IF (tIx >= before + startRepay) THEN
      histV(tIx) = linInterp(bGrid, V0srv, bNow)
      bPrime = linInterp(bGrid, EbPrSrv, bNow)
      bPrime = MAX(minB, MIN(maxB, bPrime))
      histB(tIx+1) = bPrime
      theQ = linInterp(bGrid, q0srv, bPrime)
    ELSE
      histV(tIx) = linInterp(bGrid, V0, bNow)
      bPrime = linInterp(bGrid, EbPr, bNow)
      bPrime = MAX(minB, MIN(maxB, bPrime))
      histB(tIx+1) = bPrime
      theQ = linInterp(bGrid, q0, bPrime)
    END IF

    histSp(tIx) = (1.0_wp + pay * (1.0_wp / theQ - 1.0_wp))**52_wp - 1.0_wp

    dStarHere = ((1.0_wp - kappa * theQ) * pay * bNow ) / &
      (meanZ * (1.0_wp - gamm2) * gamm0 * gamm1)

    IF (dStarHere >= 1.0_wp) THEN
      dStarHere = 1.0_wp
    ELSE IF (dStarHere <= 0_wp) THEN
      dStarHere = 0.0_wp
    ELSE
      dStarHere = dStarHere**(1.0_wp / (gamm1 - 1.0_wp))
    END IF

    CdZero = meanZ - pay * bNow + theQ * (bPrime - bNow * (1.0_wp - delta)) + histCash(tIx)
    CdInterior = meanZ * phi(dStarHere) - pay * (1.0_wp - dStarHere) * bNow + theQ &
      * (bPrime - bNow * (1.0_wp - delta + kappa * pay * dStarHere)) + histCash(tIx)

    IF (CdInterior > CdZero) THEN
      histD(tIx) = dStarHere
      histC(tIx) = CdInterior
      histY(tIx) = meanZ * phi(dStarHere)
    ELSE
      histD(tIx) = 0.0_wp
      histC(tIx) = CdZero
      histY(tIx) = meanZ
    END IF

     WRITE (OUTPUT_UNIT, "(2A,I4,3E15.7)") TAB, "S", tIx, bNow, bPrime, histCash(tIx)
     FLUSH (OUTPUT_UNIT)
  END SUBROUTINE goForwardsStat

  SUBROUTINE goForwards(tIx, inHix)
    INTEGER, INTENT(IN) :: tIx, inHix

    REAL(wp) :: cSS, cII, deaths, nSS, nII
    REAL(wp) :: bNow, bPrime, Lhere, qHere
    REAL(wp) :: dStarHere, CdZero, CdInterior

    CALL loadSirPeriod(inHix)
    cSS = histS(tIx)
    cII = histI(tIx)

    bNow = histB(tIx)
    bPrime = iB(cSS, cII, bNow)
    histB(tIx+1) = bPrime
    Lhere = iL(cSS, cII, bNow)
    histL(tIx) = Lhere

    CALL sirLOM(histWedge(tIx) * withL(Lhere), cSS, cII, deaths, nSS, nII)
    qHere = iQ(nSS, nII, bPrime)

    histSp(tIx) = (1.0_wp + pay * (1.0_wp / qHere - 1.0_wp))**52_wp - 1.0_wp
    histV(tIx) = iV(cSS, cII, bNow)

    histS(tIx+1) = nSS
    histI(tIx+1) = nII
    histDeath(tIx+1) = histDeath(tIx) + deaths
    histPop(tIx+1) = histPop(tIx) - deaths
    histR(tIx+1) = 1.0_wp - histDeath(tIx+1) - histS(tIx+1) - histI(tIx+1)

    dStarHere = ((1.0_wp - kappa * qHere) * pay * bNow ) / &
      (histZ(tIx) * (1.0_wp - thetaY * Lhere)**drs * (1.0_wp - gamm2) * gamm0 * gamm1)

    IF (dStarHere >= 1.0_wp) THEN
      dStarHere = 1.0_wp
    ELSE IF (dStarHere <= 0_wp) THEN
      dStarHere = 0.0_wp
    ELSE
      dStarHere = dStarHere**(1.0_wp / (gamm1 - 1.0_wp))
    END IF

    CdZero = histZ(tIx) * (1.0_wp - thetaY * Lhere)**drs - pay * bNow &
      + qHere * (bPrime - bNow * (1.0_wp - delta)) + histCash(tIx)
    CdInterior = histZ(tIx) * (1.0_wp - thetaY * Lhere)**drs * phi(dStarHere) &
      - pay * (1.0_wp - dStarHere) * bNow + qHere * (bPrime - bNow * (1.0_wp - delta + kappa * pay * dStarHere)) &
      + histCash(tIx)

    IF (CdInterior > CdZero) THEN
      histD(tIx) = dStarHere
      histC(tIx) = CdInterior
      histY(tIx) = histZ(tIx) * (1.0_wp - thetaY * Lhere)**drs * phi(dStarHere)
    ELSE
      histD(tIx) = 0.0_wp
      histC(tIx) = CdZero
      histY(tIx) = histZ(tIx) * (1.0_wp - thetaY * Lhere)**drs
    END IF
    
     WRITE (OUTPUT_UNIT, "(2A,I4,4E15.7)") TAB, "P", tIx, bNow, bPrime, Lhere, histCash(tIx)
     FLUSH (OUTPUT_UNIT)
  END SUBROUTINE goForwards

END MODULE VFI