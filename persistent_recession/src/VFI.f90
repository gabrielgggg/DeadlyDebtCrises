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
    REAL(wp) :: theMax, theSum, dStarHere, CdZero, CdInterior, theZ
    REAL(wp), DIMENSION(bSz) :: theExps
    REAL(wp), DIMENSION(bSz, bSz) :: uVal

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(bIx,bPrIx,dStarHere,CdZero,CdInterior,theZ)
    DO bIx = 1,bSz
    DO bPrIx = 1,bSz
      IF (inMTU) THEN
        theZ = mtuZ
      ELSE
        theZ = meanZ
      END IF

      dStarHere = ((1.0_wp - kappa * q0(bPrIx)) * pay * bGrid(bIx) ) / &
        (theZ * (1.0_wp - gamm2) * gamm0 * gamm1)

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

      CdZero = theZ - pay * bGrid(bIx) &
        + q0(bPrIx) * (bGrid(bPrIx) - bGrid(bIx) * (1.0_wp - delta))
      CdInterior = theZ * phi(dStarHere) - pay * (1.0_wp - dStarHere) * bGrid(bIx) &
        + q0(bPrIx) * (bGrid(bPrIx) - bGrid(bIx) * (1.0_wp - delta + kappa * pay * dStarHere))

      IF (CdInterior > CdZero) THEN
        d(bIx, bPrIx) = dStarHere
        C(bIx, bPrIx) = CdInterior
      ELSE
        d(bIx, bPrIx) = 0.0_wp
        C(bIx, bPrIx) = CdZero
      END IF
        IF (inMTU) C(bIx, bPrIx) = C(bIx, bPrIx) - exoLoanPay
    END DO
    END DO

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(bIx, bPrIx)
    DO bIx = 1,bSz
    DO bPrIx = 1,bSz
      IF (C(bIx, bPrIx) <= 0.0_wp) THEN
        uVal(bIx, bPrIx) = veryNeg
      ELSE
        IF (inMTU) THEN
          uVal(bIx, bPrIx) = uu(C(bIx, bPrIx)) &
            + beta * ((1.0_wp - mtuP) * V0(bPrIx) + mtuP * Vrecover(bPrIx))
        ELSE
          uVal(bIx, bPrIx) = uu(C(bIx, bPrIx)) + beta * V0(bPrIx)
        END IF
      END IF
    END DO
    END DO

    !$OMP PARALLEL DO PRIVATE(bIx,theExps,theMax,theSum,bPrIx)
    DO bIx = 1,bSz
      theMax = MAXVAL( uVal(bIx, :) ) ! , uVal(bIx, :) > veryNeg )
      theExps = EXP( (uVal(bIx, :) - theMax) / tasteB )
      theSum = SUM( theExps )

      IF (theSum > 0.0_wp) THEN
        polPr(bIx, :) = theExps / theSum
        V1(bIx) = theMax + tasteB * (LOG(theSum) - eulerMascheroni)
      ELSE
        polPr(bIx, :) = 0.0_wp
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
    REAL(wp) :: qTmp

    !$OMP PARALLEL DO PRIVATE(bIx,bPrIx,qTmp) 
    DO bIx = 1,bSz
      qTmp = 0.0_wp
      DO bPrIx = 1,bSz
        IF (inMTU) THEN
          qTmp = qTmp + &
            (1.0_wp - mtuP) * polPr(bIx, bPrIx) * &
              ( pay * (1.0_wp - d(bIx, bPrIx)) &
              + (1.0_wp - delta + kappa * pay * d(bIx, bPrIx) ) * q0(bPrIx) ) &
            + mtuP * polPrRecover(bIx, bPrIx) * &
              ( pay * (1.0_wp - dRecover(bIx, bPrIx)) &
              + (1.0_wp - delta + kappa * pay * dRecover(bIx, bPrIx) ) * qRecover(bPrIx) )
        ELSE
          qTmp = qTmp + polPr(bIx, bPrIx) * &
              ( pay * (1.0_wp - d(bIx, bPrIx)) &
              + (1.0_wp - delta + kappa * pay * d(bIx, bPrIx) ) * q0(bPrIx) )
        END IF
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
      sirV1(sIx, iIx, :) = V1mtu
      sirq1(sIx, iIx, :) = q1mtu
      sirEbPr(sIx, iIx, :) = EbPrMtu
      sird(sIx, iIx, :, :) = dMtu
      sirPolPr(sIx, iIx, :, :) = polPrMtu
    END DO
    END DO

    sirEL = 0.0_wp
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

        CdZero = histZ(hIx) * (1.0_wp - thetaY * L(Lix))**drs &
          - pay * bGrid(bIx) + qTmp * (bGrid(bPrIx) - bGrid(bIx) * (1.0_wp - delta))
          
        CdInterior = histZ(hIx) * (1.0_wp - thetaY * L(Lix))**drs * phi(dStarHere) & 
          - pay * (1.0_wp - dStarHere) * bGrid(bIx) + qTmp &
          * (bGrid(bPrIx) - bGrid(bIx) * (1.0_wp - delta + kappa * pay * dStarHere))

        IF (CdInterior > CdZero) THEN
          dvals(Lix) = dStarHere
          Cvals(Lix) = CdInterior
        ELSE
          dvals(Lix) = 0.0_wp
          Cvals(Lix) = CdZero
        END IF

        IF (hIx == 1) Cvals(Lix) = Cvals(Lix) + exoLoan

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
  PURE FUNCTION goForwardsStat(bNow, mtuFlag) RESULT(simp)
    REAL(wp), INTENT(IN) :: bNow
    LOGICAL, INTENT(IN) :: mtuFlag
    TYPE(simPoint) :: simp

    REAL(wp) :: bPrime, theQ, dStarHere, CdZero, CdInterior, theZ

    IF (.NOT. mtuFlag) THEN
      simp%V = linInterp(bGrid, V0, bNow)
      bPrime = linInterp(bGrid, EbPr, bNow)
    ELSE
      simp%V = linInterp(bGrid, V0mtu, bNow)
      bPrime = linInterp(bGrid, EbPrMtu, bNow)
    END IF

    bPrime = MAX(minB, MIN(maxB, bPrime))
    simp%bPrime = bPrime

    IF (.NOT. mtuFlag) THEN
      theQ = linInterp(bGrid, q0, bPrime)
      theZ = meanZ
    ELSE
      theQ = linInterp(bGrid, q0mtu, bPrime)
      theZ = mtuZ
    END IF

    simp%sp = (1.0_wp + pay * (1.0_wp / theQ - 1.0_wp))**52_wp - 1.0_wp

    dStarHere = ((1.0_wp - kappa * theQ) * pay * bNow ) / &
      (theZ * (1.0_wp - gamm2) * gamm0 * gamm1)

    IF (dStarHere >= 1.0_wp) THEN
      dStarHere = 1.0_wp
    ELSE IF (dStarHere <= 0_wp) THEN
      dStarHere = 0.0_wp
    ELSE
      dStarHere = dStarHere**(1.0_wp / (gamm1 - 1.0_wp))
    END IF

    CdZero = theZ - pay * bNow + theQ * (bPrime - bNow * (1.0_wp - delta))
    CdInterior = theZ * phi(dStarHere) - pay * (1.0_wp - dStarHere) * bNow + theQ &
      * (bPrime - bNow * (1.0_wp - delta + kappa * pay * dStarHere))

    IF (CdInterior > CdZero) THEN
      simp%d = dStarHere
      simp%C = CdInterior
      simp%Y = theZ * phi(dStarHere)
    ELSE
      simp%d = 0.0_wp
      simp%C = CdZero
      simp%Y = theZ
    END IF
    IF (mtuFlag) simp%C = simp%C - exoLoanPay
    simp%L = 0.0_wp
  END FUNCTION goForwardsStat

  PURE FUNCTION goForwards(cSS, cII, bNow, wwedge, useZ, hIx) RESULT(simp)
    REAL(wp), INTENT(IN) :: cSS, cII, bNow, wwedge, useZ
    INTEGER, INTENT(IN) :: hIx
    TYPE(simPoint) :: simp

    REAL(wp) :: deaths, nSS, nII
    REAL(wp) :: bPrime, Lhere, qHere
    REAL(wp) :: dStarHere, CdZero, CdInterior

    bPrime = iB(cSS, cII, bNow)
    simp%bPrime = bPrime
    Lhere = iL(cSS, cII, bNow)
    simp%L = Lhere

    CALL sirLOM(wwedge * withL(Lhere), cSS, cII, deaths, nSS, nII)
    simp%muSprime = nSS
    simp%muIprime = nII
    simp%newDeaths = deaths

    qHere = iQ(nSS, nII, bPrime)
    simp%sp = (1.0_wp + pay * (1.0_wp / qHere - 1.0_wp))**52_wp - 1.0_wp

    simp%V = iV(cSS, cII, bNow)

    dStarHere = ((1.0_wp - kappa * qHere) * pay * bNow ) / &
      (useZ * (1.0_wp - thetaY * Lhere)**drs * (1.0_wp - gamm2) * gamm0 * gamm1)

    IF (dStarHere >= 1.0_wp) THEN
      dStarHere = 1.0_wp
    ELSE IF (dStarHere <= 0_wp) THEN
      dStarHere = 0.0_wp
    ELSE
      dStarHere = dStarHere**(1.0_wp / (gamm1 - 1.0_wp))
    END IF

    CdZero = useZ * (1.0_wp - thetaY * Lhere)**drs - pay * bNow &
      + qHere * (bPrime - bNow * (1.0_wp - delta))
    CdInterior = useZ * (1.0_wp - thetaY * Lhere)**drs * phi(dStarHere) &
      - pay * (1.0_wp - dStarHere) * bNow + qHere * (bPrime - bNow * (1.0_wp - delta + kappa * pay * dStarHere))

    IF (CdInterior > CdZero) THEN
      simp%d = dStarHere
      simp%C = CdInterior
      simp%Y = useZ * (1.0_wp - thetaY * Lhere)**drs * phi(dStarHere)
    ELSE
      simp%d = 0.0_wp
      simp%C = CdZero
      simp%Y = useZ * (1.0_wp - thetaY * Lhere)**drs
    END IF
    IF (hIx == 1) simp%C = simp%C + exoLoan
  END FUNCTION goForwards

END MODULE VFI
