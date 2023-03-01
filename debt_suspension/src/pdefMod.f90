MODULE pdefMod
  USE Param
  USE NL
  USE sim
  USE MPIMod
  USE datetime_module
  IMPLICIT NONE

  LOGICAL :: withService = .FALSE.

  REAL(wp), DIMENSION(wedgeSz) :: dataWedge, dataMuD, dataL

  REAL(wp), DIMENSION(muSsz) :: muSgrid
  REAL(wp), DIMENSION(muIsz) :: muIgrid

  REAL(wp), DIMENSION(bSz) :: bGrid
  REAL(wp), DIMENSION(Lsz) :: L, thetaLsq

  REAL(wp), DIMENSION(:), ALLOCATABLE :: exoZ, exoCash
  REAL(wp), DIMENSION(:), ALLOCATABLE :: histZ, histL, histWedge, &
    histS, histI, histR, histDeath, histPop, &
    histB, histC, histD, histSp, histY, histV, histCash

  REAL(wp), DIMENSION(:, :), ALLOCATABLE :: C, d, dStar, polPr
  REAL(wp), DIMENSION(:), ALLOCATABLE :: EC, Ed, EbPr
  REAL(wp), DIMENSION(:), ALLOCATABLE :: V0, V1, q0, q1

  REAL(wp), DIMENSION(:, :), ALLOCATABLE :: Csrv, dSrv, dStarSrv, polPrSrv
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ECsrv, EdSrv, EbPrSrv
  REAL(wp), DIMENSION(:), ALLOCATABLE :: V0srv, V1srv, q0srv, q1srv

  REAL(wp), DIMENSION(:, :, :), ALLOCATABLE :: sirV0, sirV1, sirq0, sirq1
  REAL(wp), DIMENSION(:, :, :), ALLOCATABLE :: sirEL, sirEbPr, sirEd, sirEC
  REAL(wp), DIMENSION(:, :, :, :), ALLOCATABLE :: sird, sirPolPr
  INTEGER(sip), DIMENSION(:, :, :, :), ALLOCATABLE :: sirLix

CONTAINS

  SUBROUTINE initialValues()
    INTEGER :: ix

    V0 = uu(meanZ) / (1.0_wp - beta)
    q0 = 1.0_wp
    C = meanZ
    d = 0.0_wp
    polPr = 0.0_wp
    DO ix = 1,bSz
      polPr(ix, ix) = 1.0_wp
    END DO
  END SUBROUTINE initialValues

  PURE SUBROUTINE sirLOM(infect, cSS, cII, deaths, nSS, nII)
    REAL(wp), INTENT(IN) :: infect, cSS, cII
    REAL(wp), INTENT(OUT) :: deaths, nSS, nII
    REAL(wp) :: newInfect, deathRate

    deathRate = piD + piDsq * cII
    newInfect = infect * cSS * cII

    nSS = MAX(0.0_wp, cSS - newInfect)
    nII = (1.0_wp - pp) * cII + newInfect
    deaths = deathRate * cII
  END SUBROUTINE sirLOM

  ELEMENTAL FUNCTION uu(cc) RESULT(util)
    REAL(wp), INTENT(IN) :: cc
    REAL(wp) :: util

    IF (cc <= 0.0_wp) THEN
      util = veryNeg
    ELSE
!      IF (crra == 1.0_wp) THEN
!        util = LOG(cc)
!      ELSE IF (crra == 2.0_wp) THEN
         util = 1.0_wp - 1.0_wp / cc
!      ELSE IF (crra == 0.0_wp) THEN
!        util = cc
!      ELSE
!        util = (cc**(1.0_wp - crra) - 1.0_wp) / (1.0_wp - crra)
!      END IF
    END IF
  END FUNCTION uu

  ELEMENTAL PURE FUNCTION withL(lVal)
    REAL(wp), INTENT(IN) :: lVal
    REAL(wp) :: withL

    withL = (1.0_wp - theta * lVal)**2
  END FUNCTION withL

  SUBROUTINE generateGrids()
    CALL linspace(bGrid, minB, maxB, bSz)
    CALL linspace(L, minL, maxL, Lsz)
    thetaLsq = withL(L)
  END SUBROUTINE generateGrids

  SUBROUTINE generateSirGrid()
    REAL(wp), PARAMETER :: shiftS = 0.01_wp, shiftI = 0.005_wp
    INTEGER, PARAMETER :: sszz = 52*5
    INTEGER :: ix, iunit
    REAL(wp) :: minI, maxI, minS, maxS, deaths
    REAL(wp), DIMENSION(muSsz) :: workS
    REAL(wp) :: wedgePath(sszz), pathS(sszz), pathI(sszz)

    wedgePath(1:(wave1Sz+wave2Sz)) = dataWedge(startWave1:endWave2)
    DO ix = (wave1Sz+wave2Sz+1),sszz
      wedgePath(ix) = rho_psi * wedgePath(ix-1) + (1.0_wp - rho_psi) * psi_wave1_end * pp
    END DO

    pathS = 1.0_wp - defaultInitI - defaultInitR
    pathI = defaultInitI

    DO ix = 2,sszz
      CALL sirLOM(wedgePath(ix-1), pathS(ix-1), pathI(ix-1), deaths, pathS(ix), pathI(ix))
    END DO

    minI = 0.5_wp * defaultInitI
    maxI = MAXVAL(pathI) !  + 0.001_wp

    minS = MINVAL(pathS)
    maxS = MIN(1.0_wp, 1.0_wp - defaultInitI - defaultInitR + 0.005_wp)
    CALL linspace(workS, LOG(shiftS), LOG(maxS-minS+shiftS), muSsz)
    workS = maxS - (EXP(workS) - shiftS)
    DO ix = 1,muSsz
      muSgrid(ix) = workS(muSsz - ix + 1)
    END DO

    CALL linspace(muIgrid, LOG(minI+shiftI), LOG(maxI+shiftI), muIsz)
    muIgrid = EXP(muIgrid) - shiftI

    IF (workerId == 0) THEN
      OPEN(newunit=iunit, file=outDir // "muS.tab")
      WRITE (iunit, "(F23.10)") muSgrid
      CLOSE(iunit)

      OPEN(newunit=iunit, file=outDir // "muI.tab")
      WRITE (iunit, "(F23.10)") muIgrid
      CLOSE(iunit)

      OPEN(newunit=iunit, file=outDir // "noLpaths.tab")
      WRITE (iunit, "(3(A23,A))") "Wedge", TAB, "muS", TAB, "muI", TAB
      DO ix = 1,sszz
        WRITE (iunit, "(3(F23.10,A))") wedgePath(ix), TAB, pathS(ix), TAB, pathI(ix), TAB
      END DO
      CLOSE(iunit)
    END IF
  END SUBROUTINE generateSirGrid

  ELEMENTAL FUNCTION phi(dintensity) RESULT(phiVal)
    REAL(wp), INTENT(IN) :: dintensity
    REAL(wp) :: phiVal

    IF (dintensity > 0.0_wp) THEN
      phiVal = (1.0_wp - gamm0 * dintensity**gamm1) * (1.0_wp - gamm2)
    ELSE
      phiVal = 1.0_wp
    END IF
  END FUNCTION phi

  ! Grids and memory  
  SUBROUTINE allocateAll()
    ALLOCATE( exoZ(H), exoCash(H) )
    ALLOCATE( &
      histZ(Tsz), histL(Tsz), histWedge(Tsz), &
      histS(Tsz), histI(Tsz), histR(Tsz), &
      histDeath(Tsz), histPop(Tsz), histB(Tsz), &
      histC(Tsz), histD(Tsz), histSp(Tsz), &
      histY(Tsz), histV(Tsz), histCash(Tsz) )

    ALLOCATE( C(bSz, bSz) )
    ALLOCATE( d(bSz, bSz) )
    ALLOCATE( dStar(bSz, bSz) )
    ALLOCATE( polPr(bSz, bSz) )

    ALLOCATE( EC(bSz) )
    ALLOCATE( Ed(bSz) )
    ALLOCATE( EbPr(bSz) )

    ALLOCATE( V0(bSz) )
    ALLOCATE( V1(bSz) )
    ALLOCATE( q0(bSz) )
    ALLOCATE( q1(bSz) )

    ALLOCATE( Csrv(bSz, bSz) )
    ALLOCATE( dSrv(bSz, bSz) )
    ALLOCATE( dStarSrv(bSz, bSz) )
    ALLOCATE( polPrSrv(bSz, bSz) )

    ALLOCATE( ECsrv(bSz) )
    ALLOCATE( EdSrv(bSz) )
    ALLOCATE( EbPrSrv(bSz) )

    ALLOCATE( V0srv(bSz) )
    ALLOCATE( V1srv(bSz) )
    ALLOCATE( q0srv(bSz) )
    ALLOCATE( q1srv(bSz) )

    ALLOCATE( sirV0(muSsz, muIsz, bSz) )
    ALLOCATE( sirV1(muSsz, muIsz, bSz) )
    ALLOCATE( sirq0(muSsz, muIsz, bSz) )
    ALLOCATE( sirq1(muSsz, muIsz, bSz) )

    ALLOCATE( sirEL(muSsz, muIsz, bSz) )
    ALLOCATE( sirEbPr(muSsz, muIsz, bSz) )
    ALLOCATE( sirEd(muSsz, muIsz, bSz) )
    ALLOCATE( sirEC(muSsz, muIsz, bSz) )

    ALLOCATE( sird(muSsz, muIsz, bSz, bSz) )
    ALLOCATE( sirPolPr(muSsz, muIsz, bSz, bSz) )
    ALLOCATE( sirLix(muSsz, muIsz, bSz, bSz) )
  END SUBROUTINE allocateAll

  ! I/O
  SUBROUTINE saveEquilibrium()
    CHARACTER (LEN=4) :: stub
    INTEGER :: iunit

    IF (withService) THEN
      stub = "srv_"
    ELSE
      stub = "reg_"
    END IF

    OPEN(newunit=iunit, file=outDir // "parameters.tab")
    WRITE (iunit, "(I20)") muSsz
    WRITE (iunit, "(I20)") muIsz
    WRITE (iunit, "(I20)") bSz
    WRITE (iunit, "(I20)") Tsz
    WRITE (iunit, "(I20)") H
    WRITE (iunit, "(I20)") before
    WRITE (iunit, "(I20)") after
    WRITE (iunit, "(E20.10)") chi
    WRITE (iunit, "(E20.10)") crra
    WRITE (iunit, "(E20.10)") beta
    WRITE (iunit, "(E20.10)") rf
    WRITE (iunit, "(E20.10)") duration
    WRITE (iunit, "(E20.10)") delta
    WRITE (iunit, "(E20.10)") pay
    WRITE (iunit, "(E20.10)") kappa
    WRITE (iunit, "(E20.10)") gamm0
    WRITE (iunit, "(E20.10)") gamm1
    WRITE (iunit, "(E20.10)") gamm2
    WRITE (iunit, "(E20.10)") pp
    WRITE (iunit, "(E20.10)") rho_psi
    WRITE (iunit, "(E20.10)") psi_wave1_end
    WRITE (iunit, "(E20.10)") psi_wave2_end
    WRITE (iunit, "(E20.10)") piD
    WRITE (iunit, "(E20.10)") piDsq
    WRITE (iunit, "(E20.10)") theta
    WRITE (iunit, "(E20.10)") thetaY
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "bGrid.tab")
    WRITE (iunit, "(F30.20)") bGrid
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "V.tab")
    WRITE (iunit, "(F30.20)") V1
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "q.tab")
    WRITE (iunit, "(F30.20)") q1
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "polPr.tab")
    WRITE (iunit, "(F30.20)") polPr
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "C.tab")
    WRITE (iunit, "(F30.20)") C
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "d.tab")
    WRITE (iunit, "(F30.20)") d
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "dStar.tab")
    WRITE (iunit, "(F30.20)") dStar
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "EC.tab")
    WRITE (iunit, "(F30.20)") EC
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "Ed.tab")
    WRITE (iunit, "(F30.20)") Ed
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "EbPr.tab")
    WRITE (iunit, "(F30.20)") EbPr
    CLOSE(iunit)
  END SUBROUTINE saveEquilibrium

  SUBROUTINE loadEquilibrium()
    CHARACTER (LEN=4) :: stub
    INTEGER :: iunit

    IF (withService) THEN
      stub = "srv_"
    ELSE
      stub = "reg_"
    END IF

    OPEN(newunit=iunit, file=outDir // stub // "V.tab")
    IF (withService) THEN
      READ (iunit, "(F30.20)") V0srv
      V1srv = V0srv
    ELSE
      READ (iunit, "(F30.20)") V0
      V1 = V0
    END IF
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "q.tab")
    IF (withService) THEN
      READ (iunit, "(F30.20)") q0srv
      q1srv = q0srv
    ELSE
      READ (iunit, "(F30.20)") q0
      q1 = q0
    END IF
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // stub // "polPr.tab")
    IF (withService) THEN
      READ (iunit, "(F30.20)") polPrSrv
    ELSE
      READ (iunit, "(F30.20)") polPr
    END IF

    OPEN(newunit=iunit, file=outDir // stub // "C.tab")
    IF (withService) THEN
      READ (iunit, "(F30.20)") Csrv
    ELSE
      READ (iunit, "(F30.20)") C
    END IF

    OPEN(newunit=iunit, file=outDir // stub // "d.tab")
    IF (withService) THEN
      READ (iunit, "(F30.20)") dSrv
    ELSE
      READ (iunit, "(F30.20)") d
    END IF

    OPEN(newunit=iunit, file=outDir // stub // "Ed.tab")
    IF (withService) THEN
      READ (iunit, "(F30.20)") EdSrv
    ELSE
      READ (iunit, "(F30.20)") Ed
    END IF

    OPEN(newunit=iunit, file=outDir // stub // "EC.tab")
    IF (withService) THEN
      READ (iunit, "(F30.20)") ECsrv
    ELSE
      READ (iunit, "(F30.20)") EC
    END IF

    OPEN(newunit=iunit, file=outDir // stub // "EbPr.tab")
    IF (withService) THEN
      READ (iunit, "(F30.20)") EbPrSrv
    ELSE
      READ (iunit, "(F30.20)") EbPr
    END IF
  END SUBROUTINE loadEquilibrium

  SUBROUTINE saveSirPeriod(tt)
    CHARACTER (LEN=13) :: fileName
    INTEGER, INTENT(IN) :: tt
    INTEGER :: iunit

    WRITE (fileName, "(A,I3.3,A)") "sir_V_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) sirV1
    CLOSE(iunit)

    WRITE (fileName, "(A,I3.3,A)") "sir_q_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) sirq1
    CLOSE(iunit)

    WRITE (fileName, "(A,I3.3,A)") "sirEL_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) sirEL
    CLOSE(iunit)

    WRITE (fileName, "(A,I3.3,A)") "sirEC_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) sirEC
    CLOSE(iunit)

    WRITE (fileName, "(A,I3.3,A)") "sirEd_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) sirEd
    CLOSE(iunit)

    WRITE (fileName, "(A,I3.3,A)") "sirEb_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) sirEbPr
    CLOSE(iunit)
  END SUBROUTINE saveSirPeriod

  SUBROUTINE loadSirPeriod(tt)
    CHARACTER (LEN=13) :: fileName
    INTEGER, INTENT(IN) :: tt
    INTEGER :: iunit

    WRITE (fileName, "(A,I3.3,A)") "sir_V_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="old")
    READ (iunit) sirV0
    CLOSE(iunit)

    WRITE (fileName, "(A,I3.3,A)") "sir_q_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="old")
    READ (iunit) sirq0
    CLOSE(iunit)

    WRITE (fileName, "(A,I3.3,A)") "sirEL_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="old")
    READ (iunit) sirEL
    CLOSE(iunit)

    WRITE (fileName, "(A,I3.3,A)") "sirEC_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="old")
    READ (iunit) sirEC
    CLOSE(iunit)

    WRITE (fileName, "(A,I3.3,A)") "sirEd_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="old")
    READ (iunit) sirEd
    CLOSE(iunit)

    WRITE (fileName, "(A,I3.3,A)") "sirEb_", tt, ".bin"
    OPEN(newunit=iunit, file=outDir // fileName, &
      FORM="unformatted", ACCESS="stream", STATUS="old")
    READ (iunit) sirEbPr
    CLOSE(iunit)
  END SUBROUTINE loadSirPeriod

  SUBROUTINE dataImpliedWedge()
    REAL(wp), DIMENSION(wedgeSz) :: newDeaths, lockdown, muI, newI, muS, bbeta, dataR
    TYPE(datetime), DIMENSION(wedgeSz) :: calDates
    INTEGER, DIMENSION(3, wedgeSz) :: dataDates
    INTEGER :: iunit, ix
    LOGICAL, PARAMETER :: printData = .FALSE.

    OPEN(newunit=iunit, file=dataDir // "data_L.txt", status="OLD")
    READ (iunit, "(ES15.9)") lockdown
    CLOSE(iunit)

    OPEN(newunit=iunit, file=dataDir // "data_D.txt", status="OLD")
    READ (iunit, "(ES15.9)") newDeaths
    CLOSE(iunit)

    OPEN(newunit=iunit, file=dataDir // "data_dates.txt", status="OLD")
    READ (iunit, *) dataDates
    CLOSE(iunit)

    muI = (-piD + SQRT(piD * piD + 4.0_wp * piDsq * newDeaths) ) / (2.0_wp * piDsq)
    newI = 0.0_wp
    DO ix = 1,wedgeSz-1
      newI(ix) = muI(ix+1) - (1_wp - pp) * muI(ix)
    END DO
    
    muS(1) = defaultInitS
    DO ix = 2,wedgeSz
      muS(ix) = muS(ix-1) - newI(ix-1)
    END DO

    bbeta = newI / (muI * muS)
    dataR = bbeta * muS / pp

    dataWedge = bbeta / ((1.0_wp - theta * lockdown)**2)

    dataL = lockdown
    dataMuD(1) = newDeaths(1)
    DO ix = 2,wedgeSz
      dataMuD(ix) = dataMuD(ix-1) + newDeaths(ix)
    END DO

    DO ix = 1,wedgeSz
      calDates(ix) = datetime(year=dataDates(1, ix), month=dataDates(2, ix), day=dataDates(3, ix))
    END DO
    eventYear = dataDates(1, wedgeStart)
    eventMonth = dataDates(2, wedgeStart)
    eventDay =  dataDates(3, wedgeStart)

    defaultInitI = muI(wedgeStart)
    defaultInitS = muS(wedgeStart)
    defaultInitR = 1.0_wp  - defaultInitI - defaultInitS

    IF (workerId == 0) THEN
      OPEN(newunit=iunit, file=outDir // "data.tab")
      DO ix = wedgeStart,wedgeEnd
        WRITE(iunit, "(3(I4,A),3(ES12.4,A))") &
          dataDates(1, ix), TAB, dataDates(2, ix), TAB, dataDates(3, ix), TAB, &
          dataL(ix), TAB, dataMuD(ix), TAB, dataWedge(ix), TAB
      END DO 
      CLOSE(iunit)
    END IF

    IF (printData .AND. workerId == 0) THEN
      WRITE (*, "(A25,A,9(A12,A))") "week", TAB, "ix", TAB, "Lockdown", TAB, "New Deaths", TAB, &
        "newI", TAB, "muI", TAB, "muS", TAB, "beta", TAB, "R", TAB, "wedge", TAB
      DO ix = 1,wedgeSz
        WRITE (*, "(A25,A,I12,A,8(ES12.4,A))") calDates(ix)%isoformat(), TAB, ix, TAB, &
          lockdown(ix), TAB, newDeaths(ix), TAB, newI(ix), TAB, muI(ix), TAB, muS(ix), TAB, &
          bbeta(ix), TAB, dataR(ix), TAB, dataWedge(ix), TAB
      END DO
    END IF
  END SUBROUTINE

  !
  ! Interpolations
  !
  PURE FUNCTION iV(wS, wI, wB)
    REAL(wp), INTENT(IN) :: wS, wI, wB
    REAL(wp) :: iV
  
    iV = trilinear(muSgrid, muIgrid, bGrid, sirV0, wS, wI, wB)
  END FUNCTION iV

  PURE FUNCTION iVi(wS, wI, bPrIx)
    REAL(wp), INTENT(IN) :: wS, wI
    INTEGER, INTENT(IN) :: bPrIx
    REAL(wp) :: iVi
  
    iVi = bilinear(muSgrid, muIgrid, sirV0(:, :, bPrIx), wS, wI)
  END FUNCTION iVi

  PURE FUNCTION iQ(wS, wI, wB, useOne)
    REAL(wp), INTENT(IN) :: wS, wI, wB
    LOGICAL, INTENT(IN), OPTIONAL :: useOne
    REAL(wp) :: iQ
   
    IF (PRESENT(useOne) .AND. useOne) THEN
      iQ = trilinear(muSgrid, muIgrid, bGrid, sirq1, wS, wI, wB)
    ELSE
      iQ = trilinear(muSgrid, muIgrid, bGrid, sirq0, wS, wI, wB)
    END IF
    iQ = MAX(0.0_wp, MIN(1.0_wp, iQ))
  END FUNCTION iQ

  PURE FUNCTION iQi(wS, wI, bPrIx, useOne)
    REAL(wp), INTENT(IN) :: wS, wI
    INTEGER, INTENT(IN) :: bPrIx
    LOGICAL, INTENT(IN), OPTIONAL :: useOne
    REAL(wp) :: iQi
   
    IF (PRESENT(useOne) .AND. useOne) THEN
      iQi = bilinear(muSgrid, muIgrid, sirq1(:, :, bPrIx), wS, wI)
    ELSE
      iQi = bilinear(muSgrid, muIgrid, sirq0(:, :, bPrIx), wS, wI)
    END IF
    iQi = MAX(0.0_wp, MIN(1.0_wp, iQi))
  END FUNCTION iQi

  PURE FUNCTION iB(wS, wI, wB)
    REAL(wp), INTENT(IN) :: wS, wI, wB
    REAL(wp) :: iB

    iB = trilinear(muSgrid, muIgrid, bGrid, sirEbPr, wS, wI, wB)
    iB = MAX(minB, MIN(maxB, iB))
  END FUNCTION iB

  PURE FUNCTION iL(wS, wI, wB)
    REAL(wp), INTENT(IN) :: wS, wI, wB
    REAL(wp) :: iL
   
    iL = trilinear(muSgrid, muIgrid, bGrid, sirEL, wS, wI, wB)
    iL = MAX(minL, MIN(maxL, iL))
  END FUNCTION iL

  SUBROUTINE tstamp()
    TYPE(datetime) :: dt

    dt = dt%now()
    WRITE (OUTPUT_UNIT, "(A,3(I0.2,A))", ADVANCE="no") " ~> ", dt%getHour(), ":", dt%getMinute(), ":", dt%getSecond(), TAB
  END SUBROUTINE tstamp

END MODULE pdefMod
