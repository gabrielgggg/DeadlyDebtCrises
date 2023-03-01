!
! ``Deadly Debt Crises: COVID-19 in Emerging Markets''
!     Debt Suspension case
!     Arellano Bai Mihalache 2023
!
PROGRAM pDef
  USE Param
  USE NL
  USE datetime_module
  USE pdefMod
  USE MPIMod
  USE VFI
  IMPLICIT NONE

  INTEGER :: tIx, okFlag, iunit
  REAL(wp) :: deaths

  !
  ! Start MPI
  !
  CALL mpiSetup(okFlag)
  IF (okFlag == 0) THEN
    CALL mpiFinalize()
    IF (workerId == 0) THEN
      WRITE (OUTPUT_UNIT, "(A,I10,A,I3)") "MPI okFlag == 0. ", spaceSz, " ", workerNo
      FLUSH (OUTPUT_UNIT)
    END IF
    STOP
  ELSE
    IF (workerId == 0) CALL tstamp()
  END IF

  ! Init
  CALL dataImpliedWedge()
  CALL allocateAll()
  CALL generateGrids()

  !
  ! Stationary
  !
  IF (.NOT. loadResults) THEN
    CALL initialValues()
    withService = .TRUE.
    CALL compute()
    CALL expectedPolicies()
    IF (workerId == 0) CALL saveEquilibrium()

    Csrv = C
    dSrv = d
    dStarSrv = dStar
    polPrSrv = polPr
    ECsrv = EC
    EdSrv = Ed
    EbPrSrv = EbPr
    V0srv = V0
    V1srv = V1
    q0srv = q0
    q1srv = q1

    withService = .FALSE.
    CALL compute()
    CALL expectedPolicies()
    IF (workerId == 0) CALL saveEquilibrium()
  ELSE
    withService = .TRUE.
    CALL loadEquilibrium()

    withService = .FALSE.
    CALL loadEquilibrium()
  END IF
  FLUSH (OUTPUT_UNIT)

  histZ = meanZ

  histL = 0.0_wp
  histCash = 0.0_wp
  exoCash = 0.0_wp
  IF (exoLoan) THEN
    histCash(before+1) = exoLoanSz
    histCash(before+startRepay:TSz) = -exoLoanPay
    exoCash = histCash(before+1:before+H)
  END IF

  CALL generateSirGrid()

  CALL mpiBarrier()
  IF (runBack) THEN
    IF (workerId == 0) THEN
      WRITE (OUTPUT_UNIT, *) "Running backwards..."
      FLUSH (OUTPUT_UNIT)
    END IF

    ! Wave 1 wedge path
    histWedge = dataWedge(startWave1)
    histWedge((before+1):(before+wave1Sz+wave2Sz)) = dataWedge(startWave1:endWave2)
    DO tIx = (before+wave1Sz+wave2Sz+1),Tsz
      histWedge(tIx) = rho_psi * histWedge(tIx - 1) + (1.0_wp - rho_psi) * psi_wave1_end * pp
    END DO

    CALL setupTerminal()
    CALL mpiBarrier()
    DO tIx = H-1,1,-1
      CALL goBackwards(tIx)
    END DO

    IF (workerId == 0) THEN
      CALL forwardAndSave(simFileName)
      WRITE (OUTPUT_UNIT, *) "Done backwards, pass 1."
      FLUSH (OUTPUT_UNIT)
    END IF
  END IF ! runBack

  CALL mpiFinalize()
CONTAINS

  SUBROUTINE forwardAndSave(label)
    CHARACTER (LEN=*), INTENT(IN) :: label
    TYPE(datetime) :: dateObj, thisObs
    TYPE(timedelta) :: dateStep
    INTEGER :: ix, daysShift

    !
    ! Forward
    !
    histDeath = 0.0_wp
    histPop = 1.0_wp
  
    histS(1:before) = 1.0_wp - defaultInitR
    histI(1:before) = 0.0_wp
    histR(1:before) = defaultInitR
  
    histS(before+1) = defaultInitS
    histI(before+1) = defaultInitI
    histR(before+1) = defaultInitR
  
    histDeath(before+1) = (piD + piDsq * defaultInitI) * defaultInitI
    histPop(before+1) = 1.0_wp - histDeath(before+1)
  
    ! Before Pandemic
    histB = initB
    IF (useResetB) histB(1) = resetB
    DO ix = 1,Tsz-1
      IF (ix <= before) THEN
        ! Before pandemic: (..., before]
        CALL goForwardsStat(ix)
      ELSE IF (ix >= before + H) THEN
        ! After pandemic: [before+H, ...)
        IF (ix == before + H) histS(ix) = 0.0_wp ! Vaccination
  
        CALL sirLOM(histWedge(ix), &
          histS(ix), histI(ix), deaths, histS(ix+1), histI(ix+1)) 
  
        histPop(ix) = histPop(ix-1) - deaths
        histDeath(ix) = histDeath(ix-1) + deaths
        histR(ix) = 1.0_wp - histS(ix) - histI(ix) - histDeath(ix)
  
        CALL goForwardsStat(ix)
      ELSE
        ! During pandemic: [before+1, before+H-1]
        IF  (ix == before + 1 .AND. useResetB) &
          histB(ix) = resetB
        CALL goForwards(ix, ix - before)
      END IF
    END DO
  
    !
    ! Save
    !
    dateObj = datetime(year=eventYear, month=eventMonth, day=eventDay)
  
    OPEN(newunit=iunit, file=outDir // "sir_" // TRIM(label) // ".tab")
    DO ix = 2,Tsz-1
      daysShift = 7 * (ix - (before + 1))
      dateStep = timedelta(days=daysShift)
      thisObs = dateObj + dateStep
  
      WRITE (iunit, "(I8,14(A,F25.10),3(A,I4))") ix-before-1, &
        TAB, histS(ix), TAB, histI(ix), TAB, histR(ix), TAB, histDeath(ix), &
        TAB, histB(ix), TAB, histC(ix), TAB, histD(ix), TAB, histSp(ix), &
        TAB, histY(ix), TAB, histV(ix), TAB, histZ(ix), TAB, histL(ix), &
        TAB, histCash(ix), TAB, histWedge(ix), &
        TAB, thisObs%getYear(), TAB, thisObs%getMonth(), TAB, thisObs%getDay()
    END DO
    CLOSE(iunit)
  
    WRITE (OUTPUT_UNIT, "(3A)") TAB, "Done with ", TRIM(label)
    FLUSH (OUTPUT_UNIT)
  END SUBROUTINE forwardAndSave

END PROGRAM pDef
