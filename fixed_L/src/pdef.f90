!
! ``Deadly Debt Crises: COVID-19 in Emerging Markets''
!     Fixed L case for Appendix
!     Arellano Bai Mihalache 2023
!
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

  IF (workerId == 0) WRITE (OUTPUT_UNIT, "(A,2F12.4)") "Range of Ed: ", MINVAL(Ed), MAXVAL(Ed)

  histZ = meanZ
  histL = 0.0_wp
  
  fixedL = 0.0_wp
  CALL loadFixedL()

  CALL mpiBarrier()
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

  CALL mpiFinalize()
CONTAINS

  SUBROUTINE loadFixedL()
    INTEGER :: ix

    OPEN(newunit=iunit, file=ioDir // "perfectL.txt")
    DO ix = 1,H
      READ(iunit, "(ES15.9)") fixedL(ix)
      WRITE (*, *) " Period ", ix, " L of ", fixedL(ix)
    END DO
    CLOSE(iunit)
  END SUBROUTINE loadFixedL

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
        TAB, 0.0_wp, TAB, histWedge(ix), &
        TAB, thisObs%getYear(), TAB, thisObs%getMonth(), TAB, thisObs%getDay()
    END DO
    CLOSE(iunit)

    WRITE (OUTPUT_UNIT, "(3A)") TAB, "Done with ", TRIM(label)
    FLUSH (OUTPUT_UNIT)
  END SUBROUTINE forwardAndSave

END PROGRAM pDef
