!
! ``Deadly Debt Crises: COVID-19 in Emerging Markets''
!           Persistent recession model 
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
    inMTU = .FALSE.
    CALL initialValues()
    CALL compute()
    CALL expectedPolicies()
    IF (workerId == 0) CALL saveEquilibrium()
    IF (workerId == 0) WRITE (OUTPUT_UNIT, *) "Wrote down regular equilibrium."

    Vrecover = V1
    qRecover = q1
    dRecover = d
    polPrRecover = polPr

    inMTU = .TRUE.
    CALL initialValues()
    CALL compute()
    CALL expectedPolicies()
    IF (workerId == 0) CALL saveEquilibrium()
    IF (workerId == 0) WRITE (OUTPUT_UNIT, *) "Wrote down MTU equilibrium."

    Cmtu = C
    dmtu = d
    dStarMtu = dStar
    polPrMtu = polPr
    ECmtu = EC
    Edmtu = Ed
    EbPrMtu = EbPr
    V0mtu = V0
    V1mtu = V1
    q0mtu = q0
    q1mtu = q1

    IF (workerId == 0) WRITE (OUTPUT_UNIT, *) "At barrier..."
    FLUSH (OUTPUT_UNIT)
    CALL mpiBarrier()
    IF (workerId == 0) WRITE (OUTPUT_UNIT, *) "Now all load..."

    inMTU = .FALSE.
    CALL loadEquilibrium()
  ELSE
    IF (workerId == 0) WRITE (OUTPUT_UNIT, *) "Load MTU."
    inMTU = .TRUE.
    CALL loadEquilibrium()

    IF (workerId == 0) WRITE (OUTPUT_UNIT, *) "Load regular."
    inMTU = .FALSE.
    CALL loadEquilibrium()
  END IF
  FLUSH (OUTPUT_UNIT)

  histZ = meanZ
  histL = 0.0_wp

  CALL generateSirGrid()

  ! Wave 1 wedge path
  histWedge = dataWedge(startWave1)
  histWedge((before + 1):(before+wave1Sz+wave2Sz)) = dataWedge(startWave1:endWave2)
  DO tIx = (before+wave1Sz+wave2Sz+1),Tsz
    histWedge(tIx) = rho_psi * histWedge(tIx - 1) + (1.0_wp - rho_psi) * psi_wave1_end * pp
  END DO

  CALL mpiBarrier()
  IF (runBack) THEN
    IF (workerId == 0) THEN
      WRITE (OUTPUT_UNIT, *) "Running backwards..."
      FLUSH (OUTPUT_UNIT)
    END IF

    CALL setupTerminal()
    CALL mpiBarrier()
    DO tIx = H-1,1,-1
      CALL goBackwards(tIx)
    END DO
  END IF ! runBack
  DEALLOCATE( sird )
  DEALLOCATE( sirPolPr )
  DEALLOCATE( sirLix )

  IF (workerId == 0) THEN
    CALL forwardAndSave("persistent")
    WRITE (OUTPUT_UNIT, *) "Done backwards, pass 1."
    FLUSH (OUTPUT_UNIT)
  END IF
  CALL mpiFinalize()
CONTAINS

  SUBROUTINE forwardAndSave(label)
    CHARACTER (LEN=*), INTENT(IN) :: label
    TYPE(datetime) :: dateObj, thisObs
    TYPE(timedelta) :: dateStep
    INTEGER :: ix, pid, daysShift
    REAL(wp) :: deaths, draw, roll
    LOGICAL :: off

    TYPE(simPoint) :: simp

    INTEGER, PARAMETER :: simWide = 10000 
    REAL(wp), DIMENSION(after+1) :: &
      ahistZ, ahistL, ahistB, ahistC, ahistD, ahistSp, ahistY, ahistV

    !
    ! Forward
    !
    histDeath = 0.0_wp
    histPop = 1.0_wp
  
    histS(1:before) = 1.0_wp - defaultInitR
    histI(1:before) = 0.0_wp
    histR(1:before) = defaultInitR
  
    histS(before + 1) = defaultInitS
    histI(before + 1) = defaultInitI
    histR(before + 1) = defaultInitR
  
    histDeath(before + 1) = (piD + piDsq * defaultInitI) * defaultInitI
    histPop(before + 1) = 1.0_wp - histDeath(before + 1)
 
    histZ = meanZ
    histB = initB

    ! Before pandemic
    DO ix = 1,before
      simp = goForwardsStat(histB(ix), .FALSE.)

      histB(ix + 1) = simp%bPrime
      histV(ix) = simp%V
      histSp(ix) = simp%sp
      histD(ix) = simp%d
      histL(ix) = simp%L
      histC(ix) = simp%C
      histY(ix) = simp%Y

      IF (workerId == 0) WRITE (OUTPUT_UNIT, *) "Before ", ix, " B= ",  histB(ix)/52
    END DO

    ! Pandemic
    DO ix = before + 1,before+H-1
      CALL loadSirPeriod(ix - before)

      simp = goForwards(histS(ix), histI(ix), histB(ix), histWedge(ix), histZ(ix), ix-before)

      histS(ix + 1) = simp%muSprime
      histI(ix + 1) = simp%muIprime
      IF (ix > before + 1) THEN
        histPop(ix) = histPop(ix - 1) - simp%newDeaths
        histDeath(ix) = histDeath(ix - 1) + simp%newDeaths
        histR(ix) = 1.0_wp - histS(ix) - histI(ix) - histDeath(ix)
      END IF

      histB(ix + 1) = simp%bPrime
      histV(ix) = simp%V
      histSp(ix) = simp%sp
      histD(ix) = simp%d
      histL(ix) = simp%L
      histC(ix) = simp%C
      histY(ix) = simp%Y

      IF (workerId == 0) WRITE (OUTPUT_UNIT, *) "During ", ix, " B = ", histB(ix)/52
    END DO

    ! After pandemic -- SIRD
    histS(before+H:Tsz) = 0.0_wp
    DO ix = before+H,Tsz-1
      CALL sirLOM(histWedge(ix), &
        histS(ix), histI(ix), deaths, histS(ix+1), histI(ix+1)) 

      histPop(ix) = histPop(ix - 1) - deaths
      histDeath(ix) = histDeath(ix - 1) + deaths
      histR(ix) = 1.0_wp - histS(ix) - histI(ix) - histDeath(ix)
    END DO

    ! After pandemic panel
    DO pid = 1,simWide
      ahistB(1) = histB(before+H)
      off = .FALSE.
      DO ix = 1,after
        IF (ix > 1) THEN
          IF (.NOT. off) THEN
            CALL simUniform(draw)
            IF (draw <= mtuP) THEN
              off = .TRUE.
            END IF
          END IF
        END IF

        simp = goForwardsStat(ahistB(ix), .NOT. off)

        IF (off) THEN
          ahistZ(ix) = meanZ
        ELSE
          ahistZ(ix) = mtuZ
        END IF

        ahistB(ix + 1) = simp%bPrime
        ahistV(ix) = simp%V
        ahistSp(ix) = simp%sp
        ahistD(ix) = simp%d
        ahistL(ix) = simp%L
        ahistC(ix) = simp%C
        ahistY(ix) = simp%Y
      END DO

      roll = (pid - 1.0_wp) / (pid + 0.0_wp)
      
      histB(before + H + 1:Tsz) = roll * histB(before + H + 1:Tsz) + (1.0_wp - roll) * ahistB(2:after+1)

      histZ(before + H:Tsz) = roll * histZ(before + H:Tsz) + (1.0_wp - roll) * ahistZ
      histV(before + H:Tsz) = roll * histV(before + H:Tsz) + (1.0_wp - roll) * ahistV
      histSp(before + H:Tsz) = roll * histSp(before + H:Tsz) + (1.0_wp - roll) * ahistSp
      histD(before + H:Tsz) = roll * histD(before + H:Tsz) + (1.0_wp - roll) * ahistD
      histL(before + H:Tsz) = roll * histL(before + H:Tsz) + (1.0_wp - roll) * ahistL
      histC(before + H:Tsz) = roll * histC(before + H:Tsz) + (1.0_wp - roll) * ahistC
      histY(before + H:Tsz) = roll * histY(before + H:Tsz) + (1.0_wp - roll) * ahistY

      IF (workerId == 0 .AND. MOD(pid, 100) == 1) &
        WRITE (OUTPUT_UNIT, "(A,I5,A,I5,A,F10.4)") "Sim after ", pid, " of ", simWide, ", roll ", roll
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
        TAB, 0.0_wp, TAB, histWedge(ix), TAB, thisObs%getYear(), TAB, &
        thisObs%getMonth(), TAB, thisObs%getDay()
    END DO
    CLOSE(iunit)
  
    WRITE (OUTPUT_UNIT, "(3A)") TAB, "Done with ", TRIM(label)
    FLUSH (OUTPUT_UNIT)
  END SUBROUTINE forwardAndSave

END PROGRAM pDef
