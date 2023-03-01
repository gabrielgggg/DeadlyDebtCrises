MODULE perfectMod
  USE iso_Fortran_env, ONLY: wp => real64, OUTPUT_UNIT
  USE NL
  USE sim
  USE datetime_module
  USE nlopt_wrap
  USE nlopt_enum
  IMPLICIT NONE

  REAL(wp), PARAMETER :: xTol = 1.0D-9
  CHARACTER, PARAMETER:: tab = CHAR(9)
  CHARACTER (LEN=*), PARAMETER :: ioDir = "./"
  CHARACTER (LEN=*), PARAMETER :: resultsDir = ioDir // "results/"
  CHARACTER (LEN=*), PARAMETER :: dataDir = ioDir // "../data/"
  CHARACTER (LEN=*), PARAMETER :: simFileName = "perfect"

  LOGICAL, PARAMETER :: surpriseWave2 = .FALSE.

  REAL(wp), PARAMETER :: rho_psi = 0.99_wp
  REAL(wp), PARAMETER :: psi_wave1_end = 1.35_wp
  REAL(wp), PARAMETER :: psi_wave2_end = psi_wave1_end

  REAL(wp), PARAMETER :: crra = 2.0_wp
  REAL(wp), PARAMETER :: beta = 0.98_wp**(1.0_wp / 52.0_wp)
  REAL(wp), PARAMETER :: chi = 3500.0_wp
  REAL(wp), PARAMETER :: drs = 0.67_wp

  INTEGER, PARAMETER :: days = 6
  REAL(wp), PARAMETER :: pp = 1.0_wp - (1.0_wp - 1.0_wp / days)**7
  REAL(wp), PARAMETER :: piD0 = 0.0085_wp
  REAL(wp), PARAMETER :: piD1 = 0.08_wp
  REAL(wp), PARAMETER :: piD = piD0 * pp
  REAL(wp), PARAMETER :: piDsq = piD1 * pp

  REAL(wp), PARAMETER :: minL = 0.0_wp
  REAL(wp), PARAMETER :: maxL = 0.8_wp
  REAL(wp), PARAMETER :: theta = 0.5_wp
  REAL(wp), PARAMETER :: thetaY = 0.8_wp

  REAL(wp), PARAMETER :: rfRate = 1.01_wp**(1.0_wp / 52.0_wp)
  REAL(wp), PARAMETER :: defaultInitB = 0.612966_wp * 52.0_wp

  REAL(wp), PARAMETER :: shrinkRate = (beta * rfRate)**(1.0_wp / crra)
  REAL(wp), PARAMETER :: wealthShare = 1.0_wp - beta**(1.0_wp / crra) * rfRate**((1.0_wp - crra) / crra)
  !
  ! c_t^{-crra} = beta rfRate c_{t+1}^{-\crra}
  ! [ c_{t+1} / c_t ]^{crra} = beta rfRate
  ! c_{t+1} / c_t = ( beta RfRate )^{1 / crra}
  !

  REAL(wp) :: defaultInitI = 0.001_wp
  REAL(wp) :: defaultInitS = 0.965_wp
  REAL(wp) :: defaultInitR = 1.0_wp - 0.965_wp - 0.001_wp 

  ! Set from data, upon execution
  INTEGER :: eventYear = 2020, eventMonth = 3, eventDay = 1

  INTEGER, PARAMETER :: wedgeSz = 97, wedgeStart = 8, wedgeEnd = 91, dataWave2 = 45
  INTEGER, PARAMETER :: startWave1 = wedgeStart, endWave1 = dataWave2-1, wave1Sz = endWave1 - startWave1 + 1
  INTEGER, PARAMETER :: startWave2 = dataWave2, endWave2 = wedgeEnd, wave2Sz = endWave2 - startWave2 + 1
  REAL(wp), DIMENSION(wedgeSz) :: dataWedge, dataMuD, dataL

  TYPE Config
    REAL(wp) :: initS, initI, initR, initB
    REAL(wp) :: cashProgramSize, lossValL, lossValD
    INTEGER :: maxT, freeVaxAt, waveTwoAt

    REAL(wp), DIMENSION(:), ALLOCATABLE :: discountsBeta, discountsRf
    REAL(wp), DIMENSION(:), ALLOCATABLE :: pushCash
    REAL(wp), DIMENSION(:), ALLOCATABLE :: wedge
  END TYPE

  TYPE Solution
    TYPE(Config) :: configuration
    REAL(wp), DIMENSION(:), ALLOCATABLE :: muS, muI, muD, newD
    REAL(wp), DIMENSION(:), ALLOCATABLE :: Y, C, L, B, cfC, cfB
    REAL(wp) :: lifetime, lifetimeNoCOVID
  END TYPE

CONTAINS

  ELEMENTAL FUNCTION u(cons)
    REAL(wp), INTENT(IN) :: cons
    REAL(wp) :: u

    IF (crra == 1.0) THEN
      u = LOG(cons)
    ELSEIF (crra == 2.0) THEN
       u = 1.0_wp - 1.0_wp / cons
    ELSE
      u = (cons**(1.0_wp - crra) - 1.0_wp) / (1.0_wp - crra)
    END IF
  END FUNCTION u

  PURE FUNCTION ce(valFun)
    REAL(wp), INTENT(IN) :: valFun
    REAL(wp) :: ce

    ce = (1.0_wp + (1.0_wp - crra) * (1 - beta) * valFun)**(1.0_wp / (1.0_wp - crra))
  END FUNCTION ce

  PURE SUBROUTINE consDynamics(Lpath, conf, Ypath, Cpath, Bpath, valFromC)
    REAL(wp), DIMENSION(:), INTENT(IN) :: Lpath
    TYPE(Config), INTENT(IN) :: conf
    REAL(wp), DIMENSION(SIZE(Lpath)), INTENT(OUT) :: Ypath, Cpath, Bpath
    REAL(wp), INTENT(OUT) :: valFromC
    REAL(wp) :: wealth
    INTEGER :: ix, sz

    sz = SIZE(Lpath)

    Ypath = (1.0_wp - thetaY * Lpath)**drs

    wealth = -conf%initB + DOT_PRODUCT( conf%discountsRf, Ypath + conf%pushCash ) &
      + conf%discountsRf(sz) * ( 1.0_wp + conf%pushCash(sz) ) / (rfRate - 1.0_wp) / rfRate

    Cpath(1) = wealthShare * wealth 
    Bpath(1) = conf%initB
    DO ix = 2,sz
      Cpath(ix) = Cpath(ix - 1) * shrinkRate
      Bpath(ix) = rfRate * Bpath(ix - 1) - Ypath(ix - 1) + Cpath(ix - 1)
    END DO

    valFromC = wealth**(1.0_wp - crra) / ( (1.0_wp - crra) * wealthShare**crra ) &
      - 1.0_wp / ( (1.0_wp - crra) * (1.0_wp - beta) )
  END SUBROUTINE consDynamics

  PURE SUBROUTINE sirDynamics(Lpath, conf, muS, muI, muD, newD)
    REAL(wp), DIMENSION(:), INTENT(IN) :: Lpath
    TYPE(Config), INTENT(IN) :: conf
    REAL(wp), DIMENSION(SIZE(Lpath)), INTENT(OUT) :: muS, muI, muD, newD

    REAL(wp) :: newI
    INTEGER :: ix, sz
    
    muS = conf%initS
    muI = conf%initI
    muD = 0.0_wp
    newD = 0.0_wp

    sz = SIZE(Lpath)
    DO ix = 1,sz-1
      newI = conf%wedge(ix) * (1.0_wp - theta * Lpath(ix))**2 * muS(ix) * muI(ix)
      IF (ix == conf%freeVaxAt) THEN
        muS(ix+1) = 0.0_wp
      ELSE
        muS(ix+1) = MAX(0.0_wp, muS(ix) - newI)
      END IF
      muI(ix+1) = (1.0_wp - pp) * muI(ix) + newI
      newD(ix) = (piD + piDsq * muI(ix)) * muI(ix)
      IF (ix > 1) THEN
        muD(ix) = muD(ix-1) + newD(ix)
      ELSE
        muD(ix) = newD(ix) + (1.0_wp - conf%initS - conf%initI - conf%initR)
      END IF
    END DO

    muS(sz) = muS(sz-1)
    muI(sz) = muI(sz-1)
    muD(sz) = muD(sz-1)
  END SUBROUTINE sirDynamics

  FUNCTION sirObj(guess, grad, func_data) RESULT(f)
    REAL(wp), DIMENSION(:), INTENT(IN) :: guess
    REAL(wp), DIMENSION(:), INTENT(INOUT), OPTIONAL :: grad
    CLASS(*), INTENT(IN), OPTIONAL :: func_data
    REAL(wp) :: f, valFromC
    TYPE(Config) :: cfg
    INTEGER  :: h, lSz
    REAL(wp), DIMENSION(:), ALLOCATABLE :: Lpath, muS, muI, muD, newD
    REAL(wp), DIMENSION(:), ALLOCATABLE :: Ypath, Cpath, Bpath

    INTEGER, SAVE :: evalCt = 0

    SELECT TYPE(func_data)
    TYPE IS(Config)
      cfg = func_data
    END SELECT

    IF (PRESENT(grad)) THEN
      ! Why?!
    END IF

    lSz = SIZE(guess)
    h = cfg%maxT

    ALLOCATE(Lpath(h), muS(h), muI(h), muD(h), newD(h))
    ALLOCATE(Ypath(h), Cpath(h), Bpath(h))

    Lpath = 0.0_wp
    Lpath(1:lSz) = guess

    CALL sirDynamics(Lpath, cfg, muS, muI, muD, newD)
    CALL consDynamics(Lpath, cfg, Ypath, Cpath, Bpath, valFromC)

    IF (ANY(Cpath <= 0.0_wp)) THEN
      f = -1.0D+8 + MINVAL(Cpath)
    ELSE
      f = valFromC - chi * DOT_PRODUCT(cfg%discountsBeta, newD)
    END IF

    evalCt = evalCt + 1
    IF (MOD(evalCt, 500000) == 0) THEN
      WRITE (OUTPUT_UNIT, "(A)", ADVANCE="NO") "."
      FLUSH(OUTPUT_UNIT)
    END IF
  END FUNCTION sirObj

  SUBROUTINE writeSolution(sol, label)
    TYPE(Solution), INTENT(IN) :: sol
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER :: iunit, ix, sz, ixData
    REAL(wp) :: doStat, dataMuDval, dataLval, bbeta, newI
    REAL(wp), DIMENSION(:), ALLOCATABLE :: R0path

    TYPE(datetime) :: dateObj
    TYPE(timedelta) :: dateStep

    sz = sol%configuration%maxT

    OPEN(newunit=iunit, file=resultsDir // TRIM(label) // "_param.tab")
    WRITE (iunit, "(ES25.6)") sol%configuration%initS
    WRITE (iunit, "(ES25.6)") sol%configuration%initI
    WRITE (iunit, "(ES25.6)") sol%configuration%initR
    WRITE (iunit, "(ES25.6)") pp
    WRITE (iunit, "(ES25.6)") piD0
    WRITE (iunit, "(ES25.6)") piD1
    WRITE (iunit, "(ES25.6)") theta
    WRITE (iunit, "(ES25.6)") thetaY
    WRITE (iunit, "(ES25.6)") rfRate
    WRITE (iunit, "(ES25.6)") beta
    WRITE (iunit, "(ES25.6)") chi
    WRITE (iunit, "(ES25.6)") crra 
    WRITE (iunit, "(ES25.6)") drs
    WRITE (iunit, "(I25)") days 
    WRITE (iunit, "(I25)") sol%configuration%maxT
    WRITE (iunit, "(I25)") sol%configuration%freeVaxAt
    WRITE (iunit, "(I25)") sol%configuration%waveTwoAt
    WRITE (iunit, "(ES25.6)") rho_psi
    WRITE (iunit, "(ES25.6)") psi_wave1_end
    WRITE (iunit, "(ES25.6)") psi_wave2_end

    ! Welfare
    WRITE(iunit, "(ES25.6)") sol%lifetime
    WRITE(iunit, "(ES25.6)") sol%lifetimeNoCOVID

    ! CE
    doStat = ce(sol%lifetime)
    WRITE(iunit, "(ES25.6)") doStat
    doStat = ce(sol%lifetimeNoCOVID)
    WRITE(iunit, "(ES25.6)") doStat

    ! Output loss NPV, % annual output
    doStat = DOT_PRODUCT(sol%configuration%discountsRf, 1.0_wp - sol%Y) / 52.0_wp
    WRITE(iunit, "(ES25.6)") doStat

    ! Total deaths
    doStat = sol%muD(sz)
    WRITE(iunit, "(ES25.6)") doStat

    WRITE (iunit, "(ES25.6)") sol%configuration%lossValL
    WRITE (iunit, "(ES25.6)") sol%configuration%lossValD
    CLOSE(iunit)

    dateObj = datetime(year=eventYear, month=eventMonth, day=eventDay)
    dateStep = timedelta(days=7)

    ALLOCATE(R0path(sz))
    R0path = sol%configuration%wedge * (1.0-theta*sol%L)**2.0*sol%muS / pp

    OPEN(newunit=iunit, file=resultsDir // TRIM(label) // "_sim.tab")
    DO ix = 1,sz
      ixData = ix + wedgeStart - 1
      IF (ixData <= wedgeEnd) THEN
        dataMuDval = dataMuD(ixData)
        dataLval = dataL(ixData)
      ELSE
        dataMuDval = -100_wp
        dataLval = -100_wp
      END IF
      IF (ix == 1) THEN
        newI = sol%muI(1)
      ELSE
        newI = sol%muI(ix) - (1_wp - pp) * sol%muI(ix-1)
      END IF
      bbeta = newI / (sol%muS(ix) * sol%muI(ix))

      WRITE (iunit, "(I6,A,15(ES25.6,A),I25,A)") &
        ix, tab, R0path(ix), tab, sol%configuration%wedge(ix), tab, &
        sol%muS(ix), tab, sol%muI(ix), tab, sol%muD(ix), tab, &
        sol%Y(ix), tab, sol%C(ix), tab, sol%cfC(ix), tab, sol%cfB(ix), tab, &
        sol%L(ix), tab, sol%B(ix), tab, &
        sol%configuration%pushCash(ix), tab, dataMuDval, tab, dataLval, tab, &
        bbeta, tab, dateObj%secondsSinceEpoch(), tab
      dateObj = dateObj + dateStep
    END DO
    CLOSE(iunit)
  END SUBROUTINE writeSolution

  SUBROUTINE dataImpliedWedge()
    REAL(wp), DIMENSION(wedgeSz) :: newDeaths, lockdown, muI, newI, muS, bbeta, dataR
    TYPE(datetime), DIMENSION(wedgeSz) :: calDates
    INTEGER, DIMENSION(3, wedgeSz) :: dataDates
    INTEGER :: iunit, ix
    LOGICAL, PARAMETER :: printData = .FALSE.

    OPEN(newunit=iunit, file=dataDir // "data_L.txt", status="OLD")
    READ (iunit, *) lockdown
    CLOSE(iunit)

    OPEN(newunit=iunit, file=dataDir // "data_D.txt", status="OLD")
    READ (iunit, *) newDeaths
    CLOSE(iunit)

    OPEN(newunit=iunit, file=dataDir // "data_dates.txt", status="OLD")
    READ (iunit, *) dataDates
    CLOSE(iunit)

    muI = (-piD  + SQRT(piD**2 +  4_wp * piDsq * newDeaths) ) / (2_wp * piDsq)
    newI = 0_wp
    DO ix = 1,wedgeSz-1
      newI(ix) = muI(ix+1) - (1_wp - pp) * muI(ix)
    END DO
    
    muS(1) = defaultInitS
    DO ix = 2,wedgeSz
      muS(ix) = muS(ix-1) - newI(ix-1)
    END DO

    bbeta = newI / (muI * muS)
    dataR = bbeta * muS / pp

    dataWedge = bbeta / ((1_wp - theta * lockdown)**2)
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

    IF (printData) THEN
      WRITE (*, "(A25,A,9(A12,A))") "week", TAB, "ix", TAB, "Lockdown", TAB, "New Deaths", TAB, &
        "newI", TAB, "muI", TAB, "muS", TAB, "beta", TAB, "R", TAB, "wedge", TAB
      DO ix = 1,wedgeSz
        WRITE (*, "(A25,A,I12,A,8(ES12.4,A))") calDates(ix)%isoformat(), TAB, ix, TAB, &
          lockdown(ix), TAB, newDeaths(ix), TAB, newI(ix), TAB, muI(ix), TAB, muS(ix), TAB, &
          bbeta(ix), TAB, dataR(ix), TAB, dataWedge(ix), TAB
      END DO
    END IF
  END SUBROUTINE

END MODULE perfectMod
