!
! Perfect fin markets case
! Arellano Bai Mihalache 2023
!
PROGRAM perfect
  USE iso_Fortran_env, ONLY: wp => real64, OUTPUT_UNIT
  USE perfectMod
  USE NL
  IMPLICIT NONE

  INTEGER :: sz, ix, mitTime, mitSz
  TYPE(Config) :: baseConf, mitConf
  TYPE(Solution) :: baseSol, mitSol

  CALL dataImpliedWedge()

  !
  ! Baseline
  !
  WRITE (*, *) "Baseline..."
  sz = 6 * 52

  baseConf%initI = defaultInitI
  baseConf%initR = defaultInitR
  baseConf%initS = defaultInitS

  baseConf%initB = defaultInitB

  baseConf%maxT = sz
  baseConf%freeVaxAt = sz - 3 * 52
  IF (surpriseWave2) THEN
    baseConf%waveTwoAt = dataWave2 - wedgeStart + 1
  ELSE
    baseConf%waveTwoAt = sz-1
  END IF

  ALLOCATE(baseConf%discountsBeta(sz))
  ALLOCATE(baseConf%discountsRf(sz))
  DO ix = 0,sz-1
    baseConf%discountsBeta(ix+1) = beta**ix
    baseConf%discountsRf(ix+1) = (1.0_wp / rfRate)**ix
  END DO

  ALLOCATE(baseConf%wedge(sz))
  IF (surpriseWave2) THEN
    baseConf%wedge(1:wave1Sz) = dataWedge(startWave1:endWave1)
    DO ix = wave1Sz+1,sz
      baseConf%wedge(ix) = rho_psi * baseConf%wedge(ix-1) + (1.0-rho_psi) * psi_wave1_end * pp
    END DO
  ELSE 
    baseConf%wedge(1:(wave1Sz+wave2Sz)) = dataWedge(startWave1:endWave2)
    DO ix = wave1Sz+wave2Sz+1,sz
      baseConf%wedge(ix) = rho_psi * baseConf%wedge(ix-1) + (1.0-rho_psi) * psi_wave1_end * pp
    END DO
  END IF

  ALLOCATE(baseConf%pushCash(sz))
  baseConf%cashProgramSize = 0.0_wp
  baseConf%pushCash = 0.0_wp

  baseSol = computeSolution(baseConf)
  IF (surpriseWave2) THEN
    CALL writeSolution(baseSol, "perfect_wave1")
  ELSE
    ! CALL writeSolution(baseSol, "overallSol")
  END IF

  IF (surpriseWave2) THEN
    !
    ! Second wave
    !
    WRITE (*, *) "Second wave..."
    mitTime = baseConf%waveTwoAt
  
    mitSz = sz - mitTime + 1
  
    mitConf%initI = baseSol%muI(mitTime)
    mitConf%initR = 1.0_wp - baseSol%muI(mitTime) - baseSol%muS(mitTime) - baseSol%muD(mitTime)
    mitConf%initS = baseSol%muS(mitTime)

    mitConf%initB = baseSol%B(mitTime)
  
    mitConf%maxT = mitSz
    mitConf%freeVaxAt = baseConf%freeVaxAt - mitTime + 1
    mitConf%waveTwoAt = 1
  
    ALLOCATE(mitConf%discountsBeta(mitSz))
    ALLOCATE(mitConf%discountsRf(mitSz))
    DO ix = 0,mitSz-1
      mitConf%discountsBeta(ix+1) = beta**ix
      mitConf%discountsRf(ix+1) = (1.0_wp / rfRate)**ix
    END DO
  
    ALLOCATE(mitConf%wedge(mitSz))
    mitConf%wedge(1:wave2Sz) = dataWedge(startWave2:endWave2)
    DO ix = wave2Sz+1,mitSz
      mitConf%wedge(ix) = rho_psi * mitConf%wedge(ix-1) + (1.0-rho_psi) * psi_wave2_end * pp
    END DO
  
    ALLOCATE(mitConf%pushCash(mitSz))
    mitConf%cashProgramSize = 0.0_wp
    mitConf%pushCash = 0.0_wp
  
    mitSol = computeSolution(mitConf)
    CALL writeSolution(mitSol, "perfect_wave2")
  
    !
    ! Pasted
    !
    WRITE (*, *) "Pasting..."
    baseSol%configuration%wedge(mitTime:sz) = mitSol%configuration%wedge
    baseSol%configuration%pushCash(mitTime:sz) = mitSol%configuration%pushCash
    baseSol%muS(mitTime:sz) = mitSol%muS
    baseSol%muI(mitTime:sz) = mitSol%muI
    baseSol%muD(mitTime:sz) = mitSol%muD
    baseSol%newD(mitTime:sz) = mitSol%newD
    baseSol%Y(mitTime:sz) = mitSol%Y
    baseSol%C(mitTime:sz) = mitSol%C
    baseSol%L(mitTime:sz) = mitSol%L
    baseSol%B(mitTime:sz) = mitSol%B

    baseSol%lifetime = beta**mitTime * mitSol%lifetime + DOT_PRODUCT( baseSol%configuration%discountsBeta(1:mitTime-1), &
      u(baseSol%C(1:mitTime-1)) - chi * baseSol%newD(1:mitTime-1) )
  END IF ! surpriseWave2

  baseSol%configuration%lossValL = SQRT(SUM( (baseSol%L(1:(wave1Sz+wave2Sz)) - dataL(wedgeStart:wedgeEnd))**2 )) &
                    /(1.0+SQRT(SUM((dataL(wedgeStart:wedgeEnd)**2))))
  baseSol%configuration%lossValD = SQRT(SUM( (baseSol%muD(1:(wave1Sz+wave2Sz)) - dataMuD(wedgeStart:wedgeEnd))**2 )) &
                    /(1.0+SQRT(SUM((dataMuD(wedgeStart:wedgeEnd)**2))))

  CALL writeSolution(baseSol, "perfect")
  
  WRITE (*, *) ""
  WRITE (*, "(A,ES15.4,A,A,ES15.4)") "Loss L: ", &
    baseSol%configuration%lossValL, tab, &
    " Loss muD: ", baseSol%configuration%lossValD 
  WRITE (*, *) ""
  WRITE (*, *) "The End."

CONTAINS
  !
  !
  ! Contains...
  !
  !

  FUNCTION computeSolution(conf) RESULT(sol)
    TYPE(Config), INTENT(IN), TARGET :: conf
    TYPE(Solution) :: sol
    INTEGER :: h
    INTEGER :: lSz
    REAL(wp), DIMENSION(:), ALLOCATABLE :: lGuess, iniStep, Lzeros
    REAL(wp), DIMENSION(:), ALLOCATABLE :: lLB
    REAL(wp), DIMENSION(:), ALLOCATABLE :: lUB

    TYPE(nlopt_opt) :: opt
    INTEGER :: stat
    REAL(wp) :: fmax, valFromC

    h = conf%maxT

    lSz = conf%freeVaxAt - 1
    ALLOCATE(lGuess(lSz), lLB(lSz), lUB(lSz))
    lUB = maxL
    lLB = minL
    lGuess = 0_wp
    lGuess(1:52) = 0.5_wp
    ! lGuess(53:104) = 0.25_wp

    ALLOCATE(iniStep(lSz))
    iniStep = 0.1_wp

    CALL create(opt, algorithm_from_string("LN_SBPLX"), lSz)
    CALL opt%set_lower_bounds(lLB)
    CALL opt%set_upper_bounds(lUB)
    CALL opt%set_initial_step(iniStep)
    CALL opt%set_xtol_rel(xTol)
    CALL opt%set_ftol_rel(xTol / 10.0_wp)
    CALL opt%set_max_objective(nlopt_func(sirObj, conf))
    CALL opt%optimize(lGuess, fmax, stat)
    CALL destroy(opt)

    IF (stat < NLOPT_SUCCESS) THEN
      WRITE (*, "(A,I5)") "NLopt failed with code ", stat
    ELSE
      WRITE (*, "(A,I5)") "NLopt success with code ", stat
    END IF

    sol%configuration = conf
    ALLOCATE(sol%muS(h), sol%muI(h), sol%muD(h), sol%newD(h))
    ALLOCATE(sol%Y(h), sol%C(h), sol%B(h), sol%L(h), sol%cfC(h), sol%cfB(h))

    ALLOCATE(Lzeros(h))
    Lzeros = 0.0_wp

    sol%L = 0.0_wp
    sol%L(1:lSz) = lGuess

    CALL sirDynamics(sol%L, conf, sol%muS, sol%muI, sol%muD, sol%newD)
    CALL consDynamics(Lzeros, conf, sol%Y, sol%cfC, sol%cfB, valFromC)
    sol%lifetimeNoCOVID = valFromC

    CALL consDynamics(sol%L, conf, sol%Y, sol%C, sol%B, valFromC)
    sol%lifetime = valFromC - chi * DOT_PRODUCT( conf%discountsBeta, sol%newD)
  END FUNCTION computeSolution

END PROGRAM perfect
