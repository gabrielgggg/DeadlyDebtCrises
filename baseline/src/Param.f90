MODULE Param
  USE iso_Fortran_env, ONLY: wp => real64, sip => int8, OUTPUT_UNIT
  IMPLICIT NONE

  ! I/O
  CHARACTER, PARAMETER:: TAB = CHAR(9)
  CHARACTER (LEN=*), PARAMETER :: ioDir = "./"
  CHARACTER (LEN=*), PARAMETER :: outDir = ioDir // "results/"
  CHARACTER (LEN=*), PARAMETER :: dataDir = ioDir // "../data/"
  CHARACTER (LEN=*), PARAMETER :: simFileName = "baseline"

  ! Configuration
  LOGICAL, PARAMETER :: loadResults = .FALSE.
  LOGICAL, PARAMETER :: runBack = .TRUE.

  LOGICAL, PARAMETER :: useResetB = .FALSE.
  LOGICAL, PARAMETER :: exoLoan = .FALSE.

  LOGICAL, PARAMETER :: surpriseWave2 = .FALSE.

  ! Set from data, upon execution
  INTEGER :: eventYear = 2020, eventMonth = 4, eventDay = 1

  INTEGER, PARAMETER :: wedgeSz = 97, wedgeStart = 8, wedgeEnd = 91, dataWave2 = 45
  INTEGER, PARAMETER :: startWave1 = wedgeStart, endWave1 = dataWave2-1, wave1Sz = endWave1 - startWave1 + 1
  INTEGER, PARAMETER :: startWave2 = dataWave2, endWave2 = wedgeEnd, wave2Sz = endWave2 - startWave2 + 1
  
  ! Parameterization
  REAL(wp), PARAMETER :: crra = 2.0_wp
  REAL(wp), PARAMETER :: beta = 0.98_wp**(1.0_wp / 52.0_wp)
  REAL(wp), PARAMETER :: chi = 3500.0_wp

  REAL(wp), PARAMETER :: rf = 1.01_wp**(1.0_wp / 52.0_wp)
  REAL(wp), PARAMETER :: duration = 5.0_wp * 52.0_wp
  REAL(wp), PARAMETER :: delta = rf / duration - rf + 1.0_wp
  REAL(wp), PARAMETER :: pay = delta + rf - 1.0_wp
  REAL(wp), PARAMETER :: kappa = 0.54_wp

  REAL(wp), PARAMETER :: gamm0 = 0.04_wp
  REAL(wp), PARAMETER :: gamm1 = 1.62_wp
  REAL(wp), PARAMETER :: gamm2 = 1.7825D-2

  REAL(wp), PARAMETER :: meanZ = 1.0_wp
  REAL(wp), PARAMETER :: drs = 0.67_wp

  INTEGER, PARAMETER :: H = 3 * 52, before = 10 * 52, after = before - H
  INTEGER, PARAMETER :: Tsz = H  + before + after

  INTEGER, PARAMETER :: startRepay = 3 * 52
  REAL(wp), PARAMETER :: exoLoanSz = 0.0_wp * 52_wp
  REAL(wp), PARAMETER :: exoLoanPay = (rf - 1.0_wp) * exoLoanSz * rf**(startRepay)

  ! Infection
  REAL(wp), PARAMETER :: days = 6.0_wp
  REAL(wp), PARAMETER :: pp = 1.0_wp - (1.0_wp - 1.0_wp / days)**7.0_wp
  REAL(wp), PARAMETER :: piD = 0.0085_wp * pp
  REAL(wp), PARAMETER :: piDsq = 0.08_wp * pp
  REAL(wp), PARAMETER :: theta = 0.5_wp
  REAL(wp), PARAMETER :: thetaY = 0.8_wp

  REAL(wp), PARAMETER :: rho_psi = 0.99_wp
  REAL(wp), PARAMETER :: psi_wave1_end = 1.35_wp
  REAL(wp), PARAMETER :: psi_wave2_end = psi_wave1_end

  REAL(wp), PARAMETER :: minL = 0.0_wp, maxL = 0.8_wp
  REAL(wp) :: defaultInitI = 0.001_wp
  REAL(wp) :: defaultInitS = 0.965_wp
  REAL(wp) :: defaultInitR = 1.0_wp - 0.001_wp - 0.965_wp

  ! State space size
  INTEGER, PARAMETER :: muSsz = 80
  INTEGER, PARAMETER :: muIsz = 81
  INTEGER, PARAMETER :: bSz = 350
  INTEGER, PARAMETER :: Lsz = 123
  INTEGER, PARAMETER :: spaceSz = muSsz * muIsz * bSz
  
  REAL(wp), PARAMETER :: minB = 0.425_wp * 52_wp
  REAL(wp), PARAMETER :: maxB = 0.825_wp * 52_wp
  REAL(wp), PARAMETER :: initB = 0.612966_wp * 52.0_wp
  REAL(wp), PARAMETER :: resetB = initB
 
  ! Numerical
  REAL(wp), PARAMETER :: eulerMascheroni = 0.577215664901532_wp
  REAL(wp), PARAMETER :: tasteB = 1.0D-4
  REAL(wp), PARAMETER :: veryNeg = -99999_wp
  REAL(wp), PARAMETER :: epsTol = 1.0D-6

  INTEGER, PARAMETER :: maxIter = 50000
END MODULE Param
