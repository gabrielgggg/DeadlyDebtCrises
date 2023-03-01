MODULE sim
  USE iso_Fortran_env, ONLY: wp => real64
  IMPLICIT NONE

  REAL(wp), PARAMETER :: piwp = 3.14159265358979323846_wp
CONTAINS

  SUBROUTINE fixSeed()
    INTEGER, ALLOCATABLE, DIMENSION(:) :: seed
    INTEGER :: n
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    seed = 1989
    CALL RANDOM_SEED(put = seed)
  END SUBROUTINE fixSeed

  FUNCTION simDiscrete(pmf) RESULT(val)
    REAL(wp), DIMENSION(:), INTENT(IN) :: pmf
    INTEGER :: val
    INTEGER :: sz, ix
    REAL(wp) :: draw, theSum

    sz = SIZE(pmf)
    CALL simUniform(draw)
    theSum = 0_wp
    val = -1
    DO ix = 1,sz
      theSum = theSum + pmf(ix)
      IF (draw <= theSum) THEN
        val = ix
        EXIT
      END IF
    END DO
  END FUNCTION simDiscrete

  SUBROUTINE simMarkov(piMat, newIx, oldIx)
    INTEGER, INTENT(IN) :: oldIx
    REAL(wp), DIMENSION(:, :), INTENT(IN) :: piMat
    INTEGER, INTENT(OUT) :: newIx

    newIx = simDiscrete( piMat(oldIx, :) )
  END SUBROUTINE simMarkov

  SUBROUTINE simStdNormal(val)
    REAL(wp), INTENT(OUT) :: val
    REAL(wp) :: v1, v2

    CALL RANDOM_NUMBER(v1)
    CALL RANDOM_NUMBER(v2)
    val = SQRT(-2_wp * LOG(v1)) * COS(2_wp * piwp * v2)
  END SUBROUTINE simStdNormal

  SUBROUTINE simUniform(val)
    REAL(wp), INTENT(OUT) :: val

    CALL RANDOM_NUMBER(val)
  END SUBROUTINE simUniform

  SUBROUTINE simDiscreteUniform(sz, retIx)
    INTEGER, INTENT(IN) :: sz
    INTEGER, INTENT(OUT) :: retIx
    REAL(wp) :: rno

    CALL RANDOM_NUMBER(rno)
    retIx = FLOOR(sz * rno) + 1
  END SUBROUTINE simDiscreteUniform

END MODULE sim

