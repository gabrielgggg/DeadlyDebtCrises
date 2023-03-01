MODULE MPIMod
  USE iso_Fortran_env, ONLY: wp => real64, sip => int8
  USE mpi
  USE Param
  IMPLICIT NONE

  INTEGER :: workerId, workerNo, chunk
  INTEGER :: mpiErr

  INTEGER, DIMENSION(:, :), ALLOCATABLE :: longToWide
  INTEGER, DIMENSION(:), ALLOCATABLE :: wideToLong

  REAL(wp), DIMENSION(:), ALLOCATABLE :: buffq
  REAL(wp), DIMENSION(:), ALLOCATABLE :: buffv, buffBIG

  REAL(wp), DIMENSION(:), ALLOCATABLE :: buffd
  REAL(wp), DIMENSION(:), ALLOCATABLE :: buffc
  REAL(wp), DIMENSION(:), ALLOCATABLE :: buffb, buffHUGE
  INTEGER(sip), DIMENSION(:), ALLOCATABLE :: buffl, buffHUGEi

CONTAINS

  SUBROUTINE mpiSetup(okFlag)
    INTEGER, INTENT(OUT) :: okFlag

    CALL MPI_Init(mpiErr)
    CALL MPI_Comm_size(MPI_COMM_WORLD, workerNo, mpiErr)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, workerId, mpiErr)

    chunk = spaceSz / workerNo

    IF (MOD(spaceSz, workerNo) == 0) THEN
      okFlag = 1

      CALL translation()

      ALLOCATE(buffq(chunk))
      ALLOCATE(buffv(chunk))
      ALLOCATE(buffc(chunk))

      ALLOCATE(buffBIG(spaceSz))

      ALLOCATE(buffl(chunk * bSz))
      ALLOCATE(buffd(chunk * bSz))
      ALLOCATE(buffb(chunk * bSz))

      ALLOCATE(buffHUGE(spaceSz * bSz))
      ALLOCATE(buffHUGEi(spaceSz * bSz))
    ELSE
      okFlag = 0
    END IF

    IF (workerId == 0) THEN
      WRITE (*, *) "workerNo = ", workerNo, ", okFlag = ", okFlag
    END IF
  END SUBROUTINE mpiSetup

  SUBROUTINE translation()
    INTEGER :: bIx, ix

    ALLOCATE( longToWide(spaceSz, 1) )
    ALLOCATE( wideToLong(bSz) )

    ix = 1
    DO bIx = 1,bSz
      longToWide(ix, :) = [ bIx ]
      wideToLong(bIx) = ix
      ix = ix + 1
    END DO
  END SUBROUTINE translation

  SUBROUTINE mpiFinalize()
    CALL MPI_Finalize(mpiErr)
  END SUBROUTINE mpiFinalize

  SUBROUTINE mpiBarrier()
    CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
  END SUBROUTINE

END MODULE MPIMod
