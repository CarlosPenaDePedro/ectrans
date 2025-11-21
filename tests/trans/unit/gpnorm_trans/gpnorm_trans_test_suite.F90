! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE GPNORM_TRANS_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB, JPRD
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_MYRANK, MPL_ALLREDUCE, MPL_END

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "dist_grid.h"
#include "gpnorm_trans.h"
#include "trans_end.h"

! Spectral truncation used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: TRUNCATION = 79

! Number of latitudes used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: NDGL = 2 * (TRUNCATION + 1)

! Tolerance for "close to zero"
REAL(KIND=JPRB), PARAMETER :: TOLERANCE = 100.0_JPRB * EPSILON(1.0_JPRB)

! NPROMA blocking factor
INTEGER(KIND=JPIM), PARAMETER :: NPROMA = 16

CONTAINS

!---------------------------------------------------------------------------------------------------

! Approximate equality check for reals
ELEMENTAL LOGICAL FUNCTION APPROX_EQ(A, B, TOL) RESULT(RET)
  REAL(KIND=JPRB), INTENT(IN) :: A
  REAL(KIND=JPRB), INTENT(IN) :: B
  REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: TOL

  IF (PRESENT(TOL)) THEN
    RET = ABS(A - B) < TOL
  ELSE
    RET = ABS(A - B) < TOLERANCE
  END IF
END FUNCTION APPROX_EQ

! Initialise global field with all ones and distribute it
FUNCTION GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, NFIELDS) RESULT(PGP)
  INTEGER(KIND=JPIM), INTENT(IN) :: MY_PROC
  INTEGER(KIND=JPIM), INTENT(IN) :: NGPTOTG
  INTEGER(KIND=JPIM), INTENT(IN) :: NGPBLKS
  INTEGER(KIND=JPIM), INTENT(IN) :: NFIELDS

  REAL(KIND=JPRB), ALLOCATABLE :: ZGPG(:,:), PGP(:,:,:)
  INTEGER(KIND=JPIM) :: I

  ! Initialise global field
  IF (MY_PROC == 1) THEN
    ALLOCATE(ZGPG(NGPTOTG,NFIELDS))

    ! Set one element to one, all others to zero
    ZGPG(1,1) = 1.0_JPRB
    ZGPG(2:,2:) = 0.0_JPRB
  ENDIF

  ! Initialise distributed fields
  ALLOCATE(PGP(NPROMA,NFIELDS,NGPBLKS))

  ! Distribute field from first task to other tasks
  IF (MY_PROC == 1) THEN
    CALL DIST_GRID(PGPG=ZGPG, KFDISTG=NFIELDS, KFROM=(/(1, I = 1,NFIELDS)/), PGP=PGP, KPROMA=NPROMA)
  ELSE
    CALL DIST_GRID(KFDISTG=NFIELDS, KFROM=(/ (1, I = 1, NFIELDS) /), PGP=PGP, KPROMA=NPROMA)
  ENDIF

  IF (MY_PROC == 1) DEALLOCATE(ZGPG)
END FUNCTION GET_INPUT_FIELD

! Setup fixture
SUBROUTINE SETUP_TEST(KSPEC2, KGPTOTG, KGPTOT, KGPBLKS, LUSE_MPI, MY_PROC)
  USE UTIL, ONLY: DETECT_MPIRUN

  INTEGER(KIND=JPIM), INTENT(OUT) :: KSPEC2
  INTEGER(KIND=JPIM), INTENT(OUT) :: KGPTOTG
  INTEGER(KIND=JPIM), INTENT(OUT) :: KGPTOT
  INTEGER(KIND=JPIM), INTENT(OUT) :: KGPBLKS
  LOGICAL, INTENT(OUT) :: LUSE_MPI
  INTEGER(KIND=JPIM), INTENT(OUT) :: MY_PROC

  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I
  INTEGER(KIND=JPIM) :: NPROC

  ! Set up MPI
  LUSE_MPI = DETECT_MPIRUN()
  IF (LUSE_MPI) THEN
    CALL MPL_INIT
    NPROC = MPL_NPROC()
    MY_PROC = MPL_MYRANK()
  ELSE
    NPROC = 1
    MY_PROC = 1
  ENDIF

  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC, KPRGPEW=1, KPRTRW=NPROC)

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN)

  CALL TRANS_INQ(KSPEC2=KSPEC2, KGPTOTG=KGPTOTG, KGPTOT=KGPTOT)

  ! Number of NPROMA blocks
  KGPBLKS = (KGPTOT - 1) / NPROMA + 1
END SUBROUTINE SETUP_TEST

! Tear down fixture
SUBROUTINE CLEANUP_TEST(LUSE_MPI)
  LOGICAL, INTENT(IN) :: LUSE_MPI

  CALL TRANS_END

  IF (LUSE_MPI) THEN
    CALL MPL_END
  ENDIF
END SUBROUTINE CLEANUP_TEST

!---------------------------------------------------------------------------------------------------

! Test GPNORM_TRANS minimum functionality
INTEGER FUNCTION UNIT_TEST_GPNORM_TRANS_MIN() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:)
  REAL(KIND=JPRB) :: ZAVE(1), ZMIN(1), ZMAX(1)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI

  ! Set up everything
  ! ZGP will contain a single one for one task and zeros elsewhere
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)
  ZGP = GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, 1)

  ! Make the task a little harder - negate ZGP so there's one negative one and zeros elsewhere
  ZGP = -ZGP

  ! Calculate average, minimum, and maximum of field (LDAVE_ONLY = .FALSE.)
  ! Note: not possible to compute only one of average, min, max - all three must be computed
  ! (unless LDAVE_ONLY = .TRUE. in which case only average is computed)
  CALL GPNORM_TRANS(ZGP, 1, NPROMA, ZAVE, ZMIN, ZMAX, .FALSE.)

  ! If I am task 1, check that the minimum is negative one
  IF (MY_PROC == 1) THEN
    RET = MERGE(0, 1, ZMIN(1) == -1.0_JPRB)
  ELSE
    RET = 0
  ENDIF

  ! Communicate results to other tasks
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZGP)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_GPNORM_TRANS_MIN

!---------------------------------------------------------------------------------------------------

! Test GPNORM_TRANS maximum functionality
INTEGER FUNCTION UNIT_TEST_GPNORM_TRANS_MAX() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:)
  REAL(KIND=JPRB) :: ZAVE(1), ZMIN(1), ZMAX(1)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI

  ! Set up everything
  ! ZGP will contain a single one for one task and zeros elsewhere
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)
  ZGP = GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, 1)

  ! Calculate average, minimum, and maximum of field (LDAVE_ONLY = .FALSE.)
  ! Note: not possible to compute only one of average, min, max - all three must be computed
  ! (unless LDAVE_ONLY = .TRUE. in which case only average is computed)
  CALL GPNORM_TRANS(ZGP, 1, NPROMA, ZAVE, ZMIN, ZMAX, .FALSE.)

  ! If I am task 1, check that the maximum is one
  IF (MY_PROC == 1) THEN
    RET = MERGE(0, 1, ZMAX(1) == 1.0_JPRB)
  ELSE
    RET = 0
  ENDIF

  ! Communicate results to other tasks
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZGP)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_GPNORM_TRANS_MAX

!---------------------------------------------------------------------------------------------------

! Test GPNORM_TRANS average functionality
INTEGER FUNCTION UNIT_TEST_GPNORM_TRANS_AVE() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:)
  REAL(KIND=JPRB) :: ZAVE(1), ZMIN(1), ZMAX(1)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI

  ! Set up everything
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)

  ! Initialise distributed fields
  ALLOCATE(ZGP(NPROMA,1,NGPBLKS))

  ! Set all values to one
  ZGP(:,:,:) = 1.0_JPRB

  ! Calculate average, minimum, and maximum of field (LDAVE_ONLY = .FALSE.)
  ! Note: not possible to compute only one of average, min, max - all three must be computed
  ! (unless LDAVE_ONLY = .TRUE. in which case only average is computed)
  CALL GPNORM_TRANS(ZGP, 1, NPROMA, ZAVE, ZMIN, ZMAX, .FALSE.)

  ! If I am task 1, check that the average is one
  IF (MY_PROC == 1) THEN
    RET = MERGE(0, 1, APPROX_EQ(ZAVE(1), 1.0_JPRB))
  ELSE
    RET = 0
  ENDIF

  ! Communicate results to other tasks
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZGP)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_GPNORM_TRANS_AVE


INTEGER FUNCTION UNIT_TEST_GPNORM_TRANS_COMPARE() RESULT(RET) BIND(C)

  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:), ZGPG(:,:)
  REAL(KIND=JPRB) :: AVE_A(1), MIN_A(1), MAX_A(1)
  REAL(KIND=JPRB) :: AVE_B(1), MIN_B(1), MAX_B(1)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I, J, K
  INTEGER :: FAILS
  REAL(KIND=JPRB), PARAMETER :: REL_TOL = 1.0e-6_JPRB
  REAL(KIND=JPRB), PARAMETER :: ABS_TOL = 1.0e2_JPRB * EPSILON(1.0_JPRB)

  ! ----------------- Test body -----------------
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)

  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  FAILS = 0

  ! Case 1: delta
  IF (MY_PROC == 1) THEN
    ALLOCATE(ZGPG(NGPTOTG,1)); ZGPG = 0.0_JPRB; ZGPG(1,1) = 1.0_JPRB
  END IF
  ALLOCATE(ZGP(NPROMA,1,NGPBLKS))
  IF (MY_PROC == 1) THEN
    CALL DIST_GRID(PGPG=ZGPG, KFDISTG=1, KFROM=(/1/), PGP=ZGP, KPROMA=NPROMA)
  ELSE
    CALL DIST_GRID(KFDISTG=1, KFROM=(/1/), PGP=ZGP, KPROMA=NPROMA)
  END IF
  IF (MY_PROC == 1) DEALLOCATE(ZGPG)
  CALL DO_COMPARE('delta ave/min/max', LDAVE=.FALSE.)
  DEALLOCATE(ZGP)

  ! Case 2: ones
  ALLOCATE(ZGP(NPROMA,1,NGPBLKS)); ZGP = 1.0_JPRB
  CALL DO_COMPARE('ones ave/min/max', LDAVE=.FALSE.)
  DEALLOCATE(ZGP)

  ! Case 3: PW-sensitive row-constant
  IF (MY_PROC == 1) THEN
    ALLOCATE(ZGPG(NGPTOTG,1)); ZGPG = 0.0_JPRB
    K = 1
    DO I = 1, NDGL
      DO J = 1, ILOEN(I)
        ZGPG(K,1) = REAL(I, JPRB)
        K = K + 1
      END DO
    END DO
  END IF
  ALLOCATE(ZGP(NPROMA,1,NGPBLKS))
  IF (MY_PROC == 1) THEN
    CALL DIST_GRID(PGPG=ZGPG, KFDISTG=1, KFROM=(/1/), PGP=ZGP, KPROMA=NPROMA)
  ELSE
    CALL DIST_GRID(KFDISTG=1, KFROM=(/1/), PGP=ZGP, KPROMA=NPROMA)
  END IF
  IF (MY_PROC == 1) DEALLOCATE(ZGPG)

  ! 3a) Average-only path (exercise ave-only branch)
  CALL DO_COMPARE('row-const ave only', LDAVE=.TRUE.)
  ! 3b) Full path (also checks min/max)
  CALL DO_COMPARE('row-const ave+minmax', LDAVE=.FALSE.)
  DEALLOCATE(ZGP)

  CALL CLEANUP_TEST(LUSE_MPI)
  RET = MERGE(0, 1, FAILS == 0)
CONTAINS
  PURE LOGICAL FUNCTION OK_EQ(X, Y) RESULT(OK)
    REAL(KIND=JPRB), INTENT(IN) :: X, Y
    REAL(KIND=JPRB) :: REL
    REL = ABS(X - Y) / MAX(1.0_JPRB, ABS(X), ABS(Y))
    OK = (ABS(X - Y) <= ABS_TOL) .OR. (REL <= REL_TOL)
  END FUNCTION OK_EQ

  SUBROUTINE DO_COMPARE(TAG, LDAVE)
    CHARACTER(*), INTENT(IN) :: TAG
    LOGICAL,      INTENT(IN) :: LDAVE
    INTEGER :: LOC_FAIL
    REAL(KIND=JPRB) :: DAVE, DMIN, DMAX, RELAVE, RELMIN, RELMAX

    LOC_FAIL = 0

    CALL GPNORM_TRANS(ZGP, 1, NPROMA, AVE_A, MIN_A, MAX_A, LDAVE, LREPRO=.FALSE.)
    CALL GPNORM_TRANS(ZGP, 1, NPROMA, AVE_B, MIN_B, MAX_B, LDAVE, LREPRO=.TRUE.)

    IF (MY_PROC == 1) THEN
    DAVE = AVE_A(1) - AVE_B(1)
    RELAVE = ABS(DAVE) / MAX(1.0_JPRB, ABS(AVE_A(1)), ABS(AVE_B(1)))

    IF (.NOT. OK_EQ(AVE_A(1), AVE_B(1))) LOC_FAIL = 1
      IF (.NOT. LDAVE) THEN
        DMIN = MIN_A(1) - MIN_B(1)
        DMAX = MAX_A(1) - MAX_B(1)
        RELMIN = ABS(DMIN) / MAX(1.0_JPRB, ABS(MIN_A(1)), ABS(MIN_B(1)))
        RELMAX = ABS(DMAX) / MAX(1.0_JPRB, ABS(MAX_A(1)), ABS(MAX_B(1)))
        IF (.NOT. OK_EQ(MIN_A(1), MIN_B(1))) LOC_FAIL = 1
        IF (.NOT. OK_EQ(MAX_A(1), MAX_B(1))) LOC_FAIL = 1
      ELSE
        DMIN = 0.0_JPRB; DMAX = 0.0_JPRB
        RELMIN = 0.0_JPRB; RELMAX = 0.0_JPRB
      END IF

        WRITE(*,'(A)') ''
        WRITE(*,'(A)') '--- ' // TRIM(TAG) // ' ---'
        WRITE(*,'(A,L1)') ' LDAVE_ONLY = ', LDAVE
        WRITE(*,'(A,1PE16.8)') ' A_avg : ', AVE_A(1)
        WRITE(*,'(A,1PE16.8)') ' B_avg : ', AVE_B(1)
        WRITE(*,'(A,1PE16.8,2X,A,1PE10.3)') ' Delta avg : ', DAVE, ' rel=', RELAVE

        IF (.NOT. LDAVE) THEN
          WRITE(*,'(A,1PE16.8)') ' A_min : ', MIN_A(1)
          WRITE(*,'(A,1PE16.8)') ' B_min : ', MIN_B(1)
          WRITE(*,'(A,1PE16.8,2X,A,1PE10.3)') ' Delta min : ', DMIN, ' rel=', RELMIN
          WRITE(*,'(A,1PE16.8)') ' A_max : ', MAX_A(1)
          WRITE(*,'(A,1PE16.8)') ' B_max : ', MAX_B(1)
          WRITE(*,'(A,1PE16.8,2X,A,1PE10.3)') ' Delta max : ', DMAX, ' rel=', RELMAX
        END IF

        IF (LOC_FAIL /= 0) THEN
          WRITE(*,'(A)') ' WARNING: difference exceeds tolerance.'
        ELSE
          WRITE(*,'(A)') ' OK: results within tolerance.'
        END IF
    END IF

    IF (LUSE_MPI) CALL MPL_ALLREDUCE(LOC_FAIL, CDOPER="MAX")
    FAILS = FAILS + LOC_FAIL
  END SUBROUTINE DO_COMPARE
END FUNCTION UNIT_TEST_GPNORM_TRANS_COMPARE

!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------

END MODULE GPNORM_TRANS_TEST_SUITE
