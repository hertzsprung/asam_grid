MODULE BoundaryCond_Mod

  USE Kind_Mod

  IMPLICIT NONE

  INTEGER, PRIVATE :: InputUnit=10

  TYPE BoundaryCon_T
    CHARACTER(10) :: West
    CHARACTER(10) :: East
    CHARACTER(10) :: South
    CHARACTER(10) :: North
    CHARACTER(10) :: Bottom
    CHARACTER(10) :: Top
  END TYPE BoundaryCon_T
  CHARACTER(10) :: WestC
  CHARACTER(10) :: EastC
  CHARACTER(10) :: SouthC
  CHARACTER(10) :: NorthC
  CHARACTER(10) :: BottomC
  CHARACTER(10) :: TopC

  TYPE (BoundaryCon_T) :: BCVel

  NAMELIST /ModelBCVel/ BCVel

CONTAINS

SUBROUTINE InputModelBCVel(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: ic
  CHARACTER(300) :: Line

  BCVel%West  ='OutFlow'
  BCVel%East  ='OutFlow'
  BCVel%South ='OutFlow'
  BCVel%North ='OutFlow'
  BCVel%Bottom='OutFlow'
  BCVel%Top   ='OutFlow'

! Find line 'BCVel' (first appearance) and modify

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ModelBCVel')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ModelBCVel)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(InputUnit)

END SUBROUTINE InputModelBCVel


END MODULE BoundaryCond_Mod

