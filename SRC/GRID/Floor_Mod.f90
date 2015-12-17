MODULE Floor_Mod

  USE Kind_Mod
  USE Domain_Mod

  IMPLICIT NONE

  INTEGER :: nb
  TYPE (Domain_T), POINTER :: Floor(:)
  INTEGER :: nzLang 

  INTERFACE Allocate
    MODULE PROCEDURE FloorAllocate
  END INTERFACE  
  INTERFACE Deallocate
    MODULE PROCEDURE FloorDeAllocate
  END INTERFACE

CONTAINS

SUBROUTINE FloorAllocate(Floor)

  TYPE (Domain_T), POINTER :: Floor(:)

  INTEGER :: ib

  DO ib=1,nb
    CALL Allocate(Floor(ib))
  END DO

END SUBROUTINE FloorAllocate


SUBROUTINE FloorDeAllocate(Floor)

   TYPE (Domain_T), POINTER :: Floor(:)

   INTEGER :: ib

   DO ib=1,nb
     CALL DeAllocate(Floor(ib))
   END DO

END SUBROUTINE FloorDeAllocate


END MODULE Floor_Mod


