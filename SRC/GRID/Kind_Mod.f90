MODULE Kind_Mod

  INTEGER, PARAMETER :: RealKind=8
  INTEGER, PARAMETER :: Real4Kind=4
  INTEGER, PARAMETER :: IntKind=4

  ! 1. Kind-Parameters for the Program:
! -----------------------------------

  INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers

CONTAINS


  SUBROUTINE native_4byte_real( realIn, realOut )
  IMPLICIT NONE
  REAL(4), INTENT(IN)                              :: realIn
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element
  REAL(4), INTENT(OUT)                             :: realOut
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element, with
                                                   ! reverse byte order to
                                                   ! that of realIn
  realOut=realIn
  END SUBROUTINE native_4byte_real

END MODULE Kind_Mod
