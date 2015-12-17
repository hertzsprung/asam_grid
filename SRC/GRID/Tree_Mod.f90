MODULE Tree_Mod

  USE Geometry_Mod
  USE Haus_Mod
  IMPLICIT NONE

  TYPE Box_T
    CHARACTER :: TypeCut
    TYPE(Point_T) :: P0,P1
    REAL(8) :: Height
    TYPE(Box_T), POINTER :: Parent=>NULL()
    TYPE(Box_T), POINTER :: Child1=>NULL()
    TYPE(Box_T), POINTER :: Child2=>NULL()
    INTEGER, POINTER :: List(:)
    INTEGER :: Index
  END TYPE Box_T

CONTAINS 

RECURSIVE SUBROUTINE CreateTree(Root,Depth)

  TYPE(Box_T), POINTER :: Root
  INTEGER :: Depth

  IF (Depth>0) THEN
    ALLOCATE(Root%Child1)   
    Root%Child1%Parent=>Root
    NULLIFY(Root%Child1%List)
    ALLOCATE(Root%Child2)   
    Root%Child2%Parent=>Root
    NULLIFY(Root%Child2%List)
    IF (Root%TypeCut=='x') THEN
      Root%Child1%TypeCut='y'
      Root%Child1%P0=Root%P0
      Root%Child1%P1=Root%P1
      Root%Child1%P1%x=0.5d0*(Root%P0%x+Root%P1%x)
      Root%Child1%Index=2*Root%Index+1
      Root%Child2%TypeCut='y'
      Root%Child2%P0=Root%P0
      Root%Child2%P1=Root%P1
      Root%Child2%P0%x=0.5d0*(Root%P0%x+Root%P1%x)
      Root%Child2%Index=2*Root%Index+2
    ELSE
      Root%Child1%TypeCut='x'
      Root%Child1%P0=Root%P0
      Root%Child1%P1=Root%P1
      Root%Child1%P1%y=0.5d0*(Root%P0%y+Root%P1%y)
      Root%Child1%Index=2*Root%Index+1
      Root%Child2%TypeCut='x'
      Root%Child2%P0=Root%P0
      Root%Child2%P1=Root%P1
      Root%Child2%P0%y=0.5d0*(Root%P0%y+Root%P1%y)
      Root%Child2%Index=2*Root%Index+2
   END IF
   CALL CreateTree(Root%Child1,Depth-1)
   CALL CreateTree(Root%Child2,Depth-1)
 END IF

END SUBROUTINE CreateTree

RECURSIVE SUBROUTINE FindBoxPoint(Root,Box,P)

  TYPE(Box_T), TARGET :: Root
  TYPE(Box_T), POINTER :: Box
  TYPE(Point_T) :: P

  Box=>Root
  IF (ASSOCIATED(Root%Child1)) THEN
    IF (InsideBoxPoint(Root%Child1,P)) THEN
      CALL FindBoxPoint(Root%Child1,Box,P)
    END IF
  END IF
  IF (ASSOCIATED(Root%Child2)) THEN
    IF (InsideBoxPoint(Root%Child2,P)) THEN
      CALL FindBoxPoint(Root%Child2,Box,P)
    END IF
  END IF
END SUBROUTINE FindBoxPoint

FUNCTION InsideBoxPoint(Box,P)

   LOGICAL :: InsideBoxPoint
   TYPE(Box_T) :: Box
   TYPE(Point_T) :: P

   IF (Box%P0<=P.AND.P<=Box%P1) THEN
     InsideBoxPoint=.TRUE.
   ELSE
     InsideBoxPoint=.FALSE.
   END IF

END FUNCTION InsideBoxPoint

FUNCTION InsideBoxHaus(Box,Haus)

   LOGICAL :: InsideBoxHaus
   TYPE(Box_T) :: Box
   TYPE(Haus_T) :: Haus
 
   INTEGER :: i,j

   InsideBoxHaus=.TRUE.
   S1:DO i=1,Haus%NumberOfFaces
     DO j=1,Haus%Faces(i)%NumberOfPoints
       IF (.NOT.InsideBoxPoint(Box,Haus%Faces(i)%Points(j))) THEN
         InsideBoxHaus=.FALSE.
         EXIT S1
       END IF
     END DO
   END DO S1

END FUNCTION InsideBoxHaus

RECURSIVE SUBROUTINE InsertHaus(Box,Haus,Number)

  TYPE(Box_T) :: Box
  TYPE(Haus_T) :: Haus
  INTEGER :: Number

  INTEGER, ALLOCATABLE :: List(:)
  INTEGER :: LenList

  DO 
    IF (ASSOCIATED(Box%Child1)) THEN
      IF (InsideBoxHaus(Box%Child1,Haus)) THEN
        CALL InsertHaus(Box%Child1,Haus,Number)
        EXIT
      END IF
    END IF
    IF (ASSOCIATED(Box%Child2)) THEN
      IF (InsideBoxHaus(Box%Child2,Haus)) THEN
        CALL InsertHaus(Box%Child2,Haus,Number)
        EXIT
      END IF
    END IF
    IF (ASSOCIATED(Box%List)) THEN
      LenList=SIZE(Box%List)
      ALLOCATE(List(1:LenList))
      List=Box%List
      DEALLOCATE(Box%List)
      ALLOCATE(Box%List(1:LenList+1))
      Box%List(1:LenList)=List
      Box%List(LenList+1)=Number
      DEALLOCATE(List)
    ELSE
      ALLOCATE(Box%List(1:1))
      Box%List(1)=Number
    END IF
    EXIT
  END DO
END SUBROUTINE InsertHaus

END MODULE Tree_Mod
