MODULE Haus_Mod

  USE Geometry_Mod 
  IMPLICIT NONE

  TYPE GFace_T
    INTEGER :: NumberOfPoints
    TYPE(Point_T), POINTER :: Points(:)
    INTEGER, POINTER :: ListOfPoints(:)
    TYPE(Point_T) :: Normal 
    CHARACTER :: Type         ! r=roof,w=wall
  END TYPE GFace_T

  TYPE Haus_T
    INTEGER :: Number
    INTEGER :: NumberOfPoints
    INTEGER :: NumberOfFaces
    TYPE(Point_T), POINTER :: Points(:)
    TYPE(GFace_T), POINTER :: Faces(:)
    Type(Point_T) :: P0,P1
  END TYPE Haus_T
 

CONTAINS 

FUNCTION Dist(P,Haus)
  REAL(8) :: Dist
  TYPE(Point_T) :: P
  TYPE(Haus_T) :: Haus

  INTEGER :: i
  REAL(8) :: DistFace
  TYPE(Point_T) :: P1
  TYPE(GFace_T), POINTER :: Face

  IF (Haus%P0<=P.AND.P<=Haus%P1) THEN
    Dist=-1.0d99
    DO i=1,Haus%NumberOfFaces
      Face=>Haus%Faces(i)
      P1=Face%Points(1)
      DistFace=(P-P1)*Face%Normal
      Dist=MAX(Dist,DistFace)
    END DO
  ELSE
    Dist=1.0d0
  END IF 
  Dist=Dist+1.d-8 !OSSI
END FUNCTION Dist

SUBROUTINE BoundingBox(Haus)
  TYPE(Haus_T) :: Haus

  INTEGER :: i,j
  TYPE(GFace_T), POINTER :: Face
  Haus%P0%x=1.d99
  Haus%P0%y=1.d99
  Haus%P0%z=1.d99
  Haus%P1%x=-1.d99
  Haus%P1%y=-1.d99
  Haus%P1%z=-1.d99
  DO i=1,Haus%NumberOfFaces
    Face=>Haus%Faces(i)
    DO j=1,Face%NumberOfPoints
      Haus%P0%x=MIN(Haus%P0%x,Face%Points(j)%x)
      Haus%P0%y=MIN(Haus%P0%y,Face%Points(j)%y)
      Haus%P0%z=MIN(Haus%P0%z,Face%Points(j)%z)
      Haus%P1%x=MAX(Haus%P1%x,Face%Points(j)%x)
      Haus%P1%y=MAX(Haus%P1%y,Face%Points(j)%y)
      Haus%P1%z=MAX(Haus%P1%z,Face%Points(j)%z)
    END DO
  END DO
END SUBROUTINE BoundingBox

SUBROUTINE NormalForm(Haus)
  TYPE(Haus_T) :: Haus

  INTEGER :: i
  TYPE(Point_T) :: P1,P2,P3
  TYPE(Point_T) :: P3mP1,P3mP2
  TYPE(GFace_T), POINTER :: Face
  REAL(8) :: Temp

  DO i=1,Haus%NumberOfFaces
    Face=>Haus%Faces(i)
    P1=Face%Points(1)
    P2=Face%Points(2)
    P3=Face%Points(3)
    P3mP1=P3-P1
    P3mP2=P3-P2
    Face%Normal=P3mP1.CROSS.P3mP2
    Temp=Norm(Face%Normal)
    Face%Normal=Face%Normal/Temp
  END DO
END SUBROUTINE NormalForm

END MODULE Haus_Mod

