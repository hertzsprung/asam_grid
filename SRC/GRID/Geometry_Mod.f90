MODULE Geometry_Mod

  USE Parameter_Mod

  IMPLICIT NONE

  TYPE Point_T
    REAL(8) :: x
    REAL(8) :: y
    REAL(8) :: z
  END TYPE Point_T

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE Add
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE Diff 
  END INTERFACE
  INTERFACE OPERATOR(/)
    MODULE PROCEDURE Div 
  END INTERFACE
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE Mult,ScalarProduct
  END INTERFACE
  INTERFACE OPERATOR (.CROSS.)
    MODULE PROCEDURE CrossProduct
  END INTERFACE
  INTERFACE OPERATOR(<)
    MODULE PROCEDURE Lower
  END INTERFACE
  INTERFACE OPERATOR(<=)
    MODULE PROCEDURE LowerEQ
  END INTERFACE
  INTERFACE OPERATOR(>)
    MODULE PROCEDURE Greater
  END INTERFACE
  INTERFACE OPERATOR(>=)
    MODULE PROCEDURE GreaterEQ
  END INTERFACE
  INTERFACE OPERATOR(.PEQ.)
    MODULE PROCEDURE PEqual
  END INTERFACE
CONTAINS 

FUNCTION Lower(P1,P2)
  LOGICAL :: Lower
  TYPE(Point_T), INTENT(IN) :: P1,P2

  IF (P1%x<P2%x.AND.P1%y<P2%y.AND.P1%z<P2%z) THEN
    Lower=.TRUE.
  ELSE
    Lower=.FALSE.
  END IF
END FUNCTION Lower

FUNCTION LowerEQ(P1,P2)
  LOGICAL :: LowerEQ
  TYPE(Point_T), INTENT(IN) :: P1,P2
  
  IF (P1%x<=P2%x.AND.P1%y<=P2%y.AND.P1%z<=P2%z) THEN
    LowerEQ=.TRUE. 
  ELSE
    LowerEQ=.FALSE.
  END IF
END FUNCTION LowerEQ

FUNCTION PointEQ(P1,P2)
  LOGICAL :: PointEQ
  TYPE(Point_T), INTENT(IN) :: P1,P2
  
  IF (P1%x==P2%x.AND.P1%y==P2%y.AND.P1%z==P2%z) THEN
    PointEQ=.TRUE. 
  ELSE
    PointEQ=.FALSE.
  END IF
END FUNCTION PointEQ

FUNCTION Greater(P1,P2)
  LOGICAL :: Greater
  TYPE(Point_T), INTENT(IN) :: P1,P2

  IF (P1%x>P2%x.AND.P1%y>P2%y.AND.P1%z>P2%z) THEN
    Greater=.TRUE.
  ELSE
    Greater=.FALSE.
  END IF
END FUNCTION Greater

FUNCTION GreaterEQ(P1,P2)
  LOGICAL :: GreaterEQ
  TYPE(Point_T), INTENT(IN) :: P1,P2

  IF (P1%x>=P2%x.AND.P1%y>=P2%y.AND.P1%z>=P2%z) THEN
    GreaterEQ=.TRUE.
  ELSE
    GreaterEQ=.FALSE.
  END IF
END FUNCTION GreaterEQ

FUNCTION Add(P1,P2)
  TYPE(Point_T) :: Add
  TYPE(Point_T), INTENT(IN) :: P1,P2

  Add%x=P1%x+P2%x
  Add%y=P1%y+P2%y
  Add%z=P1%z+P2%z
END FUNCTION Add

FUNCTION Diff(V2,V1)
  TYPE(Point_T) :: Diff
  TYPE(Point_T), INTENT(IN)  :: V2,V1
  Diff%x=V2%x-V1%x
  Diff%y=V2%y-V1%y
  Diff%z=V2%z-V1%z
END FUNCTION Diff 

FUNCTION Mult(alpha,P)
  TYPE(Point_T) :: Mult
  REAL(8), INTENT(IN) :: alpha
  TYPE(Point_T), INTENT(IN) :: P

  Mult%x=alpha*P%x
  Mult%y=alpha*P%y
  Mult%z=alpha*P%z
END FUNCTION Mult

FUNCTION Div(P,temp)
  TYPE(Point_T) :: Div 
  REAL(8), INTENT(IN) :: temp
  TYPE(Point_T), INTENT(IN) :: P

  Div%x=P%x/temp
  Div%y=P%y/temp
  Div%z=P%z/temp
END FUNCTION Div 

FUNCTION ScalarProduct(V1,V2)
  REAL(8)  :: ScalarProduct
  TYPE(Point_T), INTENT(IN)  :: V1,V2

  ScalarProduct=V1%x*V2%x &
               +V1%y*V2%y &
               +V1%z*V2%z
END FUNCTION ScalarProduct 

FUNCTION CrossProduct(r1,r2)
  TYPE(Point_T) :: CrossProduct
  TYPE(Point_T), INTENT(IN) :: r1,r2

  CrossProduct%x=r1%y*r2%z-r1%z*r2%y
  CrossProduct%y=-r1%x*r2%z+r1%z*r2%x
  CrossProduct%z=r1%x*r2%y-r1%y*r2%x
END FUNCTION CrossProduct

FUNCTION PEqual(P1,P2)
  INTEGER :: PEqual
  TYPE(Point_T), INTENT(IN) :: P1,P2
  INTEGER :: diffx,diffy,diffz 
  diffx=1;diffy=1;diffz=1
  if(P1%x==P2%x)then
     diffx=0
  end if
  if(P1%y==P2%y)then
     diffy=0
  end if
  if(P1%z==P2%z)then
     diffz=0
  end if
  if(diffx==0.AND.diffy==0.AND.diffz==0)then
     PEqual=0
  else
     PEqual=1
  end if 
END FUNCTION PEqual 

FUNCTION Norm(r)
  REAL(8) :: Norm
  TYPE(Point_T) :: r

  Norm=SQRT(r%x*r%x &
           +r%y*r%y &
           +r%z*r%z)
END FUNCTION Norm

FUNCTION Vol4(P1,P2,P3,P4)
  REAL(8) :: Vol4
  TYPE(Point_T) :: P1,P2,P3,P4

  TYPE(Point_T) :: r1,r2,r3

  r1=P2-P1
  r2=P3-P1
  r3=P4-P1
  Vol4=1.0d0/6.0d0*(r1.CROSS.r2)*r3
END FUNCTION Vol4


FUNCTION VolTetra(P1,P2,P3,P4)
  REAL(8) :: VolTetra
  TYPE(Point_T) :: P1,P2,P3,P4
  !......................................
  !             4     !Point-Folge P1-P4
  !            /*\  
  !           / * \   
  !          /  *  \ 
  !         /  ,*,  \
  !        / ,* 3 *, \
  !       /,*       *,\
  !      1*-----------*2
  !......................................
  TYPE(Point_T) :: r1,r2,r3
  r1=P2-P1
  r2=P3-P1
  r3=P4-P1
  VolTetra=1.0d0/6.0d0*(r1.CROSS.r2)*r3
END FUNCTION VolTetra


FUNCTION VolTriPrisma(P1,P2,P3,P4,P5,P6)
  REAL(8) :: VolTriPrisma,VolTriPrisma1,VolTriPrisma2,VolTriPrisma3
  TYPE(Point_T) :: P1,P2,P3,P4,P5,P6
  !.....................................
  !           6\      !Point-Folge P1-P6
  !          /| \
  !         / |  \      
  !        5\ |   \     
  !        | \3----4
  !        | /\   /
  !        |/  \ /
  !        1----2
  !.....................................
  TYPE(Point_T) :: r1,r2,r3
  !...Test zum Vergleich 
  !r1=P2-P1
  !r2=P4-P1
  !r3=P5-P1
  !VolTriPrisma=1.0d0/2.0d0*(r1.CROSS.r2)*r3 
  !Write(*,*) "VolTriPrisma=",VolTriPrisma, "Berechnung, wie Wuerfel und 1/2"
  !...Berechnung ueber Tetraeder
  r1=P2-P1
  r2=P5-P1
  r3=P4-P1
  VolTriPrisma1=1.0d0/6.0d0*(r1.CROSS.r2)*r3
  r1=P6-P4
  r2=P3-P4
  r3=P1-P4
  VolTriPrisma2=1.0d0/6.0d0*(r1.CROSS.r2)*r3
  r1=P5-P1
  r2=P6-P1
  r3=P4-P1
  VolTriPrisma3=1.0d0/6.0d0*(r1.CROSS.r2)*r3
  VolTriPrisma=ABS(VolTriPrisma1)+ABS(VolTriPrisma2)+ABS(VolTriPrisma3)
  !...
  !Write(*,*) "VolTriPrisma1=",VolTriPrisma1
  !Write(*,*) "VolTriPrisma2=",VolTriPrisma2
  !Write(*,*) "VolTriPrisma3=",VolTriPrisma3
  !Write(*,*) "VolTriPrismaGesamt=",VolTriPrisma 
END FUNCTION VolTriPrisma

FUNCTION VolP8(P0,P1,P2,P3,P4,P5,P6,P7)
  REAL(8) :: VolP8
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(Point_T) :: r1,r2,r3,r4,r5,r6,r7,r8,r9
  REAL(8) :: VolTetra1,VolTetra2,VolTetra3
  !......................................................
  !          6--------7    !Point-Folge Vol-Berechnung
  !         /|       /|    !P0-P7
  !        4--------5 |
  !        | |      | |
  !        | 2------|-3
  !        |/       |/
  !        0--------1
  !.......................................................
  ! nach Jeffry Grandy:(Efficient Computation of Volume of
  !                    Hexahedral Cells (Gl:12)
  !.......................................................
  r1=(P7-P1)+(P6-P0)
  r2=P7-P2
  r3=P3-P0
  VolTetra1=(r1.CROSS.r2)*r3
  r4=(P6-P0)
  r5=(P7-P2)+(P5-P0)
  r6=P7-P4
  VolTetra2=(r4.CROSS.r5)*r6
  r7=P7-P1
  r8=P5-P0
  r9=(P7-P4)+(P3-P0)
  VolTetra3=(r7.CROSS.r8)*r9
  VolP8=1.0d0/12.0d0*(VolTetra1+VolTetra2+VolTetra3)

END FUNCTION VolP8

FUNCTION FaceTriangle(P1,P2,P3)
  REAL(8) :: FaceTriangle
  TYPE(Point_T) :: P1,P2,P3

  TYPE(Point_T) :: r1,r2,r3

  r1=P1-P2
  r2=P1-P3
  r3=r1.CROSS.r2
  FaceTriangle=0.5d0*Norm(r3)
END FUNCTION FaceTriangle


FUNCTION PcFace(P1,P2,P3,P4)
  TYPE(Point_T) :: PcFace
  TYPE(Point_T) :: P1,P2,P3,P4

  TYPE(Point_T) :: PSum,r1,r2,r3
  PcFace=0.25d0*(P1+P2+P3+P4)
END FUNCTION PcFace

FUNCTION VolQuadrilateral(P1,P2,P3,P4)
  REAL(8) :: VolQuadrilateral
  TYPE(Point_T) :: P1,P2,P3,P4
  TYPE(Point_T) :: Pc
  REAL(8) :: Vol1,Vol2,Vol3,Vol4
  Pc=PcFace(P1,P2,P3,P4)
  Vol1=FaceTriangle(P1,P2,Pc)
  Vol2=FaceTriangle(P2,P3,Pc)
  Vol3=FaceTriangle(P3,P4,Pc)
  Vol4=FaceTriangle(P4,P1,Pc)
  VolQuadrilateral=Vol1+Vol2+Vol3+Vol4
END FUNCTION VolQuadrilateral


FUNCTION FaceP4(P1,P2,P3,P4)
  REAL(8) :: FaceP4 
  TYPE(Point_T) :: P1,P2,P3,P4
!    4-----3   !Point-Folge Face-Berechnung
!    | \   |
!    |   \ |
!    1-----2
  FaceP4=FaceTriangle(P1,P2,P4)
  FaceP4=FaceP4+FaceTriangle(P2,P3,P4)
!    4-----3   !Point-Folge Face-Berechnung
!    |  /  |
!    |/    |
!    1-----2
  FaceP4=FaceP4+FaceTriangle(P1,P2,P3)
  FaceP4=FaceP4+FaceTriangle(P1,P3,P4)
  FaceP4=0.5d0*FaceP4
END FUNCTION FaceP4

FUNCTION FaceMidP3(P1,P2,P3)
  TYPE(Point_T) :: FaceMidP3
  TYPE(Point_T) :: P1,P2,P3
  FaceMidP3=P1+P2+P3
  FaceMidP3=FaceMidP3/3.0d0
END FUNCTION FaceMidP3

FUNCTION FaceMidP4(P1,P2,P3,P4)
  TYPE(Point_T) :: FaceMidP4
  TYPE(Point_T) :: P1,P2,P3,P4
  !FaceMidP4=0.25d0*(P1+P2+P3+P4)
  FaceMidP4=P1+P2+P3+P4
  FaceMidP4=FaceMidP4/4.0d0
END FUNCTION FaceMidP4


FUNCTION FaceMidP5(P1,P2,P3,P4,P5)
  TYPE(Point_T) :: FaceMidP5
  TYPE(Point_T) :: P1,P2,P3,P4,P5
  FaceMidP5=P1+P2+P3+P4+P5
  FaceMidP5=FaceMidP5/5.0d0
END FUNCTION FaceMidP5

FUNCTION FaceMidP6(P1,P2,P3,P4,P5,P6)
  TYPE(Point_T) :: FaceMidP6
  TYPE(Point_T) :: P1,P2,P3,P4,P5,P6
  FaceMidP6=P1+P2+P3+P4+P5+P6
  FaceMidP6=FaceMidP6/6.0d0
END FUNCTION FaceMidP6


FUNCTION VolMidP4(P0,P1,P2,P3)
  TYPE(Point_T) :: VolMidP4
  TYPE(Point_T) :: P0,P1,P2,P3

  VolMidP4=P0+P1+P2+P3
  VolMidP4=VolMidP4/4.0d0
END FUNCTION VolMidP4

FUNCTION VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
  TYPE(Point_T) :: VolMidP8
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  VolMidP8=P0+P1+P2+P3+P4+P5+P6+P7
  VolMidP8=VolMidP8/8.0d0
END FUNCTION VolMidP8

FUNCTION FaceGlobeXY(r,phi1,phi2,lam1,lam2)
  REAL(8) :: FaceGlobeXY
  REAL(8) :: r
  REAL(8) :: phi1,phi2,lam1,lam2
  !Wertebereich: phi [-Pi/2;Pi/2]
  !Wertebereich: lambda [0;2*Pi]
  FaceGlobeXY=r**2*(SIN(phi2)-SIN(phi1))*(lam2-lam1)
END FUNCTION FaceGlobeXY

FUNCTION FaceGlobeYZ(r1,r2,phi1,phi2)
  REAL(8) :: FaceGlobeYZ
  REAL(8) :: r1,r2,phi1,phi2
  FaceGlobeYZ=0.5*(phi1-phi2)*(r2**2-r1**2)
END FUNCTION FaceGlobeYZ

FUNCTION FaceGlobeZX(r1,r2,phi,lam1,lam2)
  REAL(8) :: FaceGlobeZX
  REAL(8) :: r1,r2,phi,lam1,lam2
  FaceGlobeZX=0.5*COS(phi)*(lam2-lam1)*(r2**2-r1**2)
END FUNCTION FaceGlobeZX

FUNCTION VolCellGlobe(r1,r2,phi1,phi2,lam1,lam2)
  REAL(8) :: VolCellGlobe
  REAL(8) :: r1,r2,phi1,phi2,lam1,lam2
  VolCellGlobe=1/3*((lam2-lam1)*(r2**3-r1**3)*(SIN(phi1)-SIN(phi2)))
END FUNCTION VolCellGlobe


END MODULE Geometry_Mod

