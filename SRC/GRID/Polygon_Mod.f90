MODULE Polygon_Mod

  USE Kind_Mod
  USE Geometry_Mod

  IMPLICIT NONE

  TYPE VertexS_T
    TYPE (Point_T) :: Point
    INTEGER :: ix,iy,iz
    INTEGER :: nrP=0
    INTEGER :: nrInP=0
    INTEGER :: nrCutP=0
    INTEGER :: nrI=0  !for ParaView
  END TYPE

  TYPE Vertex_T
    TYPE (Point_T) :: Point
    INTEGER :: ix=-1000
    INTEGER :: iy=-1000
    INTEGER :: iz=-1000
    INTEGER :: in_out=1
    INTEGER :: nrP=0
    INTEGER :: nrInP=0
    INTEGER :: nrCutP=0
    INTEGER :: Shift=1
    INTEGER :: iType
    REAL(RealKind) :: Dist
    INTEGER :: Number=0  !for ParaView
    TYPE(Vertex_T), POINTER :: VertX=>NULL()
    TYPE(Vertex_T), POINTER :: VertY=>NULL()
    TYPE(Vertex_T), POINTER :: VertZ=>NULL()
  END TYPE

  TYPE VertexP_T
    TYPE(Vertex_T), POINTER :: Vertex=>NULL()
    INTEGER :: Number=0  !for ParaView
  END TYPE

  TYPE Edge_T
    TYPE (Vertex_T), POINTER :: Vert1=>NULL()
    TYPE (Vertex_T), POINTER :: Vert2=>NULL()
    INTEGER :: eg_nr=0
    INTEGER :: in_out
    INTEGER :: yes_sp=-1
    CHARACTER*1 :: EdgeType=''
    INTEGER :: Number=0  !for ParaView
    INTEGER :: NumberNode=0  !for ParaView
  END TYPE

  TYPE EdgeP_T
    TYPE(Edge_T), POINTER :: Edge=>NULL()
    INTEGER :: Number=0  !for ParaView
  END TYPE

  TYPE Face_T
    TYPE (Edge_T), POINTER :: Edge1=>NULL()
    TYPE (Edge_T), POINTER :: Edge2=>NULL()
    TYPE (Edge_T), POINTER :: Edge3=>NULL()
    TYPE (Edge_T), POINTER :: Edge4=>NULL()
    INTEGER :: in_out
    INTEGER :: NumberVert
    INTEGER :: VertexList(1:8)
    INTEGER :: ec
    INTEGER :: ec_vs=0
    INTEGER :: ec_vg=0 
    INTEGER :: mp=0
    INTEGER :: egc_nr=0
    INTEGER :: EdgeCut(1:4)     !EdgeCut(1:2)!EdgeCut(1:6)
    REAL(8) :: Vol=-1111111.111
    TYPE (Point_T) :: MidPoint
    CHARACTER*2 :: FaceType=''
    INTEGER :: Number=0  !for ParaView
    INTEGER :: NumberEdge=0
  END TYPE

  TYPE FaceP_T
    TYPE(Face_T), POINTER :: Face=>NULL()
    INTEGER :: Number=0  !for ParaView
  END TYPE

  TYPE Cell_T
    TYPE (Face_T), POINTER :: Face1=>NULL()
    TYPE (Face_T), POINTER :: Face2=>NULL()
    TYPE (Face_T), POINTER :: Face3=>NULL()
    TYPE (Face_T), POINTER :: Face4=>NULL()
    TYPE (Face_T), POINTER :: Face5=>NULL()
    TYPE (Face_T), POINTER :: Face6=>NULL()
    INTEGER :: in_out
    INTEGER :: vc
    INTEGER :: mp=0
    INTEGER, POINTER :: VertCut(:)=>NULL()
    INTEGER, POINTER :: VertCutOnlyP(:)=>NULL()
    REAL(8) :: Vol=-1111111.111
    TYPE (Point_T) :: MidPoint
    TYPE (Point_T) :: CutF_MidP
    INTEGER, POINTER :: layer_soil(:)=>NULL()
    INTEGER :: LandClass=0
    INTEGER :: Number=0  !for ParaView
    INTEGER :: NumberFace=0
  END TYPE

  TYPE CellP_T
    TYPE(Cell_T), POINTER :: Cell=>NULL()
    INTEGER :: Number=0  !for ParaView
  END TYPE

  TYPE FaceP4_T
    TYPE (Vertex_T) :: V0,V1,V2,V3
    TYPE (Edge_T) :: Edge1
    TYPE (Edge_T) :: Edge2
    TYPE (Edge_T) :: Edge3
    TYPE (Edge_T) :: Edge4
    INTEGER :: in_out
    INTEGER :: NumberVert
    INTEGER :: VertexList(1:8)
    INTEGER :: ec
    INTEGER :: mp=0
    INTEGER :: EdgeCut(1:2)       !EdgeCut(1:6)
    REAL(8) :: Vol=-1111111.111
    TYPE (Point_T) :: MidPoint
  END TYPE

  TYPE CellP8_T
    TYPE (Vertex_T) :: V0,V1,V2,V3,V4,V5,V6,V7
    TYPE (Face_T) :: Face1
    TYPE (Face_T) :: Face2
    TYPE (Face_T) :: Face3
    TYPE (Face_T) :: Face4
    TYPE (Face_T) :: Face5
    TYPE (Face_T) :: Face6
    INTEGER :: in_out    
    INTEGER :: vc
    INTEGER :: mp=0
    REAL(8) :: Vol=-1111111.111
    TYPE (Point_T) :: MidPoint
  END TYPE

CONTAINS

SUBROUTINE WriteFace(Face)
  TYPE(Face_T), POINTER :: Face
  WRITE(*,*) ' Face%--->'
  WRITE(*,*) "    !(TYPE Vertex_T:", "  Point(x,y,z)    in_out   nrP   nrInP   nrCutP   Shift)"
  WRITE(*,*) "    Edge1%Vert1:",Face%Edge1%Vert1%Point 
  WRITE(*,*) "    Edge1%Vert2:",Face%Edge1%Vert2%Point
  WRITE(*,*) "    Edge3%Vert1:",Face%Edge3%Vert1%Point
  WRITE(*,*) "    Edge3%Vert2:",Face%Edge3%Vert2%Point
  WRITE(*,*) "    ..........."
  WRITE(*,*) '    Vol        =',Face%Vol
  WRITE(*,*) '    NumberVert =',Face%NumberVert
  WRITE(*,*) '    VertexList =',Face%VertexList(1:Face%NumberVert)
  WRITE(*,*) '    in_out     =',Face%in_out
  WRITE(*,*) '    ec         =',Face%ec
  WRITE(*,*) '    EdgeCut    =',Face%EdgeCut 
END SUBROUTINE WriteFace

SUBROUTINE WriteCellDebug(Cell)
  TYPE(Cell_T), POINTER :: Cell

  WRITE(*,*) 'Cell'
  CALL WriteFaceDebug(Cell%Face1)
  CALL WriteFaceDebug(Cell%Face2)
  CALL WriteFaceDebug(Cell%Face3)
  CALL WriteFaceDebug(Cell%Face4)
  CALL WriteFaceDebug(Cell%Face5)
  CALL WriteFaceDebug(Cell%Face6)
END SUBROUTINE WriteCellDebug
SUBROUTINE WriteFaceDebug(Face)
  TYPE(Face_T), POINTER :: Face

  WRITE(*,*) '--------------------------------------'
  WRITE(*,*) 'Face Type ',Face%FaceType
  CALL WriteEdgeDebug(Face%Edge1,Face%FaceType(1:1))
  CALL WriteEdgeDebug(Face%Edge2,Face%FaceType(1:1))
  CALL WriteEdgeDebug(Face%Edge3,Face%FaceType(1:1))
  CALL WriteEdgeDebug(Face%Edge4,Face%FaceType(1:1))
END SUBROUTINE WriteFaceDebug

SUBROUTINE WriteEdgeDebug(Edge,EdgeType)
  TYPE(Edge_T), POINTER :: Edge
  CHARACTER :: EdgeType

  IF (Edge%EdgeType==EdgeType) THEN
    WRITE(*,*) '--------------------------------------'
    WRITE(*,*) 'Edge Type ',EdgeType
    WRITE(*,*) 'yes_sp',Edge%yes_sp
    CALL WriteVertDebug(Edge%Vert1)
    CALL WriteVertDebug(Edge%Vert2)
  END IF  

END SUBROUTINE WriteEdgeDebug

SUBROUTINE WriteVertDebug(Vert)
  TYPE(Vertex_T), POINTER :: Vert

  WRITE(*,*) '--------------------------------------'
   WRITE(*,*) 'Vert ',Vert%Point
   WRITE(*,*) 'ix,iy,iz',Vert%ix,Vert%iy,Vert%iz
   WRITE(*,*) 'nrp   ',Vert%nrp
   WRITE(*,*) 'in_out',Vert%in_out
   IF (ASSOCIATED(Vert%VertX)) THEN
     WRITE(*,*) 'VertX P  ',Vert%VertX%Point
     WRITE(*,*) 'VertX nrp',Vert%VertX%nrp
   END IF  
   IF (ASSOCIATED(Vert%VertY)) THEN
     WRITE(*,*) 'VertY P  ',Vert%VertY%Point
     WRITE(*,*) 'VertY nrp',Vert%VertY%nrp
   END IF  
   IF (ASSOCIATED(Vert%VertZ)) THEN
     WRITE(*,*) 'VertZ P  ',Vert%VertZ%Point
     WRITE(*,*) 'VertZ nrp',Vert%VertZ%nrp
   END IF  

END SUBROUTINE WriteVertDebug

SUBROUTINE WriteCell(Cell,i,j,k)

  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: F1V1,F1V2,F1V3,F1V4
  INTEGER :: F2V1,F2V2,F2V3,F2V4
  !.............................. 
  IF (ASSOCIATED(Cell)) THEN
    !.................................
    F1V1=Cell%Face1%Edge1%Vert1%in_out
    F1V2=Cell%Face1%Edge1%Vert2%in_out
    F1V3=Cell%Face1%Edge3%Vert2%in_out
    F1V4=Cell%Face1%Edge3%Vert1%in_out 
    !.................................
    F2V1=Cell%Face2%Edge1%Vert1%in_out
    F2V2=Cell%Face2%Edge1%Vert2%in_out
    F2V3=Cell%Face2%Edge3%Vert2%in_out
    F2V4=Cell%Face2%Edge3%Vert1%in_out
    !................................. 

    WRITE(*,*) '-------------- WriteCell() ---------------------------------'
    WRITE(*,*) '............'
    WRITE(*,*) "        ",F2V4,"-----------",F2V3
    WRITE(*,*) "        /|  F2      /|" 
    WRITE(*,*) "      ",F2V1,"-----------",F2V2," |"  
    WRITE(*,*) "       | ",F1V4,"---------|-",F1V3
    WRITE(*,*) "       |/   F1     |/" 
    WRITE(*,*) "      ",F1V1,"-----------",F1V2   
    WRITE(*,*) '* Cell--> ',' \/',i,'\/',j,'\/',k,'\/'
    WRITE(*,*) '    Vol       =',Cell%Vol
    WRITE(*,*) '    MidPoint  =',Cell%MidPoint
    WRITE(*,*) '    vc        =',Cell%vc
    WRITE(*,*) '    in_out    =',Cell%in_out
    WRITE(*,*) '    EdgeCut   =',Cell%VertCut
    WRITE(*,*) '    CutF_MidP =',Cell%CutF_MidP
    WRITE(*,*) '    Face1%Vol =',Cell%Face1%Vol 
    WRITE(*,*) '    Face2%Vol =',Cell%Face2%Vol 
    WRITE(*,*) '    Face3%Vol =',Cell%Face3%Vol 
    WRITE(*,*) '    Face4%Vol =',Cell%Face4%Vol 
    WRITE(*,*) '    Face5%Vol =',Cell%Face5%Vol 
    WRITE(*,*) '    Face6%Vol =',Cell%Face6%Vol 
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceXY','  (F1)'
    CALL WriteFace(Cell%Face1)
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceXY','  (F2)'
    CALL WriteFace(Cell%Face2)
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceXZ','  (F3)'
    CALL WriteFace(Cell%Face3)
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceXZ','  (F4)'
    CALL WriteFace(Cell%Face4)
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceYZ','  (F5)'
    CALL WriteFace(Cell%Face5)
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceYZ','  (F6)'
    CALL WriteFace(Cell%Face6)
    WRITE(*,*) 'Cell','(',i,',',j,',',k,')'
    WRITE(*,*) '-------------- WriteCell() Ende ----------------------------'
  END IF

END SUBROUTINE WriteCell

SUBROUTINE Display_CellVol(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  Write(*,*)
  Write(*,*) "Check :"
  Write(*,*) "Cell%Vol "," ix=",i," iy=",j," iz=",k  
  IF (ASSOCIATED(Cell)) THEN
     Write(*,*) "Cell%Vol=",Cell%Vol
  ELSE
    Write(*,*) "Celle not associated"
  END IF
END SUBROUTINE Display_CellVol


SUBROUTINE Display_Cell_vc(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  Write(*,*)
  Write(*,*) "Check :"
  IF (ASSOCIATED(Cell)) THEN
    Write(*,*) "Cell", "(", i, ",", j, ",", k, ")", &
             & " Cell%vc=",Cell%vc," Cell%in_out=",Cell%in_out  
  ELSE
    Write(*,*) "Cell", "(", i, ",", j, ",", k, ")", & 
             & " Celle not associated "
  END IF
END SUBROUTINE Display_Cell_vc


SUBROUTINE DEALLOC_VertP(Edge)
  TYPE (Edge_T) Edge

END SUBROUTINE DEALLOC_VertP

SUBROUTINE DEALLOC_EdgeP(Face)
  TYPE (Face_T) :: Face

   IF (ASSOCIATED(Face%Edge1)) THEN
     CALL DEALLOC_VertP(Face%Edge1) 
     DEALLOCATE(Face%Edge1)
   END IF
   IF (ASSOCIATED(Face%Edge2)) THEN
     CALL DEALLOC_VertP(Face%Edge2) 
     DEALLOCATE(Face%Edge2)
   END IF
   IF (ASSOCIATED(Face%Edge3)) THEN
     CALL DEALLOC_VertP(Face%Edge3) 
     DEALLOCATE(Face%Edge3)
   END IF
   IF (ASSOCIATED(Face%Edge4)) THEN
     CALL DEALLOC_VertP(Face%Edge4)
     DEALLOCATE(Face%Edge4)
   END IF
END SUBROUTINE DEALLOC_EdgeP

SUBROUTINE DEALLOC_Cell_FaceP(Cell)
   TYPE (Cell_T) ::Cell

   IF (ASSOCIATED(Cell%Face1)) THEN
     CALL DEALLOC_EdgeP(Cell%Face1)
     DEALLOCATE(Cell%Face1)
   END IF
   IF (ASSOCIATED(Cell%Face2)) THEN
     CALL DEALLOC_EdgeP(Cell%Face2)
     DEALLOCATE(Cell%Face2)
   END IF
   IF (ASSOCIATED(Cell%Face3)) THEN
     CALL DEALLOC_EdgeP(Cell%Face3)
     DEALLOCATE(Cell%Face3)
   END IF
   IF (ASSOCIATED(Cell%Face4)) THEN
     CALL DEALLOC_EdgeP(Cell%Face4)
     DEALLOCATE(Cell%Face4)
   END IF
   IF (ASSOCIATED(Cell%Face5)) THEN
     CALL DEALLOC_EdgeP(Cell%Face5)
     DEALLOCATE(Cell%Face5)
   END IF
   IF (ASSOCIATED(Cell%Face6)) THEN
     CALL DEALLOC_EdgeP(Cell%Face6)
     DEALLOCATE(Cell%Face6)
   END IF
END SUBROUTINE DEALLOC_Cell_FaceP

SUBROUTINE DEALLOC_Cell_CutsListP(Cell)
   TYPE (Cell_T) ::Cell

   IF (ASSOCIATED(Cell%VertCut)) THEN
      DEALLOCATE(Cell%VertCut)
   END IF
   IF (ASSOCIATED(Cell%VertCutOnlyP)) THEN
      DEALLOCATE(Cell%VertCutOnlyP)
   END IF
END SUBROUTINE DEALLOC_Cell_CutsListP


SUBROUTINE DEALLOC_Cell_Parts(Cell)
   TYPE (Cell_T) ::Cell
   
   CALL DEALLOC_Cell_CutsListP(Cell)
   CALL DEALLOC_Cell_FaceP(Cell)
END SUBROUTINE DEALLOC_Cell_Parts

END MODULE Polygon_Mod

