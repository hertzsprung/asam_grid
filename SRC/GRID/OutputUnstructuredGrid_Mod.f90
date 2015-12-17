MODULE OutputUnstructuredGrid_Mod

  USE Domain_Mod
  USE Floor_Mod
  IMPLICIT NONE

  INTEGER, PRIVATE :: OutputUnit=10

CONTAINS
SUBROUTINE WriteUnstructuredGrid(FileName)
  CHARACTER*50 :: FileName

  INTEGER :: ib
  INTEGER :: ix,iy,iz
  INTEGER :: NumberNodes
  INTEGER :: NumberEdges
  INTEGER :: NumberFaces
  INTEGER :: NumberCells
  INTEGER :: NumberMin
  TYPE(Vertex_T), POINTER :: Vertex

  NumberNodes=0
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0,iz1
      DO iy=iy0,iy1
        DO ix=ix0,ix1
          Vertex=>Vertices(ix,iy,iz)%Vertex
          IF (ASSOCIATED(Vertex)) THEN
            NumberNodes=NumberNodes+1
            Vertices(ix,iy,iz)%Number=NumberNodes
          END IF  
        END DO  
      END DO  
    END DO  
  END DO  
  NumberEdges=0
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0,iz1
      DO iy=iy0,iy1
        DO ix=ix0+1,ix1
          CALL CountEdge(Edges_X(ix,iy,iz),NumberEdges,NumberNodes)
        END DO  
      END DO  
    END DO  
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0,ix1
          CALL CountEdge(Edges_Y(ix,iy,iz),NumberEdges,NumberNodes)
        END DO  
      END DO  
    END DO  
    DO iz=iz0+1,iz1
      DO iy=iy0,iy1
        DO ix=ix0,ix1
          CALL CountEdge(Edges_Z(ix,iy,iz),NumberEdges,NumberNodes)
        END DO  
      END DO  
    END DO  
  END DO  
  NumberFaces=0
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          CALL CountFace(Faces_XY(ix,iy,iz),NumberFaces,NumberEdges)
        END DO  
      END DO  
    END DO  
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0,ix1
          CALL CountFace(Faces_YZ(ix,iy,iz),NumberFaces,NumberEdges)
        END DO  
      END DO  
    END DO  
    DO iz=iz0+1,iz1
      DO iy=iy0,iy1
        DO ix=ix0+1,ix1
          CALL CountFace(Faces_ZX(ix,iy,iz),NumberFaces,NumberEdges)
        END DO  
      END DO  
    END DO  
  END DO  
  NumberCells=0
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          CALL CountCell(Cells(ix,iy,iz),NumberCells,NumberFaces)
        END DO  
      END DO  
    END DO  
  END DO  
  OPEN(UNIT=OutputUnit,FILE=TRIM(FileName)//'.unstr',STATUS='UNKNOWN')
  WRITE(OutputUnit,*) TRIM(FileName)
  WRITE(OutputUnit,*) 'NumberNodes  ',NumberNodes
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0,iz1
      DO iy=iy0,iy1
        DO ix=ix0,ix1
          Vertex=>Vertices(ix,iy,iz)%Vertex
          IF (Vertices(ix,iy,iz)%Number>0) THEN
            WRITE(OutputUnit,*) 'n   ',Vertices(ix,iy,iz)%Number,Vertex%Point
          END IF  
        END DO  
      END DO  
    END DO  
  END DO  
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0,iz1
      DO iy=iy0,iy1
        DO ix=ix0+1,ix1
          CALL WriteVertexCut(Edges_X(ix,iy,iz),OutputUnit)
        END DO  
      END DO  
    END DO  
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0,ix1
          CALL WriteVertexCut(Edges_Y(ix,iy,iz),OutputUnit)
        END DO  
      END DO  
    END DO  
    DO iz=iz0+1,iz1
      DO iy=iy0,iy1
        DO ix=ix0,ix1
          CALL WriteVertexCut(Edges_Z(ix,iy,iz),OutputUnit)
        END DO  
      END DO  
    END DO  
  END DO  

  WRITE(OutputUnit,*) 'NumberEdges  ',NumberEdges
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0,iz1
      DO iy=iy0,iy1
        DO ix=ix0+1,ix1
          IF (Edges_X(ix,iy,iz)%Number>0) THEN
            WRITE(OutputUnit,*) 'ex  ',Edges_X(ix,iy,iz)%Number,Vertices(ix-1,iy,iz)%Number &
                                      ,Vertices(ix,iy,iz)%Number
          ELSE IF (Edges_X(ix,iy,iz)%Edge%NumberNode==0) THEN
            WRITE(OutputUnit,*) 'ex  ',Edges_X(ix,iy,iz)%Edge%Number,Vertices(ix-1,iy,iz)%Number &
                                      ,Vertices(ix,iy,iz)%Number
          ELSE
            IF (Edges_X(ix,iy,iz)%Edge%Vert1%in_out>=0) THEN
              WRITE(OutputUnit,*) 'exc ',Edges_X(ix,iy,iz)%Edge%Number,Vertices(ix-1,iy,iz)%Number &
                                        ,Edges_X(ix,iy,iz)%Edge%NumberNode
              WRITE(OutputUnit,*) 'exc ',Edges_X(ix,iy,iz)%Edge%Number+1,Edges_X(ix,iy,iz)%Edge%NumberNode &
                                        ,Vertices(ix,iy,iz)%Number 
            ELSE
              WRITE(OutputUnit,*) 'exc ',Edges_X(ix,iy,iz)%Edge%Number,Edges_X(ix,iy,iz)%Edge%NumberNode &
                                        ,Vertices(ix,iy,iz)%Number 
              WRITE(OutputUnit,*) 'exc ',Edges_X(ix,iy,iz)%Edge%Number+1,Vertices(ix-1,iy,iz)%Number &
                                        ,Edges_X(ix,iy,iz)%Edge%NumberNode
            END IF
          END IF  
        END DO  
      END DO  
    END DO  
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0,ix1
          IF (Edges_Y(ix,iy,iz)%Number>0) THEN
            WRITE(OutputUnit,*) 'ey  ',Edges_Y(ix,iy,iz)%Number,Vertices(ix,iy-1,iz)%Number &
                                      ,Vertices(ix,iy,iz)%Number
          ELSE IF (Edges_Y(ix,iy,iz)%Edge%NumberNode==0) THEN
            WRITE(OutputUnit,*) 'ey  ',Edges_Y(ix,iy,iz)%Edge%Number,Vertices(ix,iy-1,iz)%Number &
                                      ,Vertices(ix,iy,iz)%Number
          ELSE
            IF (Edges_Y(ix,iy,iz)%Edge%Vert1%in_out>=0) THEN
              WRITE(OutputUnit,*) 'eyc ',Edges_Y(ix,iy,iz)%Edge%Number,Vertices(ix,iy-1,iz)%Number &
                                        ,Edges_Y(ix,iy,iz)%Edge%NumberNode
              WRITE(OutputUnit,*) 'eyc ',Edges_Y(ix,iy,iz)%Edge%Number+1,Edges_Y(ix,iy,iz)%Edge%NumberNode &
                                        ,Vertices(ix,iy,iz)%Number 
            ELSE
              WRITE(OutputUnit,*) 'eyc ',Edges_Y(ix,iy,iz)%Edge%Number,Edges_Y(ix,iy,iz)%Edge%NumberNode &
                                        ,Vertices(ix,iy,iz)%Number 
              WRITE(OutputUnit,*) 'eyc ',Edges_Y(ix,iy,iz)%Edge%Number+1,Vertices(ix,iy-1,iz)%Number &
                                        ,Edges_Y(ix,iy,iz)%Edge%NumberNode
            END IF
          END IF  
        END DO  
      END DO  
    END DO  
    DO iz=iz0+1,iz1
      DO iy=iy0,iy1
        DO ix=ix0,ix1
          IF (Edges_Z(ix,iy,iz)%Number>0) THEN
            WRITE(OutputUnit,*) 'ez  ',Edges_Z(ix,iy,iz)%Number,Vertices(ix,iy,iz-1)%Number &
                                      ,Vertices(ix,iy,iz)%Number
          ELSE IF (Edges_Z(ix,iy,iz)%Edge%NumberNode==0) THEN
            WRITE(OutputUnit,*) 'ez  ',Edges_Z(ix,iy,iz)%Edge%Number,Vertices(ix,iy,iz-1)%Number &
                                      ,Vertices(ix,iy,iz)%Number
          ELSE
            IF (Edges_Z(ix,iy,iz)%Edge%Vert1%in_out>=0) THEN
              WRITE(OutputUnit,*) 'ezc ',Edges_Z(ix,iy,iz)%Edge%Number,Vertices(ix,iy,iz-1)%Number &
                                        ,Edges_Z(ix,iy,iz)%Edge%NumberNode
              WRITE(OutputUnit,*) 'ezc ',Edges_Z(ix,iy,iz)%Edge%Number+1,Edges_Z(ix,iy,iz)%Edge%NumberNode &
                                        ,Vertices(ix,iy,iz)%Number 
            ELSE
              WRITE(OutputUnit,*) 'ezc ',Edges_Z(ix,iy,iz)%Edge%Number,Edges_Z(ix,iy,iz)%Edge%NumberNode &
                                        ,Vertices(ix,iy,iz)%Number 
              WRITE(OutputUnit,*) 'ezc ',Edges_Z(ix,iy,iz)%Edge%Number+1,Vertices(ix,iy,iz-1)%Number &
                                        ,Edges_Z(ix,iy,iz)%Edge%NumberNode
            END IF
          END IF  
        END DO  
      END DO  
    END DO  
  END DO  
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          CALL WriteEdgeCut(Faces_XY(ix,iy,iz),OutputUnit)
        END DO  
      END DO  
    END DO  
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0,ix1
          CALL WriteEdgeCut(Faces_YZ(ix,iy,iz),OutputUnit)
        END DO  
      END DO  
    END DO  
    DO iz=iz0+1,iz1
      DO iy=iy0,iy1
        DO ix=ix0+1,ix1
          CALL WriteEdgeCut(Faces_ZX(ix,iy,iz),OutputUnit)
        END DO  
      END DO  
    END DO  
  END DO  

  WRITE(OutputUnit,*) 'NumberFaces  ',NumberFaces
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          IF (Faces_XY(ix,iy,iz)%Number>0) THEN
            WRITE(OutputUnit,*) 'fxy ',Faces_XY(ix,iy,iz)%Number &
                                      ,4,1                         
            WRITE(OutputUnit,*)        EdgeNumber(Edges_X(ix,iy-1,iz)) &
                                      ,EdgeNumber(Edges_Y(ix,iy,iz)) &
                                      ,EdgeNumber(Edges_X(ix,iy,iz)) &
                                      ,EdgeNumber(Edges_Y(ix-1,iy,iz)) 
          ELSE IF (Faces_XY(ix,iy,iz)%Face%NumberEdge==0) THEN                            
            WRITE(OutputUnit,*) 'fxy ',Faces_XY(ix,iy,iz)%Face%Number &
                                      ,4,1                         
            WRITE(OutputUnit,*)        Edges_X(ix,iy-1,iz)%Edge%Number &
                                      ,Edges_Y(ix,iy,iz)%Edge%Number &
                                      ,Edges_X(ix,iy,iz)%Edge%Number &
                                      ,Edges_Y(ix-1,iy,iz)%Edge%Number 
          ELSE
            CALL WriteCutFace(Faces_XY(ix,iy,iz),OutputUnit,'xy')
          END IF                            
        END DO  
      END DO  
    END DO  
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0,ix1
          IF (Faces_YZ(ix,iy,iz)%Number>0) THEN
            WRITE(OutputUnit,*) 'fyz ',Faces_YZ(ix,iy,iz)%Number &
                                      ,4,1                         
            WRITE(OutputUnit,*)        EdgeNumber(Edges_Y(ix,iy,iz-1)) &
                                      ,EdgeNumber(Edges_Z(ix,iy,iz)) &
                                      ,EdgeNumber(Edges_Y(ix,iy,iz)) &
                                      ,EdgeNumber(Edges_Z(ix,iy-1,iz)) 
          ELSE IF (Faces_YZ(ix,iy,iz)%Face%NumberEdge==0) THEN
            WRITE(OutputUnit,*) 'fyz ',Faces_YZ(ix,iy,iz)%Face%Number &
                                      ,4,1                         
            WRITE(OutputUnit,*)        Edges_Y(ix,iy,iz-1)%Edge%Number &
                                      ,Edges_Z(ix,iy,iz)%Edge%Number &
                                      ,Edges_Y(ix,iy,iz)%Edge%Number &
                                      ,Edges_Z(ix,iy-1,iz)%Edge%Number 
          ELSE
            CALL WriteCutFace(Faces_YZ(ix,iy,iz),OutputUnit,'yz')
          END IF                            
        END DO  
      END DO  
    END DO  
    DO iz=iz0+1,iz1
      DO iy=iy0,iy1
        DO ix=ix0+1,ix1
          IF (Faces_ZX(ix,iy,iz)%Number>0) THEN
            WRITE(OutputUnit,*) 'fzx ',Faces_ZX(ix,iy,iz)%Number &
                                      ,4,1                         
            WRITE(OutputUnit,*)        EdgeNumber(Edges_Z(ix-1,iy,iz)) &
                                      ,EdgeNumber(Edges_X(ix,iy,iz)) &
                                      ,EdgeNumber(Edges_Z(ix,iy,iz)) &
                                      ,EdgeNumber(Edges_X(ix,iy,iz-1)) 
          ELSE IF (Faces_ZX(ix,iy,iz)%Face%NumberEdge==0) THEN
            WRITE(OutputUnit,*) 'fzx ',Faces_ZX(ix,iy,iz)%Face%Number &
                                      ,4,1                         
            WRITE(OutputUnit,*)        Edges_Z(ix-1,iy,iz)%Edge%Number &
                                      ,Edges_X(ix,iy,iz)%Edge%Number &
                                      ,Edges_Z(ix,iy,iz)%Edge%Number &
                                      ,Edges_X(ix,iy,iz-1)%Edge%Number 
          ELSE
            CALL WriteCutFace(Faces_ZX(ix,iy,iz),OutputUnit,'zx')
          END IF                            
        END DO  
      END DO  
    END DO  
  END DO  
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          CALL WriteFaceCut(Cells(ix,iy,iz),OutputUnit)
        END DO  
      END DO  
    END DO  
  END DO  
  WRITE(OutputUnit,*) 'NumberCells  ',NumberCells
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          IF (Cells(ix,iy,iz)%Number>0) THEN
            WRITE(OutputUnit,*) 'c   ',Cells(ix,iy,iz)%Number &
                                      ,6,0                      
            WRITE(OutputUnit,*)        FaceNumber(Faces_XY(ix,iy,iz-1)) &
                                      ,FaceNumber(Faces_XY(ix,iy,iz)) &
                                      ,FaceNumber(Faces_YZ(ix-1,iy,iz)) &
                                      ,FaceNumber(Faces_YZ(ix,iy,iz)) &
                                      ,FaceNumber(Faces_ZX(ix,iy-1,iz)) &
                                      ,FaceNumber(Faces_ZX(ix,iy,iz))  
            WRITE(100,*) 'c   ',Cells(ix,iy,iz)%Number &
                                      ,6,0                      
            WRITE(100,*)        FaceNumber(Faces_XY(ix,iy,iz-1)) &
                                      ,FaceNumber(Faces_XY(ix,iy,iz)) &
                                      ,FaceNumber(Faces_YZ(ix-1,iy,iz)) &
                                      ,FaceNumber(Faces_YZ(ix,iy,iz)) &
                                      ,FaceNumber(Faces_ZX(ix,iy-1,iz)) &
                                      ,FaceNumber(Faces_ZX(ix,iy,iz))  
            WRITE(100,*) Vertices(ix-1,iy-1,iz-1)%Vertex%Point                          
            WRITE(100,*) Vertices(ix  ,iy-1,iz-1)%Vertex%Point                          
            WRITE(100,*) Vertices(ix  ,iy  ,iz-1)%Vertex%Point                          
            WRITE(100,*) Vertices(ix-1,iy  ,iz-1)%Vertex%Point                          
            WRITE(100,*) Vertices(ix-1,iy-1,iz  )%Vertex%Point                          
            WRITE(100,*) Vertices(ix  ,iy-1,iz  )%Vertex%Point                          
            WRITE(100,*) Vertices(ix  ,iy  ,iz  )%Vertex%Point                          
            WRITE(100,*) Vertices(ix-1,iy  ,iz  )%Vertex%Point                          
          ELSE 
            CALL WriteCellCut(Cells(ix,iy,iz),OutputUnit)
          END IF  
        END DO  
      END DO  
    END DO  
  END DO  
  CLOSE(OutputUnit)
END SUBROUTINE WriteUnstructuredGrid

FUNCTION FaceNumber(FaceP)
  INTEGER :: FaceNumber
  TYPE(FaceP_T) :: FaceP

  IF (FaceP%Number>0) THEN
    FaceNumber=FaceP%Number
  ELSE  
    FaceNumber=FaceP%Face%Number
  END IF  
END FUNCTION FaceNumber

FUNCTION EdgeNumber(EdgeP)
  INTEGER :: EdgeNumber
  TYPE(EdgeP_T) :: EdgeP

  IF (EdgeP%Number>0) THEN
    EdgeNumber=EdgeP%Number
  ELSE  
    EdgeNumber=EdgeP%Edge%Number
  END IF  
END FUNCTION EdgeNumber

SUBROUTINE WriteCellCut(CellP,OutputUnit)
  TYPE(CellP_T) :: CellP
  INTEGER :: OutputUnit

  TYPE(Cell_T), POINTER :: Cell
  TYPE(Face_T), POINTER :: Face
  INTEGER :: iList1,List1Face(8)
  INTEGER :: iList2,List2Face(8)

  Cell=>CellP%Cell
  IF (ASSOCIATED(Cell).AND.Cell%NumberFace>0) THEN
    iList1=1
    iList2=1
    List1Face(iList1)=Cell%NumberFace
    List2Face(iList2)=Cell%NumberFace
    CALL TestFace(Cell%Face1)
    CALL TestFace(Cell%Face2)
    CALL TestFace(Cell%Face3)
    CALL TestFace(Cell%Face4)
    CALL TestFace(Cell%Face5)
    CALL TestFace(Cell%Face6)
    WRITE(OutputUnit,*) 'cc  ',Cell%Number &
                              ,iList1,0
    WRITE(OutputUnit,*) List1Face(1:iList1)                           
    WRITE(OutputUnit,*) 'cc  ',Cell%Number+1 &
                              ,iList2,-1
    WRITE(OutputUnit,*) List2Face(1:iList2)                           
  ELSE  
    WRITE(OutputUnit,*) 'ci  ',Cell%Number &
                              ,6,-1
    WRITE(OutputUnit,*)        Cell%Face1%Number      &
                              ,Cell%Face2%Number      &
                              ,Cell%Face3%Number      &
                              ,Cell%Face4%Number      &
                              ,Cell%Face5%Number      &
                              ,Cell%Face6%Number      
   END IF   
CONTAINS
SUBROUTINE TestFace(Face)
  TYPE(Face_T) :: Face

  IF (Face%NumberEdge>0) THEN
    iList1=iList1+1
    List1Face(iList1)=Face%Number
    iList2=iList2+1
    List2Face(iList2)=Face%Number+1
  ELSE  
    IF (Face%Edge1%Vert1%in_out>=0) THEN
      iList1=iList1+1
      List1Face(iList1)=Face%Number
    ELSE  
      iList2=iList2+1
      List2Face(iList2)=Face%Number
    END IF  
  END IF  
END SUBROUTINE TestFace  
END SUBROUTINE WriteCellCut

SUBROUTINE WriteVertexCut(EdgeP,OutputUnit)
  TYPE(EdgeP_T) :: EdgeP
  INTEGER :: OutputUnit

  TYPE(Edge_T), POINTER :: Edge

  Edge=>EdgeP%Edge
  IF (ASSOCIATED(Edge)) THEN
    SELECT CASE (Edge%EdgeType)
      CASE ('X')
        IF (ASSOCIATED(Edge%Vert1%VertX)) THEN
          WRITE(OutputUnit,*) 'nXc ',Edge%NumberNode,Edge%Vert1%VertX%Point
        END IF  
      CASE ('Y')
        IF (ASSOCIATED(Edge%Vert1%VertY)) THEN
          WRITE(OutputUnit,*) 'nYc ',Edge%NumberNode,Edge%Vert1%VertY%Point
        END IF  
      CASE ('Z')
        IF (ASSOCIATED(Edge%Vert1%VertZ)) THEN
          WRITE(OutputUnit,*) 'nZc ',Edge%NumberNode,Edge%Vert1%VertZ%Point
        END IF  
    END SELECT  
  END IF
END SUBROUTINE WriteVertexCut

SUBROUTINE WriteCutFace(FaceP,OutputUnit,Type)
  TYPE(FaceP_T) :: FaceP
  INTEGER :: OutputUnit
  CHARACTER*2 :: Type

  TYPE(Face_T), POINTER :: Face
  INTEGER :: iList1,List1Edge(5)
  INTEGER :: iList2,List2Edge(5)

  Face=>FaceP%Face
  IF (ASSOCIATED(Face).AND.Face%NumberEdge>0) THEN
    iList1=1
    iList2=1
    List1Edge(iList1)=Face%NumberEdge
    List2Edge(iList2)=Face%NumberEdge
    IF (Face%Edge1%NumberNode>0) THEN
      iList1=iList1+1
      iList2=iList2+1
      List1Edge(iList1)=Face%Edge1%Number 
      List2Edge(iList2)=Face%Edge1%Number+1 
    ELSE  
      IF (Face%Edge1%Vert1%in_out>=0) THEN
        iList1=iList1+1
        List1Edge(iList1)=Face%Edge1%Number 
      ELSE  
        iList2=iList2+1
        List2Edge(iList2)=Face%Edge1%Number 
      END IF  
    END IF  
    IF (Face%Edge2%NumberNode>0) THEN
      iList1=iList1+1
      iList2=iList2+1
      List1Edge(iList1)=Face%Edge2%Number 
      List2Edge(iList2)=Face%Edge2%Number+1 
    ELSE  
      IF (Face%Edge2%Vert1%in_out>=0) THEN
        iList1=iList1+1
        List1Edge(iList1)=Face%Edge2%Number 
      ELSE  
        iList2=iList2+1
        List2Edge(iList2)=Face%Edge2%Number 
      END IF  
    END IF  
    IF (Face%Edge3%NumberNode>0) THEN
      iList1=iList1+1
      iList2=iList2+1
      List1Edge(iList1)=Face%Edge3%Number 
      List2Edge(iList2)=Face%Edge3%Number+1 
    ELSE  
      IF (Face%Edge3%Vert1%in_out>=0) THEN
        iList1=iList1+1
        List1Edge(iList1)=Face%Edge3%Number 
      ELSE  
        iList2=iList2+1
        List2Edge(iList2)=Face%Edge3%Number 
      END IF  
    END IF  
    IF (Face%Edge4%NumberNode>0) THEN
      iList1=iList1+1
      iList2=iList2+1
      List1Edge(iList1)=Face%Edge4%Number 
      List2Edge(iList2)=Face%Edge4%Number+1 
    ELSE  
      IF (Face%Edge4%Vert1%in_out>=0) THEN
        iList1=iList1+1
        List1Edge(iList1)=Face%Edge4%Number 
      ELSE  
        iList2=iList2+1
        List2Edge(iList2)=Face%Edge4%Number 
      END IF  
    END IF  
    WRITE(OutputUnit,*) 'f'//Type//'c',Face%Number,iList1,1
    WRITE(OutputUnit,*) List1Edge(1:iList1)
    WRITE(OutputUnit,*) 'f'//Type//'c',Face%Number+1,iList2,1
    WRITE(OutputUnit,*) List2Edge(1:iList2)
  ELSE
    WRITE(OutputUnit,*) 'f'//Type//' ',Face%Number &
                              ,4,1
    WRITE(OutputUnit,*)        Face%Edge1%Number &
                              ,Face%Edge2%Number &
                              ,Face%Edge3%Number &
                              ,Face%Edge4%Number 
  END IF
END SUBROUTINE WriteCutFace

SUBROUTINE WriteFaceCut(CellP,OutputUnit)
  TYPE(CellP_T) :: CellP
  INTEGER :: OutputUnit

  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: iList,ListEdge(6)

  Cell=>CellP%Cell
  IF (ASSOCIATED(Cell).AND.Cell%NumberFace>0) THEN
    iList=0
    IF (Cell%Face1%NumberEdge>0) THEN
      iList=iList+1
      ListEdge(iList)=Cell%Face1%NumberEdge
    END IF  
    IF (Cell%Face2%NumberEdge>0) THEN
      iList=iList+1
      ListEdge(iList)=Cell%Face2%NumberEdge
    END IF  
    IF (Cell%Face3%NumberEdge>0) THEN
      iList=iList+1
      ListEdge(iList)=Cell%Face3%NumberEdge
    END IF  
    IF (Cell%Face4%NumberEdge>0) THEN
      iList=iList+1
      ListEdge(iList)=Cell%Face4%NumberEdge
    END IF  
    IF (Cell%Face5%NumberEdge>0) THEN
      iList=iList+1
      ListEdge(iList)=Cell%Face5%NumberEdge
    END IF  
    IF (Cell%Face6%NumberEdge>0) THEN
      iList=iList+1
      ListEdge(iList)=Cell%Face6%NumberEdge
    END IF  
    IF (iList>=3) THEN
      WRITE(OutputUnit,*) 'fc  ',Cell%NumberFace,iList,1
      WRITE(OutputUnit,*) ListEdge(1:iList)
    ELSE
      WRITE(*,*) 'ListEdge',ListEdge(1:iList)
      STOP
    END IF  
  END IF  
END SUBROUTINE WriteFaceCut

SUBROUTINE WriteEdgeCut(FaceP,OutputUnit)
  TYPE(FaceP_T) :: FaceP
  INTEGER :: OutputUnit

  TYPE(Face_T), POINTER :: Face
  INTEGER :: iList,ListNode(4)
  Face=>FaceP%Face
  IF (ASSOCIATED(Face).AND.Face%NumberEdge>0) THEN
    iList=0
    IF (Face%Edge1%NumberNode>0) THEN
      iList=iList+1
      ListNode(iList)=Face%Edge1%NumberNode 
    END IF  
    IF (Face%Edge2%NumberNode>0) THEN
      iList=iList+1
      ListNode(iList)=Face%Edge2%NumberNode 
    END IF  
    IF (Face%Edge3%NumberNode>0) THEN
      iList=iList+1
      ListNode(iList)=Face%Edge3%NumberNode 
    END IF  
    IF (Face%Edge4%NumberNode>0) THEN
      iList=iList+1
      ListNode(iList)=Face%Edge4%NumberNode 
    END IF  
    IF (iList==2) THEN
      WRITE(OutputUnit,*) 'ec  ',Face%NumberEdge,ListNode(1:2)
    ELSE
      WRITE(OutputUnit,*) 'ec0 ',Face%NumberEdge,ListNode(1:2)
      WRITE(*,*) 'ListNode',iList,ListNode(1:iList)
      WRITE(*,*) 'Face%Number',Face%Number,'Face%NumberEdge',Face%NumberEdge
      WRITE(*,*) 'Face%Edge1%NumberNode',Face%Edge1%NumberNode,Face%Edge1%Vert1%in_out,Face%Edge1%Vert2%in_out
      WRITE(*,*) Face%Edge1%Vert1%ix,Face%Edge1%Vert1%iy,Face%Edge1%Vert1%iz
      WRITE(*,*) Face%Edge1%Vert2%ix,Face%Edge1%Vert2%iy,Face%Edge1%Vert2%iz
      WRITE(*,*) 'Face%Edge2%NumberNode',Face%Edge2%NumberNode,Face%Edge2%Vert1%in_out,Face%Edge2%Vert2%in_out
      WRITE(*,*) Face%Edge2%Vert1%ix,Face%Edge2%Vert1%iy,Face%Edge2%Vert1%iz
      WRITE(*,*) Face%Edge2%Vert2%ix,Face%Edge2%Vert2%iy,Face%Edge2%Vert2%iz
      WRITE(*,*) 'Face%Edge3%NumberNode',Face%Edge3%NumberNode,Face%Edge3%Vert1%in_out,Face%Edge3%Vert2%in_out
      WRITE(*,*) Face%Edge3%Vert1%ix,Face%Edge3%Vert1%iy,Face%Edge3%Vert1%iz
      WRITE(*,*) Face%Edge3%Vert2%ix,Face%Edge3%Vert2%iy,Face%Edge3%Vert2%iz
      WRITE(*,*) 'Face%Edge4%NumberNode',Face%Edge4%NumberNode,Face%Edge4%Vert1%in_out,Face%Edge4%Vert2%in_out
      WRITE(*,*) Face%Edge4%Vert1%ix,Face%Edge4%Vert1%iy,Face%Edge4%Vert1%iz
      WRITE(*,*) Face%Edge4%Vert2%ix,Face%Edge4%Vert2%iy,Face%Edge4%Vert2%iz
      STOP
    END IF  
  END IF
END SUBROUTINE WriteEdgeCut


SUBROUTINE CountEdge(EdgeP,NumberEdges,NumberNodes)
  TYPE(EdgeP_T) :: EdgeP
  INTEGER :: NumberEdges
  INTEGER :: NumberNodes

  TYPE(Edge_T), POINTER :: Edge

  Edge=>EdgeP%Edge
  IF (ASSOCIATED(Edge)) THEN
    SELECT CASE (Edge%EdgeType)
      CASE ('X')
        IF (ASSOCIATED(Edge%Vert1%VertX)) THEN
          NumberEdges=NumberEdges+1
          Edge%Number=NumberEdges
          NumberNodes=NumberNodes+1
          Edge%NumberNode=NumberNodes
          NumberEdges=NumberEdges+1
        ELSE
          NumberEdges=NumberEdges+1
          Edge%Number=NumberEdges
        END IF  
      CASE ('Y')
        IF (ASSOCIATED(Edge%Vert1%VertY)) THEN
          NumberEdges=NumberEdges+1
          Edge%Number=NumberEdges
          NumberNodes=NumberNodes+1
          Edge%NumberNode=NumberNodes
          NumberEdges=NumberEdges+1
        ELSE
          NumberEdges=NumberEdges+1
          Edge%Number=NumberEdges
        END IF  
      CASE ('Z')
        IF (ASSOCIATED(Edge%Vert1%VertZ)) THEN
          NumberEdges=NumberEdges+1
          Edge%Number=NumberEdges
          NumberNodes=NumberNodes+1
          Edge%NumberNode=NumberNodes
          NumberEdges=NumberEdges+1
        ELSE
          NumberEdges=NumberEdges+1
          Edge%Number=NumberEdges
        END IF  
    END SELECT  
  ELSE
    NumberEdges=NumberEdges+1
    EdgeP%Number=NumberEdges
  END IF
END SUBROUTINE CountEdge

SUBROUTINE CountFace(FaceP,NumberFaces,NumberEdges)
  TYPE(FaceP_T) :: FaceP
  INTEGER :: NumberFaces
  INTEGER :: NumberEdges

  INTEGER :: NumberEdgeMax
  TYPE(Face_T), POINTER :: Face

  Face=>FaceP%Face
  IF (ASSOCIATED(Face)) THEN
    NumberEdgeMax=MAX(Face%Edge1%NumberNode &
                     ,Face%Edge2%NumberNode &
                     ,Face%Edge3%NumberNode &
                     ,Face%Edge4%NumberNode)
    IF (NumberEdgeMax>0) THEN
      NumberFaces=NumberFaces+1
      Face%Number=NumberFaces
      NumberEdges=NumberEdges+1
      Face%NumberEdge=NumberEdges
      NumberFaces=NumberFaces+1
    ELSE
      NumberFaces=NumberFaces+1
      Face%Number=NumberFaces
    END IF  
  ELSE
    NumberFaces=NumberFaces+1
    FaceP%Number=NumberFaces
  END IF
END SUBROUTINE CountFace

SUBROUTINE CountCell(CellP,NumberCells,NumberFaces)
  TYPE(CellP_T) :: CellP
  INTEGER :: NumberCells
  INTEGER :: NumberFaces

  INTEGER :: NumberFaceMax
  TYPE(Cell_T), POINTER :: Cell

  Cell=>CellP%Cell
  IF (ASSOCIATED(Cell)) THEN
    NumberFaceMax=MAX(Cell%Face1%NumberEdge &
                     ,Cell%Face2%NumberEdge &
                     ,Cell%Face3%NumberEdge &
                     ,Cell%Face4%NumberEdge &
                     ,Cell%Face5%NumberEdge &
                     ,Cell%Face6%NumberEdge)
    IF (NumberFaceMax>0) THEN
      NumberCells=NumberCells+1
      Cell%Number=NumberCells
      NumberFaces=NumberFaces+1
      Cell%NumberFace=NumberFaces
      NumberCells=NumberCells+1
    ELSE
      NumberCells=NumberCells+1
      Cell%Number=NumberCells
    END IF  
  ELSE
    NumberCells=NumberCells+1
    CellP%Number=NumberCells
  END IF
END SUBROUTINE CountCell

END MODULE OutputUnstructuredGrid_Mod

