MODULE OutputOutGmvGNeu_Mod

  USE Function_Mod
  USE Floor_Mod
  USE IOControl_Mod
  USE Parametric_Mod
  USE GridNeu_Mod

  IMPLICIT NONE

CONTAINS

SUBROUTINE CountAllCellsGMV(NumberCells)
  INTEGER :: NumberCells

  INTEGER :: ib,ix,iy,iz

  NumberCells=0
  IF(out_area=='y') THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO iz=v_z0+1,v_z1
          DO iy=v_y0+1,v_y1
            DO ix=v_x0+1,v_x1
              CALL CountCellGMV(Cells(ix,iy,iz)%Cell,ix,iy,iz,NumberCells)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            CALL CountCellGMV(Cells(ix,iy,iz)%Cell,ix,iy,iz,NumberCells)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
END SUBROUTINE CountAllCellsGMV

SUBROUTINE CountCellGMV(Cell,i,j,k,NumberCells)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k
  INTEGER :: NumberCells

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out

  IF (ASSOCIATED(Cell)) THEN
    IF ((Cell%in_out>-8.AND.Cell%Vol>0.0d0) &
         .OR.Cell%in_out==8) THEN
      NumberCells=NumberCells+1
    END IF
  ELSE
    in_out=Vertices(i-1,j-1,k-1)%Vertex%in_out &
          +Vertices(i,j-1,k-1)%Vertex%in_out &
          +Vertices(i-1,j,k-1)%Vertex%in_out &
          +Vertices(i-1,j-1,k)%Vertex%in_out &
          +Vertices(i-1,j,k)%Vertex%in_out &
          +Vertices(i,j-1,k)%Vertex%in_out &
          +Vertices(i,j,k-1)%Vertex%in_out &
          +Vertices(i,j,k)%Vertex%in_out
    IF (in_out>=0) THEN
      NumberCells=NumberCells+1
    END IF
  END IF
END SUBROUTINE CountCellGMV

SUBROUTINE CountAllCellsOroGMV(NumberCells)
  INTEGER :: NumberCells

  INTEGER :: ib,ix,iy,iz

  NumberCells=0
  IF(out_area=='y') THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO iz=v_z0+1,v_z1
          DO iy=v_y0+1,v_y1
            DO ix=v_x0+1,v_x1
              CALL CountCellOroGMV(Cells(ix,iy,iz)%Cell,ix,iy,iz,NumberCells)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            CALL CountCellOroGMV(Cells(ix,iy,iz)%Cell,ix,iy,iz,NumberCells)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
END SUBROUTINE CountAllCellsOroGMV

SUBROUTINE CountCellOroGMV(Cell,i,j,k,NumberCells)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k
  INTEGER :: NumberCells

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out

  IF (ASSOCIATED(Cell)) THEN
    NumberCells=NumberCells+1
  END IF
END SUBROUTINE CountCellOroGMV

SUBROUTINE WriteAllCellsBinGMV(FileName)
  CHARACTER*50 :: FileName

  INTEGER :: i,ib,ix,iy,iz,lvertout
  TYPE (Vertex_T) :: V(1:2)
  REAL(8) :: sfac
  INTEGER :: NumberCellsGMV
  REAL(8) :: tempw1,tempw2

  CALL CountAllCellsGMV(NumberCellsGMV)
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.out.gmvG',STATUS='UNKNOWN'&
    &         ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)
  CALL Display_OutGMVBlk(FileName) 
  nRec=0
  !.. WRITE(10) gmvinput,ascii/iecxi4r4
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(5:8)
  !.. WRITE(10) nodes,nr_view_out
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) nr_view_out
  !WRITE(10,*) (xParametricOut(VertOutView(i)),i=1,nr_view_out))
  DO i=1,nr_view_out
    nRec=nRec+1
    tempw=xParametricOut(VertOutView(i)%Vertex)
!   tempw1=yParametricOut(VertOutView(i)%Vertex)
!   tempw2=zParametricOut(VertOutView(i)%Vertex)
    WRITE(10,REC=nRec) tempw
!   WRITE(200,*) i,tempw,tempw1,tempw2
  END DO
  !.. WRITE(10,*) VertOutView%Point%y
  DO i=1,nr_view_out
    nRec=nRec+1
    tempw=yParametricOut(VertOutView(i)%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10,*) VertOutView%Point%z
  DO i=1,nr_view_out
    nRec=nRec+1
    tempw=zParametricOut(VertOutView(i)%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10) cells,nr_viewcells
  nRec=nRec+1
  WRITE(10,REC=nRec) cellsGMV(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) cellsGMV(5:8)
  nRec=nRec+1

  IF(out_area=='y') THEN
    WRITE(10,REC=nRec) NumberCellsGMV
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO iz=v_z0+1,v_z1
          DO iy=v_y0+1,v_y1
            DO ix=v_x0+1,v_x1
              CALL WriteCellGMVBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    WRITE(10,REC=nRec) NumberCellsGMV
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            CALL WriteCellGMVBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
  ! --------------------------
  ! ----- Polygone -----------
  ! CutPlane,Bloecke,Haeuser
  !---------------------------
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(5:8)
  ! CutPlane:
  !----------
  IF(out_area=='y') THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO iz=v_z0+1,v_z1
          DO iy=v_y0+1,v_y1
            DO ix=v_x0+1,v_x1
              CALL WriteCutPlanePolyGMVBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            CALL WriteCutPlanePolyGMVBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
  !LandClass-Plane
  !---------------
  IF(nr_landdef>0) THEN
    IF(out_area=='y') THEN
      DO ib=1,nb
        CALL Set(Floor(ib))
        CALL Search_OutArea
        IF(view_cell=='y') THEN
          DO iz=v_z0+1,v_z1
            DO iy=v_y0+1,v_y1
              DO ix=v_x0+1,v_x1
! OSSI          CALL WriteLandKlPlanePolyGMVBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
              END DO
            END DO
          END DO
       END IF
      END DO  !ib
    ELSE
      DO ib=1,nb
        CALL Set(Floor(ib))
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1
!  OSSI       CALL WriteLandKlPlanePolyGMVBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
            END DO
          END DO
        END DO
      END DO  !ib
    END IF
  END IF  !if(nr_landdef)
  ! Soil-layering:
  !-----------
  !dzi_soil auf Grundgitter bezogen eingestellt
  IF(nr_sb>0) THEN
    V(1:2)%Point%x=Domain%xP(0)
    V(1:2)%Point%y=Domain%yP(0)
    V(1)%Point%z=Domain%zP(0)
    V(2)%Point%z=Domain%zP(1)
    dzi_soil=zParametricOut(V(2))-zParametricOut(V(1))
    sfac=(Domain%zP(nz)-Domain%zP(0))/(dzi_soil*100)
    dzi_soil=dzi_soil*sfac !skaliere dzi_soil auf 1% von Domain%zP(nz)
    dzi_soil=dzi_soil*ScaleSoil   !extern hinzu, 1% Standard
    IF(out_wahlgrid=='C') THEN
      IF(out_area=='y') THEN
        DO ib=1,nb
          CALL Set(Floor(ib))
          CALL Search_OutArea
          IF(view_cell=='y') THEN
            DO iz=v_z0+1,v_z1
              DO iy=v_y0+1,v_y1
                DO ix=v_x0+1,v_x1
!  OSSI            CALL WritePolySoilLayerGMVAsciiBin(Cells(ix,iy,iz)%Cell,ix,iy,iz)
                END DO
              END DO
            END DO
          END IF
        END DO  !ib
      ELSE
        DO ib=1,nb
          CALL Set(Floor(ib))
          DO iz=iz0+1,iz1
            DO iy=iy0+1,iy1
              DO ix=ix0+1,ix1
!  OSSI         CALL WritePolySoilLayerGMVAsciiBin(Cells(ix,iy,iz)%Cell,ix,iy,iz)
              END DO
            END DO
          END DO
        END DO  !ib
      END IF
    END IF
  END IF
  ! Bloecke:
  !----------
  DO ib=1,nb
     CALL Set(Floor(ib))
     IF(out_area=='y') THEN
       ! Output-View
       CALL Search_OutArea
       IF(view_cell=='y') THEN
          CALL WriteBlockGMVBinary
       END IF
     ELSE  
       ! Output-Global
       CALL WriteBlockGMVBinary
     END IF
   END DO  !ib

  ! Haeuser:
  !----------
  CALL WriteAllHausMultiColorBinGMV
  !.. WRITE(10) endpoly
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(5:8)
 
  !.. WRITE(10) endgmv
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(5:8)

  CLOSE(10)

END SUBROUTINE WriteAllCellsBinGMV

SUBROUTINE WriteAllCellsBinary(FileName)
  CHARACTER*50 :: FileName

  INTEGER :: i,ib,ix,iy,iz,lvertout
  TYPE (Vertex_T) :: V(1:2)
  REAL(8) :: sfac
  INTEGER :: NumberCellsGMV
  REAL(8) :: Temp

  CALL CountAllCellsGMV(NumberCellsGMV)
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.gmv',STATUS='UNKNOWN'&
    &         ,FORM='UNFORMATTED')
  WRITE(10) gmvinput(1:4)
  WRITE(10) gmvinput(5:8)
  WRITE(10) iecxi4r4(1:4)
  WRITE(10) iecxi4r4(5:8)
  WRITE(10) nodes(1:4)
  WRITE(10) nodes(5:8)
  WRITE(10) nr_view_out
  DO i=1,nr_view_out
    temp=xParametricOut(VertOutView(i)%Vertex)
    WRITE(10) temp
!   WRITE(200,*) i,temp
  END DO
  DO i=1,nr_view_out
    temp=yParametricOut(VertOutView(i)%Vertex)
    WRITE(10) temp
  END DO
  DO i=1,nr_view_out
    temp=zParametricOut(VertOutView(i)%Vertex)
    WRITE(10) temp
  END DO
  WRITE(10) cellsGMV(1:4)
  WRITE(10) cellsGMV(5:8)
  IF(out_area=='y') THEN
    WRITE(10) NumberCellsGMV
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO iz=v_z0+1,v_z1
          DO iy=v_y0+1,v_y1
            DO ix=v_x0+1,v_x1
              CALL WriteCellBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    WRITE(10) NumberCellsGMV
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            CALL WriteCellBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
  WRITE(10) endgmv(1:4)
  WRITE(10) endgmv(5:8)
  CLOSE(10)

END SUBROUTINE WriteAllCellsBinary

SUBROUTINE WriteCellGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out

  IF (ASSOCIATED(Cell)) THEN
    IF ((Cell%in_out>-8.AND.Cell%Vol>0.0d0) &
         .OR.Cell%in_out==8) THEN
      IF (Cell%vc==0) THEN
       ! hex
        nRec=nRec+1
        WRITE(10,REC=nRec) hex(1:4)
        nRec=nRec+1
        WRITE(10,REC=nRec) hex(5:8)
        nRec=nRec+1
        WRITE(10,REC=nRec) 8
        DO iVert=1,4
          nRec=nRec+1
          WRITE(10,REC=nRec) VertViewScale(Cell%Face2%VertexList(iVert))
        END DO
        DO iVert=1,4
          nRec=nRec+1
          WRITE(10,REC=nRec) VertViewScale(Cell%Face1%VertexList(iVert))
        END DO
      ELSE
        ! general
        nFaces=1
        IF (Cell%Face1%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face1%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face2%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face2%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face3%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face3%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face4%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face4%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face5%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face5%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face6%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face6%NumberVert
          nFaces=nFaces+1
        END IF
        nVerts(nFaces)=Cell%vc

        !.. WRITE(10) general,nfaces
        nRec=nRec+1
        WRITE(10,REC=nRec) general(1:4)
        nRec=nRec+1
        WRITE(10,REC=nRec) general(5:8)
        nRec=nRec+1
!       WRITE(200,*) general(1:8),i,j,k
!       WRITE(200,*) nfaces
        WRITE(10,REC=nRec) nfaces
        !.. WRITE(10) (nVerts(i),i=1,nFaces)
!       WRITE(200,*) nVerts(1:nFaces)
        DO iVert=1,nFaces
          nRec=nRec+1
          WRITE(10,REC=nRec) nVerts(iVert)
        END DO
        IF (Cell%Face1%NumberVert>2) THEN
!         WRITE(200,*) 'Face1:',Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
          DO iVert=1,Cell%Face1%NumberVert
            nRec=nRec+1
            WRITE(10,REC=nRec)  VertViewScale(Cell%Face1%VertexList(iVert))
            ListVert(iVert)=VertViewScale(Cell%Face1%VertexList(iVert))
          END DO
!         WRITE(200,*) 'Face1:',ListVert(1:Cell%Face1%NumberVert)
        END IF
        IF (Cell%Face2%NumberVert>2) THEN
!         WRITE(200,*) 'Face2:',Cell%Face2%VertexList(1:Cell%Face2%NumberVert)
          DO iVert=1,Cell%Face2%NumberVert
            nRec=nRec+1
            WRITE(10,REC=nRec) VertViewScale(Cell%Face2%VertexList(iVert))
            ListVert(iVert)=VertViewScale(Cell%Face2%VertexList(iVert))
          END DO
!         WRITE(200,*) 'Face2:',ListVert(1:Cell%Face2%NumberVert)
        END IF
        IF (Cell%Face3%NumberVert>2) THEN
!         WRITE(200,*) 'Face3:',Cell%Face3%VertexList(1:Cell%Face3%NumberVert)
          DO iVert=1,Cell%Face3%NumberVert
!         DO iVert=Cell%Face3%NumberVert,1,-1
            nRec=nRec+1
            WRITE(10,REC=nRec) VertViewScale(Cell%Face3%VertexList(iVert))
            ListVert(iVert)=VertViewScale(Cell%Face3%VertexList(iVert))
          END DO
!         WRITE(200,*) 'Face3:',ListVert(1:Cell%Face3%NumberVert)
        END IF
        IF (Cell%Face4%NumberVert>2) THEN
!         WRITE(200,*) 'Face4:',Cell%Face4%VertexList(1:Cell%Face4%NumberVert)
          DO iVert=1,Cell%Face4%NumberVert
            nRec=nRec+1
            WRITE(10,REC=nRec) VertViewScale(Cell%Face4%VertexList(iVert))
            ListVert(iVert)=VertViewScale(Cell%Face4%VertexList(iVert))
          END DO
!         WRITE(200,*) 'Face4:',ListVert(1:Cell%Face4%NumberVert)
        END IF
        IF (Cell%Face5%NumberVert>2) THEN
!         WRITE(200,*) 'Face5:',Cell%Face5%VertexList(1:Cell%Face5%NumberVert)
          DO iVert=1,Cell%Face5%NumberVert
!         DO iVert=Cell%Face5%NumberVert,1,-1
            nRec=nRec+1
            WRITE(10,REC=nRec) VertViewScale(Cell%Face5%VertexList(iVert))
            ListVert(iVert)=VertViewScale(Cell%Face5%VertexList(iVert))
          END DO
!         WRITE(200,*) 'Face5:',ListVert(1:Cell%Face5%NumberVert)
        END IF
        IF (Cell%Face6%NumberVert>2) THEN
!         WRITE(200,*) 'Face6:',Cell%Face6%VertexList(1:Cell%Face6%NumberVert)
          DO iVert=1,Cell%Face6%NumberVert
            nRec=nRec+1
            WRITE(10,REC=nRec) VertViewScale(Cell%Face6%VertexList(iVert))
            ListVert(iVert)=VertViewScale(Cell%Face6%VertexList(iVert))
          END DO
!         WRITE(200,*) 'Face6:',ListVert(1:Cell%Face6%NumberVert)
        END IF
!       WRITE(200,*) 'Cut 6:',Cell%VertCut(1:Cell%vc)
        DO iVert=1,Cell%vc
          nRec=nRec+1
          WRITE(10,REC=nRec) VertViewScale(Cell%VertCut(iVert))
          ListVert(iVert)=VertViewScale(Cell%VertCut(iVert))
        END DO
!       WRITE(200,*) 'Cut 6:',ListVert(1:Cell%vc)
      END IF
    END IF
  ELSE
    in_out=Vertices(i-1,j-1,k-1)%Vertex%in_out &
          +Vertices(i,j-1,k-1)%Vertex%in_out &
          +Vertices(i-1,j,k-1)%Vertex%in_out &
          +Vertices(i-1,j-1,k)%Vertex%in_out &
          +Vertices(i-1,j,k)%Vertex%in_out &
          +Vertices(i,j-1,k)%Vertex%in_out &
          +Vertices(i,j,k-1)%Vertex%in_out &
          +Vertices(i,j,k)%Vertex%in_out
    IF (in_out>=0) THEN
     ! hex
      !          4--------3    !Point-Folge Output-Hex-GMV-Vorlage-Def
      !         /|       /|
      !        1--------2 |
      !        | 8------|-7
      !        |/       |/
      !        5--------6
      nRec=nRec+1
      WRITE(10,REC=nRec) hex(1:4)
      nRec=nRec+1
      WRITE(10,REC=nRec) hex(5:8)
      nRec=nRec+1
      WRITE(10,REC=nRec) 8
      nRec=nRec+1
      WRITE(10,REC=nRec) VertViewScale(Vertices(i-1,j-1,k)%Vertex%nrP)
      nRec=nRec+1
      WRITE(10,REC=nRec) VertViewScale(Vertices(i,j-1,k)%Vertex%nrP)
      nRec=nRec+1
      WRITE(10,REC=nRec) VertViewScale(Vertices(i,j,k)%Vertex%nrP)
      nRec=nRec+1
      WRITE(10,REC=nRec) VertViewScale(Vertices(i-1,j,k)%Vertex%nrP)
      nRec=nRec+1
      WRITE(10,REC=nRec) VertViewScale(Vertices(i-1,j-1,k-1)%Vertex%nrP)
      nRec=nRec+1
      WRITE(10,REC=nRec) VertViewScale(Vertices(i,j-1,k-1)%Vertex%nrP)
      nRec=nRec+1
      WRITE(10,REC=nRec) VertViewScale(Vertices(i,j,k-1)%Vertex%nrP)
      nRec=nRec+1
      WRITE(10,REC=nRec) VertViewScale(Vertices(i-1,j,k-1)%Vertex%nrP)
    END IF
  END IF
END SUBROUTINE WriteCellGMVBinary

SUBROUTINE WriteCellBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out

  IF (ASSOCIATED(Cell)) THEN
    IF ((Cell%in_out>-8.AND.Cell%Vol>0.0d0) &
         .OR.Cell%in_out==8) THEN
      IF (Cell%vc==0) THEN
       ! hex
        WRITE(10) hex(1:4)
        WRITE(10) hex(5:8)
        WRITE(10) 8
        DO iVert=1,4
          WRITE(10) VertViewScale(Cell%Face2%VertexList(iVert))
        END DO
        DO iVert=1,4
          WRITE(10) VertViewScale(Cell%Face1%VertexList(iVert))
        END DO
      ELSE
        IF (Cell%Vol<1.d-10) THEN
          WRITE(*,*) Cell%Face1%NumberVert
          WRITE(*,*) Cell%Face1%VertexList
          WRITE(*,*) Cell%Face2%NumberVert
          WRITE(*,*) Cell%Face2%VertexList
          WRITE(*,*) Cell%Face3%NumberVert
          WRITE(*,*) Cell%Face3%VertexList
          WRITE(*,*) Cell%Face4%NumberVert
          WRITE(*,*) Cell%Face4%VertexList
          WRITE(*,*) Cell%Face5%NumberVert
          WRITE(*,*) Cell%Face5%VertexList
          WRITE(*,*) Cell%Face6%NumberVert
          WRITE(*,*) Cell%Face6%VertexList
        END IF  
        ! general
        nFaces=1
        IF (Cell%Face1%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face1%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face2%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face2%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face3%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face3%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face4%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face4%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face5%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face5%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face6%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face6%NumberVert
          nFaces=nFaces+1
        END IF
        nVerts(nFaces)=Cell%vc

        WRITE(10) general(1:4)
        WRITE(10) general(5:8)
        WRITE(10) nfaces
        DO iVert=1,nFaces
          WRITE(10) nVerts(iVert)
        END DO
        IF (Cell%Face1%NumberVert>2) THEN
          DO iVert=1,Cell%Face1%NumberVert
            WRITE(10)  VertViewScale(Cell%Face1%VertexList(iVert))
          END DO
        END IF
        IF (Cell%Face2%NumberVert>2) THEN
          DO iVert=1,Cell%Face2%NumberVert
            WRITE(10) VertViewScale(Cell%Face2%VertexList(iVert))
          END DO
        END IF
        IF (Cell%Face3%NumberVert>2) THEN
          DO iVert=1,Cell%Face3%NumberVert
             WRITE(10) VertViewScale(Cell%Face3%VertexList(iVert))
          END DO
        END IF
        IF (Cell%Face4%NumberVert>2) THEN
          DO iVert=1,Cell%Face4%NumberVert
             WRITE(10) VertViewScale(Cell%Face4%VertexList(iVert))
          END DO
        END IF
        IF (Cell%Face5%NumberVert>2) THEN
          DO iVert=1,Cell%Face5%NumberVert
             WRITE(10) VertViewScale(Cell%Face5%VertexList(iVert))
          END DO
        END IF
        IF (Cell%Face6%NumberVert>2) THEN
          DO iVert=1,Cell%Face6%NumberVert
             WRITE(10) VertViewScale(Cell%Face6%VertexList(iVert))
          END DO
        END IF
        DO iVert=1,Cell%vc
          WRITE(10) VertViewScale(Cell%VertCut(iVert))
        END DO
      END IF
    END IF
  ELSE
    in_out=Vertices(i-1,j-1,k-1)%Vertex%in_out &
          +Vertices(i,j-1,k-1)%Vertex%in_out &
          +Vertices(i-1,j,k-1)%Vertex%in_out &
          +Vertices(i-1,j-1,k)%Vertex%in_out &
          +Vertices(i-1,j,k)%Vertex%in_out &
          +Vertices(i,j-1,k)%Vertex%in_out &
          +Vertices(i,j,k-1)%Vertex%in_out &
          +Vertices(i,j,k)%Vertex%in_out
    IF (in_out>=0) THEN
     ! hex
      !          4--------3    !Point-Folge Output-Hex-GMV-Vorlage-Def
      !         /|       /|
      !        1--------2 |
      !        | 8------|-7
      !        |/       |/
      !        5--------6
      WRITE(10) hex(1:4)
      WRITE(10) hex(5:8)
      WRITE(10) 8
      WRITE(10) VertViewScale(Vertices(i-1,j-1,k)%Vertex%nrP)
      WRITE(10) VertViewScale(Vertices(i,j-1,k)%Vertex%nrP)
      WRITE(10) VertViewScale(Vertices(i,j,k)%Vertex%nrP)
      WRITE(10) VertViewScale(Vertices(i-1,j,k)%Vertex%nrP)
      WRITE(10) VertViewScale(Vertices(i-1,j-1,k-1)%Vertex%nrP)
      WRITE(10) VertViewScale(Vertices(i,j-1,k-1)%Vertex%nrP)
      WRITE(10) VertViewScale(Vertices(i,j,k-1)%Vertex%nrP)
      WRITE(10) VertViewScale(Vertices(i-1,j,k-1)%Vertex%nrP)
    END IF
  END IF
END SUBROUTINE WriteCellBinary

SUBROUTINE WriteAllCellsOroBinGMV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: i,ib,ix,iy,iz
  INTEGER :: NumberCellsGMV
  nr_cell_vc_oro=0
  nr_cell_novc_oro=0
  nr_insidehex=0

  CALL CountAllCellsOroGMV(NumberCellsGMV)
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.Oro.out.gmvG',STATUS='UNKNOWN'&
    &         ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)

  nRec=0
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) nr_inside
  DO i=1,nr_inside
    nRec=nRec+1
    tempw=xParametricOut(VertIn(i)%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO
  DO i=1,nr_inside
    nRec=nRec+1
    tempw=yParametricOut(VertIn(i)%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO
  DO i=1,nr_inside
    nRec=nRec+1
    tempw=zParametricOut(VertIn(i)%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO
  nRec=nRec+1
  WRITE(10,REC=nRec) cellsGMV(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) cellsGMV(5:8)
  nRec=nRec+1

  IF(out_area=='y') THEN
    nr_vieworocells=nr_viewincells+nr_viewcutcells+nr_viewcutf1
    WRITE(10,REC=nRec) NumberCellsGMV
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO iz=v_z0+1,v_z1
          DO iy=v_y0+1,v_y1
            DO ix=v_x0+1,v_x1
              CALL WriteCellOroGMVBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
            END DO
          END DO
        END DO
      END IF
    END DO
  ELSE
    nr_orocells=nr_incells+nr_cutcells+nr_cutf1
    WRITE(11,REC=nRec) NumberCellsGMV
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            CALL WriteCellOroGMVBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
          END DO
        END DO
      END DO
    END DO
  END IF
  nr_ges_oro_out=nr_gen+nr_grenzeF1+nr_grenzeF2+nr_grenzeSp3+nr_insidehex
  
  ! --------------------------
  ! ----- Polygone -----------
  ! CutPlane,Bloecke,Haeuser
  !---------------------------
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(5:8)
  ! CutPlane:
  !----------
  IF(out_area=='y') THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO iz=v_z0+1,v_z1
          DO iy=v_y0+1,v_y1
            DO ix=v_x0+1,v_x1
              CALL WriteOroCutPlanePolyGMVBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
            END DO
          END DO
        END DO
      END IF
    END DO
  ELSE
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            CALL WriteOroCutPlanePolyGMVBinary(Cells(ix,iy,iz)%Cell,ix,iy,iz)
          END DO
        END DO
      END DO
    END DO
  END IF

  ! Bloecke:
  !----------
  DO ib=1,nb
   CALL Set(Floor(ib))
   CALL WriteBlockGMVBinary
  END DO  !ib
  ! Haeuser:
  !----------
  CALL WriteAllHausMultiColorBinGMV
  !.. WRITE(10) endpoly
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(5:8)

  !.. WRITE(10) endgmv
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(5:8)

  CLOSE(10)

END SUBROUTINE WriteAllCellsOroBinGMV

SUBROUTINE WriteCellOroGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7),yFace(6)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out
 
  yFace(1:6)=0 
  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN
      ! general
      !--------
      nr_gen=nr_gen+1
      nr_cell_vc_oro=nr_cell_vc_oro+1
      nFaces=1
      IF ((MIN(Cell%Face1%Edge1%Vert1%in_out,Cell%Face1%Edge1%Vert2%in_out, &
               Cell%Face1%Edge3%Vert1%in_out,Cell%Face1%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face1%NumberVert>2) THEN
        nVerts(nFaces)=Cell%Face1%NumberVert
        nFaces=nFaces+1
        yFace(1)=1
      END IF
      IF ((MIN(Cell%Face2%Edge1%Vert1%in_out,Cell%Face2%Edge1%Vert2%in_out, &
               Cell%Face2%Edge3%Vert1%in_out,Cell%Face2%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face2%NumberVert>2) THEN
        nVerts(nFaces)=Cell%Face2%NumberVert
        nFaces=nFaces+1
        yFace(2)=1
      END IF
      IF ((MIN(Cell%Face3%Edge1%Vert1%in_out,Cell%Face3%Edge1%Vert2%in_out, &
               Cell%Face3%Edge3%Vert1%in_out,Cell%Face3%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face3%NumberVert>2) THEN
        nVerts(nFaces)=Cell%Face3%NumberVert
        nFaces=nFaces+1
        yFace(3)=1
      END IF
      IF ((MIN(Cell%Face4%Edge1%Vert1%in_out,Cell%Face4%Edge1%Vert2%in_out, &
               Cell%Face4%Edge3%Vert1%in_out,Cell%Face4%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face4%NumberVert>2) THEN
        nVerts(nFaces)=Cell%Face4%NumberVert
        nFaces=nFaces+1
        yFace(4)=1
      END IF
      IF ((MIN(Cell%Face5%Edge1%Vert1%in_out,Cell%Face5%Edge1%Vert2%in_out, &
               Cell%Face5%Edge3%Vert1%in_out,Cell%Face5%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face5%NumberVert>2) THEN
        nVerts(nFaces)=Cell%Face5%NumberVert
        nFaces=nFaces+1
        yFace(5)=1
      END IF
      IF ((MIN(Cell%Face6%Edge1%Vert1%in_out,Cell%Face6%Edge1%Vert2%in_out, &
               Cell%Face6%Edge3%Vert1%in_out,Cell%Face6%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face6%NumberVert>2) THEN
        nVerts(nFaces)=Cell%Face6%NumberVert
        nFaces=nFaces+1
        yFace(6)=1
      END IF
      nVerts(nFaces)=Cell%vc
      nRec=nRec+1
      WRITE(10,REC=nRec) general(1:4)
      nRec=nRec+1
      WRITE(10,REC=nRec) general(5:8)
      nRec=nRec+1
      WRITE(10,REC=nRec) nfaces
!     WRITE(100,*) general(1:8),i,j,k
!     WRITE(100,*) nfaces
!     WRITE(100,*) nVerts(1:nfaces)
      DO iVert=1,nFaces
        nRec=nRec+1
        WRITE(10,REC=nRec) nVerts(iVert)
      END DO
      IF (yFace(1)==1) THEN
!       WRITE(100,*) 'Face1:',Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
        DO iVert=1,Cell%Face1%NumberVert
          nRec=nRec+1
          WRITE(10,REC=nRec)  Cell%Face1%VertexList(iVert)
        END DO
      END IF
      IF (yFace(2)==1) THEN
        WRITE(100,*) 'Face2:',Cell%Face2%VertexList(1:Cell%Face2%NumberVert)
        DO iVert=1,Cell%Face2%NumberVert
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%Face2%VertexList(iVert)
        END DO
      END IF
      IF (yFace(3)==1) THEN
!       WRITE(100,*) 'Face3:',Cell%Face3%VertexList(1:Cell%Face3%NumberVert)
        DO iVert=1,Cell%Face3%NumberVert
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%Face3%VertexList(iVert)
        END DO
      END IF
      IF (yFace(4)==1) THEN
!       WRITE(100,*) 'Face4:',Cell%Face4%VertexList(1:Cell%Face4%NumberVert)
        DO iVert=1,Cell%Face4%NumberVert
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%Face4%VertexList(iVert)
        END DO
      END IF
      IF (yFace(5)==1) THEN
!       WRITE(100,*) 'Face5:',Cell%Face5%VertexList(1:Cell%Face5%NumberVert)
        DO iVert=1,Cell%Face5%NumberVert
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%Face5%VertexList(iVert)
        END DO
      END IF
      IF (yFace(6)==1) THEN
!       WRITE(100,*) 'Face6:',Cell%Face6%VertexList(1:Cell%Face6%NumberVert)
        DO iVert=1,Cell%Face6%NumberVert
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%Face6%VertexList(iVert)
        END DO
      END IF
      IF (nfaces/=1) THEN
!       WRITE(100,*) 'Cut  :',Cell%VertCut(1:Cell%vc)
        DO iVert=1,Cell%vc
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%VertCut(iVert)
        END DO
      END IF
    ELSE
      ! hex
      !------
      !          4--------3    !Point-Folge Output-Hex-GMV-Vorlage-Def
      !         /|       /|
      !        1--------2 |
      !        | 8------|-7
      !        |/       |/
      !        5--------6
      nr_insidehex=nr_insidehex+1
      nRec=nRec+1
      WRITE(10,REC=nRec) hex(1:4)
      nRec=nRec+1
      WRITE(10,REC=nRec) hex(5:8)
      nRec=nRec+1
      WRITE(10,REC=nRec) 8
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j-1,k)%Vertex%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j-1,k)%Vertex%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j,k)%Vertex%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j,k)%Vertex%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j-1,k-1)%Vertex%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j-1,k-1)%Vertex%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j,k-1)%Vertex%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j,k-1)%Vertex%nrInP
    END IF
  END IF
END SUBROUTINE WriteCellOroGMVBinary

SUBROUTINE WriteOroCutPlanePolyGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  ! local
  INTEGER :: mat,nvert
  INTEGER :: li

  mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN
      ! polygon
      !---------
      mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
      nvert=Cell%vc
      nRec=nRec+1
      WRITE(10,REC=nRec) mat
      nRec=nRec+1
      WRITE(10,REC=nRec) nvert
      DO li=1,nvert
        nRec=nRec+1
        tempw=xParametricOut(VertIn(Cell%VertCut(li))%Vertex)
        WRITE(10,REC=nRec) tempw
      END DO
      DO li=1,nvert
        nRec=nRec+1
        tempw=yParametricOut(VertIn(Cell%VertCut(li))%Vertex)
        WRITE(10,REC=nRec) tempw
      END DO
      DO li=1,nvert
        nRec=nRec+1
        tempw=zParametricOut(VertIn(Cell%VertCut(li))%Vertex)
        WRITE(10,REC=nRec) tempw
      END DO
    ELSE IF (Cell%vc==0) THEN
      IF (Cell%Face1%Edge1%Vert1%in_out==0.AND. &
          Cell%Face1%Edge1%Vert2%in_out==0.AND. &
          Cell%Face1%Edge3%Vert1%in_out==0.AND. &
          Cell%Face1%Edge3%Vert2%in_out==0 ) THEN
        mat=25   !mat=26 !zum check
        nvert=Cell%Face1%NumberVert !4
        nRec=nRec+1
        WRITE(10,REC=nRec) mat
        nRec=nRec+1
        WRITE(10,REC=nRec) nvert
        DO li=1,nvert
          nRec=nRec+1
          tempw=xParametricOut(VertIn(Cell%Face1%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=yParametricOut(VertIn(Cell%Face1%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=zParametricOut(VertIn(Cell%Face1%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
      END IF
      IF (Cell%Face2%Edge1%Vert1%in_out==0.AND. &
          Cell%Face2%Edge1%Vert2%in_out==0.AND. &
          Cell%Face2%Edge3%Vert1%in_out==0.AND. &
          Cell%Face2%Edge3%Vert2%in_out==0 ) THEN
        mat=25 !mat=27 ! zum check
        nvert=Cell%Face2%NumberVert !4
        nRec=nRec+1
        WRITE(10,REC=nRec) mat
        nRec=nRec+1
        WRITE(10,REC=nRec) nvert
        DO li=1,nvert
          nRec=nRec+1
          tempw=xParametricOut(VertIn(Cell%Face2%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=yParametricOut(VertIn(Cell%Face2%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=zParametricOut(VertIn(Cell%Face2%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
      END IF
      IF (Cell%Face3%Edge1%Vert1%in_out==0.AND. &
          Cell%Face3%Edge1%Vert2%in_out==0.AND. &
          Cell%Face3%Edge3%Vert1%in_out==0.AND. &
          Cell%Face3%Edge3%Vert2%in_out==0 ) THEN
        nvert=Cell%Face3%NumberVert !4
        nRec=nRec+1
        WRITE(10,REC=nRec) mat
        nRec=nRec+1
        WRITE(10,REC=nRec) nvert
        DO li=1,nvert
          nRec=nRec+1
          tempw=xParametricOut(VertIn(Cell%Face3%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=yParametricOut(VertIn(Cell%Face3%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=zParametricOut(VertIn(Cell%Face3%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
      END IF
      IF (Cell%Face4%Edge1%Vert1%in_out==0.AND. &
          Cell%Face4%Edge1%Vert2%in_out==0.AND. &
          Cell%Face4%Edge3%Vert1%in_out==0.AND. &
          Cell%Face4%Edge3%Vert2%in_out==0 ) THEN
        nvert=Cell%Face4%NumberVert !4
        nRec=nRec+1
        WRITE(10,REC=nRec) mat
        nRec=nRec+1
        WRITE(10,REC=nRec) nvert
        DO li=1,nvert
          nRec=nRec+1
          tempw=xParametricOut(VertIn(Cell%Face4%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=yParametricOut(VertIn(Cell%Face4%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=zParametricOut(VertIn(Cell%Face4%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
      END IF
      IF (Cell%Face5%Edge1%Vert1%in_out==0.AND. &
          Cell%Face5%Edge1%Vert2%in_out==0.AND. &
          Cell%Face5%Edge3%Vert1%in_out==0.AND. &
          Cell%Face5%Edge3%Vert2%in_out==0 ) THEN
        nvert=Cell%Face5%NumberVert !4
        nRec=nRec+1
        WRITE(10,REC=nRec) mat
        nRec=nRec+1
        WRITE(10,REC=nRec) nvert
        DO li=1,nvert
          nRec=nRec+1
          tempw=xParametricOut(VertIn(Cell%Face5%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=yParametricOut(VertIn(Cell%Face5%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=zParametricOut(VertIn(Cell%Face5%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
      END IF
      IF (Cell%Face6%Edge1%Vert1%in_out==0.AND. &
          Cell%Face6%Edge1%Vert2%in_out==0.AND. &
          Cell%Face6%Edge3%Vert1%in_out==0.AND. &
          Cell%Face6%Edge3%Vert2%in_out==0 ) THEN
        nvert=Cell%Face6%NumberVert !4
        nRec=nRec+1
        WRITE(10,REC=nRec) mat
        nRec=nRec+1
        WRITE(10,REC=nRec) nvert
        DO li=1,nvert
          nRec=nRec+1
          tempw=xParametricOut(VertIn(Cell%Face6%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=yParametricOut(VertIn(Cell%Face6%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          tempw=zParametricOut(VertIn(Cell%Face6%VertexList(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
      END IF
    END IF
  END IF
END SUBROUTINE WriteOroCutPlanePolyGMVBinary

SUBROUTINE Search_OutArea()

  v_x0=Domain%view_ixe
  v_x1=Domain%view_ixa
  v_y0=Domain%view_iye
  v_y1=Domain%view_iya
  v_z0=Domain%view_ize
  v_z1=Domain%view_iza
  !...........................................................
  IF (igx0>=Domain%view_ixa.AND.igx0<Domain%view_ixe) THEN
    v_x0=igx0
  ELSE IF (igx0<Domain%view_ixa.AND.igx1>Domain%view_ixa) THEN
    v_x0=Domain%view_ixa
  END IF
  IF (igx1<=Domain%view_ixe.AND.igx1>Domain%view_ixa) THEN
    v_x1=igx1
  ELSE IF (igx1>Domain%view_ixe.AND.igx0<Domain%view_ixe) THEN
    v_x1=Domain%view_ixe
  END IF
  !............................................................
  IF (igy0>=Domain%view_iya.AND.igy0<Domain%view_iye) THEN
    v_y0=igy0
  ELSE IF (igy0<Domain%view_iya.AND.igy1>Domain%view_iya) THEN
    v_y0=Domain%view_iya
  END IF
  IF (igy1<=Domain%view_iye.AND.igy1>Domain%view_iya) THEN
    v_y1=igy1
  ELSE IF (igy1>Domain%view_iye.AND.igy0<Domain%view_iye) THEN
    v_y1=Domain%view_iye
  END IF
  !.............................................................
  IF (igz0>=Domain%view_iza.AND.igz0<Domain%view_ize) THEN
    v_z0=igz0
  ELSE IF (igz0<Domain%view_iza.AND.igz1>Domain%view_iza) THEN
    v_z0=Domain%view_iza
  END IF
  IF (igz1<=Domain%view_ize.AND.igz1>Domain%view_iza) THEN
    v_z1=igz1
  ELSE IF (igz1>Domain%view_ize.AND.igz0<Domain%view_ize) THEN
    v_z1=Domain%view_ize
  END IF
  !.............................................................
  IF (v_x0.NE.Domain%view_ixe .AND. v_x1.NE.Domain%view_ixa .AND. &
      v_y0.NE.Domain%view_iye .AND. v_y1.NE.Domain%view_iya .AND. &
      v_z0.NE.Domain%view_ize .AND. v_z1.NE.Domain%view_iza ) THEN
    IF (RefineX<0) THEN
      v_zw=v_x0
      v_x0=v_zw/2**(-RefineX)
      v_zw=v_x1
      v_x1=v_zw/2**(-RefineX)
    END IF
    IF (RefineY<0) THEN
      v_zw=v_y0
      v_y0=v_zw/2**(-RefineY)
      v_zw=v_y1
      v_y1=v_zw/2**(-RefineY)
    END IF
    IF (RefineZ<0) THEN
      v_zw=v_z0
      v_z0=v_zw/2**(-RefineZ)
      v_zw=v_z1
      v_z1=v_zw/2**(-RefineZ)
    END IF
    view_cell='y'
  ELSE
    view_cell='n'
  END IF
END SUBROUTINE Search_OutArea

SUBROUTINE WriteBlockGMVBinary
  INTEGER :: ix,iy,iz
  INTEGER :: mat,nvert

  ! polygons
  mat=6
  nvert=4
  DO iy=iy0,iy1,iy1-iy0
    nRec=nRec+1
    WRITE(10,REC=nRec) mat
    nRec=nRec+1
    WRITE(10,REC=nRec) nvert
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix0,iy,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix1,iy,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix1,iy,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix0,iy,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix0,iy,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix1,iy,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix1,iy,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix0,iy,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix0,iy,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix1,iy,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix1,iy,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix0,iy,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO

  DO iz=iz0,iz1,iz1-iz0
    nRec=nRec+1
    WRITE(10,REC=nRec) mat
    nRec=nRec+1
    WRITE(10,REC=nRec) nvert
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix0,iy0,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix1,iy0,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix1,iy1,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix0,iy1,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix0,iy0,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix1,iy0,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix1,iy1,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix0,iy1,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix0,iy0,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix1,iy0,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix1,iy1,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix0,iy1,iz)%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO

  DO ix=ix0,ix1,ix1-ix0
    nRec=nRec+1
    WRITE(10,REC=nRec) mat
    nRec=nRec+1
    WRITE(10,REC=nRec) nvert
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix,iy0,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix,iy0,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix,iy1,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix,iy1,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix,iy0,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix,iy0,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix,iy1,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix,iy1,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix,iy0,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix,iy0,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix,iy1,iz1)%Vertex)
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix,iy1,iz0)%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO

END SUBROUTINE WriteBlockGMVBinary

SUBROUTINE WriteAllHausMultiColorBinGMV
  INTEGER   :: h,f,i,nr_p_faces
  INTEGER   :: nr_color_w,nr_farbe,nr_fvor
  CHARACTER :: vor_type

  nr_color_w=1
  DO h=1,NumberHaus
    DO f=1,Haus(h)%NumberOfFaces
      nr_p_faces=Haus(h)%Faces(f)%NumberOfPoints
      !Farbe setzen, Standard 1-5 routieren, 6 für Blöcke
      IF (Haus(h)%Faces(f)%type=='r') THEN
        nr_farbe=5
      ELSE !Haus(h)%Faces(f)%type=='w'
        IF (nr_color_w<4) THEN
          nr_color_w=nr_color_w+1
          nr_farbe=nr_color_w
        ELSE
          nr_color_w=2
          nr_farbe=nr_color_w
        END IF
      END IF ! r,w
     !Output
      nRec=nRec+1
      WRITE(10,REC=nRec) nr_farbe
      nRec=nRec+1
      WRITE(10,REC=nRec) nr_p_faces
      DO i=1,nr_p_faces
        nRec=nRec+1
        tempw=xPointParametricOut(Haus(h)%Faces(f)%Points(i))
        WRITE(10,REC=nRec) tempw
      END DO
      DO i=1,nr_p_faces
        nRec=nRec+1
        tempw=yPointParametricOut(Haus(h)%Faces(f)%Points(i))
        WRITE(10,REC=nRec) tempw
      END DO
      DO i=1,nr_p_faces
        nRec=nRec+1
        tempw=zPointParametricOut(Haus(h)%Faces(f)%Points(i))
        WRITE(10,REC=nRec) tempw
      END DO
      vor_type=Haus(h)%Faces(f)%type
    END DO
  END DO  !h

END SUBROUTINE WriteAllHausMultiColorBinGMV

SUBROUTINE WriteCutPlanePolyGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  ! loacal
  INTEGER :: mat1,mat2,mat3,nvert
  INTEGER :: li,ix,iy,liv1,liv2,gl

  ! Für Specialsuche:  Farbe mat1-3 -> für Var1-3
  mat1=25 !(rotbraun); Cell-Schnitte Standard
  mat2=26 !(grün)
  mat3=27 !(blau); z-Ebene-0,auch Wasser
  ! Standard-View
  !mat1=25;mat2=25;mat3=25 

  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN
    !~~~~~~~~~~~~~~~~~~
      ! polygon_cutface 
      ! --------------
      nvert=Cell%vc
      !............
      nRec=nRec+1
      WRITE(10,REC=nRec) mat1
      nRec=nRec+1
      WRITE(10,REC=nRec) nvert
      DO li=1,nvert
        nRec=nRec+1
        tempw=xParametricOut(VertOut(Cell%VertCut(li))%Vertex)
        WRITE(10,REC=nRec) tempw
      END DO
      DO li=1,nvert
        nRec=nRec+1
        tempw=yParametricOut(VertOut(Cell%VertCut(li))%Vertex)
        WRITE(10,REC=nRec) tempw
      END DO
      DO li=1,nvert
        nRec=nRec+1
        tempw=zParametricOut(VertOut(Cell%VertCut(li))%Vertex)
        WRITE(10,REC=nRec) tempw
      END DO
      IF(k==1) THEN
         gl=0
         IF(Cell%vc==Cell%Face1%NumberVert) THEN
            DO liv1=1,Cell%vc
              DO liv2=1,Cell%vc
               IF(Cell%VertCut(liv1)==Cell%Face1%VertexList(liv2)) THEN
                 gl=gl+1
               END IF
             END DO
           END DO
         END IF
         IF(gl/=Cell%vc) THEN
           IF((Cell%Face1%Edge1%yes_sp==1.OR.Cell%Face1%Edge2%yes_sp==1.OR. &
               Cell%Face1%Edge3%yes_sp==1.OR.Cell%Face1%Edge4%yes_sp==1) .OR. &
              (Cell%Face1%ec==-1.and.Cell%Face1%NumberVert>2)) THEN
              ! polygon: Special-Hügel, ContainerT.grid v8931
              ! --------- 
              CALL WriteCutFacePolyGMVBinary(Cell%Face1,mat2)

           ELSE IF(Cell%Face1%in_out==-1) THEN
              ! polygon: (nur Face1-3points)
              ! ---------
              CALL WriteCutFacePolyGMVBinary(Cell%Face1,mat3)
           END IF
         END IF  ! Ungleich Check
      END IF  ! cell%vc und k==1
    ELSE  !Cell%vc==0
    !~~~~~~~~~~~~~~~~~~~~
      IF (k==1) THEN
        IF(Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
           Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0 ) THEN
          ! polygon: (nur Face1)
          ! ---------
           CALL WriteCutFacePolyGMVBinary(Cell%Face1,mat3)
        ELSE IF(Cell%Face1%in_out==-1) THEN
          ! polygon: (nur Face1-3points)
          ! ---------
           CALL WriteCutFacePolyGMVBinary(Cell%Face1,mat3)
        END IF    ! IF Face1 Grenze ((V0-V4)%in_out==0)
      END IF ! (k==1)
    END IF  ! if-else cell%vc
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END IF  ! if associated celle

END SUBROUTINE WriteCutPlanePolyGMVBinary

SUBROUTINE WriteCutFacePolyGMVBinary(Face1,matv)
  TYPE(Face_T), POINTER :: Face1
  INTEGER :: matv
  !..............................
  INTEGER :: nvert
  INTEGER :: li,ix,iy
  !---------
  ! polygon: (nur Face1)
  !---------
  nvert=Face1%NumberVert
  !................................
  nRec=nRec+1
  WRITE(10,REC=nRec) matv
  nRec=nRec+1
  WRITE(10,REC=nRec) nvert
  DO li=1,nvert
    nRec=nRec+1
    tempw=xParametricOut(VertOut(Face1%VertexList(li))%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO
  DO li=1,nvert
    nRec=nRec+1
    tempw=yParametricOut(VertOut(Face1%VertexList(li))%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO
  DO li=1,nvert
    nRec=nRec+1
    tempw=zParametricOut(VertOut(Face1%VertexList(li))%Vertex)
    WRITE(10,REC=nRec) tempw
  END DO
END SUBROUTINE WriteCutFacePolyGMVBinary

SUBROUTINE WriteAllCellsCutPlaneBinGMV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k

  INTEGER :: NumberCells


  OPEN(UNIT=10,FILE=TRIM(FileName)//'.Cut.out.gmvG',STATUS='UNKNOWN'&
    &         ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)

  CALL CountAllCellsCutPlaneGMV(NumberCells)
  nRec=0
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) nr_cutplane
  DO i=1,nr_out
    IF (VertOut(i)%Vertex%nrCutP>0) THEN
      nRec=nRec+1
      tempw=xParametricOut(VertOut(i)%Vertex)
      WRITE(10,REC=nRec) tempw
    END IF
  END DO
  DO i=1,nr_out
    IF (VertOut(i)%Vertex%nrCutP>0) THEN
      nRec=nRec+1
      tempw=yParametricOut(VertOut(i)%Vertex)
      WRITE(10,REC=nRec) tempw
    END IF
  END DO
  DO i=1,nr_out
    IF (VertOut(i)%Vertex%nrCutP>0) THEN
      nRec=nRec+1
      tempw=zParametricOut(VertOut(i)%Vertex)
      WRITE(10,REC=nRec) tempw
    END IF
  END DO
  nRec=nRec+1
  WRITE(10,REC=nRec) cellsGMV(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) cellsGMV(5:8)
  nRec=nRec+1

  IF(out_area=='y') THEN
    WRITE(10,REC=nRec) NumberCells
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
              CALL WriteCellCutPlaneGMVBinary(Cells(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    WRITE(10,REC=nRec) NumberCells
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteCellCutPlaneGMVBinary(Cells(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF

  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(5:8)
  DO ib=1,nb
    CALL Set(Floor(ib))
    CALL WriteBlockGMVBinary
  END DO  !ib
  CALL WriteAllHausMultiColorBinGMV
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(5:8)

  CLOSE(10)

END SUBROUTINE WriteAllCellsCutPlaneBinGMV

SUBROUTINE CountAllCellsCutPlaneGMV(NumberCells)
  INTEGER :: NumberCells

  INTEGER :: ib,ix,iy,iz

  NumberCells=0
  IF (out_area=='y') THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO iz=v_z0+1,v_z1
          DO iy=v_y0+1,v_y1
            DO ix=v_x0+1,v_x1
              CALL CountCellCutPlaneGMV(Cells(ix,iy,iz)%Cell,ix,iy,iz,NumberCells)
            END DO
          END DO
        END DO
      END IF
    END DO  
  ELSE
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            CALL CountCellCutPlaneGMV(Cells(ix,iy,iz)%Cell,ix,iy,iz,NumberCells)
          END DO
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE CountAllCellsCutPlaneGMV

SUBROUTINE CountCellCutPlaneGMV(Cell,i,j,k,NumberCells)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k
  INTEGER :: NumberCells
  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN
      NumberCells=NumberCells+1
    END IF  
  END IF  
END SUBROUTINE CountCellCutPlaneGMV

SUBROUTINE WriteCellCutPlaneGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out
  INTEGER :: liv1,liv2,gl

  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN
      ! general
      nFaces=1
      nVerts(nFaces)=Cell%vc
      !.. WRITE(10) general,nfaces
      nRec=nRec+1
      WRITE(10,REC=nRec) general(1:4)
      nRec=nRec+1
      WRITE(10,REC=nRec) general(5:8)
      nRec=nRec+1
      WRITE(10,REC=nRec) nfaces
      !.. WRITE(10) (nVerts(i),i=1,nFaces)
      nRec=nRec+1
      WRITE(10,REC=nRec) nVerts(1)
      DO iVert=1,Cell%vc
        nRec=nRec+1
        WRITE(10,REC=nRec) VertOut(Cell%VertCut(iVert))%Vertex%nrCutP
      END DO
    END IF
  END IF
END SUBROUTINE WriteCellCutPlaneGMVBinary

SUBROUTINE  WriteAllCellsSoilBinGMV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k

  CHARACTER*120 :: workstr
  CHARACTER(4)  :: tempfrom
  CHARACTER, POINTER ::p
  INTEGER       :: lNameOut,lerg
  INTEGER :: NumberCellsSoil
  TYPE (Vertex_T) :: V(1:2)

  CALL CountAllCellsCutPlaneGMV(NumberCellsSoil)
  IF (nr_sb>0) THEN
    OPEN(UNIT=10,FILE=TRIM(FileName)//'.Soil.out.gmvG',STATUS='UNKNOWN'&
       &         ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)

     nRec=0
     !.. WRITE(10) gmvinput,ascii/iecxi4r4
     nRec=nRec+1
     WRITE(10,REC=nRec) gmvinput(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) gmvinput(5:8)
     nRec=nRec+1
     WRITE(10,REC=nRec) iecxi4r4(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) iecxi4r4(5:8)
     nRec=nRec+1
     WRITE(10,REC=nRec) nodes(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) nodes(5:8)
     nRec=nRec+1
     WRITE(10,REC=nRec) nr_out_soil
     DO i=1,nr_out_soil
       nRec=nRec+1
       tempw=xParametricOut(VertSoilOut(i)%Vertex)
       WRITE(10,REC=nRec) tempw
     END DO
     DO i=1,nr_out_soil
       nRec=nRec+1
       tempw=yParametricOut(VertSoilOut(i)%Vertex)
       WRITE(10,REC=nRec) tempw
     END DO
     DO i=1,nr_out_soil
       nRec=nRec+1
       tempw=zParametricOut(VertSoilOut(i)%Vertex)
       WRITE(10,REC=nRec) tempw
     END DO
     nRec=nRec+1
     WRITE(10,REC=nRec) cellsGMV(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) cellsGMV(5:8)
     nRec=nRec+1
    
     IF(out_area=='y') THEN
       !.. WRITE(10) nr_cells for soil
       WRITE(10,REC=nRec) NumberCellsSoil*nr_lsoil
       DO ib=1,nb
         CALL Set(Floor(ib))
         CALL Search_OutArea
         IF(view_cell=='y') THEN
           DO k=v_z0+1,v_z1
             DO j=v_y0+1,v_y1
               DO i=v_x0+1,v_x1
                 CALL WriteCellSoilLayerGMVBin(Cells(i,j,k)%Cell,i,j,k)
               END DO 
             END DO
           END DO
         END IF
       END DO  !ib
     ELSE
       WRITE(10,REC=nRec) NumberCellsSoil*nr_lsoil
       DO ib=1,nb
         CALL Set(Floor(ib))
         DO k=iz0+1,iz1
           DO j=iy0+1,iy1
             DO i=ix0+1,ix1
               CALL WriteCellSoilLayerGMVBin(Cells(i,j,k)%Cell,i,j,k)
             END DO
           END DO
         END DO
       END DO  !ib
     END IF

     !.. WRITE(10) polygons
     nRec=nRec+1
     WRITE(10,REC=nRec) polygons(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) polygons(5:8)
   
     !LandClass-Plane
     !---------------
     IF(nr_landdef>0) THEN
       IF(out_area=='y') THEN
         DO ib=1,nb
           CALL Set(Floor(ib))
           CALL Search_OutArea
           IF(view_cell=='y') THEN
             DO k=v_z0+1,v_z1
               DO j=v_y0+1,v_y1
                 DO i=v_x0+1,v_x1
                   CALL WriteLandKlPlanePolyGMVBinary(Cells(i,j,k)%Cell,i,j,k)
                 END DO
               END DO
             END DO
          END IF
         END DO  !ib
       ELSE
         DO ib=1,nb
           CALL Set(Floor(ib))
           DO k=iz0+1,iz1
             DO j=iy0+1,iy1
               DO i=ix0+1,ix1
                 CALL WriteLandKlPlanePolyGMVBinary(Cells(i,j,k)%Cell,i,j,k)
               END DO
             END DO
           END DO
         END DO  !ib
       END IF
     END IF  !if(nr_landdef)

     ! Soil-layering:
     !-----------
     IF(nr_sb>0) THEN
       V(1:2)%Point%x=Domain%xP(0)
       V(1:2)%Point%y=Domain%yP(0)
       V(1)%Point%z=Domain%zP(0)
       V(2)%Point%z=Domain%zP(1)
       dzi_soil=zParametricOut(V(2))-zParametricOut(V(1))
       IF(out_wahlgrid=='C') THEN
         IF(out_area=='y') THEN
           DO ib=1,nb
             CALL Set(Floor(ib))
             CALL Search_OutArea
             IF(view_cell=='y') THEN
               DO k=v_z0+1,v_z1
                 DO j=v_y0+1,v_y1
                   DO i=v_x0+1,v_x1
                     CALL WritePolySoilLayerGMVAsciiBin(Cells(i,j,k)%Cell,i,j,k)
                   END DO
                 END DO
               END DO
             END IF
           END DO  !ib
         ELSE
           DO ib=1,nb
             CALL Set(Floor(ib))
             DO k=iz0+1,iz1
               DO j=iy0+1,iy1
                 DO i=ix0+1,ix1
                   CALL WritePolySoilLayerGMVAsciiBin(Cells(i,j,k)%Cell,i,j,k)
                 END DO
               END DO
             END DO
           END DO  !ib
         END IF
       END IF
     END IF
   
     !.. WRITE(10) endpoly
     nRec=nRec+1
     WRITE(10,REC=nRec) endpoly(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) endpoly(5:8)
     !.. WRITE(10) endgmv
     nRec=nRec+1
     WRITE(10,REC=nRec) endgmv(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) endgmv(5:8)
   
     CLOSE(10)
   END IF

END SUBROUTINE  WriteAllCellsSoilBinGMV

SUBROUTINE WritePolySoilLayerGMVAsciiBin(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k
  !...............................
  
  INTEGER :: li,la,lix,liy,liz,ix,iy,iz
  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7),yFace(6)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out,nvert,posvert(1:8)
  INTEGER :: soiltype(1:nr_lsoil),i_st,mat,fmat(1:7)
  INTEGER :: st(1:9)
  INTEGER :: basis,allface,nrs,is
  INTEGER :: fv,nrf,aktfx,aktfy,aktfz,fakt,f1,f2,n,l
  INTEGER :: ifvc,ifgf,akt_vc,pos
  INTEGER, POINTER :: fnvert(:) 
  REAL(8), POINTER :: xc(:),cfxc(:,:)
  REAL(8), POINTER :: yc(:),cfyc(:,:)
  REAL(8), POINTER :: zc(:),cfzc(:,:)
  REAL(8) :: dzi
  TYPE (Vertex_T) :: V1,VCell(1:8),Fvert(1:4),V(1:2)


  !soiltype: 1=ice   2=rock  3=sand  4=sandy loam  5=loam  6=clay loam
  !          7=clay  8=peat  9=water
  !          clay=Ton; peat=Torf
  !z.Zt: Zuordnung BodenType-->Farbe-View-GMV
  !      1=ice=hellblau(41)  2=rock=grau(49)  3=sand=gelb(36) 4=sandy loam=khaki(55)
  !      5=loam=rot(53)47 6=clay loam(57) 7=clay=dunkelrot(39)  8=peat=dunkel_lila(43)
  !      9=water=blue(44)          39                 57
  mat=16
  yFace(1:6)=0 
  fmat(:)   =(/8,9,10,11,12,13,14/)
  st(:)     =(/41,49,36,55,47,39,57,43,44/)

  !Testvarianten
  !soiltype(1:8)=25   ! 14=grau,16=lila-braun,25=rot-braun,33=gruen

  soiltype(:)=0
  ifvc=0
  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN  
      nvert=Cell%vc    ! Nummer Point Cell%vc, resultiert anz. Flächen (F3-nFace) 
      nrf=nvert+2      ! nFace = anz.Cell%vc + (face1 und face2)
      ALLOCATE(fnvert(1:nrf))
      ifvc=1
    ELSE 
      IF (k==1 .and. &  !(Cell%in_out==4.AND. &
        (Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
         Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0) )THEN
        nvert=4          ! Nummer Point Face1, resultiert anz. Flächen (F3-nFace)
        nrf=nvert+2      ! nFace = anz. Punkte(aus Celle-Face1=4) + (face1 und face2)
        ALLOCATE(fnvert(1:nrf))
        IF (.NOT.ASSOCIATED(Cell%VertCut)) THEN
          ALLOCATE(Cell%VertCut(nvert))
        END IF
         Cell%VertCut(1:nvert)=Cell%Face1%VertexList(1:nvert)
         !Cell%vc   !!!!nicht belegen da Orography-Abfrage noch folgt!!!!
         ifvc=1
      ELSE
      END IF
    END IF
!...........................................................
! 2.Variante:
! Beschreibt Bodenschichten entlang Orography mit Polygone
!...........................................................
!-Soile-Boden unterhalb Cut-line
!-Face2 entspricht Cut-Fläche, Face2 kopiert auf Face1,
!-Face1-zebene mit minus dzi_soil verschoben
!-Face3-Face(nvert) entgegen Uhrzeigersinn aufgestellt
!-Ermitteln der jeweiligen Faces, Checken gl. Points x-,bzw y-Richtung sind
!-alle Polygone je Face erhalten Points mit Anzahl Face(nvert+1) 
!-je Soil-Boden werden alle Faces je Ebene Soil(1-7)
! in z-Ebene komplett mit minus dzi_soil verschoben
!...........................................................

    IF (ifvc==1) THEN            
      !Belege Soil-Liste mit akt.Soiltypen der akt Grund-Celle 
      !-----------------------------
      DO nrs=1,nr_soildef
        IF (i>(s_ixa(nrs)*2.e0**RefineX).AND.i<=(s_ixe(nrs)*2.e0**RefineX) .AND. &
          j>(s_iya(nrs)*2.e0**RefineY).AND.j<=(s_iye(nrs)*2.e0**RefineY)) THEN
          DO is=1,nr_lsoil !Std. 7 mit #NrSoilLayers- dynamisch 
            soiltype(is)=st(soil_type(nrs,is))
          END DO
          EXIT
        END IF
      END DO
      dzi=dzi_soil
  
      !Points Face1-Face(nvert+2) zusammenstellen
      !-------------------------------------------
      ALLOCATE(cfxc(1:nrf,1:8)) ! max. 8 Points möglich
      ALLOCATE(cfyc(1:nrf,1:8))
      ALLOCATE(cfzc(1:nrf,1:8))
      cfxc(:,:)=11111.11d0 !bessere Sicht Debugger
      cfyc(:,:)=11111.11d0
      cfzc(:,:)=11111.11d0
      !-----
      !Face2
      !-----
      DO li=1,nvert
        cfxc(2,li)=xParametricOut(VertOut(Cell%VertCut(li))%Vertex)
        cfyc(2,li)=yParametricOut(VertOut(Cell%VertCut(li))%Vertex)
        cfzc(2,li)=zParametricOut(VertOut(Cell%VertCut(li))%Vertex)
      END DO
      !für Polygone Anfangswert==Endwert; fnvert um 1 akumulieren
      cfxc(2,li)=cfxc(2,1)
      cfyc(2,li)=cfyc(2,1)
      cfzc(2,li)=cfzc(2,1)
      !-----
      !Face1
      !-----
      cfxc(1,:)=cfxc(2,:)
      cfyc(1,:)=cfyc(2,:)
      cfzc(1,:)=cfzc(2,:)-dzi_soil
      fnvert(1:2)=nvert+1  ! +1 für anhängen Anfangswert
  
      fakt=3
      f1=1;f2=2;pos=1
      DO While (fakt<=nrf)
        !.................................
        !immer 2Point(Vert) Face2/1 eine Face(fakt)
        cfxc(fakt,1) = cfxc(f1,pos) !AnfangPoint
        cfxc(fakt,2) = cfxc(f1,pos+1)
        cfxc(fakt,3) = cfxc(f2,pos+1)
        cfxc(fakt,4) = cfxc(f2,pos)
        cfxc(fakt,5) = cfxc(f1,pos) !EndPoint=AnfangPoint
        !.................................
        cfyc(fakt,1) = cfyc(f1,pos)
        cfyc(fakt,2) = cfyc(f1,pos+1)
        cfyc(fakt,3) = cfyc(f2,pos+1)
        cfyc(fakt,4) = cfyc(f2,pos)
        cfyc(fakt,5) = cfyc(f1,pos)
        !.................................
        cfzc(fakt,1) = cfzc(f1,pos)
        cfzc(fakt,2) = cfzc(f1,pos+1)
        cfzc(fakt,3) = cfzc(f2,pos+1)
        cfzc(fakt,4) = cfzc(f2,pos)
        cfzc(fakt,5) = cfzc(f1,pos)
        !.................................
        fnvert(fakt)=4+1 ! +1 für anhängen Anfangswert
        pos=pos+1;fakt=fakt+1
        !.................................
      END DO !DO While(fakt<=nrf)


      !Output Soil-Layer to file *out.gmvG for ascii
      IF (out_type=='a') then
        DO i_st=1,nr_lsoil
          DO fv=1,nrf
            Write(10,*) soiltype(i_st), fnvert(fv) &
                   ,(cfxc(fv,li),li=1,fnvert(fv)) &
                   ,(cfyc(fv,li),li=1,fnvert(fv)) &
                  ,(cfzc(fv,li),li=1,fnvert(fv))
          END DO
          DO fv=1,nrf
            cfzc(fv,:)=cfzc(fv,:)-dzi_soil
          END DO
        END DO
      ELSE
        !Output Soil-Layer to file *out.gmvG for binary
        DO i_st=1,nr_lsoil
          DO fv=1,nrf
            nRec=nRec+1
            WRITE(10,REC=nRec) soiltype(i_st)
            nRec=nRec+1
            WRITE(10,REC=nRec) fnvert(fv)
            DO li=1,fnvert(fv)
               nRec=nRec+1
               tempw=cfxc(fv,li)
               WRITE(10,REC=nRec) tempw
            END DO
            DO li=1,fnvert(fv)
              nRec=nRec+1
              tempw=cfyc(fv,li)
              WRITE(10,REC=nRec) tempw
            END DO
            DO li=1,fnvert(fv)
              nRec=nRec+1
              tempw=cfzc(fv,li)
              WRITE(10,REC=nRec) tempw
            END DO
          END DO
          DO fv=1,nrf
            cfzc(fv,:)=cfzc(fv,:)-dzi_soil
          END DO
        END DO
      END IF    

      DEALLOCATE(fnvert)
      DEALLOCATE(cfxc)
      DEALLOCATE(cfyc)
      DEALLOCATE(cfzc)
  
    END IF  ! IF(ifvc)
  END IF  ! If(ASSOCIATED(Cell))

END SUBROUTINE WritePolySoilLayerGMVAsciiBin

SUBROUTINE WriteLandKlPlanePolyGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  ! loacal
  INTEGER :: mat,nvert
  INTEGER :: li,ix,iy
  INTEGER :: Lk(1:9)
  !LandKlassen: 
  !Land use:  1       2        3          4          5       6    
  !         urban savannah  deciduous  coniferous  mixed  shrubland 
  !         area            forest     forest      forest          
  !           7       8         9
  !         annual grassland  water
  !         crops
  !Zuordnung Lk:(1) (2) (3) (4) (5) (6)  (7)   (8) (9)
  !             79  85  92  90  77  82  81(89)  96  76 
  Lk(:)=(/79,85,92,90,77,82,89,96,76/)
  !Lk(:)=(/74,77,82,87,90,92,96,97/)  ! Testfarben: Grün-Nuancen 
                                      ! für Wald-,Wiesen-,Sträucherdef.
  ! Im Nachtrag (7) 81 dark gray auf ( 89) braun 
  IF (ASSOCIATED(Cell)) THEN
   IF (Cell%LandClass /=0 ) THEN
    IF (Cell%vc>0) THEN
        ! polygons 
        mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
        IF(Cell%LandClass/=0) mat=Lk(Cell%LandClass)
        nvert=Cell%vc
        nRec=nRec+1
        WRITE(10,REC=nRec) mat
        nRec=nRec+1
        WRITE(10,REC=nRec) nvert
        DO li=1,nvert
          nRec=nRec+1
          !tempw=VertOut(Cell%VertCut(li))%Point%x
          tempw=xParametricOut(VertOut(Cell%VertCut(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          !tempw=VertOut(Cell%VertCut(li))%Point%y
          tempw=yParametricOut(VertOut(Cell%VertCut(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          !tempw=VertOut(Cell%VertCut(li))%Point%z
          tempw=zParametricOut(VertOut(Cell%VertCut(li))%Vertex)
          WRITE(10,REC=nRec) tempw
        END DO
        IF(k==1) THEN
          IF(Cell%Face1%Edge1%yes_sp==1.OR.Cell%Face1%Edge2%yes_sp==1.OR. &
             Cell%Face1%Edge3%yes_sp==1.OR.Cell%Face1%Edge4%yes_sp==1) THEN
             mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
             IF(Cell%LandClass/=0) mat=Lk(Cell%LandClass)
             nvert=Cell%Face1%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) mat
             nRec=nRec+1
             WRITE(10,REC=nRec) nvert
             DO li=1,nvert
               nRec=nRec+1
               tempw=xParametricOut(VertOut(Cell%Face1%VertexList(li))%Vertex)
               WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
               nRec=nRec+1
               tempw=yParametricOut(VertOut(Cell%Face1%VertexList(li))%Vertex)
               WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
               nRec=nRec+1
               tempw=zParametricOut(VertOut(Cell%Face1%VertexList(li))%Vertex)
               WRITE(10,REC=nRec) tempw
             END DO
          END IF
        END IF

    ELSE IF  & !(Cell%in_out==4) THEN !Cell%Vol==maxVoloder >0.0d0 ausreichend
      (k==1) THEN
      IF (Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
          Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0 ) THEN
             ! polygon: (Face1)
             !---------
         mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
         IF(Cell%LandClass/=0) mat=Lk(Cell%LandClass)
         nvert=4
         nRec=nRec+1
         WRITE(10,REC=nRec) mat
         nRec=nRec+1
         WRITE(10,REC=nRec) nvert
         DO li=1,nvert
           SELECT CASE (li)
                CASE (1)        !.. Vertices(i-1,j-1,k)
                    ix=i-1;iy=j-1
                CASE (2)        !.. Vertices(i,j-1,k
                    ix=i
                CASE (3)        !..Vertices(i,j,k)
                    iy=j
                CASE (4)        !..Vertices(i-1,j,k)
                    ix=i-1
               END SELECT
           nRec=nRec+1
           tempw=xParametricOut(Vertices(ix,iy,k-1)%Vertex)
           WRITE(10,REC=nRec) tempw
         END DO
         DO li=1,nvert
           SELECT CASE (li)
                CASE (1)        !.. Vertices(i-1,j-1,k)
                    ix=i-1;iy=j-1
                CASE (2)        !.. Vertices(i,j-1,k
                    ix=i
                CASE (3)        !..Vertices(i,j,k)
                    iy=j
                CASE (4)        !..Vertices(i-1,j,k)
                    ix=i-1
               END SELECT
           nRec=nRec+1
           tempw=yParametricOut(Vertices(ix,iy,k-1)%Vertex)
           WRITE(10,REC=nRec) tempw
         END DO
         DO li=1,nvert
           SELECT CASE (li)
                CASE (1)        !.. Vertices(i-1,j-1,k)
                    ix=i-1;iy=j-1
                CASE (2)        !.. Vertices(i,j-1,k
                    ix=i
                CASE (3)        !..Vertices(i,j,k)
                    iy=j
                CASE (4)        !..Vertices(i-1,j,k)
                    ix=i-1
               END SELECT
           nRec=nRec+1
           tempw=zParametricOut(Vertices(ix,iy,k-1)%Vertex)
           WRITE(10,REC=nRec) tempw
         END DO
      END IF   ! IF Face1 Grenze ((V0-V4)%in_out==0)
    END IF   ! IF(Cell%vc>0) ELSE IF (k==1)
   END IF  !  IF (Cell%LandClass/=0)
  END IF  ! IF (ASSOCIATED(Cell)) 
END SUBROUTINE WriteLandKlPlanePolyGMVBinary

SUBROUTINE WriteCellSoilLayerGMVBin(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k
  !..............................
  INTEGER :: iv,isol,nFaces
  
  IF (ASSOCIATED(Cell)) THEN   
    IF (Cell%vc>0) THEN
!     WRITE(*,*) 'Cell',i,j,k,Cell%vc,nRec
      nFaces=Cell%vc+2
      DO isol=1,nr_lsoil
        nRec=nRec+1
        WRITE(10,REC=nRec) general(1:4)
        nRec=nRec+1
        WRITE(10,REC=nRec) general(5:8)
        nRec=nRec+1
!       WRITE(*,*) 'Face',isol
        WRITE(10,REC=nRec) nFaces
!       WRITE(*,*) 'nFaces',nFaces
        DO iv=1,2
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%vc
!         WRITE(*,*) Cell%vc
        END DO
        DO iv=3,nFaces
          nRec=nRec+1
          WRITE(10,REC=nRec) 4
!         WRITE(*,*) 4
        END DO
        ! Face 1
        DO iv=1,Cell%vc
          nRec=nRec+1
          WRITE(10,REC=nRec) VertOut(Cell%VertCut(iv))%Vertex%nrCutP+(iSol-1)*nr_cutplane
!         WRITE(*,*) VertOut(Cell%VertCut(i))%Vertex%nrCutP+(iSol-1)*nr_cutplane
        END DO  
        ! Face 2
        DO iv=1,Cell%vc
          nRec=nRec+1
          WRITE(10,REC=nRec) VertOut(Cell%VertCut(iv))%Vertex%nrCutP+iSol*nr_cutplane
!         WRITE(*,*) VertOut(Cell%VertCut(i))%Vertex%nrCutP+iSol*nr_cutplane
        END DO  
        ! Face 3...Cell%vc+2
        DO iv=1,Cell%vc-1
          nRec=nRec+1
          WRITE(10,REC=nRec) VertOut(Cell%VertCut(iv+1))%Vertex%nrCutP+(iSol-1)*nr_cutplane
!         WRITE(*,*) VertOut(Cell%VertCut(i+1))%Vertex%nrCutP+(iSol-1)*nr_cutplane
          nRec=nRec+1
          WRITE(10,REC=nRec) VertOut(Cell%VertCut(iv+1))%Vertex%nrCutP+iSol*nr_cutplane
!         WRITE(*,*) VertOut(Cell%VertCut(i+1))%Vertex%nrCutP+iSol*nr_cutplane
          nRec=nRec+1
          WRITE(10,REC=nRec) VertOut(Cell%VertCut(iv))%Vertex%nrCutP+iSol*nr_cutplane
!         WRITE(*,*) VertOut(Cell%VertCut(i))%Vertex%nrCutP+iSol*nr_cutplane
          nRec=nRec+1
          WRITE(10,REC=nRec) VertOut(Cell%VertCut(iv))%Vertex%nrCutP+(iSol-1)*nr_cutplane
!         WRITE(*,*) VertOut(Cell%VertCut(i))%Vertex%nrCutP+(iSol-1)*nr_cutplane
        END DO  
        nRec=nRec+1
        WRITE(10,REC=nRec) VertOut(Cell%VertCut(1))%Vertex%nrCutP+(iSol-1)*nr_cutplane
!       WRITE(*,*) VertOut(Cell%VertCut(1))%Vertex%nrCutP+(iSol-1)*nr_cutplane
        nRec=nRec+1
        WRITE(10,REC=nRec) VertOut(Cell%VertCut(1))%Vertex%nrCutP+iSol*nr_cutplane
!       WRITE(*,*) VertOut(Cell%VertCut(1))%Vertex%nrCutP+iSol*nr_cutplane
        nRec=nRec+1
        WRITE(10,REC=nRec) VertOut(Cell%VertCut(Cell%vc))%Vertex%nrCutP+iSol*nr_cutplane
!       WRITE(*,*) VertOut(Cell%VertCut(Cell%vc))%Vertex%nrCutP+iSol*nr_cutplane
        nRec=nRec+1
        WRITE(10,REC=nRec) VertOut(Cell%VertCut(Cell%vc))%Vertex%nrCutP+(iSol-1)*nr_cutplane
!       WRITE(*,*) VertOut(Cell%VertCut(Cell%vc))%Vertex%nrCutP+(iSol-1)*nr_cutplane
      END DO
    END IF  
  END IF  


!   ! according to the selected NrB_Cells for Weight2
!   IF (((Cell%in_out<5.AND.Cell%Vol>0.0d0) .OR. &
!        (Cell%in_out==6.AND.Cell%vc>0) .OR. &
!        (Cell%in_out==5.AND.Cell%vc>0)) &
!        .and. &
!       (Cell%Face1%Vol/ VolFace_XY(Cell%Face1)<1.d0-dist_scMaxCell .OR. &
!        Cell%Face2%Vol/ VolFace_XY(Cell%Face2)<1.d0-dist_scMaxCell .OR. &
!        Cell%Face3%Vol/ VolFace_ZX(Cell%Face3)<1.d0-dist_scMaxCell .OR. &
!        Cell%Face4%Vol/ VolFace_ZX(Cell%Face4)<1.d0-dist_scMaxCell .OR. &
!        Cell%Face5%Vol/ VolFace_YZ(Cell%Face5)<1.d0-dist_scMaxCell .OR. &
!        Cell%Face6%Vol/ VolFace_YZ(Cell%Face6)<1.d0-dist_scMaxCell) ) THEN
!
!       !-Soile-Boden unterhalb Cut-line
!       IF(Cell%Face1%in_out==-1.and.Cell%Face1%ec==-1.and.Cell%Face1%NumberVert==3) THEN
!          ! Cut_Face (Face1 und FCut)neu zusammenstellen, sollte aber schon in AnalyzeCell erfolgen
!          ! Face1 und FaceCut einzeln im einem 'general' listen 
!          nrPCut=Cell%vc
!          nrPFace=Cell%Face1%NumberVert
!          nFaces=nrPCut+2+nrPFace+2
!          ALLOCATE(fnvert(1:nFaces))
!          fnvert(1:2)=NrPCut
!          fnvert(nrPCut+3:nrPCut+4)=nrPFace
!          ifvc=1
!        ELSE IF(Cell%Face1%numberVert==4.AND.Cell%Face2%numberVert==4.AND. &
!                Cell%Face1%in_out==0 .AND. Cell%vc==0) THEN
!          nrPCut=4      ! Nummer Point Face1, resultiert anz. Flächen (F3-nFace) 
!          nFaces=nrPCut+2   ! nFaces = anz. Punkte(aus Celle-Face1=4) + (face1 und face2)
!          ALLOCATE(fnvert(1:nFaces))
!          IF(.NOT.ASSOCIATED(Cell%VertCut)) THEN
!            ALLOCATE(Cell%VertCut(nrPCut))
!          END IF
!          IF(.NOT.ASSOCIATED(Cell%VertCutOnlyP)) THEN
!            ALLOCATE(Cell%VertCutOnlyP(nrPCut))
!          END IF
!          fnvert(1:2)=nrPCut
!          Cell%VertCut(1:nrPCut)=Cell%Face1%VertexList(1:nrPCut)
!          ifvc=1
!       ELSE !Standard
!         nrPCut=Cell%vc     ! Nummer Point Cell%vc, resultiert anz. Flächen (F3-nFace) 
!         nFaces=nrPCut+2    ! nFaces = anz. Cell%vc + (face1 und face2)
!         ALLOCATE(fnvert(1:nFaces))
!         fnvert(1:2)=nrPCut
!         ifvc=1
!       END IF
!   ELSE  ! (vc==0 .OR. (vc>0.AND.Vol=0.0d0))
!       IF((k==1) .AND. & 
!          (Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
!           Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0) )THEN
!           !Face1, Grenzflaeche-Gelaende
!          nrPCut=4      ! Nummer Point Face1, resultiert anz. Flächen (F3-nFace) 
!          nFaces=nrPCut+2   ! nFaces = anz. Punkte(aus Celle-Face1=4) + (face1 und face2)
!          ALLOCATE(fnvert(1:nFaces))
!          IF(.NOT.ASSOCIATED(Cell%VertCut)) THEN
!            ALLOCATE(Cell%VertCut(nrPCut))
!          END IF
!          IF(.NOT.ASSOCIATED(Cell%VertCutOnlyP)) THEN
!            ALLOCATE(Cell%VertCutOnlyP(nrPCut))
!          END IF
!          fnvert(1:2)=nrPCut
!          Cell%VertCut(1:nrPCut)=Cell%Face1%VertexList(1:nrPCut)
!          ifvc=1
!       ELSE
!       END IF
!   END IF  ! IF (Cell%vc>0.AND.Cell%Vol>0.0d0.AND...<MAXVol... ELSE
!   !------------------------------------------------------------------
!   !------------------------------------------------------------------
!   IF(ifvc==1) THEN
!       ALLOCATE(LocCVertCut(nrPCut+1))      ! Local VertCut-Liste
!       IF(nrPCut>nrPFace) THEN
!            ALLOCATE(LocCVertCutOnlyP(nrPCut+1))  ! Local VertCutOnlyP
!       ELSE 
!            ALLOCATE(LocCVertCutOnlyP(nrPFace+1)) ! Local VertCutOnlyP
!       END IF
!       IF(nrPFace>0) THEN
!            ALLOCATE(LocCVertFace(nrPFace+1))
!       END IF
!       ALLOCATE(L_NrPCell(1:nFaces,1:10)) ! max. 8 Points möglich
!       ALLOCATE(NrPFaces(1:nFaces,1:10))  ! nrCutP je Face
!       ALLOCATE(cfxc(1:nFaces,1:10))      ! max. 8 Points möglich
!       ALLOCATE(cfyc(1:nFaces,1:10))
!       ALLOCATE(cfzc(1:nFaces,1:10))
!       L_NrPCell(:,:)=0
!       NrPFaces(:,:)=0
!       cfxc(:,:)=11111.11d0 !bessere Sicht Debugger
!       cfyc(:,:)=11111.11d0
!       cfzc(:,:)=11111.11d0
!
!       IF(nrPFace>0) THEN
!         ! NummerVerts aus aus Face1 und CutFace in Liste
!         LocCVertCut(1:nrPCut)=Cell%VertCut(1:nrPCut)
!         LocCVertFace(1:nrPFace)=Cell%Face1%VertexList(1:nrPFace)
!       ELSE
!         LocCVertCut(1:nrPCut)=Cell%VertCut(1:nrPCut)
!       END IF
!       DO isol=1,nr_lsoil
!         DO ivc=1,nrPCut
!           nrP=LocCVertCut(ivc)      !nrP=Cell%VertCut(ivc)
!           LocCVertCutOnlyP(ivc)=VertOut(nrP)%Vertex%nrCutP
!           L_NrPCell(2,ivc)=LocCVertCutOnlyP(ivc)+(isol-1)*nr_cutplane
!           L_NrPCell(1,ivc)=LocCVertCutOnlyP(ivc)+isol*nr_cutplane
!           NrPFaces(2,ivc)=L_NrPCell(2,ivc)
!           NrPFaces(1,ivc)=L_NrPCell(1,ivc)
!         END DO
!         NrPFaces(2,ivc)=L_NrPCell(2,1)  !ivc = nrPCut+1
!         NrPFaces(1,ivc)=L_NrPCell(1,1)  !ivc = nrPCut+1
!
!         IF(nrPFace>0) THEN
!           DO ivc=1,nrPFace
!             nrP=LocCVertFace(ivc)
!             LocCVertCutOnlyP(ivc)=VertOut(nrP)%Vertex%nrCutP
!             L_NrPCell(2,ivc)=LocCVertCutOnlyP(ivc)+(isol-1)*nr_cutplane
!             L_NrPCell(1,ivc)=LocCVertCutOnlyP(ivc)+isol*nr_cutplane
!             NrPFaces(nrPCut+2+2,ivc)=L_NrPCell(2,ivc)
!             NrPFaces(nrPCut+2+1,ivc)=L_NrPCell(1,ivc)
!           END DO
!           NrPFaces(nrPCut+2+2,ivc)=L_NrPCell(2,1)
!           NrPFaces(nrPCut+2+1,ivc)=L_NrPCell(1,1)
!         END IF
!
!         !Face3-Face(nFaces)
!         li=1
!         DO fakt=3,nrPCut+2
!           cfxc(fakt,1)=cfxc(1,li)
!           cfyc(fakt,1)=cfyc(1,li)
!           cfzc(fakt,1)=cfzc(1,li)
!           cfxc(fakt,2)=cfxc(1,li+1)
!           cfyc(fakt,2)=cfyc(1,li+1)
!           cfzc(fakt,2)=cfzc(1,li+1)
!           cfxc(fakt,3)=cfxc(2,li+1)
!           cfyc(fakt,3)=cfyc(2,li+1)
!           cfzc(fakt,3)=cfzc(2,li+1)
!           cfxc(fakt,4)=cfxc(2,li)
!           cfyc(fakt,4)=cfyc(2,li)
!           cfzc(fakt,4)=cfzc(2,li)
!           NrPFaces(fakt,1)=NrPFaces(1,li)
!           NrPFaces(fakt,2)=NrPFaces(1,li+1)
!           NrPFaces(fakt,3)=NrPFaces(2,li+1)
!           NrPFaces(fakt,4)=NrPFaces(2,li)
!           fnvert(fakt)=4 !+1 !nur symbolisch; Anfangs-Endwert gleich
!           li=li+1
!         END DO
!       
!         IF(nrPFace>0) THEN
!           li=1
!           DO fakt=nrPCut+2+3,nFaces
!             NrPFaces(fakt,1)=NrPFaces(nrPCut+2+1,li)
!             NrPFaces(fakt,2)=NrPFaces(nrPCut+2+1,li+1)
!             NrPFaces(fakt,3)=NrPFaces(nrPCut+2+2,li+1)
!             NrPFaces(fakt,4)=NrPFaces(nrPCut+2+2,li)
!             fnvert(fakt)=4 !+1 !nur symbolisch; Anfangs-Endwert gleich
!             li=li+1
!           END DO
!         END IF 
!
!         if(out_type=='a') then
!            !WRITE(10,*) i,j,k, "Celle"
!            WRITE(10,'(a8,I8)') general,nFaces
!            WRITE(10,*) (fnvert(iFace),iFace=1,nFaces)
!            DO iFace=1,nFaces
!               WRITE(10,*) (NrPFaces(iFace,iNrP),iNrP=1,fnvert(iFace)) !identisch, Ergebnis aus Analyze
!            END DO
!         else   !Output CellSoil-Layer to file *.Soil.out.gmvG for binary
!            !.. WRITE(10) general,nfaces
!            nRec=nRec+1
!            WRITE(10,REC=nRec) general(1:4)
!            nRec=nRec+1
!            WRITE(10,REC=nRec) general(5:8)
!            nRec=nRec+1
!            WRITE(10,REC=nRec) nFaces
!            !.. WRITE(10) (fnvert(iFace),iFace=1,nFaces)
!            DO iFace=1,nFaces
!              nRec=nRec+1
!              WRITE(10,REC=nRec) fnvert(iFace)
!            END DO
!            !.. WRITE(10,*) (L_NrPCell(iFace,iNrP),iNrP=1,fnvert(iFace)) 
!            DO iFace=1,nFaces
!              DO iNrP=1,fnvert(iFace)
!                 nRec=nRec+1
!                 !WRITE(10,REC=nRec) L_NrPCell(iFace,iNrP)
!                 WRITE(10,REC=nRec) NrPFaces(iFace,iNrP)
!              END DO
!            END DO
!         end if
!
!       END DO  !isol
!       IF (ASSOCIATED(LocCVertCutOnlyP)) THEN
!         DEALLOCATE(LocCVertCutOnlyP)
!       END IF
!       IF (ASSOCIATED(LocCVertCut)) THEN
!         DEALLOCATE(LocCVertCut)
!       END IF
!       IF (ASSOCIATED(LocCVertFace)) THEN
!         DEALLOCATE(LocCVertFace)
!       END IF
!       DEALLOCATE(fnvert)
!       DEALLOCATE(L_NrPCell)
!       DEALLOCATE(NrPFaces)
!       DEALLOCATE(cfxc)
!       DEALLOCATE(cfyc)
!       DEALLOCATE(cfzc)
!
!   END IF !IF (ifvc)
! END IF  ! (ASSOCIATED(Cell)) 
END SUBROUTINE WriteCellSoilLayerGMVBin
END MODULE OutputOutGmvGNeu_Mod
