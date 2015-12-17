MODULE GridNeu_Mod

  USE Function_Mod
  USE Floor_Mod
  USE Parametric_Mod
  USE IOControl_Mod
  USE BoundaryCond_Mod
  USE GridInput_Mod
  IMPLICIT NONE

  REAL(RealKind) :: dxLoc           ! (Domain%x1-Domain%x0)/Domain%nx*distx_coeff 
  REAL(RealKind) :: dyLoc           ! (Domain%y1-Domain%y0)/Domain%ny*disty_coeff
  REAL(RealKind) :: dzLoc           ! (Domain%z1-Domain%z0)/Domain%nz*distz_coeff

  TYPE (VertexP_T), POINTER :: VertOut(:)=>NULL()
  TYPE (VertexP_T), POINTER :: VertSoilOut(:)=>NULL()
  TYPE (VertexP_T), POINTER :: VertOutView(:)=>NULL()
  INTEGER :: nr_inside
  TYPE (VertexP_T), POINTER :: VertIn(:)=>NULL()
  INTEGER, PARAMETER :: LenListP=10000

  INTERFACE Allocate
    MODULE PROCEDURE AllocateVertex
  END INTERFACE

CONTAINS

SUBROUTINE AllocateVertex(Vertex,ix,iy,iz,Domain)
  TYPE(Vertex_T), POINTER :: Vertex
  INTEGER :: ix,iy,iz
  TYPE(Domain_T) :: Domain

  IF (.NOT.ASSOCIATED(Vertex)) THEN
    ALLOCATE(Vertex)
    Vertex%ix=ix
    Vertex%iy=iy
    Vertex%iz=iz
    Vertex%Point%x=Domain%xP(ix)
    Vertex%Point%y=Domain%yP(iy)
    Vertex%Point%z=Domain%zP(iz)
    Vertex%in_out=iz*1000
    Vertex%nrP=-1
    Vertex%nrInP=-1
    Vertex%nrCutP=-1
!   Vertex%nrI=-1
  END IF      
END SUBROUTINE AllocateVertex

SUBROUTINE InitAllVertices
  INTEGER :: ix,iy,iz
  INTEGER :: igx,igy,igz
  INTEGER :: jx,jy,jz
  INTEGER :: ib,in
  REAL(RealKind) :: iLoc
  REAL(RealKind) :: dxLoc,dyLoc,dzLoc
  
  ALLOCATE(Domain%Vertices(Domain%ix0:Domain%ix1 &
                          ,Domain%iy0:Domain%iy1 &
                          ,Domain%iz0:Domain%iz1)) 

  DO ib=1,nb
    CALL Set(Floor(ib))
    DO ix=ix0,ix1
      igx=ix*2**(-RefineX)
      DO iy=iy0,iy1
        igy=iy*2**(-RefineY)
        DO iz=iz0,iz1
          igz=iz*2**(-RefineZ)
          CALL Allocate(Domain%Vertices(igx,igy,igz)%Vertex,igx,igy,igz,Domain)
          Vertices(ix,iy,iz)%Vertex=>Domain%Vertices(igx,igy,igz)%Vertex
        END DO
      END DO
    END DO
    DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      IF (nType=='iw'.OR.nType=='ie') THEN
        jx=jx0
        igx=jNxO*2**(-RefineNachbarX)
        DO jy=jy0,jy1
          igy=jy*2**(-RefineY)
          DO jz=jz0,jz1
            igz=jz*2**(-RefineZ)
            CALL Allocate(Domain%Vertices(igx,igy,igz)%Vertex,igx,igy,igz,Domain)
            VerticesNachbar(jx,jy,jz)%Vertex=>Domain%Vertices(igx,igy,igz)%Vertex
          END DO
        END DO
        jx=jx1
        igx=jNxI*2**(-RefineNachbarX)
        DO jy=jy0,jy1
          igy=jy*2**(-RefineY)
          DO jz=jz0,jz1
            igz=jz*2**(-RefineZ)
            CALL Allocate(Domain%Vertices(igx,igy,igz)%Vertex,igx,igy,igz,Domain)
            VerticesNachbar(jx,jy,jz)%Vertex=>Domain%Vertices(igx,igy,igz)%Vertex
          END DO
        END DO
      END IF
      IF (nType=='is'.OR.nType=='in') THEN
        jy=jy0
        igy=jNyO*2**(-RefineNachbarY)
        DO jx=jx0,jx1
          igx=jx*2**(-RefineX)
          DO jz=jz0,jz1
            igz=jz*2**(-RefineZ)
            CALL Allocate(Domain%Vertices(igx,igy,igz)%Vertex,igx,igy,igz,Domain)
            VerticesNachbar(jx,jy,jz)%Vertex=>Domain%Vertices(igx,igy,igz)%Vertex
          END DO
        END DO
        jy=jy1
        igy=jNyI*2**(-RefineNachbarY)
        DO jx=jx0,jx1
          igx=jx*2**(-RefineX)
          DO jz=jz0,jz1
            igz=jz*2**(-RefineZ)
            CALL Allocate(Domain%Vertices(igx,igy,igz)%Vertex,igx,igy,igz,Domain)
            VerticesNachbar(jx,jy,jz)%Vertex=>Domain%Vertices(igx,igy,igz)%Vertex
          END DO
        END DO
      END IF
      IF (nType=='ib'.OR.nType=='it') THEN
        jz=jz0
        igz=jNzO*2**(-RefineNachbarZ)
        DO jx=jx0,jx1
          igx=jx*2**(-RefineX)
          DO jy=jy0,jy1
            igy=jy*2**(-RefineY)
            CALL Allocate(Domain%Vertices(igx,igy,igz)%Vertex,igx,igy,igz,Domain)
            VerticesNachbar(jx,jy,jz)%Vertex=>Domain%Vertices(igx,igy,igz)%Vertex
          END DO
        END DO
        jz=jz1
        igz=jNzI*2**(-RefineNachbarZ)
        DO jx=jx0,jx1
          igx=jx*2**(-RefineX)
          DO jy=jy0,jy1
            igy=jy*2**(-RefineY)
            CALL Allocate(Domain%Vertices(igx,igy,igz)%Vertex,igx,igy,igz,Domain)
            VerticesNachbar(jx,jy,jz)%Vertex=>Domain%Vertices(igx,igy,igz)%Vertex
          END DO
        END DO
      END IF
    END DO
  END DO

END SUBROUTINE InitAllVertices

SUBROUTINE AnalyzeAllVertices
  INTEGER :: ix,iy,iz,jx,jy,jz
  INTEGER :: ib,in
  INTEGER :: i1,i2,iNum
  INTEGER :: ListP(4,LenListP),NumList

  dxLoc=(Domain%x1-Domain%x0)/Domain%nx*distx_coeff
  dyLoc=(Domain%y1-Domain%y0)/Domain%ny*disty_coeff
  dzLoc=(Domain%z1-Domain%z0)/Domain%nz*distz_coeff
  DO iz=Domain%iz0,Domain%iz1
    DO iy=Domain%iy0,Domain%iy1
      DO ix=Domain%ix0,Domain%ix1
        IF (ASSOCIATED(Domain%Vertices(ix,iy,iz)%Vertex).AND.Domain%Vertices(ix,iy,iz)%Vertex%nrP==-1) THEN
          CALL AnalyzeVertex(Domain%Vertices(ix,iy,iz)%Vertex)
        END IF  
      END DO
    END DO
  END DO
         
  NumList=0       
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO ix=ix0,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
!         YZ Face
          CALL CheckVertexFace(Vertices(ix  ,iy-1,iz-1)%Vertex, &
                               Vertices(ix  ,iy  ,iz-1)%Vertex, &
                               Vertices(ix  ,iy  ,iz  )%Vertex, &
                               Vertices(ix  ,iy-1,iz  )%Vertex, &
                               ListP,NumList)
        END DO
      END DO
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0,iy1
        DO iz=iz0+1,iz1
!         ZX Face
          CALL CheckVertexFace(Vertices(ix-1,iy  ,iz-1)%Vertex, &
                               Vertices(ix  ,iy  ,iz-1)%Vertex, &
                               Vertices(ix  ,iy  ,iz  )%Vertex, &
                               Vertices(ix-1,iy  ,iz  )%Vertex, &
                               ListP,NumList)
        END DO
      END DO
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0,iz1
!         XY Face
          CALL CheckVertexFace(Vertices(ix-1,iy-1,iz  )%Vertex, &
                               Vertices(ix  ,iy-1,iz  )%Vertex, &
                               Vertices(ix  ,iy  ,iz  )%Vertex, &
                               Vertices(ix-1,iy  ,iz  )%Vertex, &
                               ListP,NumList)
        END DO
      END DO
    END DO
  END DO  
!   DO in=1,AnzahlNachbar
!     CALL Set(Nachbars(in))
!     IF (nType=='iw'.OR.nType=='ie') THEN
!       DO jy=jy0+1,jy1
!         DO jz=jz0+1,jz1
!           YZ Face
!           CALL CheckVertexFace(VerticesNachbar(jy-1,jz-1)%Vertex, &
!                                VerticesNachbar(jy  ,jz-1)%Vertex, &
!                                VerticesNachbar(jy  ,jz  )%Vertex, &
!                                VerticesNachbar(jy-1,jz  )%Vertex, &
!                                ListP,NumList)
!         END DO
!       END DO
!       DO jy=jy0,jy1
!         DO jz=jz0+1,jz1
!           ZX Face
!           CALL CheckVertexFace(VerticesNachbar(jy  ,jz-1)%Vertex, &
!                                VerticesNachbar(jy  ,jz  )%Vertex, &
!                                Vertices(jx0 ,jy  ,jz  )%Vertex, &
!                                Vertices(jx0 ,jy  ,jz-1)%Vertex, &
!                                ListP,NumList)
!         END DO
!       END DO
!       DO jy=jy0+1,jy1
!         DO jz=jz0,jz1
!           XY Face
!           CALL CheckVertexFace(VerticesNachbar(jy-1,jz  )%Vertex, &
!                                VerticesNachbar(jy  ,jz  )%Vertex, &
!                                Vertices(jx0 ,jy  ,jz  )%Vertex, &
!                                Vertices(jx0 ,jy-1,jz  )%Vertex, &
!                                ListP,NumList)
!         END DO
!       END DO
!     ELSE IF (nType=='is'.OR.nType=='in') THEN
!       DO jx=jx0,jx1
!         DO jz=jz0+1,jz1
!           YZ Face
!           CALL CheckVertexFace(VerticesNachbar(jx  ,jz-1)%Vertex, &
!                                VerticesNachbar(jx  ,jz  )%Vertex, &
!                                Vertices(jx  ,jy0 ,jz  )%Vertex, &
!                                Vertices(jx  ,jy0 ,jz-1)%Vertex, &
!                                ListP,NumList)
!         END DO
!       END DO
!       DO jx=jx0+1,jx1
!         DO jz=jz0+1,jz1
!           ZX Face
!           CALL CheckVertexFace(VerticesNachbar(jx-1,jz-1)%Vertex, &
!                                VerticesNachbar(jx  ,jz-1)%Vertex, &
!                                VerticesNachbar(jx  ,jz  )%Vertex, &
!                                VerticesNachbar(jx-1,jz  )%Vertex, &
!                                ListP,NumList)
!         END DO
!       END DO
!       DO jx=jx0+1,jx1
!         DO jz=jz0,jz1
!           XY Face
!           CALL CheckVertexFace(VerticesNachbar(jx-1,jz  )%Vertex, &
!                                VerticesNachbar(jx  ,jz  )%Vertex, &
!                                Vertices(jx  ,jy0 ,jz  )%Vertex, &
!                                Vertices(jx-1,jy0 ,jz  )%Vertex, &
!                                ListP,NumList)
!         END DO
!       END DO
!     ELSE IF (nType=='ib'.OR.nType=='it') THEN
!       DO jx=jx0,jx1
!         DO jz=jy0+1,jy1
!           YZ Face
!           CALL CheckVertexFace(VerticesNachbar(jx  ,jy-1)%Vertex, &
!                                VerticesNachbar(jx  ,jy  )%Vertex, &
!                                Vertices(jx  ,jy  ,jz0 )%Vertex, &
!                                Vertices(jx  ,jy-1,jz0 )%Vertex, &
!                                ListP,NumList)
!         END DO
!       END DO
!       DO jx=jx0+1,jx1
!         DO jz=jy0,jy1
!           ZX Face
!           CALL CheckVertexFace(VerticesNachbar(jx-1,jy  )%Vertex, &
!                                VerticesNachbar(jx  ,jy  )%Vertex, &
!                                Vertices(jx  ,jy  ,jz0 )%Vertex, &
!                                Vertices(jx-1,jy  ,jz0 )%Vertex, &
!                                ListP,NumList)
!         END DO
!       END DO
!       DO jx=jx0+1,jx1
!         DO jz=jy0+1,jy1
!           XY Face
!           CALL CheckVertexFace(VerticesNachbar(jx-1,jy-1)%Vertex, &
!                                VerticesNachbar(jx  ,jy-1)%Vertex, &
!                                VerticesNachbar(jx  ,jy  )%Vertex, &
!                                VerticesNachbar(jx-1,jy  )%Vertex, &
!                                ListP,NumList)
!         END DO
!       END DO
!     END IF  
!   END DO
! END DO  
  DO iNum=1,NumList
    ix=ListP(1,iNum)
    iy=ListP(2,iNum)
    iz=ListP(3,iNum)
    Domain%Vertices(ix,iy,iz)%Vertex%in_out=0
    Domain%Vertices(ix,iy,iz)%Vertex%Shift=0
    IF (Domain%Vertices(ix,iy,iz)%Vertex%nrP<=0) THEN
      nr_out=nr_out+1
      Domain%Vertices(ix,iy,iz)%Vertex%nrP=nr_out
    END IF
    IF (Domain%Vertices(ix,iy,iz)%Vertex%nrCutP<=0) THEN
      nr_cutplane=nr_cutplane+1
      Domain%Vertices(ix,iy,iz)%Vertex%nrCutP=nr_cutplane
    END IF
    IF (Domain%Vertices(ix,iy,iz)%Vertex%nrInP<=0) THEN
      nr_inside=nr_inside+1
      Domain%Vertices(ix,iy,iz)%Vertex%nrInP=nr_inside
    END IF
  END DO
END SUBROUTINE AnalyzeAllVertices

SUBROUTINE AnalyzeVertex(Vertex)
  TYPE (Vertex_T) :: Vertex
  REAL(8)  :: dist
  dist=Level(Vertex%Point%x,Vertex%Point%y,Vertex%Point%z)
  Vertex%dist=dist
  IF (dist>0.0d0) THEN
    dist=MIN(dist,Level(Vertex%Point%x+dxLoc,Vertex%Point%y,Vertex%Point%z))
    dist=MIN(dist,Level(Vertex%Point%x-dxLoc,Vertex%Point%y,Vertex%Point%z))
    dist=MIN(dist,Level(Vertex%Point%x,Vertex%Point%y+dyLoc,Vertex%Point%z))
    dist=MIN(dist,Level(Vertex%Point%x,Vertex%Point%y-dyLoc,Vertex%Point%z))
    dist=MIN(dist,Level(Vertex%Point%x,Vertex%Point%y,Vertex%Point%z+dzLoc))
    dist=MIN(dist,Level(Vertex%Point%x,Vertex%Point%y,Vertex%Point%z-dzLoc))
    IF (dist>0) THEN
      Vertex%in_out=1
    ELSE
      Vertex%in_out=0
      Vertex%Shift=0
    END IF
    IF (Vertex%Point%x>=domain%x0View-dxViewLoc.AND. &
        Vertex%Point%x<=domain%x1View+dxViewLoc.AND. &
        Vertex%Point%y>=domain%y0View-dyViewLoc.AND. &
        Vertex%Point%y<=domain%y1View+dyViewLoc.AND. &
        Vertex%Point%z>=domain%z0View-dzViewLoc.AND. &
        Vertex%Point%z<=domain%z1View+dzViewLoc) THEN
      nr_out=nr_out+1
      Vertex%nrP=nr_out
      IF(Vertex%in_out==0) THEN
         nr_cutplane=nr_cutplane+1
         Vertex%nrCutP=nr_cutplane
         nr_inside=nr_inside+1
         Vertex%nrInP=nr_inside
      END IF
    ELSE
      Vertex%nrP=-2
      Vertex%nrInP=-2
      Vertex%nrCutP=-2
    END IF
  ELSE IF (dist>-dist_fscv) THEN  ! Standard -1.0d-12
    Vertex%in_out=0
    IF (Vertex%Point%x>=domain%x0View-dxViewLoc.AND. &
        Vertex%Point%x<=domain%x1View+dxViewLoc.AND. &
        Vertex%Point%y>=domain%y0View-dyViewLoc.AND. &
        Vertex%Point%y<=domain%y1View+dyViewLoc.AND. &
        Vertex%Point%z>=domain%z0View-dzViewLoc.AND. &
        Vertex%Point%z<=domain%z1View+dzViewLoc) THEN
      nr_out=nr_out+1
      Vertex%nrP=nr_out
      nr_cutplane=nr_cutplane+1

      Vertex%nrCutP=nr_cutplane
      nr_inside=nr_inside+1
      Vertex%nrInP=nr_inside
    ELSE
      Vertex%nrP=-2
      Vertex%nrInP=-2
      Vertex%nrCutP=-2
    END IF
  ELSE
    Vertex%in_out=-1
    Vertex%nrP=0
    Vertex%nrCutP=0
    IF (Vertex%Point%x>=domain%x0View-dxViewLoc.AND. &
        Vertex%Point%x<=domain%x1View+dxViewLoc.AND. &
        Vertex%Point%y>=domain%y0View-dyViewLoc.AND. &
        Vertex%Point%y<=domain%y1View+dyViewLoc.AND. &
        Vertex%Point%z>=domain%z0View-dzViewLoc.AND. &
        Vertex%Point%z<=domain%z1View+dzViewLoc) THEN
      nr_inside=nr_inside+1
      Vertex%nrInP=nr_inside
    ELSE
      Vertex%nrInP=-2
    END IF
  END IF

END SUBROUTINE AnalyzeVertex

SUBROUTINE CheckVertexFace(V1,V2,V3,V4,List,ListPos)
  TYPE(Vertex_T) :: V1,V2,V3,V4
  INTEGER :: List(:,:)
  INTEGER :: ListPos

  IF (V1%in_out==-1.AND.V3%in_out==-1.AND. &
      V2%in_out>=0.AND.V4%in_out>=0) THEN
    ListPos=ListPos+1
    IF (ListPos>LenListP) THEN
      WRITE(*,*) 'ListPos ',ListPos
      WRITE(*,*) 'LenListP',LenListP 
      STOP 'ListPos>LenListP'
    END IF  
    List(1,ListPos)=V1%ix
    List(2,ListPos)=V1%iy
    List(3,ListPos)=V1%iz
    List(4,ListPos)=V1%iType
    ListPos=ListPos+1
    List(1,ListPos)=V3%ix
    List(2,ListPos)=V3%iy
    List(3,ListPos)=V3%iz
    List(4,ListPos)=V3%iType
  ELSE IF (V2%in_out==-1.AND.V4%in_out==-1.AND. &
      V1%in_out>=0.AND.V3%in_out>=0) THEN
    ListPos=ListPos+1
    IF (ListPos>LenListP) THEN
      WRITE(*,*) 'ListPos ',ListPos
      WRITE(*,*) 'LenListP',LenListP 
      STOP 'ListPos>LenListP'
    END IF  
    List(1,ListPos)=V2%ix
    List(2,ListPos)=V2%iy
    List(3,ListPos)=V2%iz
    List(4,ListPos)=V2%iType
    ListPos=ListPos+1
    List(1,ListPos)=V4%ix
    List(2,ListPos)=V4%iy
    List(3,ListPos)=V4%iz
    List(4,ListPos)=V4%iType
  END IF  
END SUBROUTINE CheckVertexFace

SUBROUTINE InitAll
  INTEGER :: ib,in

  DO ib=1,nb
    CALL Set(Floor(ib))
    CALL InitBlock(ix0,ix1,iy0,iy1,iz0,iz1          &
                       ,Vertices,Edges_X,Edges_Y,Edges_Z &
                       ,Faces_XY,Faces_YZ,Faces_ZX       &
                       ,Cells)
    DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      IF (nType(1:1)=='i') THEN
        CALL InitBlock(jx0,jx1,jy0,jy1,jz0,jz1          &
                      ,VerticesNachbar,Edges_XNachbar,Edges_YNachbar,Edges_ZNachbar &
                      ,Faces_XYNachbar,Faces_YZNachbar,Faces_ZXNachbar       &
                      ,CellsNachbar)
      END IF  
    END DO  
  END DO  
END SUBROUTINE InitAll

SUBROUTINE InitAllEdges
  INTEGER :: ib,in,Iter
  INTEGER :: is,js,ks

  DO ib=1,nb
    CALL Set(Floor(ib))
    DO Iter=1,3
      IF (Iter==1) THEN
        is=1
        js=0
        ks=0
      ELSE IF (Iter==2) THEN
        is=0
        js=1
        ks=0
      ELSE
        is=0
        js=0
        ks=1
      END IF
      CALL InitEdgesBlock(ix0,ix1,iy0,iy1,iz0,iz1          &
                         ,is,js,ks                         &
                         ,Vertices,Edges_X,Edges_Y,Edges_Z &
                         ,Faces_XY,Faces_YZ,Faces_ZX       &
                         ,Cells)
      DO in=1,AnzahlNachbar
        CALL Set(Nachbars(in))
        IF (nType(1:1)=='i') THEN
          CALL InitEdgesBlock(jx0,jx1,jy0,jy1,jz0,jz1          &
                           ,is,js,ks                         &
                           ,VerticesNachbar,Edges_XNachbar,Edges_YNachbar,Edges_ZNachbar &
                           ,Faces_XYNachbar,Faces_YZNachbar,Faces_ZXNachbar       &
                           ,CellsNachbar)
        END IF  
      END DO  
    END DO  
  END DO  
END SUBROUTINE InitAllEdges

SUBROUTINE InitBlock(ix0,ix1,iy0,iy1,iz0,iz1          &
                         ,Vertices,Edges_X,Edges_Y,Edges_Z &
                         ,Faces_XY,Faces_YZ,Faces_ZX       &
                         ,Cells,Print)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  TYPE(VertexP_T) :: Vertices(ix0:ix1,iy0:iy1,iz0:iz1)
  TYPE(EdgeP_T) :: Edges_X(ix0+1:ix1,iy0:iy1,iz0:iz1)
  TYPE(EdgeP_T) :: Edges_Y(ix0:ix1,iy0+1:iy1,iz0:iz1)
  TYPE(EdgeP_T) :: Edges_Z(ix0:ix1,iy0:iy1,iz0+1:iz1)
  TYPE(FaceP_T) :: Faces_XY(ix0+1:ix1,iy0+1:iy1,iz0:iz1)
  TYPE(FaceP_T) :: Faces_YZ(ix0:ix1,iy0+1:iy1,iz0+1:iz1)
  TYPE(FaceP_T) :: Faces_ZX(ix0+1:ix1,iy0:iy1,iz0+1:iz1)
  TYPE(CellP_T) :: Cells(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)
  LOGICAL, OPTIONAL :: Print

  INTEGER :: ix,iy,iz
  INTEGER :: MinIn_Out

! Init Cells
  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        MinIn_Out=MIN(Vertices(ix-1,iy-1,iz-1)%Vertex%in_out &
                     ,Vertices(ix  ,iy-1,iz-1)%Vertex%in_out &
                     ,Vertices(ix-1,iy  ,iz-1)%Vertex%in_out &
                     ,Vertices(ix-1,iy-1,iz  )%Vertex%in_out &
                     ,Vertices(ix-1,iy  ,iz  )%Vertex%in_out &
                     ,Vertices(ix  ,iy-1,iz  )%Vertex%in_out &
                     ,Vertices(ix  ,iy  ,iz-1)%Vertex%in_out &
                     ,Vertices(ix  ,iy  ,iz  )%Vertex%in_out)
        IF (MinIn_Out<0) THEN 
          IF (PRESENT(Print)) THEN
            WRITE(*,*) 'Cell MinIn_Out',ix,iy,iz
          END IF  
          CALL AllocateEdge(Edges_X(ix,iy-1,iz-1)%Edge &
                           ,Vertices(ix-1,iy-1,iz-1)%Vertex &
                           ,Vertices(ix  ,iy-1,iz-1)%Vertex,'X')
          CALL AllocateEdge(Edges_X(ix,iy  ,iz-1)%Edge &
                           ,Vertices(ix-1,iy  ,iz-1)%Vertex &
                           ,Vertices(ix  ,iy  ,iz-1)%Vertex,'X')
          CALL AllocateEdge(Edges_X(ix,iy-1,iz  )%Edge &
                           ,Vertices(ix-1,iy-1,iz  )%Vertex &
                           ,Vertices(ix  ,iy-1,iz  )%Vertex,'X')
          CALL AllocateEdge(Edges_X(ix,iy  ,iz  )%Edge &
                           ,Vertices(ix-1,iy  ,iz  )%Vertex &
                           ,Vertices(ix  ,iy  ,iz  )%Vertex,'X')
          CALL AllocateEdge(Edges_Y(ix-1,iy,iz-1)%Edge &
                           ,Vertices(ix-1,iy-1,iz-1)%Vertex &
                           ,Vertices(ix-1,iy  ,iz-1)%Vertex,'Y')
          CALL AllocateEdge(Edges_Y(ix  ,iy,iz-1)%Edge &
                           ,Vertices(ix  ,iy-1,iz-1)%Vertex &
                           ,Vertices(ix  ,iy  ,iz-1)%Vertex,'Y')
          CALL AllocateEdge(Edges_Y(ix-1,iy,iz  )%Edge &
                           ,Vertices(ix-1,iy-1,iz  )%Vertex &
                           ,Vertices(ix-1,iy  ,iz  )%Vertex,'Y')
          CALL AllocateEdge(Edges_Y(ix  ,iy,iz  )%Edge &
                           ,Vertices(ix  ,iy-1,iz  )%Vertex &
                           ,Vertices(ix  ,iy  ,iz  )%Vertex,'Y')
          CALL AllocateEdge(Edges_Z(ix-1,iy-1,iz)%Edge &
                           ,Vertices(ix-1,iy-1,iz-1)%Vertex &
                           ,Vertices(ix-1,iy-1,iz  )%Vertex,'Z')
          CALL AllocateEdge(Edges_Z(ix  ,iy-1,iz)%Edge &
                           ,Vertices(ix  ,iy-1,iz-1)%Vertex &
                           ,Vertices(ix  ,iy-1,iz  )%Vertex,'Z')
          CALL AllocateEdge(Edges_Z(ix-1,iy  ,iz)%Edge &
                           ,Vertices(ix-1,iy  ,iz-1)%Vertex &
                           ,Vertices(ix-1,iy  ,iz  )%Vertex,'Z')
          CALL AllocateEdge(Edges_Z(ix  ,iy  ,iz)%Edge &
                           ,Vertices(ix  ,iy  ,iz-1)%Vertex &
                           ,Vertices(ix  ,iy  ,iz  )%Vertex,'Z')
          CALL AllocateFace(Faces_XY(ix,iy,iz-1)%Face &
                           ,Edges_X(ix  ,iy-1,iz-1)%Edge   &
                           ,Edges_Y(ix  ,iy  ,iz-1)%Edge &
                           ,Edges_X(ix  ,iy  ,iz-1)%Edge   &
                           ,Edges_Y(ix-1,iy  ,iz-1)%Edge &
                           ,VolFace_XY,'XY')
          CALL AllocateFace(Faces_XY(ix,iy,iz  )%Face &
                           ,Edges_X(ix  ,iy-1,iz  )%Edge   &
                           ,Edges_Y(ix  ,iy  ,iz  )%Edge &
                           ,Edges_X(ix  ,iy  ,iz  )%Edge   &
                           ,Edges_Y(ix-1,iy  ,iz  )%Edge &
                           ,VolFace_XY,'XY')
          CALL AllocateFace(Faces_YZ(ix-1,iy,iz)%Face &
                           ,Edges_Y(ix-1,iy  ,iz-1)%Edge   &
                           ,Edges_Z(ix-1,iy  ,iz  )%Edge &
                           ,Edges_Y(ix-1,iy  ,iz  )%Edge   &
                           ,Edges_Z(ix-1,iy-1,iz  )%Edge &
                           ,VolFace_YZ,'YZ')
          CALL AllocateFace(Faces_YZ(ix  ,iy,iz)%Face &
                           ,Edges_Y(ix  ,iy  ,iz-1)%Edge   &
                           ,Edges_Z(ix  ,iy  ,iz  )%Edge &
                           ,Edges_Y(ix  ,iy  ,iz  )%Edge   &
                           ,Edges_Z(ix  ,iy-1,iz  )%Edge &
                           ,VolFace_YZ,'YZ')
          CALL AllocateFace(Faces_ZX(ix,iy-1,iz)%Face &
                           ,Edges_Z(ix-1,iy-1,iz  )%Edge   &
                           ,Edges_X(ix  ,iy-1,iz  )%Edge &
                           ,Edges_Z(ix  ,iy-1,iz  )%Edge   &
                           ,Edges_X(ix  ,iy-1,iz-1)%Edge &
                           ,VolFace_ZX,'ZX')
          CALL AllocateFace(Faces_ZX(ix,iy  ,iz)%Face &
                           ,Edges_Z(ix-1,iy  ,iz  )%Edge   &
                           ,Edges_X(ix  ,iy  ,iz  )%Edge &
                           ,Edges_Z(ix  ,iy  ,iz  )%Edge   &
                           ,Edges_X(ix  ,iy  ,iz-1)%Edge &
                           ,VolFace_ZX,'ZX')
          CALL AllocateCell(Cells(ix,iy,iz)%Cell &
                           ,Faces_XY(ix,iy,iz-1)%Face &
                           ,Faces_XY(ix,iy,iz)%Face &
                           ,Faces_ZX(ix,iy-1,iz)%Face &
                           ,Faces_ZX(ix,iy,iz)%Face &
                           ,Faces_YZ(ix-1,iy,iz)%Face &
                           ,Faces_YZ(ix,iy,iz)%Face &
                           ,VolCell)
        END IF
      END DO
    END DO
  END DO
  IF (PRESENT(PRINT)) THEN
    WRITE(*,*) 'ASSOC IN',ASSOCIATED(Faces_YZ(17,13,6)%Face)
    WRITE(*,*) 'ASSOC IN',LBOUND(Faces_YZ,1),UBOUND(Faces_YZ,1)
    WRITE(*,*) 'ASSOC IN',LBOUND(Faces_YZ,2),UBOUND(Faces_YZ,2)
    WRITE(*,*) 'ASSOC IN',LBOUND(Faces_YZ,3),UBOUND(Faces_YZ,3)
  END IF  
END SUBROUTINE InitBlock



SUBROUTINE InitEdgesBlock(ix0,ix1,iy0,iy1,iz0,iz1          &
                         ,is,js,ks                         &
                         ,Vertices,Edges_X,Edges_Y,Edges_Z &
                         ,Faces_XY,Faces_YZ,Faces_ZX       &
                         ,Cells,Print)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  INTEGER :: is,js,ks
  TYPE(VertexP_T) :: Vertices(ix0:ix1,iy0:iy1,iz0:iz1)
  TYPE(EdgeP_T) :: Edges_X(ix0+1:ix1,iy0:iy1,iz0:iz1)
  TYPE(EdgeP_T) :: Edges_Y(ix0:ix1,iy0+1:iy1,iz0:iz1)
  TYPE(EdgeP_T) :: Edges_Z(ix0:ix1,iy0:iy1,iz0+1:iz1)
  TYPE(FaceP_T) :: Faces_XY(ix0+1:ix1,iy0+1:iy1,iz0:iz1)
  TYPE(FaceP_T) :: Faces_YZ(ix0:ix1,iy0+1:iy1,iz0+1:iz1)
  TYPE(FaceP_T) :: Faces_ZX(ix0+1:ix1,iy0:iy1,iz0+1:iz1)
  TYPE(CellP_T) :: Cells(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)
  LOGICAL, OPTIONAL :: Print

  INTEGER :: i,j,k
  INTEGER :: ii,jj,kk

  DO i=ix0+is,ix1
    DO j=iy0+js,iy1
      DO k=iz0+ks,iz1
        IF ((Vertices(i-is,j-js,k-ks)%Vertex%in_out*Vertices(i,j,k)%Vertex%in_out<=0 .AND.&
             MIN(Vertices(i-is,j-js,k-ks)%Vertex%in_out,Vertices(i,j,k)%Vertex%in_out)<0 .AND. &
             MIN(Vertices(i-is,j-js,k-ks)%Vertex%nrP,Vertices(i,j,k)%Vertex%nrP)>-2)) THEN
          DO ii=MAX(i,ix0+1),MIN(i+1-is,ix1)
            DO jj=MAX(j-1,iy0),MIN(j+is+ks,iy1)
              DO kk=MAX(k-1,iz0),MIN(k+is+js,iz1)
                CALL AllocateEdge(Edges_X(ii,jj,kk)%Edge &
                                 ,Vertices(ii-1,jj,kk)%Vertex,Vertices(ii,jj,kk)%Vertex,'X')
              END DO
            END DO
          END DO
          DO ii=MAX(i-1,ix0),MIN(i+js+ks,ix1)
            DO jj=MAX(j,iy0+1),MIN(j+1-js,iy1)
              DO kk=MAX(k-1,iz0),MIN(k+js+is,iz1)
               CALL AllocateEdge(Edges_Y(ii,jj,kk)%Edge &
                                ,Vertices(ii,jj-1,kk)%Vertex,Vertices(ii,jj,kk)%Vertex,'Y')
              END DO
            END DO
          END DO
          DO ii=MAX(i-1,ix0),MIN(i+ks+js,ix1)
            DO jj=MAX(j-1,iy0),MIN(j+ks+is,iy1)
              DO kk=MAX(k,iz0+1),MIN(k+1-ks,iz1)
                CALL AllocateEdge(Edges_Z(ii,jj,kk)%Edge &
                                 ,Vertices(ii,jj,kk-1)%Vertex,Vertices(ii,jj,kk)%Vertex,'Z')
             END DO
            END DO
          END DO
          DO ii=MAX(i-1,ix0),MIN(i+js+ks,ix1)
            DO jj=MAX(j,iy0+1),MIN(j+ks+is,iy1)
              DO kk=MAX(k,iz0+1),MIN(k+js+is,iz1)
                CALL AllocateFace(Faces_YZ(ii,jj,kk)%Face &
                                 ,Edges_Y(ii,jj,kk-1)%Edge   &
                                 ,Edges_Z(ii,jj,kk)%Edge &
                                 ,Edges_Y(ii,jj,kk)%Edge   &
                                 ,Edges_Z(ii,jj-1,kk)%Edge &
                                 ,VolFace_YZ,'YZ')
              END DO
            END DO
          END DO
          DO ii=MAX(i,ix0+1),MIN(i+js+ks,ix1)
            DO jj=MAX(j,iy0+1),MIN(j+is+ks,iy1)
              DO kk=MAX(k-1,iz0),MIN(k+is+js,iz1)
                CALL AllocateFace(Faces_XY(ii,jj,kk)%Face &
                                 ,Edges_X(ii,jj-1,kk)%Edge   &
                                 ,Edges_Y(ii,jj,kk)%Edge &
                                 ,Edges_X(ii,jj,kk)%Edge   &
                                 ,Edges_Y(ii-1,jj,kk)%Edge &
                                 ,VolFace_XY,'XY')
              END DO
            END DO
          END DO
          DO ii=MAX(i,ix0+1),MIN(i+ks+js,ix1)
            DO jj=MAX(j-1,iy0),MIN(j+is+ks,iy1)
              DO kk=MAX(k,iz0+1),MIN(k+is+js,iz1)
                IF (PRESENT(Print)) WRITE(*,*) 'AllocateFace_ZX',ii,jj,kk
                CALL AllocateFace(Faces_ZX(ii,jj,kk)%Face &
                                 ,Edges_Z(ii-1,jj,kk)%Edge   &
                                 ,Edges_X(ii,jj,kk)%Edge &
                                 ,Edges_Z(ii,jj,kk)%Edge   &
                                 ,Edges_X(ii,jj,kk-1)%Edge &
                                 ,VolFace_ZX,'ZX')
              END DO
            END DO
          END DO
          DO ii=MAX(i,ix0+1),MIN(i+js+ks,ix1)
            DO jj=MAX(j,iy0+1),MIN(j+is+ks,iy1)
              DO kk=MAX(k,iz0+1),MIN(k+is+js,iz1)
                CALL AllocateCell(Cells(ii,jj,kk)%Cell &
                                 ,Faces_XY(ii,jj,kk-1)%Face &
                                 ,Faces_XY(ii,jj,kk)%Face &
                                 ,Faces_ZX(ii,jj-1,kk)%Face &
                                 ,Faces_ZX(ii,jj,kk)%Face &
                                 ,Faces_YZ(ii-1,jj,kk)%Face &
                                 ,Faces_YZ(ii,jj,kk)%Face &
                                 ,VolCell)
              END DO
            END DO
          END DO
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE InitEdgesBlock

SUBROUTINE AnalyzeAllEdges

  INTEGER :: ix,iy,iz,jx,jy,jz
  INTEGER :: ib,in

! Suche Schnittpunkte
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO iy=iy0,iy1
        DO iz=iz0,iz1
          IF (ASSOCIATED(Edges_X(ix,iy,iz)%Edge)) THEN
            CALL AnalyzeEdge(Edges_X(ix,iy,iz)%Edge)
          END IF
        END DO
      END DO
    END DO
    DO ix=ix0,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0,iz1
          IF (ASSOCIATED(Edges_Y(ix,iy,iz)%Edge)) THEN
            CALL AnalyzeEdge(Edges_Y(ix,iy,iz)%Edge)
          END IF
        END DO
      END DO
    END DO
    DO ix=ix0,ix1
      DO iy=iy0,iy1
        DO iz=iz0+1,iz1
          IF (ASSOCIATED(Edges_Z(ix,iy,iz)%Edge)) THEN
            CALL AnalyzeEdge(Edges_Z(ix,iy,iz)%Edge)
          END IF
        END DO
      END DO
    END DO
    DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      IF (nType(1:1)=='i') THEN
        DO jx=jx0+1,jx1
          DO jy=jy0,jy1
            DO jz=jz0,jz1
              IF (ASSOCIATED(Edges_XNachbar(jx,jy,jz)%Edge)) THEN
                CALL AnalyzeEdge(Edges_XNachbar(jx,jy,jz)%Edge)
              END IF
            END DO
          END DO
        END DO
        DO jx=jx0,jx1
          DO jy=jy0+1,jy1
            DO jz=jz0,jz1
              IF (ASSOCIATED(Edges_YNachbar(jx,jy,jz)%Edge)) THEN
                CALL AnalyzeEdge(Edges_YNachbar(jx,jy,jz)%Edge)
              END IF
            END DO
          END DO
        END DO
        DO jx=jx0,jx1
          DO jy=jy0,jy1
            DO jz=jz0+1,jz1
              IF (ASSOCIATED(Edges_ZNachbar(jx,jy,jz)%Edge)) THEN
                CALL AnalyzeEdge(Edges_ZNachbar(jx,jy,jz)%Edge)
              END IF
            END DO
          END DO
        END DO
      END IF  
    END DO  
  END DO  
END SUBROUTINE AnalyzeAllEdges

SUBROUTINE AnalyzeEdge(Edge,Print)
  TYPE(Edge_T) :: Edge
  LOGICAL, OPTIONAL :: Print

  !local variable:
  TYPE (Vertex_T) :: VertS
  TYPE (Vertex_T), POINTER :: VertSP

  Edge%in_out=Edge%Vert1%in_out+Edge%Vert2%in_out
  IF (PRESENT(Print)) WRITE(*,*) ASSOCIATED(Edge%Vert1%VertZ)
  IF ((Edge%Vert1%ix<Edge%Vert2%ix.AND..NOT.ASSOCIATED(Edge%Vert1%VertX)).OR. &
      (Edge%Vert1%iy<Edge%Vert2%iy.AND..NOT.ASSOCIATED(Edge%Vert1%VertY)).OR. &
      (Edge%Vert1%iz<Edge%Vert2%iz.AND..NOT.ASSOCIATED(Edge%Vert1%VertZ))) THEN
  IF (PRESENT(Print)) WRITE(*,*) 'IF 1'
    IF (Edge%Vert1%in_out+Edge%Vert2%in_out==-1.OR.Edge%Vert1%in_out*Edge%Vert2%in_out==-1) THEN
      IF (Edge%Vert1%in_out*Edge%Vert2%in_out==-1) THEN
        IF (PRESENT(Print)) WRITE(*,*) 'CALL SearchPoint'
        CALL SearchPoint(VertS,Edge%Vert1,Edge%Vert2)
      ELSE  
        IF (Edge%Vert1%in_out==0) THEN
          VertS%Point=ShiftPoint*Edge%Vert1%Point+(1.0d0-ShiftPoint)*Edge%Vert2%Point
        ELSE  
          VertS%Point=ShiftPoint*Edge%Vert2%Point+(1.0d0-ShiftPoint)*Edge%Vert1%Point
        END IF  
      END IF  
      IF (PRESENT(Print)) WRITE(*,*) 'SearchPoint',VertS%Point  
      IF (Edge%Vert1%Point%x==Edge%Vert2%Point%x) THEN
        VertS%Point%x=Edge%Vert1%Point%x
      END IF
      IF (Edge%Vert1%Point%y==Edge%Vert2%Point%y) THEN
        VertS%Point%y=Edge%Vert1%Point%y
      END IF
      IF (Edge%Vert1%Point%z==Edge%Vert2%Point%z) THEN
        VertS%Point%z=Edge%Vert1%Point%z
      END IF
      IF (VertS%Point%x>=domain%x0View-dxViewLoc .AND. &
          VertS%Point%x<=domain%x1View+dxViewLoc .AND. &
          VertS%Point%y>=domain%y0View-dyViewLoc .AND. &
          VertS%Point%y<=domain%y1View+dyViewLoc .AND. &
          VertS%Point%z>=domain%z0View-dzViewLoc .AND. &
          VertS%Point%z<=domain%z1View+dzViewLoc) THEN
        nr_out=nr_out+1
        VertS%nrP=nr_out
        nr_cutplane=nr_cutplane+1
        VertS%nrCutP=nr_cutplane
        nr_inside=nr_inside+1
        VertS%nrInP=nr_inside
      ELSE
        VertS%nrP=-2
        VertS%nrInP=-2
        VertS%nrCutP=-2
      END IF
      ALLOCATE(VertSP)
      VertSP%Point=VertS%Point
      VertSP%nrP=VertS%nrP
      VertSP%nrInP=VertS%nrInp
      VertSP%nrCutP=VertS%nrCutP
!     VertSP%nrI=VertS%nrI
      VertSP%ix=Edge%Vert1%ix
      VertSP%iy=Edge%Vert1%iy
      VertSP%iz=Edge%Vert1%iz
      IF (Edge%Vert1%ix<Edge%Vert2%ix) THEN
        Edge%Vert1%VertX=>VertSP
      ELSE IF (Edge%Vert1%iy<Edge%Vert2%iy) THEN
        Edge%Vert1%VertY=>VertSP
      ELSE IF (Edge%Vert1%iz<Edge%Vert2%iz) THEN
        Edge%Vert1%VertZ=>VertSP
      END IF  
      Edge%yes_sp=1
    ELSE
      Edge%yes_sp=0
    END IF
  ELSE IF ((Edge%Vert1%ix<Edge%Vert2%ix.AND.ASSOCIATED(Edge%Vert1%VertX)).OR. &
           (Edge%Vert1%iy<Edge%Vert2%iy.AND.ASSOCIATED(Edge%Vert1%VertY)).OR. &
           (Edge%Vert1%iz<Edge%Vert2%iz.AND.ASSOCIATED(Edge%Vert1%VertZ))) THEN
    Edge%yes_sp=1
  END IF
END SUBROUTINE AnalyzeEdge

SUBROUTINE AnalyzeAllFaces

  INTEGER :: ib,in
  INTEGER :: in_out

  DO ib=1,nb
    CALL Set(Floor(ib))
    CALL AnalyzeFacesBlock(ix0,ix1,iy0,iy1,iz0,iz1          &
                          ,Faces_XY,Faces_YZ,Faces_ZX)
    DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      IF (nType(1:1)=='i') THEN
        CALL AnalyzeFacesBlock(jx0,jx1,jy0,jy1,jz0,jz1          &
                              ,Faces_XYNachbar,Faces_YZNachbar,Faces_ZXNachbar)
      END IF  
    END DO
  END DO
END SUBROUTINE AnalyzeAllFaces

SUBROUTINE AnalyzeFacesBlock(ix0,ix1,iy0,iy1,iz0,iz1        &
                          ,Faces_XY,Faces_YZ,Faces_ZX,Print)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  TYPE(FaceP_T) :: Faces_XY(ix0+1:ix1,iy0+1:iy1,iz0:iz1)
  TYPE(FaceP_T) :: Faces_YZ(ix0:ix1,iy0+1:iy1,iz0+1:iz1)
  TYPE(FaceP_T) :: Faces_ZX(ix0+1:ix1,iy0:iy1,iz0+1:iz1)
  LOGICAL, OPTIONAL :: Print

  INTEGER :: ib,ix,iy,iz
  INTEGER :: in_out

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0,iz1
        IF (ASSOCIATED(Faces_XY(ix,iy,iz)%Face)) THEN
          IF (PRESENT(Print)) WRITE(*,*) 'ASSOC Faces_XY',ix,iy,iz
!         CALL AnalyzeFace(Faces_XY(ix,iy,iz)%Face,.TRUE.)
          CALL AnalyzeFace(Faces_XY(ix,iy,iz)%Face)
          IF (PRESENT(Print)) WRITE(*,*) 'nach ASSOC Faces_XY',ix,iy,iz
        END IF
      END DO
    END DO
  END DO

  DO ix=ix0,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        IF (ASSOCIATED(Faces_YZ(ix,iy,iz)%Face)) THEN
        IF (PRESENT(Print)) WRITE(*,*) 'ASSOC Faces_YZ',ix,iy,iz
!         CALL AnalyzeFace(Faces_YZ(ix,iy,iz)%Face,.TRUE.)
          CALL AnalyzeFace(Faces_YZ(ix,iy,iz)%Face)
        IF (PRESENT(Print)) WRITE(*,*) 'nach ASSOC Faces_YZ',ix,iy,iz
        END IF
      END DO
    END DO
  END DO

  DO ix=ix0+1,ix1
    DO iy=iy0,iy1
      DO iz=iz0+1,iz1
        IF (ASSOCIATED(Faces_ZX(ix,iy,iz)%Face)) THEN
          IF (PRESENT(Print)) WRITE(*,*) 'ASSOC Faces_ZX',ix,iy,iz
!         CALL AnalyzeFace(Faces_ZX(ix,iy,iz)%Face,.TRUE.)
          CALL AnalyzeFace(Faces_ZX(ix,iy,iz)%Face)
          IF (PRESENT(Print)) WRITE(*,*) 'nach ASSOC Faces_ZX',ix,iy,iz
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE AnalyzeFacesBlock

SUBROUTINE AnalyzeFace(Face,Print)
  TYPE (Face_T) :: Face
  LOGICAL, OPTIONAL :: Print

  INTEGER :: NrW,NrMP

  INTEGER :: i,j,v,w,ecut,iTemp,tdiff
  INTEGER :: ic,jc,kc,nxc,nyc,nzc
  REAL(8) :: xc0,xc1,yc0,yc1,zc0,zc1,dx,dy,dz
  TYPE (Vertex_T) :: InterVert
  TYPE (Vertex_T) :: Vertex1,Vertex2
  TYPE(Point_T) :: EdgP0,EdgP1,EdgP2,EdgP3
  TYPE(Point_T) :: P0,P1,P2,P3
  TYPE(Point_T) :: PMin,PMax,SumMidP,MidPoint
  REAL(8) :: xvc0,xvc1,xvc0u,xvc0o,xvc1u,xvc1o
  REAL(8) :: yvc0,yvc1,yvc0u,yvc0o,yvc1u,yvc1o
  REAL(8) :: zvc0,zvc1
  REAL(8) :: Vol,Len
  REAL(8) :: xsi,ys
  INTEGER :: n_egl
  INTEGER :: EdgeList(1:2,1:6) 

  NrMP=0
  Face%ec=-1
  Face%NumberVert=0
  Face%VertexList=0
  Face%EdgeCut(:)=0
  i=0
  ecut=0
  Face%in_out=Face%Edge1%in_out+Face%Edge3%in_out
  !...............................................
  IF (Face%Edge1%Vert1%in_out>=0) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge1 In'
    i=i+1
    Face%VertexList(i)=Face%Edge1%Vert1%nrP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge1%yes_sp==1) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge1 Sp'
    i=i+1
    Face%VertexList(i)=GetNrP(Face%Edge1)
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%VertexList(i)
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%VertexList(i)
  ELSE IF (Face%Edge1%Vert1%in_out==0) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge1 NoSp'
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge1%Vert1%nrP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge1%Vert1%nrP
  END IF
  !............................................
  IF (Face%Edge2%Vert1%in_out>=0) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge2 In'
    i=i+1
    Face%VertexList(i)=Face%Edge2%Vert1%nrP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge2%yes_sp==1) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge2 Sp'
    i=i+1
    Face%VertexList(i)=GetNrP(Face%Edge2)
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%VertexList(i)
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%VertexList(i)
  ELSE IF (Face%Edge2%Vert1%in_out==0) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge1 NoSp'
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge2%Vert1%nrP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge2%Vert1%nrP
  END IF
  !..............................................
  IF (Face%Edge3%Vert2%in_out>=0) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge3 In'
    i=i+1
    Face%VertexList(i)=Face%Edge3%Vert2%nrP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge3%yes_sp==1) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge3 Sp'
    i=i+1
    Face%VertexList(i)=GetNrP(Face%Edge3)
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%VertexList(i)
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%VertexList(i)
  ELSE IF (Face%Edge3%Vert2%in_out==0) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge2 NoSp'
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge3%Vert2%nrP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge3%Vert2%nrP
  END IF
  !................................................
  IF (Face%Edge4%Vert2%in_out>=0) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge4 In'
    i=i+1
    Face%VertexList(i)=Face%Edge4%Vert2%nrP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge4%yes_sp==1) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge3 Sp'
    i=i+1
    Face%VertexList(i)=GetNrP(Face%Edge4)
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%VertexList(i)
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%VertexList(i)
  ELSE IF (Face%Edge4%Vert2%in_out==0) THEN
    IF(PRESENT(Print)) WRITE(*,*) 'Edge4 NoSp'
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge4%Vert2%nrP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge4%Vert2%nrP
  END IF
  IF (PRESENT(Print)) WRITE(*,*) 'ecut ',ecut
  IF (PRESENT(Print)) WRITE(*,*) 'Face%NumberVert ',Face%NumberVert
  IF(ecut>2) THEN
     j=0
     DO i=1,ecut
       IF(EdgeList(1,i)==1) THEN
         j=j+1
       END IF
     END DO 
     IF(j==2) THEN
       v=0
       Face%ec=1
       DO i=1,ecut
         IF(EdgeList(1,i)==1) THEN
           v=v+1
           Face%EdgeCut(v)=EdgeList(2,i)
         END IF
       END DO
     END IF
  END IF
  IF (Face%Edge1%Vert1%in_out<0 .AND. Face%Edge3%Vert2%in_out<0 .AND. &
     Face%Edge1%Vert2%in_out>0 .AND. Face%Edge3%Vert1%in_out>0) THEN
     Face%ec=-1
     Face%EdgeCut(:)=0
  ELSE IF (Face%Edge1%Vert1%in_out>0 .AND. Face%Edge3%Vert2%in_out>0 .AND. &
     Face%Edge1%Vert2%in_out<0 .AND. Face%Edge3%Vert1%in_out<0) THEN
     Face%ec=-1
     Face%EdgeCut(:)=0
  END IF
  IF (ecut==1) THEN
    Face%ec=-1
  END IF

  !................................................
  IF (Face%Edge1%Vert1%in_out==0.AND.Face%Edge1%Vert2%in_out==0.AND. &
      Face%Edge3%Vert1%in_out==0.AND.Face%Edge3%Vert2%in_out==0) THEN
        Face%Vol=0.0d0
        EdgP0=Face%Edge1%Vert1%Point !Ossi
        EdgP1=Face%Edge1%Vert2%Point
        EdgP2=Face%Edge3%Vert2%Point
        EdgP3=Face%Edge3%Vert1%Point
        P0=PointParametricEarth(Face%Edge1%Vert1%Point)
        P1=PointParametricEarth(Face%Edge1%Vert2%Point)
        P2=PointParametricEarth(Face%Edge3%Vert2%Point)
        P3=PointParametricEarth(Face%Edge3%Vert1%Point)
        Face%MidPoint=FaceMidP4(P0,P1,P2,P3)
        Face%MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3) ! Ossi
        NrMP=NrMP+1
        Face%mp=NrMP
  ELSE IF (Face%in_out==-4.OR. &
      MAX(Face%Edge1%Vert1%in_out,Face%Edge1%Vert2%in_out, &
          Face%Edge3%Vert1%in_out,Face%Edge3%Vert2%in_out)<=0 &
     ) THEN
    Face%Vol=0.0d0
  ELSE IF (Face%in_out==2.AND. &
           MIN(Face%Edge1%Vert1%in_out,Face%Edge1%Vert2%in_out, &
               Face%Edge3%Vert1%in_out,Face%Edge3%Vert2%in_out)==0 &
          ) THEN
    Face%ec=0   !! 1Edge Grenze dann Face=maxVol
  ELSE IF (Face%in_out==1.AND. &
           MIN(Face%Edge1%Vert1%in_out,Face%Edge1%Vert2%in_out, &
               Face%Edge3%Vert1%in_out,Face%Edge3%Vert2%in_out)==0 &
          ) THEN
    Face%ec=0   !! 2Edge Grenze dann Face=maxVol
  ELSE IF (Face%ec>0) THEN
    PMin=Face%Edge2%Vert2%Point
    PMax=Face%Edge1%Vert1%Point
    CALL SearchMinMaxFace(Face,PMin,PMax)
    Face%Vol=0.0d0
    nxc=IncrVol
    nyc=IncrVol
    nzc=IncrVol
    xc0=Face%Edge1%Vert1%Point%x
    yc0=Face%Edge1%Vert1%Point%y
    zc0=Face%Edge1%Vert1%Point%z
    xc1=Face%Edge2%Vert2%Point%x
    yc1=Face%Edge2%Vert2%Point%y
    zc1=Face%Edge2%Vert2%Point%z
    !        3-----2   ! generell Point-Folge Face-Vol-Berechnung
    !        | \   |
    !        |   \ |
    !        0-----1
    !---Face_YZ------------------------------------------------------------
    IF (xc0==xc1) THEN
      IF ((Face%Edge4%in_out==2).AND.(Face%Edge2%in_out==-2) .OR. &
          (Face%Edge4%in_out==-2).AND.(Face%Edge2%in_out==2))THEN
         !.........y-direction.............................................
         !     0 ---------3   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !     |          |   
         !     |          |
         !     1----------2
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc0
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc1
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%z-PMin%z
         zvc0=PMin%z
         zvc1=PMax%z
         dz=(zvc1-zvc0)/nzc
         DO kc=1,nzc
           Vertex1%Point%z=zvc0
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeYaxis(EdgP0,EdgP3,Vertex1,Vertex2,yc0,yc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3)  !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3)  !nur je Teilflaeche
           SumMidP=SumMidP+Vol*MidPoint
           zvc0=MIN(zvc0+dz,zvc1)
         END DO
      ELSE
         !...........z-direction................................................
         !        3-----2   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !        | \   |
         !        |   \ |
         !        0-----1
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc0
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc1
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%y-PMin%y
         yvc0=PMin%y
         yvc1=PMax%y
         dy=(yvc1-yvc0)/nyc
         DO jc=1,nyc
           Vertex1%Point%y=yvc0
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP0,EdgP3,Vertex1,Vertex2,zc0,zc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3)  !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3)  !nur je Teilflaeche
           SumMidP=SumMidP+Vol*MidPoint
           yvc0=MIN(yvc0+dy,yvc1)
         END DO
      END IF
      Face%MidPoint=SumMidP/Face%Vol
      NrMP=NrMP+1
      Face%mp=NrMP
    !---Face_XZ------------------------------------------------------------
    ELSE IF (yc0==yc1) THEN
      IF ((Face%Edge1%in_out==2).AND.(Face%Edge3%in_out==-2) .OR. &
          (Face%Edge1%in_out==-2).AND.(Face%Edge3%in_out==2))THEN
         !........x-direction.......................................... 
         !     0----------3   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !     |          |
         !     |          |
         !     1----------2
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc1
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc0
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%x-PMin%x
         zvc0=PMin%z
         zvc1=PMax%z
         dz=(zvc1-zvc0)/nzc
         DO kc=1,nzc
           Vertex1%Point%z=zvc0
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeXaxis(EdgP0,EdgP3,Vertex1,Vertex2,xc0,xc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           zvc0=zvc0+dz
         END DO
      ELSE
         !..........z-direction..............................................
         !        3-----2   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !        | \   |
         !        |   \ |
         !        0-----1
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc0
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc1
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%x-PMin%x
         xvc0=PMin%x
         xvc1=PMax%x
         dx=(xvc1-xvc0)/nxc
         DO ic=1,nxc
           Vertex1%Point%x=xvc0
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeZaxis(EdgP0,EdgP3,Vertex1,Vertex2,zc0,zc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           xvc0=xvc0+dx
         END DO
      END IF
      Face%MidPoint=SumMidP/Face%Vol
      NrMP=NrMP+1
      Face%mp=NrMP
    !---Face_XY------------------------------------------------------------
    ELSE IF (zc0==zc1) THEN
      IF ((Face%Edge4%in_out==2).AND.(Face%Edge2%in_out==-2) .OR. &
          (Face%Edge4%in_out==-2).AND.(Face%Edge2%in_out==2))THEN
         !....x-direction..................................................
         !     0----------3   !Point-Folge Cut-Plane Face-Vol-Berechnung
         ! dy  |          |
         !     |          |
         !     1----------2
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc1
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc0
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%y-PMin%y
         yvc0=PMin%y
         yvc1=PMax%y
         dy=(yvc1-yvc0)/nyc
         DO jc=1,nyc
           Vertex1%Point%y=yvc0
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP0,EdgP3,Vertex1,Vertex2,xc0,xc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           yvc0=yvc0+dy
         END DO
      ELSE
         !....y-direction...................................................
         !        3-----2   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !        | \   |
         !        |   \ |
         !        0-----1
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc1
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc0
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%x-PMin%x
         xvc0=PMin%x
         xvc1=PMax%x
         dx=(xvc1-xvc0)/nxc
         DO ic=1,nxc
           Vertex1%Point%x=xvc0
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP0,EdgP3,Vertex1,Vertex2,yc0,yc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           xvc0=xvc0+dx
         END DO
      END IF
      Face%MidPoint=SumMidP/Face%Vol
      NrMP=NrMP+1
      Face%mp=NrMP
    END IF
  END IF

  IF ((Face%in_out<4.AND.Face%ec>0).OR.(Face%in_out<0.AND.Face%ec==-1).OR.Face%Vol==0.0d0) THEN
    NrW=NrW+1
  END IF
END SUBROUTINE AnalyzeFace

FUNCTION GetNrP(Edge)
  INTEGER :: GetNrP
  TYPE(Edge_T) :: Edge
  CHARACTER*2 :: FaceType

  SELECT CASE (Edge%EdgeType)
    CASE ('X')
      GetNrP=Edge%Vert1%VertX%NrP
    CASE ('Y')
      GetNrP=Edge%Vert1%VertY%NrP
    CASE ('Z')
      GetNrP=Edge%Vert1%VertZ%NrP
  END SELECT    
        
END FUNCTION GetNrP

FUNCTION GetNrInP(Edge)
  INTEGER :: GetNrInP
  TYPE(Edge_T) :: Edge
  CHARACTER*2 :: FaceType

  SELECT CASE (Edge%EdgeType)
    CASE ('X')
      GetNrInP=Edge%Vert1%VertX%NrInP
    CASE ('Y')
      GetNrInP=Edge%Vert1%VertY%NrInP
    CASE ('Z')
      GetNrInP=Edge%Vert1%VertZ%NrInP
  END SELECT    
        
END FUNCTION GetNrInP

SUBROUTINE InitAllCells
  INTEGER :: ib,in

  DO ib=1,nb
    CALL Set(Floor(ib))
    CALL InitCellsBlock(ix0,ix1,iy0,iy1,iz0,iz1          &
                       ,Faces_XY,Faces_YZ,Faces_ZX       &
                       ,Cells)
    DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      IF (nType(1:1)=='i') THEN
        CALL InitCellsBlock(jx0,jx1,jy0,jy1,jz0,jz1          &
                         ,Faces_XYNachbar,Faces_YZNachbar,Faces_ZXNachbar       &
                         ,CellsNachbar)
      END IF  
    END DO  
  END DO  
END SUBROUTINE InitAllCells

SUBROUTINE InitCellsBlock(ix0,ix1,iy0,iy1,iz0,iz1          &
                         ,Faces_XY,Faces_YZ,Faces_ZX       &
                         ,Cells)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1 
  TYPE(FaceP_T) :: Faces_XY(ix0+1:ix1,iy0+1:iy1,iz0:iz1)
  TYPE(FaceP_T) :: Faces_YZ(ix0:ix1,iy0+1:iy1,iz0+1:iz1)
  TYPE(FaceP_T) :: Faces_ZX(ix0+1:ix1,iy0:iy1,iz0+1:iz1)
  TYPE(CellP_T) :: Cells(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)

  INTEGER :: ib,ix,iy,iz

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        IF (ASSOCIATED(Faces_XY(ix,iy,iz-1)%Face) .AND. ASSOCIATED(Cells(ix,iy,iz)%Cell) ) THEN
          Cells(ix,iy,iz)%Cell%Face1=>Faces_XY(ix,iy,iz-1)%Face
        END IF
        IF (ASSOCIATED(Faces_XY(ix,iy,iz)%Face) .AND. ASSOCIATED(Cells(ix,iy,iz)%Cell) ) THEN
          Cells(ix,iy,iz)%Cell%Face2=>Faces_XY(ix,iy,iz)%Face
        END IF
        IF (ASSOCIATED(Faces_ZX(ix,iy-1,iz)%Face) .AND. ASSOCIATED(Cells(ix,iy,iz)%Cell) ) THEN
          Cells(ix,iy,iz)%Cell%Face3=>Faces_ZX(ix,iy-1,iz)%Face
        END IF
        IF (ASSOCIATED(Faces_ZX(ix,iy,iz)%Face) .AND. ASSOCIATED(Cells(ix,iy,iz)%Cell) ) THEN
          Cells(ix,iy,iz)%Cell%Face4=>Faces_ZX(ix,iy,iz)%Face
        END IF
        IF (ASSOCIATED(Faces_YZ(ix-1,iy,iz)%Face) .AND. ASSOCIATED(Cells(ix,iy,iz)%Cell) ) THEN
          Cells(ix,iy,iz)%Cell%Face5=>Faces_YZ(ix-1,iy,iz)%Face
        END IF
        IF (ASSOCIATED(Faces_YZ(ix,iy,iz)%Face) .AND. ASSOCIATED(Cells(ix,iy,iz)%Cell) ) THEN
          Cells(ix,iy,iz)%Cell%Face6=>Faces_YZ(ix,iy,iz)%Face
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE InitCellsBlock

SUBROUTINE AnalyzeAllCells

  INTEGER :: ib,in

  INTEGER :: liv1,liv2,gl
  INTEGER :: v_x,v_y,v_z,pix
  INTEGER :: in_out,in_out_view
  INTEGER :: NrW_Cells1U

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz

  nr_cells=0
  nr_cutplanecells=0
  nr_cutplanecells1=0
  nr_cutplanecells2=0
  nr_cutcells=0
  nr_cutf1=0
  nr_soilplanecells=0
  nr_soilplanecells1=0
  nr_soilplanecells2=0
  nr_incells=0
  nr_viewcells=0
  nr_viewcutplanecells=0
  nr_viewcutplanecells1=0
  nr_viewcutplanecells2=0
  nr_viewcutcells=0
  nr_viewcutf1=0
  nr_viewsoilplanecells=0
  nr_viewsoilplanecells1=0
  nr_viewsoilplanecells2=0
  nr_viewincells=0
  NrAll_Cells=0
  NrW_All_Cells=0
  NrB_All_Cells=0
  NrMP_All_Cells=0

  !.....................................................................
  DO ib=1,nb
    CALL Set(Floor(ib))
    CALL AnalyzeCellsBlock(ix0,ix1,iy0,iy1,iz0,iz1,Cells)
    DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      IF (nType(1:1)=='i') THEN
        CALL AnalyzeCellsBlock(jx0,jx1,jy0,jy1,jz0,jz1,CellsNachbar)
      END IF  
    END DO  
    Floor(ib)%NrW_Cells=NrW_Cells
    Floor(ib)%NrMP_Cells=NrB_Cells
    Floor(ib)%NrB_Cells=NrB_Cells
    NrW_All_Cells=NrW_All_Cells+NrW_Cells
    NrMP_All_Cells=NrMP_All_Cells+NrB_Cells
    NrB_All_Cells=NrB_All_Cells+NrB_Cells
  END DO
! Overwrite cell volumes
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      IF (nType=='ie'.OR.nType=='iw') THEN
        IF (nType=='ie') THEN
          ix=Floor(ibn)%ix0+1
        ELSE
          ix=Floor(ibn)%ix1
        END IF  
        IF (MIN(RefineNachbarY,RefineNachbarZ)>MIN(RefineY,RefineZ)) THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              IF (ASSOCIATED(CellsNachbar(jx0+1,jy,jz)%Cell)) THEN
                CellsNachbar(jx0+1,jy,jz)%Cell%Vol=0.0d0
                DO iy=IncrY*(jy-1)+1,IncrY*jy
                  DO iz=IncrY*(jz-1)+1,IncrZ*jz
                    IF (ASSOCIATED(Floor(ibn)%Cells(ix,iy,iz)%Cell)) THEN
                      CellsNachbar(jx0+1,jy,jz)%Cell%Vol=CellsNachbar(jx0+1,jy,jz)%Cell%Vol &
                                                        +Floor(ibn)%Cells(ix,iy,iz)%Cell%Vol
                    ELSE  
                      CellsNachbar(jx0+1,jy,jz)%Cell%Vol=CellsNachbar(jx0+1,jy,jz)%Cell%Vol &
                                                        +Floor(ibn)%dx(ix)*Floor(ibn)%dy(iy)*Floor(ibn)%dz(iz)
                    END IF
                  END DO  
                END DO  
              END IF
            END DO  
          END DO  
        ELSE IF (MIN(RefineNachbarY,RefineNachbarZ)<MIN(RefineY,RefineZ)) THEN
          DO iy=Floor(ibn)%iy0+1,Floor(ibn)%iy1
            DO iz=Floor(ibn)%iz0+1,Floor(ibn)%iz1
              IF (ASSOCIATED(Floor(ibn)%Cells(ix,iy,iz)%Cell)) THEN
                Floor(ibn)%Cells(ix,iy,iz)%Cell%Vol=0.0d0
                DO jy=IncrY*(iy-1)+1,IncrY*iy
                  DO jz=IncrY*(iz-1)+1,IncrZ*iz
                    IF (ASSOCIATED(CellsNachbar(jx0+1,jy,jz)%Cell)) THEN
                      Floor(ibn)%Cells(ix,iy,iz)%Cell%Vol=Floor(ibn)%Cells(ix,iy,iz)%Cell%Vol &
                                                         +CellsNachbar(jx0+1,jy,jz)%Cell%Vol
                    ELSE
                      Floor(ibn)%Cells(ix,iy,iz)%Cell%Vol=Floor(ibn)%Cells(ix,iy,iz)%Cell%Vol &
                                                         +Floor(ibn)%dx(ix)*dy(jy)*dz(jz)
                    END IF
                  END DO  
                END DO  
              END IF  
            END DO  
          END DO  
        END IF  
      END IF  
    END DO  
  END DO

END SUBROUTINE AnalyzeAllCells

SUBROUTINE AnalyzeCellsBlock(ix0,ix1,iy0,iy1,iz0,iz1,Cells)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1 
  TYPE(CellP_T) :: Cells(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)

  INTEGER :: ix,iy,iz

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        IF (ASSOCIATED(Cells(ix,iy,iz)%Cell)) THEN
          CALL AnalyzeCell(Cells(ix,iy,iz)%Cell)
          CALL SetLandClass(Cells(ix,iy,iz)%Cell,ix,iy)
        END IF         
      END DO
    END DO
  END DO
END SUBROUTINE AnalyzeCellsBlock

SUBROUTINE SetLandClass(Cell,ix,iy)
  TYPE(Cell_T) :: Cell
  INTEGER :: ix,iy

  INTEGER :: igxLoc,igyLoc
  INTEGER :: nl
! IF(ifcorine) THEN
!  igxLoc=ix*2.e0**RefineX
!  igyLoc=iy*2.e0**RefineY
!  Cell%LandClass=v_clc(LandCover(igxLoc,igyLoc))%mat_c
! ELSE
    IF(ASSOCIATED(LandDef)) THEN
      DO nl=1,nr_landdef
        IF (ix>LandDef(1,1,nl)*2.e0**RefineX .AND. &
            ix<=LandDef(2,1,nl)*2.e0**RefineX  .AND. &
            iy>LandDef(1,2,nl)*2.e0**RefineY .AND. &
            iy<=LandDef(2,2,nl)*2.e0**RefineY ) THEN
          Cell%LandClass=LandDef(1,3,nl)
          EXIT
        END IF
      END DO
    END IF
! END IF
END SUBROUTINE SetLandClass

SUBROUTINE AnalyzeCell(Cell,Print)
  TYPE(Cell_T) :: Cell
  LOGICAL, OPTIONAL :: Print
  TYPE(Cell_T), POINTER :: pCell=>NULL()
  INTEGER :: ListCut(1:2,1:8)
  INTEGER :: nCut,nCutNeu,i,j,iv
  INTEGER :: ic,jc,kc,nxc,nyc,nzc
  INTEGER :: iTemp
  REAL(8) :: xc0,xc1,yc0,yc1,zc0,zc1,dx,dy,dz
  TYPE (Vertex_T) :: Vertex1,Vertex2
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(Point_T) :: EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7
  TYPE(Point_T) :: PMin,PMax,SumMidP,MidPoint
  TYPE(Point_T) :: PCut(1:8)
  REAL(8) :: xvc0,xvc1,xvu0,xvu1,xvo0,xvo1
  REAL(8) :: yvc0,yvc1,yvu0,yvu1,yvo0,yvo1
  REAL(8) :: zvc0,zvc1
  REAL(8) :: Vol,Len
  INTEGER :: VertexFace(1:6,1:9)
  LOGICAL :: Cut
  INTEGER :: iP1,iP2,jP1,jP2,iVert,jVert

  EdgP0=Cell%Face1%Edge1%Vert1%Point
  EdgP1=Cell%Face1%Edge1%Vert2%Point
  EdgP2=Cell%Face1%Edge3%Vert1%Point
  EdgP3=Cell%Face1%Edge3%Vert2%Point
  EdgP4=Cell%Face2%Edge1%Vert1%Point
  EdgP5=Cell%Face2%Edge1%Vert2%Point
  EdgP6=Cell%Face2%Edge3%Vert1%Point
  EdgP7=Cell%Face2%Edge3%Vert2%Point
  P0=PointParametricEarth(EdgP0)
  P1=PointParametricEarth(EdgP1)
  P2=PointParametricEarth(EdgP2)
  P3=PointParametricEarth(EdgP3)
  P4=PointParametricEarth(EdgP4)
  P5=PointParametricEarth(EdgP5)
  P6=PointParametricEarth(EdgP6)
  P7=PointParametricEarth(EdgP7)
  Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
  Cell%MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
  !----------------------------------------------------------------------
  VertexFace=-1
  nCut=0
  ListCut=0
  Cell%in_out=Cell%Face1%in_out+Cell%Face2%in_out
  
  IF (Cell%Face1%NumberVert>2) THEN
    VertexFace(1,1:Cell%Face1%NumberVert)=Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
    VertexFace(1,Cell%Face1%NumberVert+1)=VertexFace(1,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 1'
        WRITE(*,*) VertexFace(1,1:Cell%Face6%NumberVert+1)
      END IF  
  END IF
  IF (Cell%Face2%NumberVert>2) THEN
    VertexFace(2,1:Cell%Face2%NumberVert)=Cell%Face2%VertexList(1:Cell%Face2%NumberVert)
    VertexFace(2,Cell%Face2%NumberVert+1)=VertexFace(2,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 2'
        WRITE(*,*) VertexFace(2,1:Cell%Face6%NumberVert+1)
      END IF  
  END IF
  IF (Cell%Face3%NumberVert>2) THEN
    VertexFace(3,1:Cell%Face3%NumberVert)=Cell%Face3%VertexList(1:Cell%Face3%NumberVert)
    VertexFace(3,Cell%Face3%NumberVert+1)=VertexFace(3,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 3'
        WRITE(*,*) VertexFace(3,1:Cell%Face6%NumberVert+1)
      END IF  
  END IF
  IF (Cell%Face4%NumberVert>2) THEN
    VertexFace(4,1:Cell%Face4%NumberVert)=Cell%Face4%VertexList(1:Cell%Face4%NumberVert)
    VertexFace(4,Cell%Face4%NumberVert+1)=VertexFace(4,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 4'
        WRITE(*,*) VertexFace(4,1:Cell%Face6%NumberVert+1)
      END IF  
  END IF
  IF (Cell%Face5%NumberVert>2) THEN
    VertexFace(5,1:Cell%Face5%NumberVert)=Cell%Face5%VertexList(1:Cell%Face5%NumberVert)
    VertexFace(5,Cell%Face5%NumberVert+1)=VertexFace(5,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 5'
        WRITE(*,*) VertexFace(5,1:Cell%Face6%NumberVert+1)
      END IF  
  END IF
  IF (Cell%Face6%NumberVert>2) THEN
    VertexFace(6,1:Cell%Face6%NumberVert)=Cell%Face6%VertexList(1:Cell%Face6%NumberVert)
    VertexFace(6,Cell%Face6%NumberVert+1)=VertexFace(6,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 6'
        WRITE(*,*) VertexFace(6,1:Cell%Face6%NumberVert+1)
      END IF  
  END IF

  C1:DO i=1,6
    C2:DO iVert=1,8
      iP1=VertexFace(i,iVert)
      iP2=VertexFace(i,iVert+1)
      IF (ip2==-1) THEN
        EXIT C2
      END IF
      Cut=.TRUE.
      C3:DO j=1,6
        IF (j/=i) THEN
          C4:DO jVert=1,8
            jP1=VertexFace(j,jVert)
            jP2=VertexFace(j,jVert+1)
            IF (jp2==-1) THEN
              EXIT C4
            END IF
            IF ((iP1==jP1.AND.iP2==jP2).OR.  &
                (iP1==jP2.AND.iP2==jP1)) THEN
              Cut=.FALSE.
              EXIT C3
            END IF
          END DO C4
        END IF
      END DO C3
      IF (Cut) THEN
        ncut=ncut+1
        ListCut(1,nCut)=iP1
        ListCut(2,nCut)=iP2
      END IF
    END DO C2
  END DO C1
  IF (nCut>0) THEN
    nCutNeu=nCut
    DO i=1,ncut-1
      DO j=i+1,ncut
        IF (( (ListCut(1,i)==ListCut(1,j).AND.ListCut(2,i)==ListCut(2,j)).OR. &
              (ListCut(1,i)==ListCut(2,j).AND.ListCut(2,i)==ListCut(1,j)) ) &
           .AND.ListCut(1,i)>0) THEN
          ListCut(:,j)=0
          nCutNeu=nCutNeu-1
        END IF
      END DO
    END DO
    Cell%vc=nCutNeu
    ALLOCATE(Cell%VertCut(nCut))
    ALLOCATE(Cell%VertCutOnlyP(nCut))
    Cell%VertCut(1)= ListCut(1,1)
    Cell%VertCut(2)= ListCut(2,1)
    ListCut(1,1)=0
    ListCut(2,1)=0
    iv=2
    S2:DO
      S1:DO i=1,nCut
        IF (Cell%VertCut(iv)==ListCut(1,i)) THEN
          iv=iv+1
          Cell%VertCut(iv)=ListCut(2,i)
          ListCut(:,i)=0
          EXIT S1
        ELSE IF (Cell%VertCut(iv)==ListCut(2,i)) THEN
          iv=iv+1
          Cell%VertCut(iv)=ListCut(1,i)
          ListCut(:,i)=0
          EXIT S1
        END IF
      END DO S1
      IF (Cell%VertCut(1)==Cell%VertCut(iv).AND.iv<nCutNeu) THEN
        WRITE(*,*) 'More than one Loop'
        DO i=1,nCut
          WRITE(*,*) ListCut(:,i)
        END DO  
        WRITE(*,*) 'Face1',VertexFace(1,1:Cell%Face1%NumberVert)
        WRITE(*,*) 'Face2',VertexFace(2,1:Cell%Face2%NumberVert)
        WRITE(*,*) 'Face3',VertexFace(3,1:Cell%Face3%NumberVert)
        WRITE(*,*) 'Face4',VertexFace(4,1:Cell%Face4%NumberVert)
        WRITE(*,*) 'Face5',VertexFace(5,1:Cell%Face5%NumberVert)
        WRITE(*,*) 'Face6',VertexFace(6,1:Cell%Face6%NumberVert)
  WRITE(*,*) 'EdgP0',Cell%Face1%Edge1%Vert1%Point,Cell%Face1%Edge1%Vert1%Number
  WRITE(*,*) 'EdgP1',Cell%Face1%Edge1%Vert2%Point,Cell%Face1%Edge1%Vert2%Number
  WRITE(*,*) 'EdgP2',Cell%Face1%Edge3%Vert1%Point,Cell%Face1%Edge3%Vert1%Number
  WRITE(*,*) 'EdgP3',Cell%Face1%Edge3%Vert2%Point,Cell%Face1%Edge3%Vert2%Number
  WRITE(*,*) 'EdgP4',Cell%Face2%Edge1%Vert1%Point,Cell%Face2%Edge1%Vert1%Number
  WRITE(*,*) 'EdgP5',Cell%Face2%Edge1%Vert2%Point,Cell%Face2%Edge1%Vert2%Number
  WRITE(*,*) 'EdgP6',Cell%Face2%Edge3%Vert1%Point,Cell%Face2%Edge3%Vert1%Number
  WRITE(*,*) 'EdgP7',Cell%Face2%Edge3%Vert2%Point,Cell%Face2%Edge3%Vert2%Number
        STOP 'More than one Loop'
      END IF  
      IF (iv>=nCutNeu) THEN
        EXIT S2
      END IF
    END DO S2
  ELSE
    Cell%vc=0
    IF (Cell%in_out==-8) THEN
      Cell%Vol=0.0d0
    END IF
  END IF

  ! Volumen-Berechnung
  IF (Cell%vc>0) THEN
    Cell%Vol=0.0d0
    xc0=Cell%Face1%Edge1%Vert1%Point%x
    yc0=Cell%Face1%Edge1%Vert1%Point%y
    zc0=Cell%Face1%Edge1%Vert1%Point%z
    xc1=Cell%Face2%Edge2%Vert2%Point%x
    yc1=Cell%Face2%Edge2%Vert2%Point%y
    zc1=Cell%Face2%Edge2%Vert2%Point%z
    PMin=Cell%Face2%Edge3%Vert2%Point
    PMax=Cell%Face1%Edge1%Vert1%Point
    CALL SearchMinMaxCell(Cell,PMin,PMax)

    IF((Cell%face5%in_out>=0.AND.Cell%face6%in_out<=0).OR. &
       (Cell%face5%in_out<=0.AND.Cell%face6%in_out>=0)) THEN
       !...Plane x-direction.................................
       IF (PRESENT(Print)) WRITE(*,*) 'Volume Berechnung in x-Dir'
       zvc0=PMin%z
       zvc1=PMax%z
       yvc0=PMin%y
       yvc1=PMax%y
       yc0=yvc0
       nxc=IncrVol
       nyc=IncrVol
       nzc=IncrVol
       dy=(yvc1-yvc0)/nyc
       dz=(zvc1-zvc0)/nzc
       Vertex1%Point%x=xc0
       Vertex1%in_out=-1
       Vertex2%Point%x=xc1
       Vertex2%in_out=1
       SumMidP%x=0.0d0
       SumMidP%y=0.0d0
       SumMidP%z=0.0d0
       DO kc=1,nzc
         yvc0=yc0
         DO jc=1,nyc
           !      2--------------6    !Point-Folge Vol-Berechnung
           !     /|             /|
           !    0--------------4 |
           !    | 3------------|-7
           !    |/             |/
           !    1--------------5
           Vertex1%Point%z=zvc0
           Vertex1%Point%y=yvc0
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP1,EdgP5,Vertex1,Vertex2,xc0,xc1)
           P1=PointParametricEarth(EdgP1)
           P5=PointParametricEarth(EdgP5)
   
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex1%Point%y=yvc0
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP0,EdgP4,Vertex1,Vertex2,xc0,xc1)
           P0=PointParametricEarth(EdgP0)
           P4=PointParametricEarth(EdgP4)
   
           Vertex1%Point%z=zvc0
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP3,EdgP7,Vertex1,Vertex2,xc0,xc1)
           P3=PointParametricEarth(EdgP3)
           P7=PointParametricEarth(EdgP7)
   
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP2,EdgP6,Vertex1,Vertex2,xc0,xc1)
           P2=PointParametricEarth(EdgP2)
           P6=PointParametricEarth(EdgP6)
   
           Vol=ABS(VolP8(P0,P1,P2,P3,P4,P5,P6,P7))
           Cell%Vol=Cell%Vol+Vol
           MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
           MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
           SumMidP=Vol*MidPoint+SumMidP
         yvc0=yvc0+dy
         END DO
       zvc0=zvc0+dz
       END DO
    ELSE IF((Cell%face3%in_out>=0.AND.Cell%face4%in_out<=0).OR. &
            (Cell%face3%in_out<=0.AND.Cell%face4%in_out>=0)) THEN
       !...Plane y-direction..................................
       IF (PRESENT(Print)) WRITE(*,*) 'Volume Berechnung in y-Dir'
       zvc0=PMin%z
       zvc1=PMax%z
       xvc0=PMin%x
       xvc1=PMax%x
       xc0=xvc0
       nxc=IncrVol
       nyc=IncrVol
       nzc=IncrVol
       dx=(xvc1-xvc0)/nxc
       dz=(zvc1-zvc0)/nzc
       Vertex1%Point%y=yc0
       Vertex1%in_out=-1
       Vertex2%Point%y=yc1
       Vertex2%in_out=1
       SumMidP%x=0.0d0
       SumMidP%y=0.0d0
       SumMidP%z=0.0d0
       DO kc=1,nzc
         xvc0=xc0
         DO ic=1,nxc
           !          6----7    !Point-Folge Vol-Berechnung
           !         /|   /|
           !        / |  / |
           !       /  2 /  3 
           !      4----5  /
           !      | /  | / 
           !      |/   |/ 
           !      0----1
           Vertex1%Point%z=zvc0
           Vertex1%Point%x=xvc0
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP0,EdgP2,Vertex1,Vertex2,yc0,yc1)
           P0=PointParametricEarth(EdgP0)
           P2=PointParametricEarth(EdgP2)
   
           Vertex1%Point%z=zvc0
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP1,EdgP3,Vertex1,Vertex2,yc0,yc1)
           P1=PointParametricEarth(EdgP1)
           P3=PointParametricEarth(EdgP3)
   
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex1%Point%x=xvc0
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP4,EdgP6,Vertex1,Vertex2,yc0,yc1)
           P4=PointParametricEarth(EdgP4)
           P6=PointParametricEarth(EdgP6)
   
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP5,EdgP7,Vertex1,Vertex2,yc0,yc1)
           P5=PointParametricEarth(EdgP5)
           P7=PointParametricEarth(EdgP7)
   
           Vol=ABS(VolP8(P0,P1,P2,P3,P4,P5,P6,P7))
           Cell%Vol=Cell%Vol+Vol
           MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
           MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
           SumMidP=Vol*MidPoint+SumMidP
         xvc0=xvc0+dx
         END DO
       zvc0=zvc0+dz
       END DO
    ELSE 
       !...Plane z-direction..................................
       IF (PRESENT(Print)) WRITE(*,*) 'Volume Berechnung in z-Dir'
       xvc0=PMin%x
       xvc1=PMax%x
       yvc0=PMin%y
       yvc1=PMax%y
       yc0=yvc0
       nxc=IncrVol
       nyc=IncrVol
       nzc=IncrVol
       dx=(xvc1-xvc0)/nxc
       dy=(yvc1-yvc0)/nyc
       Vertex1%Point%z=zc0
       Vertex1%in_out=-1
       Vertex2%Point%z=zc1
       Vertex2%in_out=1
       SumMidP%x=0.0d0
       SumMidP%y=0.0d0
       SumMidP%z=0.0d0
       DO ic=1,nxc
         yvc0=yc0
         DO jc=1,nyc
           !      6--------7    !Point-Folge Vol-Berechnung
           !     /|       /|
           !    4--------5 |
           !    | |      | |
           !    | 2------|-3
           !    |/       |/
           !    0--------1
           Vertex1%Point%x=xvc0
           Vertex1%Point%y=yvc0
           Vertex2%Point%x=Vertex1%Point%x
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP0,EdgP4,Vertex1,Vertex2,zc0,zc1)
           IF (PRESENT(Print)) WRITE(*,*) 'E1',EdgP0,EdgP4
           P0=PointParametricEarth(EdgP0)
           P4=PointParametricEarth(EdgP4)
   
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex1%Point%y=yvc0
           Vertex2%Point%x=Vertex1%Point%x
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP1,EdgP5,Vertex1,Vertex2,zc0,zc1)
           IF (PRESENT(Print)) WRITE(*,*) 'E2',EdgP1,EdgP5
           P1=PointParametricEarth(EdgP1)
           P5=PointParametricEarth(EdgP5)
   
           Vertex1%Point%x=xvc0
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%x=Vertex1%Point%x
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP2,EdgP6,Vertex1,Vertex2,zc0,zc1)
           IF (PRESENT(Print)) WRITE(*,*) 'E3',EdgP2,EdgP6
           P2=PointParametricEarth(EdgP2)
           P6=PointParametricEarth(EdgP6)
   
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%x=Vertex1%Point%x
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP3,EdgP7,Vertex1,Vertex2,zc0,zc1)
           IF (PRESENT(Print)) WRITE(*,*) 'E4',EdgP3,EdgP7
           P3=PointParametricEarth(EdgP3)
           P7=PointParametricEarth(EdgP7)
   
           Vol=ABS(VolP8(P0,P1,P2,P3,P4,P5,P6,P7))
           Cell%Vol=Cell%Vol+Vol
           IF (PRESENT(Print)) WRITE(*,*) 'Vol',Vol,'SumVol',Cell%Vol
           MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
           MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
           SumMidP=Vol*MidPoint+SumMidP
         yvc0=yvc0+dy
         END DO
       xvc0=xvc0+dx
       END DO
    END IF
    IF (PRESENT(Print)) WRITE(*,*) 'vor Cell%MidPoint 0'
    Cell%MidPoint=SumMidP/Cell%Vol 
    NrMP_Cells=NrMP_Cells+1
    Cell%mp=NrMP_Cells
    IF (ncut<=2) THEN
      IF (Cell%in_out<=0) THEN
        Cell%Vol=0.0d0
      ELSE
        Cell%in_out=8
      END IF
    END IF
    IF (PRESENT(Print)) WRITE(*,*) 'Cell%vc',Cell%vc,OutUnitProt
    IF (Cell%vc>2) THEN
      i=1
      !write(*,*) "-------"
      IF(Cell%vc>6) THEN
        IF (PRESENT(Print)) WRITE(*,*) 'Write'
        Write(OutUnitProt,*) "Cell%vc>6"
      IF (PRESENT(Print)) WRITE(*,*) 'nach Write'
        !pCell=Cell
      END IF     
      IF (PRESENT(Print)) WRITE(*,*) 'nach Write 1'
      DO i=1,Cell%vc
        IF (PRESENT(Print)) WRITE(*,*) 'Pch Cut',i,SIZE(Cell%VertCut)
        IF (PRESENT(Print)) WRITE(*,*) 'PCut',i,Cell%VertCut(i)
        PCut(i)=VertOut(Cell%VertCut(i))%Vertex%Point
       ! write(*,*) PCut(i)
      END DO
      IF (PRESENT(Print)) WRITE(*,*) 'SELECT CASE'
      !write(*,*) (PCut(i),i=1,Cell%vc)
      SELECT CASE (Cell%vc)
        CASE(3)
          Cell%CutF_MidP=FaceMidP3(PCut(1),PCut(2),PCut(3))
        CASE(4)
          Cell%CutF_MidP=FaceMidP4(PCut(1),PCut(2),PCut(3),PCut(4))
        CASE(5)
          Cell%CutF_MidP=FaceMidP5(PCut(1),PCut(2),PCut(3),PCut(4),PCut(5))
        CASE(6)
          Cell%CutF_MidP=FaceMidP6(PCut(1),PCut(2),PCut(3),PCut(4),PCut(5),PCut(6))
      END SELECT
    END IF
  ELSE  ! Cell%vc==0
    IF ((Cell%Face1%Edge1%Vert1%in_out==0 .AND. &
         Cell%Face1%Edge1%Vert2%in_out==0 .AND. &
         Cell%Face1%Edge3%Vert1%in_out==0 .AND. &
         Cell%Face1%Edge3%Vert2%in_out==0) &
        .AND. Cell%Face2%in_out==4) THEN 
       !Grenzflaeche Face1
       EdgP0=Cell%Face1%Edge1%Vert1%Point
       EdgP1=Cell%Face1%Edge1%Vert2%Point
       EdgP2=Cell%Face1%Edge3%Vert1%Point
       EdgP3=Cell%Face1%Edge3%Vert2%Point
       EdgP4=Cell%Face2%Edge1%Vert1%Point
       EdgP5=Cell%Face2%Edge1%Vert2%Point
       EdgP6=Cell%Face2%Edge3%Vert1%Point
       EdgP7=Cell%Face2%Edge3%Vert2%Point
       P0=PointParametricEarth(EdgP0)
       P1=PointParametricEarth(EdgP1)
       P2=PointParametricEarth(EdgP2)
       P3=PointParametricEarth(EdgP3)
       P4=PointParametricEarth(EdgP4)
       P5=PointParametricEarth(EdgP5)
       P6=PointParametricEarth(EdgP6)
       P7=PointParametricEarth(EdgP7)
       Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
       Cell%MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
       Cell%CutF_MidP=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3)
    ELSE
       EdgP0=Cell%Face1%Edge1%Vert1%Point
       EdgP1=Cell%Face1%Edge1%Vert2%Point
       EdgP2=Cell%Face1%Edge3%Vert1%Point
       EdgP3=Cell%Face1%Edge3%Vert2%Point
       EdgP4=Cell%Face2%Edge1%Vert1%Point
       EdgP5=Cell%Face2%Edge1%Vert2%Point
       EdgP6=Cell%Face2%Edge3%Vert1%Point
       EdgP7=Cell%Face2%Edge3%Vert2%Point
       P0=PointParametricEarth(EdgP0)
       P1=PointParametricEarth(EdgP1)
       P2=PointParametricEarth(EdgP2)
       P3=PointParametricEarth(EdgP3)
       P4=PointParametricEarth(EdgP4)
       P5=PointParametricEarth(EdgP5)
       P6=PointParametricEarth(EdgP6)
       P7=PointParametricEarth(EdgP7)
       Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
       Cell%MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
       Cell%CutF_MidP=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3)
    END IF
  END IF
 
  IF (PRESENT(Print)) WRITE(*,*) 'vor Cell%MidPoint'
  IF (Cell%Vol==0.0d0) THEN
    Cell%Face1%Vol=0.0d0
    Cell%Face2%Vol=0.0d0
    Cell%Face3%Vol=0.0d0
    Cell%Face4%Vol=0.0d0
    Cell%Face5%Vol=0.0d0
    Cell%Face6%Vol=0.0d0
    IF(Cell%mp>0)THEN
      NrMP_Cells=NrMP_Cells-1
      Cell%mp=0
    END IF
  END IF
  IF (PRESENT(Print)) WRITE(*,*) 'vor Cell%MidPoint 1'
  IF (MAX(Cell%Face1%Vol,Cell%Face2%Vol,Cell%Face3%Vol,Cell%Face4%Vol,Cell%Face5%Vol,Cell%Face6%Vol)<=0.0d0) THEN
    Cell%Vol=0.0d0  !Ossi
    Cell%mp=0  ! Celle unterhalb Berg angrenzt 
  END IF
  IF (Cell%Vol<1.d-12) THEN
    Cell%vc=0
    Cell%Vol=0.0d0
  END IF  
 
END SUBROUTINE AnalyzeCell

SUBROUTINE SearchPointsEdgeXaxis(P0,P1,Vertex1,Vertex2,xc0,xc1)
  TYPE (Vertex_T) :: Vertex1,Vertex2
  REAL(8) :: xc0,xc1
  TYPE(Point_T)   :: P0,P1
  REAL(8) :: xp,dist1,dist2
  TYPE (Vertex_T) :: InterVert

        dist1=Level(Vertex1%Point%x,Vertex1%Point%y,Vertex1%Point%z)
        dist2=Level(Vertex2%Point%x,Vertex2%Point%y,Vertex2%Point%z)
        IF (dist2>=0.0d0) THEN
          xp=xc0
          IF (dist1<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            xp=Intervert%Point%x
          END IF
          P0%x=xp
          P0%y=Vertex1%Point%y
          P0%z=Vertex1%Point%z
          P1%x=xc1
          P1%y=Vertex2%Point%y
          P1%z=Vertex2%Point%z
        ELSE IF (dist1>=0.0d0) THEN
          xp=xc1
          IF (dist2<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            xp=Intervert%Point%x
          END IF
          P0%x=xc0
          P0%y=Vertex1%Point%y
          P0%z=Vertex1%Point%z
          P1%x=xp
          P1%y=Vertex2%Point%y
          P1%z=Vertex2%Point%z
        ELSE
          IF (dist1<=dist2) THEN
            P0%x=Vertex2%Point%x
            P0%y=Vertex2%Point%y
            P0%z=Vertex2%Point%z
            P1%x=Vertex2%Point%x
            P1%y=Vertex2%Point%y
            P1%z=Vertex2%Point%z
          ELSE IF (dist2<dist1) THEN
            P0%x=Vertex1%Point%x
            P0%y=Vertex1%Point%y
            P0%z=Vertex1%Point%z
            P1%x=Vertex1%Point%x
            P1%y=Vertex1%Point%y
            P1%z=Vertex1%Point%z
          END IF
        END IF
END SUBROUTINE SearchPointsEdgeXaxis

SUBROUTINE SearchPointsEdgeYaxis(P0,P1,Vertex1,Vertex2,yc0,yc1)
  TYPE (Vertex_T) :: Vertex1,Vertex2
  REAL(8) :: yc0,yc1
  TYPE(Point_T)   :: P0,P1
  REAL(8) :: yp,dist1,dist2
  TYPE (Vertex_T) :: InterVert

        dist1=Level(Vertex1%Point%x,Vertex1%Point%y,Vertex1%Point%z)
        dist2=Level(Vertex2%Point%x,Vertex2%Point%y,Vertex2%Point%z)
        IF (dist2>=0.0d0) THEN
          yp=yc0
          IF (dist1<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            yp=Intervert%Point%y
          END IF
          P0%x=Vertex1%Point%x
          P0%y=yp
          P0%z=Vertex1%Point%z
          P1%x=Vertex2%Point%x
          P1%y=yc1
          P1%z=Vertex2%Point%z
        ELSE IF (dist1>=0.0d0) THEN
          yp=yc1
          IF (dist2<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            yp=Intervert%Point%y
          END IF
          P0%x=Vertex1%Point%x
          P0%y=yc0
          P0%z=Vertex1%Point%z
          P1%x=Vertex2%Point%x
          P1%y=yp
          P1%z=Vertex2%Point%z
        ELSE
          IF (dist1<=dist2) THEN
            P0%x=Vertex2%Point%x
            P0%y=Vertex2%Point%y
            P0%z=Vertex2%Point%z
            P1%x=Vertex2%Point%x
            P1%y=Vertex2%Point%y
            P1%z=Vertex2%Point%z
          ELSE IF (dist2<dist1) THEN
            P0%x=Vertex1%Point%x
            P0%y=Vertex1%Point%y
            P0%z=Vertex1%Point%z
            P1%x=Vertex1%Point%x
            P1%y=Vertex1%Point%y
            P1%z=Vertex1%Point%z
          END IF
        END IF
END SUBROUTINE SearchPointsEdgeYaxis

SUBROUTINE SearchPointsEdgeZaxis(P0,P1,Vertex1,Vertex2,zc0,zc1)
  TYPE (Vertex_T) :: Vertex1,Vertex2
  REAL(8) :: zc0,zc1
  TYPE(Point_T)   :: P0,P1
  REAL(8) :: zp,dist1,dist2
  TYPE (Vertex_T) :: InterVert

        dist1=Level(Vertex1%Point%x,Vertex1%Point%y,Vertex1%Point%z)
        dist2=Level(Vertex2%Point%x,Vertex2%Point%y,Vertex2%Point%z)
        IF (dist2>=0.0d0) THEN
          zp=zc0
          IF (dist1<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            zp=Intervert%Point%z
          END IF
          P0%x=Vertex1%Point%x
          P0%y=Vertex1%Point%y
          P0%z=zp
          P1%x=Vertex2%Point%x
          P1%y=Vertex2%Point%y
          P1%z=zc1
        ELSE IF (dist1>=0.0d0) THEN
          zp=zc1
          IF (dist2<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            zp=Intervert%Point%z
          END IF
          P0%x=Vertex1%Point%x
          P0%y=Vertex1%Point%y
          P0%z=zc0
          P1%x=Vertex2%Point%x
          P1%y=Vertex2%Point%y
          P1%z=zp
        ELSE
          IF (dist1<=dist2) THEN
            P0%x=Vertex2%Point%x
            P0%y=Vertex2%Point%y
            P0%z=Vertex2%Point%z
            P1%x=Vertex2%Point%x
            P1%y=Vertex2%Point%y
            P1%z=Vertex2%Point%z
          ELSE IF (dist2<dist1) THEN
            P0%x=Vertex1%Point%x
            P0%y=Vertex1%Point%y
            P0%z=Vertex1%Point%z
            P1%x=Vertex1%Point%x
            P1%y=Vertex1%Point%y
            P1%z=Vertex1%Point%z
          END IF
        END IF
END SUBROUTINE SearchPointsEdgeZaxis

SUBROUTINE Set_Domain_OutView()
INTEGER :: ix,iy,iz

TYPE OutMv_XYZ
  CHARACTER(3) :: who           ! ixa,ixe,iya,iye,iza,ize
  CHARACTER(5) :: scaling       ! upper,lower
  INTEGER      :: yes           ! 1 .or. 0
END TYPE
TYPE(OutMv_XYZ) :: zg(1:6)

  IF(out_area=='y') THEN 
  IF (TRIM(Domain%def_out_domain)==def_out_coord) THEN

  ! search between Domain%ix0,Domain%ix1
  !.....................................
  ix=Domain%ix0    ! lower ix
  DO WHILE (Domain%view_xa>Domain%xP(ix).and.ix<Domain%ix1)
     ix=ix+1
  END DO
  IF(Domain%view_xa==Domain%xP(ix)) THEN
     Domain%view_ixa=ix
  ELSE
     !Differenz-Check xa/xP(ix) und skalieren
     IF((Domain%xP(ix)-Domain%view_xa)<=(Domain%view_xa-Domain%xP(ix-1))) THEN
        Domain%view_ixa=ix
        !zg(1)=(/'ixa','upper',1/)
        !Compiler : All elements in an array constructor must have the same type and type parameters
        zg(1)%who="ixa";zg(1)%scaling="upper";zg(1)%yes=1
     ELSE
        Domain%view_ixa=ix-1
        !zg(1)=(/"ixa","lower",1/)
        zg(1)%who="ixa";zg(1)%scaling="lower";zg(1)%yes=1
     END IF
     Write(*,*) "                                     ",&
                &"View-OutGmV a must resized:"
     Write(*,*) "                                     ",&
                &"Domain%view_ixa next",zg(1)%scaling,"cell border!"
  END IF
  ix=Domain%ix1  !upper ix
  DO WHILE (Domain%view_xe<Domain%xP(ix).and.ix>Domain%ix0)
     ix=ix-1
  END DO
  IF(Domain%view_xe==Domain%xP(ix)) THEN
     Domain%view_ixe=ix
  ELSE
     !Differenz-Check xe/xP(ix) und skalieren
     IF((Domain%view_xe-Domain%xP(ix)) <= (Domain%xP(ix+1)-Domain%view_xe)) THEN
        Domain%view_ixe=ix
        !zg(2)=(/"ixe","lower",1/)
        zg(2)%who="ixe";zg(2)%scaling="lower";zg(2)%yes=1
     ELSE
        Domain%view_ixe=ix+1
        !zg(2)=(/"ixe","upper",1/)
        zg(2)%who="ixe";zg(2)%scaling="upper";zg(2)%yes=1
     END IF
     if(zg(1)%yes==0) THEN
       Write(*,*) "                                     ",&
                  &"View-OutGmV a must resized:"
     END IF
     Write(*,*) "                                     ",&
                &"Domain%view_ixe next",zg(2)%scaling,"cell border!"
  END IF

  ! search between Domain%iy0,Domain%iy1
  !.....................................
  iy=Domain%iy0    ! lower iy
  DO WHILE (Domain%view_ya>Domain%yP(iy).and.iy<Domain%iy1)
     iy=iy+1
  END DO
  IF(Domain%view_ya==Domain%yP(iy)) THEN
     Domain%view_iya=iy
  ELSE
     !Differenz-Check ya/yP(iy) und skalieren
     IF((Domain%yP(iy)-Domain%view_ya) <=(Domain%view_ya-Domain%yP(iy-1))) THEN
        Domain%view_iya=iy
        !zg(3)=(/"iya","upper",1/)
        zg(3)%who="iya";zg(3)%scaling="upper";zg(3)%yes=1
     ELSE
        Domain%view_iya=iy-1
        !zg(3)=(/"iya","lower",1/)
        zg(3)%who="iya";zg(3)%scaling="lower";zg(3)%yes=1
     END IF
     IF(zg(1)%yes==0.eqv.zg(2)%yes==0) THEN
       Write(*,*) "                                     ",&
                  &"View-OutGmV a must resized:"
     END IF
     Write(*,*) "                                     ",&
                &"Domain%view_iya next",zg(3)%scaling,"cell border!"
  END IF
  iy=Domain%iy1   ! upper iy
  DO WHILE (Domain%view_ye<Domain%yP(iy).and.iy>Domain%iy0)
     iy=iy-1
  END DO
  IF(Domain%view_ye==Domain%yP(iy)) THEN
     Domain%view_iye=iy
  ELSE
     !Differenz-Check ye/yP(iy) und skalieren
     IF((Domain%view_ye-Domain%yP(iy)) <= (Domain%yP(iy+1)-Domain%view_ye)) THEN
        Domain%view_iye=iy
        !zg(4)=(/"iye","lower",1/)
        zg(4)%who="iye";zg(4)%scaling="lower";zg(4)%yes=1
     ELSE
        Domain%view_iye=iy+1
        !zg(4)=(/"iye","upper,1"/)
        zg(4)%who="iye";zg(4)%scaling="upper";zg(4)%yes=1
     END IF
     if(zg(1)%yes==0.eqv.zg(2)%yes==0.eqv.zg(3)%yes==0) THEN
       Write(*,*) "                                     ",&
                  &"View-OutGmV a must resized:"
     END IF
     Write(*,*) "                                     ",&
                &"Domain%view_iye next",zg(4)%scaling,"cell border!"
  END IF

  ! search between Domain%iz0,Domain%iz1
  !.....................................
  iz=Domain%iz0   ! lower iz
  DO WHILE (Domain%view_za>Domain%zP(iz).and.iz<Domain%iz1)
     iz=iz+1
  END DO
  IF(Domain%view_za==Domain%zP(iz)) THEN
     Domain%view_iza=iz
  ELSE
     !Differenz-Check za/zP(zy) und skalieren
     IF((Domain%zP(iz)-Domain%view_za) <=(Domain%view_za-Domain%zP(iz-1))) THEN
        Domain%view_iza=iz
        !zg(5)=(/"iza","upper",1/)
        zg(5)%who="iza";zg(5)%scaling="upper";zg(5)%yes=1
     ELSE
        Domain%view_iza=iz-1
        !zg(5)=(/"iza","lower",1/)
        zg(5)%who="iza";zg(5)%scaling="lower";zg(5)%yes=1
     END IF
     if(zg(1)%yes==0.eqv.zg(2)%yes==0.eqv.zg(3)%yes==0.eqv.zg(4)%yes==0) THEN
       Write(*,*) "                                     ",&
                  &"View-OutGmV a must resized:"
     END IF
     Write(*,*) "                                     ",&
                &"Domain%view_iza next",zg(5)%scaling,"cell border!"
  END IF
  iz=Domain%iz1    ! upper iz
  DO WHILE (Domain%view_ze<Domain%zP(iz).and.iz>Domain%iz0)
     iz=iz-1
  END DO
  IF(Domain%view_ze==Domain%zP(iz)) THEN
     Domain%view_ize=iz
  ELSE
     !Differenz-Check ze/zP(iz) und skalieren
     IF((Domain%view_ze-Domain%zP(iz)) <= (Domain%zP(iz+1)-Domain%view_ze)) THEN
        Domain%view_ize=iz
        !zg(6)=(/"ize","lower",1/)
        zg(6)%who="ize";zg(6)%scaling="lower";zg(6)%yes=1
     ELSE
        Domain%view_ize=iz+1
        !zg(6)=(/"ize","upper",1/)
        zg(6)%who="ize";zg(6)%scaling="upper";zg(6)%yes=1
     END IF
     if(zg(1)%yes==0.eqv.zg(2)%yes==0.eqv.zg(3)%yes==0.eqv.zg(4)%yes==0.eqv.zg(5)%yes==0) THEN
       Write(*,*) "                                     ",&
                  &"View-OutGmV a must resized:"
     END IF
     Write(*,*) "                                     ",&
                &"Domain%view_iza next",zg(6)%scaling,"cell border!"
  END IF

  END IF  !IF (...def_out_coord)
  END IF  !IF(out_area=='y')
END SUBROUTINE Set_Domain_OutView

FUNCTION VolCell(Cell)
   REAL(8) :: VolCell
   TYPE (Cell_T) :: Cell
   REAL(8) :: lam1,lam2,r1,r2,phi1,phi2,z2,z1
   !phi-> [-Pi/2,Pi/2]; lam ->[0,2Pi]

   !('Globe')
   IF((indata_type=='gk'.AND.conv_ingrid=='spherical') .OR. & 
      (indata_type=='gk'.AND.out_wahlgrid=='G') .OR. &
       indata_type=='geo' .OR. indata_type=='rad' .OR. &
       indata_type=='wgs84') THEN
     phi1=Cell%Face1%Edge2%Vert1%Point%y
     phi2=Cell%Face1%Edge2%Vert2%Point%y
     lam1=Cell%Face1%Edge1%Vert1%Point%x
     lam2=Cell%Face1%Edge1%Vert2%Point%x
     r1  =RadEarth+Cell%Face1%Edge1%Vert1%Point%z
     r2  =RadEarth+Cell%Face2%Edge1%Vert1%Point%z
     VolCell=1.0d0/3.0d0*((lam2-lam1)*(r2**3-r1**3)*(SIN(phi2)-SIN(phi1)))
   ELSE If (indata_type=='cyl') THEN
     lam1=Cell%Face1%Edge1%Vert1%Point%x
     lam2=Cell%Face1%Edge1%Vert2%Point%x
     r1=Cell%Face1%Edge2%Vert1%Point%y
     r2=Cell%Face1%Edge2%Vert2%Point%y
     z1=Cell%Face1%Edge1%Vert1%Point%z
     z2=Cell%Face2%Edge1%Vert1%Point%z
     VolCell=0.5d0*ABS((lam2-lam1) &
                      *(r2**2-r1**2) &
                      *(z2-z1))
   ELSE !('Cart',DEFAULT)
       VolCell=ABS((Cell%Face1%Edge1%Vert2%Point%x-Cell%Face1%Edge1%Vert1%Point%x) &
                  *(Cell%Face1%Edge2%Vert2%Point%y-Cell%Face1%Edge2%Vert1%Point%y) &
                  *(Cell%Face2%Edge1%Vert1%Point%z-Cell%Face1%Edge1%Vert1%Point%z))
   END IF
END FUNCTION VolCell

FUNCTION VolFace_XY(Face)
   REAL(8) :: VolFace_XY
   TYPE (Face_T) :: Face
   REAL(8) :: r,phi1,phi2,lam1,lam2,r1,r2
   !phi-> [-Pi/2,Pi/2]; lam ->[0,2Pi]

   !('Globe')
   IF((indata_type=='gk'.AND.conv_ingrid=='spherical') .OR. &
      (indata_type=='gk'.AND.out_wahlgrid=='G') .OR. &
       indata_type=='geo' .OR. indata_type=='rad' .OR. &
       indata_type=='wgs84') THEN
     phi1=Face%Edge2%Vert1%Point%y
     phi2=Face%Edge2%Vert2%Point%y
     lam1=Face%Edge1%Vert1%Point%x
     lam2=Face%Edge1%Vert2%Point%x
     r=RadEarth+Face%Edge1%Vert1%Point%z
     VolFace_XY=r**2*(SIN(phi2)-SIN(phi1))*(lam2-lam1)
   ELSE If (indata_type=='cyl') THEN
     lam1=Face%Edge1%Vert1%Point%x
     lam2=Face%Edge1%Vert2%Point%x
     r1=Face%Edge2%Vert1%Point%y
     r2=Face%Edge2%Vert2%Point%y
     VolFace_XY=ABS((r2+r1)*(r2-r1)*(lam2-lam1)/2.0d0)
   ELSE  ! ('Cart',DEFAULT)
     VolFace_XY=ABS((Face%Edge1%Vert2%Point%x-Face%Edge1%Vert1%Point%x) &
                       *(Face%Edge2%Vert2%Point%y-Face%Edge2%Vert1%Point%y))
   END IF
END FUNCTION VolFace_XY

FUNCTION VolFace_YZ(Face)
   REAL(8) :: VolFace_YZ
   TYPE (Face_T) :: Face
   REAL(8) :: phi1,phi2,r1,r2
   !phi-> [-Pi/2;Pi/2]; lam ->[0,2Pi]

   !('Globe')
   IF((indata_type=='gk'.AND.conv_ingrid=='spherical') .OR. & 
      (indata_type=='gk'.AND.out_wahlgrid=='G') .OR. &
       indata_type=='geo' .OR. indata_type=='rad' .OR. &
       indata_type=='wgs84') THEN
     phi1=Face%Edge1%Vert1%Point%y
     phi2=Face%Edge1%Vert2%Point%y
     r1  =RadEarth+Face%Edge2%Vert1%Point%z
     r2  =RadEarth+Face%Edge2%Vert2%Point%z 
     VolFace_YZ=0.5*(phi2-phi1)*(r2**2-r1**2)
   ELSE If (indata_type=='cyl') THEN
     VolFace_YZ=ABS((Face%Edge1%Vert2%Point%y-Face%Edge1%Vert1%Point%y) &
                   *(Face%Edge2%Vert2%Point%z-Face%Edge2%Vert1%Point%z))
   ELSE !('Cart',DEFAULT)
     VolFace_YZ=ABS((Face%Edge1%Vert2%Point%y-Face%Edge1%Vert1%Point%y) &
                   *(Face%Edge2%Vert2%Point%z-Face%Edge2%Vert1%Point%z))
   END IF
END FUNCTION VolFace_YZ

FUNCTION VolFace_ZX(Face)
   REAL(8) :: VolFace_ZX
   TYPE (Face_T) :: Face
   REAL(8) :: phi,lam1,lam2,r1,r2
   !phi-> [-Pi/2,Pi/2]; lam ->[0,2Pi]

   !('Globe')
   IF((indata_type=='gk'.AND.conv_ingrid=='spherical') .OR. &
      (indata_type=='gk'.AND.out_wahlgrid=='G') .OR. &
      indata_type=='geo' .OR. indata_type=='rad' .OR. &
      indata_type=='wgs84') THEN
     phi =Face%Edge1%Vert1%Point%y
     lam1=Face%Edge2%Vert1%Point%x
     lam2=Face%Edge2%Vert2%Point%x
     r1  =RadEarth+Face%Edge1%Vert1%Point%z
     r2  =RadEarth+Face%Edge1%Vert2%Point%z
     VolFace_ZX=0.5*COS(phi)*(lam2-lam1)*(r2**2-r1**2)
   ELSE If (indata_type=='cyl') THEN
     r2 =Face%Edge1%Vert2%Point%y
     VolFace_ZX=ABS((Face%Edge1%Vert2%Point%z-Face%Edge1%Vert1%Point%z) &
                     *(Face%Edge2%Vert2%Point%x-Face%Edge2%Vert1%Point%x)&
                     *r2)
   ELSE !('Cart',DEFAULT)
     VolFace_ZX=ABS((Face%Edge1%Vert2%Point%z-Face%Edge1%Vert1%Point%z) &
                     *(Face%Edge2%Vert2%Point%x-Face%Edge2%Vert1%Point%x))
   END IF

END FUNCTION VolFace_ZX

SUBROUTINE AllocateEdge(Edge,Vert1,Vert2,EdgeType)

  TYPE(Edge_T), POINTER :: Edge
  TYPE(Vertex_T), TARGET :: Vert1,Vert2
  CHARACTER*1 :: EdgeType

  IF (.NOT.ASSOCIATED(Edge)) THEN
    ALLOCATE(Edge)
    Edge%Vert1=>Vert1
    Edge%Vert2=>Vert2
    Edge%yes_sp=-1
    Edge%EdgeType=EdgeType
  END IF
END SUBROUTINE AllocateEdge

SUBROUTINE AllocateFace(Face,Edge1,Edge2,Edge3,Edge4,Vol,FaceType)

  TYPE(Face_T), POINTER :: Face
  TYPE(Edge_T), TARGET :: Edge1,Edge2,Edge3,Edge4
  CHARACTER*2 :: FaceType
  INTERFACE
    FUNCTION Vol(FaceV)
      USE Domain_Mod
      REAL(8) :: Vol
      TYPE(Face_T) :: FaceV
    END FUNCTION Vol
  END INTERFACE

  IF (.NOT.ASSOCIATED(Face)) THEN
    ALLOCATE(Face)
    Face%Edge1=>Edge1
    Face%Edge2=>Edge2
    Face%Edge3=>Edge3
    Face%Edge4=>Edge4
    Face%FaceType=FaceType
  END IF
  Face%Vol=Vol(Face)
END SUBROUTINE AllocateFace

SUBROUTINE AllocateCell(Cell,Face1,Face2,Face3,Face4,Face5,Face6,Vol)

  TYPE(Cell_T), POINTER :: Cell
  TYPE(Face_T), TARGET :: Face1,Face2,Face3,Face4,Face5,Face6
  INTERFACE
    FUNCTION Vol(CellV)
      USE Domain_Mod
      REAL(8) :: Vol
      TYPE(Cell_T) :: CellV
    END FUNCTION Vol
  END INTERFACE

  IF (.NOT.ASSOCIATED(Cell)) THEN
    ALLOCATE(Cell)
    Cell%Face1=>Face1
    Cell%Face2=>Face2
    Cell%Face3=>Face3
    Cell%Face4=>Face4
    Cell%Face5=>Face5
    Cell%Face6=>Face6
  END IF
  Cell%Vol=Vol(Cell)
END SUBROUTINE AllocateCell

SUBROUTINE SearchPoint(InterVert,Vertex1,Vertex2)
  TYPE (Vertex_T) :: InterVert
  TYPE (Vertex_T) :: Vertex1,Vertex2

  REAL(8) :: tL,xL,yL,zL,DistL,lHang
  REAL(8) :: tR,xR,yR,zR,DistR,rHang
  REAL(8) :: t,x,y,z,Dist

  REAL(8), PARAMETER :: EpsDist=1.d-12

  xL=Vertex1%Point%x
  yL=Vertex1%Point%y
  zL=Vertex1%Point%z
  xR=Vertex2%Point%x
  yR=Vertex2%Point%y
  zR=Vertex2%Point%z
  DistL=Level(xL,yL,zL)
  DistR=Level(xR,yR,zR)

  IF (DistL>=0.0d0) THEN
    xL=Vertex2%Point%x
    yL=Vertex2%Point%y
    zL=Vertex2%Point%z
    xR=Vertex1%Point%x
    yR=Vertex1%Point%y
    zR=Vertex1%Point%z
    Dist=DistL
    DistL=DistR
    DistR=Dist
  END IF
  tL=0.0d0
  tR=1.0d0
  IF (ABS(DistL)<ABS(DistR)) THEN
    t=tL
    Dist=DistL
  ELSE
    t=tR
    Dist=DistR
  END IF
  lHang=1.0d0
  rHang=1.0d0
  IF (DistR>=0.0d0) THEN
    DO
      t=t-(tR-tL)/(DistR-DistL)*Dist
      t=MIN(t,tR-1.0d0/3.0d0*(tR-tL))
      t=MAX(t,tL+1.0d0/3.0d0*(tR-tL))
      x=t*xR+(1.0d0-t)*xL
      y=t*yR+(1.0d0-t)*yL
      z=t*zR+(1.0d0-t)*zL
      Dist=Level(x,y,z)
      IF (Dist>=0.0d0) THEN
        tR=t
        DistR=Dist
        DistL=DistL/lHang
        lHang=lHang+1.0d0
        rHang=1.0d0
      ELSE IF (Dist<0.0d0) THEN
        tL=t
        DistL=Dist
        DistR=DistR/rHang
        rHang=rHang+1.0d0
        lHang=1.0d0
      ELSE
        IF (ABS(t-tL)>ABS(t-tR)) THEN
          IF (tL<t) THEN
            t=t-0.5d0*EpsDist
          ELSE
            t=t+0.5d0*EpsDist
          END IF
          tL=t
          x=t*xR+(1.0d0-t)*xL
          y=t*yR+(1.0d0-t)*yL
          z=t*zR+(1.0d0-t)*zL
          DistL=Level(x,y,z)
        ELSE
          IF (tR<t) THEN
            t=t-0.5d0*EpsDist
          ELSE
            t=t+0.5d0*EpsDist
          END IF
          tR=t
          x=t*xR+(1.0d0-t)*xL
          y=t*yR+(1.0d0-t)*yL
          z=t*zR+(1.0d0-t)*zL
          DistR=Level(x,y,z)
        END IF
      END IF
      IF (ABS(tR-tL)<=EpsDist) THEN
        EXIT
      END IF
    END DO
  END IF
  InterVert%Point%x=xL+(xR-xL)*t
  InterVert%Point%y=yL+(yR-yL)*t  
  InterVert%Point%z=zL+(zR-zL)*t 

END SUBROUTINE SearchPoint

SUBROUTINE SortVertex
  INTEGER :: ib,ix,iy,iz
  INTEGER :: nrP
  TYPE(VertexP_T), POINTER :: Vertices(:,:,:)

  Vertices=>Domain%Vertices
  ALLOCATE(VertOut(1:nr_out))
  DO ix=Domain%ix0,Domain%ix1
    DO iy=Domain%iy0,Domain%iy1
      DO iz=Domain%iz0,Domain%iz1
        IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex)) THEN
          nrP=Vertices(ix,iy,iz)%Vertex%nrP
          IF (nrP>0) THEN
            VertOut(nrP)%Vertex=>Vertices(ix,iy,iz)%Vertex
          END IF
          IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex%VertX)) THEN
            nrP=Vertices(ix,iy,iz)%Vertex%VertX%nrP
            VertOut(nrP)%Vertex=>Vertices(ix,iy,iz)%Vertex%VertX
          END IF  
          IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex%VertY)) THEN
            nrP=Vertices(ix,iy,iz)%Vertex%VertY%nrP
            VertOut(nrP)%Vertex=>Vertices(ix,iy,iz)%Vertex%VertY
          END IF  
          IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex%VertZ)) THEN
            nrP=Vertices(ix,iy,iz)%Vertex%VertZ%nrP
            VertOut(nrP)%Vertex=>Vertices(ix,iy,iz)%Vertex%VertZ
          END IF  
        END IF  
      END DO  
    END DO  
  END DO  

END SUBROUTINE SortVertex

SUBROUTINE SortVertexOro
  INTEGER :: ib,ix,iy,iz
  INTEGER :: nrInP
  TYPE(VertexP_T), POINTER :: Vertices(:,:,:)

  Vertices=>Domain%Vertices
  ALLOCATE(VertIn(1:nr_inside))
  DO ix=Domain%ix0,Domain%ix1
    DO iy=Domain%iy0,Domain%iy1
      DO iz=Domain%iz0,Domain%iz1
        IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex)) THEN
          nrInP=Vertices(ix,iy,iz)%Vertex%nrInP
          IF (nrInP>0) THEN
            VertIn(nrInP)%Vertex=>Vertices(ix,iy,iz)%Vertex
          END IF
          IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex%VertX)) THEN
            nrInP=Vertices(ix,iy,iz)%Vertex%VertX%nrInP
            VertIn(nrInP)%Vertex=>Vertices(ix,iy,iz)%Vertex%VertX
          END IF  
          IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex%VertY)) THEN
            nrInP=Vertices(ix,iy,iz)%Vertex%VertY%nrInP
            VertIn(nrInP)%Vertex=>Vertices(ix,iy,iz)%Vertex%VertY
          END IF  
          IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex%VertZ)) THEN
            nrInP=Vertices(ix,iy,iz)%Vertex%VertZ%nrInP
            VertIn(nrInP)%Vertex=>Vertices(ix,iy,iz)%Vertex%VertZ
          END IF  
        END IF  
      END DO  
    END DO  
  END DO  

END SUBROUTINE SortVertexOro

SUBROUTINE SortVertexView
  INTEGER :: ib,i,j,k
  INTEGER :: nrP,iout,nr_set

  Vertices=>Domain%Vertices
  ALLOCATE(VertViewScale(1:nr_out))
  VertViewScale(1:nr_out)=0
  nr_view_out=0
  nr_viewEgVerts=0

  DO ib=1,nb
    CALL Set(Floor(ib))

    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          IF ((i>=Domain%view_ixa.AND.i<=Domain%view_ixe) .AND. &
              (j>=Domain%view_iya.AND.j<=Domain%view_iye) .AND. &
              (k>=Domain%view_iza.AND.k<=Domain%view_ize)) THEN
            nrP=Vertices(i,j,k)%Vertex%nrP
            IF (nrP>0) THEN
              IF (VertViewScale(nrP)==0) THEN
                nr_view_out=nr_view_out+1        ! noch Doppelbelegung Zylinder
                VertViewScale(nrP)=nr_view_out
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          IF ((i>=Domain%view_ixa+1.AND.i<=Domain%view_ixe) .AND. &
              (j>=Domain%view_iya.AND.j<=Domain%view_iye) .AND. &
              (k>=Domain%view_iza.AND.k<=Domain%view_ize)) THEN
            IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
              IF (Edges_X(i,j,k)%Edge%yes_sp==1) THEN
!               nrP=Edges_X(i,j,k)%Edge%VertS%nrP
                nrP=GetNrP(Edges_X(i,j,k)%Edge)
                IF (nrP>0) THEN
                  IF (VertViewScale(nrP)==0) THEN
                    nr_view_out=nr_view_out+1
                    VertViewScale(nrP)=nr_view_out
                    nr_viewEgVerts=nr_viewEgVerts+1
                  END IF
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF ((i>=Domain%view_ixa.AND.i<=Domain%view_ixe) .AND. &
              (j>=Domain%view_iya+1.AND.j<=Domain%view_iye) .AND. &
              (k>=Domain%view_iza.AND.k<=Domain%view_ize)) THEN
            IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
              IF (Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
!               nrP=Edges_Y(i,j,k)%Edge%VertS%nrP
                nrP=GetNrP(Edges_Y(i,j,k)%Edge)
                IF (nrP>0)  THEN
                   IF (VertViewScale(nrP)==0) THEN
                      nr_view_out=nr_view_out+1
                      VertViewScale(nrP)=nr_view_out
                      nr_viewEgVerts=nr_viewEgVerts+1
                   END IF
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF ((i>=Domain%view_ixa.AND.i<=Domain%view_ixe) .AND. &
              (j>=Domain%view_iya.AND.j<=Domain%view_iye) .AND. &
              (k>=Domain%view_iza+1.AND.k<=Domain%view_ize)) THEN
            IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
              IF (Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
!               nrP=Edges_Z(i,j,k)%Edge%VertS%nrP
                nrP=GetNrP(Edges_Z(i,j,k)%Edge)
                IF (nrP>0) THEN
                  IF (VertViewScale(nrP)==0) THEN
                    nr_view_out=nr_view_out+1
                    VertViewScale(nrP)=nr_view_out
                    nr_viewEgVerts=nr_viewEgVerts+1
                  END IF
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO

  END DO   ! ib

  ALLOCATE(VertOutView(1:nr_view_out))
  nr_set=0
  DO iout=1,nr_out
   IF(VertViewScale(iout)>0) THEN
     nr_set=nr_set+1
     VertOutView(VertViewScale(iout))=VertOut(iout)
   END IF
  END DO
  IF(nr_set/=nr_view_out) THEN
     Write(OutUnitProt,*) "---SortVertexView---,nr_set/=nr_view_out"
  END IF
END SUBROUTINE SortVertexView

SUBROUTINE SearchMinMaxEdge(Edge,PMin,PMax)
  TYPE(Edge_T) :: Edge
  TYPE(Point_T) :: PMin,PMax

  TYPE(Vertex_T), POINTER :: VertS=>NULL()

  IF (Edge%Vert1%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Edge%Vert1%Point%x)
    PMin%y=MIN(PMin%y,Edge%Vert1%Point%y)
    PMin%z=MIN(PMin%z,Edge%Vert1%Point%z)
    PMax%x=MAX(PMax%x,Edge%Vert1%Point%x)
    PMax%y=MAX(PMax%y,Edge%Vert1%Point%y)
    PMax%z=MAX(PMax%z,Edge%Vert1%Point%z)
  END IF
  IF (Edge%Vert2%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Edge%Vert2%Point%x)
    PMin%y=MIN(PMin%y,Edge%Vert2%Point%y)
    PMin%z=MIN(PMin%z,Edge%Vert2%Point%z)
    PMax%x=MAX(PMax%x,Edge%Vert2%Point%x)
    PMax%y=MAX(PMax%y,Edge%Vert2%Point%y)
    PMax%z=MAX(PMax%z,Edge%Vert2%Point%z)
  END IF
  SELECT CASE (Edge%EdgeType)
    CASE ('X')
      VertS=>Edge%Vert1%VertX
    CASE ('Y')
      VertS=>Edge%Vert1%VertY
    CASE ('Z')
      VertS=>Edge%Vert1%VertZ
  END SELECT    
  IF (ASSOCIATED(VertS)) THEN
    PMin%x=MIN(PMin%x,VertS%Point%x)
    PMin%y=MIN(PMin%y,VertS%Point%y)
    PMin%z=MIN(PMin%z,VertS%Point%z)
    PMax%x=MAX(PMax%x,VertS%Point%x)
    PMax%y=MAX(PMax%y,VertS%Point%y)
    PMax%z=MAX(PMax%z,VertS%Point%z)
  END IF

END SUBROUTINE SearchMinMaxEdge

SUBROUTINE SearchMinMaxFace(Face,PMin,PMax)
  TYPE(Face_T) :: Face
  TYPE(Point_T) :: PMin,PMax
  CALL SearchMinMaxEdge(Face%Edge1,PMin,PMax)
  CALL SearchMinMaxEdge(Face%Edge2,PMin,PMax)
  CALL SearchMinMaxEdge(Face%Edge3,PMin,PMax)
  CALL SearchMinMaxEdge(Face%Edge4,PMin,PMax)

END SUBROUTINE SearchMinMaxFace

SUBROUTINE SearchMinMaxCell(Cell,PMin,PMax)
  TYPE(Cell_T) :: Cell
  TYPE(Point_T) :: PMin,PMax
  CALL SearchMinMaxFace(Cell%Face1,PMin,PMax)
  CALL SearchMinMaxFace(Cell%Face2,PMin,PMax)
  CALL SearchMinMaxFace(Cell%Face3,PMin,PMax)
  CALL SearchMinMaxFace(Cell%Face4,PMin,PMax)
  CALL SearchMinMaxFace(Cell%Face5,PMin,PMax)
  CALL SearchMinMaxFace(Cell%Face6,PMin,PMax)

END SUBROUTINE SearchMinMaxCell

SUBROUTINE SortVertAllFacesIn

  INTEGER :: ib,ix,iy,iz

  DO ib=1,nb
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0,iz1
          IF (ASSOCIATED(Faces_XY(ix,iy,iz)%Face)) THEN
            CALL SortVertFaceIn(Faces_XY(ix,iy,iz)%Face)
          END IF
        END DO
      END DO
    END DO
    DO ix=ix0,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          IF (ASSOCIATED(Faces_YZ(ix,iy,iz)%Face)) THEN
            CALL SortVertFaceIn(Faces_YZ(ix,iy,iz)%Face)
          END IF
        END DO
      END DO
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0,iy1
        DO iz=iz0+1,iz1
          IF (ASSOCIATED(Faces_ZX(ix,iy,iz)%Face)) THEN
            CALL SortVertFaceIn(Faces_ZX(ix,iy,iz)%Face)
          END IF
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE SortVertAllFacesIn

SUBROUTINE SortVertFaceIn(Face)
  TYPE (Face_T) :: Face

  INTEGER :: i,j,v,ecut
  INTEGER :: EdgeList(1:2,1:6)

  Face%ec=-1
  Face%NumberVert=0
  Face%VertexList=0
  Face%EdgeCut(:)=0
  Face%in_out=Face%Edge1%in_out+Face%Edge3%in_out
  i=0
  ecut=0
  !...............................................
  IF (Face%Edge1%Vert1%in_out<=0) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge1%Vert1%nrInP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge1%yes_sp==1) THEN
    i=i+1
    Face%VertexList(i)=GetNrInP(Face%Edge1)
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%VertexList(i)
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%VertexList(i)
  ELSE IF (Face%Edge1%Vert1%in_out==0) THEN
!   i=i+1
!   Face%VertexList(i)=Face%Edge1%Vert1%nrInP
!   Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge1%Vert1%nrInP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge1%Vert1%nrInP
  END IF
  !............................................
  IF (Face%Edge2%Vert1%in_out<=0) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge2%Vert1%nrInP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge2%yes_sp==1) THEN
    i=i+1
    Face%VertexList(i)=GetNrInP(Face%Edge2)
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%VertexList(i)
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%VertexList(i)
  ELSE IF (Face%Edge2%Vert1%in_out==0) THEN
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge2%Vert1%nrInP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge2%Vert1%nrInP
  END IF
  !..............................................
  IF (Face%Edge3%Vert2%in_out<=0) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge3%Vert2%nrInP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge3%yes_sp==1) THEN
    i=i+1
    Face%VertexList(i)=GetNrInP(Face%Edge3)
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%VertexList(i)
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%VertexList(i)
  ELSE IF (Face%Edge3%Vert2%in_out==0) THEN
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge3%Vert2%nrInP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge3%Vert2%nrInP
  END IF
  !................................................
  IF (Face%Edge4%Vert2%in_out<=0) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge4%Vert2%nrInP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge4%yes_sp==1) THEN
    i=i+1
    Face%VertexList(i)=GetNrInP(Face%Edge4)
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%VertexList(i)
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%VertexList(i)
  ELSE IF (Face%Edge4%Vert2%in_out==0) THEN
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge4%Vert2%nrInP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge4%Vert2%nrInP
  END IF
  IF(ecut>2) THEN
     j=0
     DO i=1,ecut
       IF(EdgeList(1,i)==1) THEN
         j=j+1
       END IF
     END DO 
     IF(j==2) THEN
       v=0
       Face%ec=1
       DO i=1,ecut
         IF(EdgeList(1,i)==1) THEN
           v=v+1
           Face%EdgeCut(v)=EdgeList(2,i)
         END IF
       END DO
     END IF
  END IF
 
  ! -------
  ! Special
  ! -------
  IF (Face%Edge1%Vert1%in_out<0 .AND. Face%Edge3%Vert2%in_out<0 .AND. &
     Face%Edge1%Vert2%in_out>0 .AND. Face%Edge3%Vert1%in_out>0) THEN
     Face%ec=-1
     Face%EdgeCut(:)=0
  ELSE IF (Face%Edge1%Vert1%in_out>0 .AND. Face%Edge3%Vert2%in_out>0 .AND. &
     Face%Edge1%Vert2%in_out<0 .AND. Face%Edge3%Vert1%in_out<0) THEN
     Face%ec=-1
     Face%EdgeCut(:)=0
  END IF
  IF (ecut==1) THEN
    Face%ec=-1
  END IF

END SUBROUTINE SortVertFaceIn 

SUBROUTINE SortVertCutAllCells

  INTEGER :: ib,ix,iy,iz
  nr_cutinside=0
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
           IF (ASSOCIATED(Cells(ix,iy,iz)%Cell)) THEN
             CALL SortVertCutCell(Cells(ix,iy,iz)%Cell,ix,iy,iz)
           END IF
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE SortVertCutAllCells

SUBROUTINE SortVertCutCell(Cell,ix,iy,iz,Print)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: ix,iy,iz
  INTEGER :: ListCut(1:2,1:8) !ListCut(1:2,1:6)
  INTEGER :: nCut,nCutNeu,i,j,iv
  INTEGER :: VertexFace(1:6,1:9)
  LOGICAL :: Cut
  INTEGER :: iP1,iP2,jP1,jP2,iVert,jVert
  INTEGER :: vc_out1
  TYPE(Face_T) :: Face
  LOGICAL, OPTIONAL :: Print

  IF (Cell%vc>0) THEN
    VertexFace=-1
    nCut=0
    ListCut=0
    Cell%in_out=Cell%Face1%in_out+Cell%Face2%in_out

    IF (Present(Print)) THEN
      CALL WriteCellDebug(Cell)
    END IF  
    IF (Cell%Face1%NumberVert>2) THEN
      VertexFace(1,1:Cell%Face1%NumberVert)=Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
      VertexFace(1,Cell%Face1%NumberVert+1)=VertexFace(1,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 1'
        WRITE(*,*) VertexFace(1,1:Cell%Face1%NumberVert+1)
      END IF  
    END IF
    IF (Cell%Face2%NumberVert>2) THEN
      VertexFace(2,1:Cell%Face2%NumberVert)=Cell%Face2%VertexList(1:Cell%Face2%NumberVert)
      VertexFace(2,Cell%Face2%NumberVert+1)=VertexFace(2,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 2'
        WRITE(*,*) VertexFace(2,1:Cell%Face2%NumberVert+1)
      END IF  
    END IF
    IF (Cell%Face3%NumberVert>2) THEN
      VertexFace(3,1:Cell%Face3%NumberVert)=Cell%Face3%VertexList(1:Cell%Face3%NumberVert)
      VertexFace(3,Cell%Face3%NumberVert+1)=VertexFace(3,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 3'
        WRITE(*,*) VertexFace(3,1:Cell%Face3%NumberVert+1)
      END IF  
    END IF
    IF (Cell%Face4%NumberVert>2) THEN
      VertexFace(4,1:Cell%Face4%NumberVert)=Cell%Face4%VertexList(1:Cell%Face4%NumberVert)
      VertexFace(4,Cell%Face4%NumberVert+1)=VertexFace(4,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 4'
        WRITE(*,*) VertexFace(4,1:Cell%Face4%NumberVert+1)
      END IF  
    END IF
    IF (Cell%Face5%NumberVert>2) THEN
      VertexFace(5,1:Cell%Face5%NumberVert)=Cell%Face5%VertexList(1:Cell%Face5%NumberVert)
      VertexFace(5,Cell%Face5%NumberVert+1)=VertexFace(5,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 5'
        WRITE(*,*) VertexFace(5,1:Cell%Face5%NumberVert+1)
      END IF  
    END IF
    IF (Cell%Face6%NumberVert>2) THEN
      VertexFace(6,1:Cell%Face6%NumberVert)=Cell%Face6%VertexList(1:Cell%Face6%NumberVert)
      VertexFace(6,Cell%Face6%NumberVert+1)=VertexFace(6,1)
      IF (PRESENT(Print)) THEN
        WRITE(*,*) 'Face 6'
        WRITE(*,*) VertexFace(6,1:Cell%Face6%NumberVert+1)
      END IF  
    END IF

    C1:DO i=1,6
      C2:DO iVert=1,8
        iP1=VertexFace(i,iVert)
        iP2=VertexFace(i,iVert+1)
        IF (ip2==-1) THEN
          EXIT C2
        END IF
        Cut=.TRUE.
        C3:DO j=1,6
          IF (j/=i) THEN
            C4:DO jVert=1,8
              jP1=VertexFace(j,jVert)
              jP2=VertexFace(j,jVert+1)
              IF (jp2==-1) THEN
                EXIT C4
              END IF
              IF ((iP1==jP1.AND.iP2==jP2).OR.  &
                  (iP1==jP2.AND.iP2==jP1)) THEN
                Cut=.FALSE.
                EXIT C3
              END IF
            END DO C4
          END IF
        END DO C3
        IF (Cut) THEN
          ncut=ncut+1
          ListCut(1,nCut)=iP1
          ListCut(2,nCut)=iP2
        END IF
      END DO C2
    END DO C1

    IF (nCut>0) THEN
      nCutNeu=nCut
      DO i=1,ncut-1
        DO j=i+1,ncut
          IF (( (ListCut(1,i)==ListCut(1,j).AND.ListCut(2,i)==ListCut(2,j)).OR. &
                (ListCut(1,i)==ListCut(2,j).AND.ListCut(2,i)==ListCut(1,j)) ) &
             .AND.ListCut(1,i)>0) THEN
            ListCut(:,j)=0
            nCutNeu=nCutNeu-1
          END IF
        END DO
      END DO
      Cell%vc=nCutNeu
      IF (.NOT.ASSOCIATED(Cell%VertCut)) THEN
         ALLOCATE(Cell%VertCut(nCut))
      ELSE
         DEALLOCATE(Cell%VertCut)
         ALLOCATE(Cell%VertCut(nCut))
      END IF
      Cell%VertCut(1)=ListCut(1,1)
      Cell%VertCut(2)=ListCut(2,1)
      ListCut(1,1)=0
      ListCut(2,1)=0
      iv=2
      S2:DO
        S1:DO i=1,nCut
          IF (Cell%VertCut(iv)==ListCut(1,i)) THEN
            iv=iv+1
            Cell%VertCut(iv)=ListCut(2,i)
            ListCut(:,i)=0
            EXIT S1
          ELSE IF (Cell%VertCut(iv)==ListCut(2,i)) THEN
            iv=iv+1
            Cell%VertCut(iv)=ListCut(1,i)
            ListCut(:,i)=0
            EXIT S1
          END IF
        END DO S1
        IF (iv>=nCutNeu) THEN
          EXIT S2
        END IF
      END DO S2
    ELSE
      Cell%vc=0
    END IF
  END IF

  IF (Cell%vc==1) THEN
    Cell%vc=0
  END IF
END SUBROUTINE SortVertCutCell

SUBROUTINE SortVertexSoil
  INTEGER :: i,j
  INTEGER :: ib,ix,iy,iz
  INTEGER :: nrP
  INTEGER :: nrCutP
  TYPE (Vertex_T) :: V(1:2)
 
  IF(nr_sb>0) THEN
     nr_out_soil=(1+nr_lsoil)*nr_cutplane
     ALLOCATE(VertSoilOut(1:nr_out_soil)) 
     DO ib=1,nb
       CALL Set(Floor(ib))
       DO ix=ix0,ix1
         DO iy=iy0,iy1
           DO iz=iz0,iz1
             nrP=Vertices(ix,iy,iz)%Vertex%nrP
             nrCutP=Vertices(ix,iy,iz)%Vertex%nrCutP
             IF (nrCutP>0) THEN
               VertSoilOut(nrCutP)%Vertex=>Vertices(ix,iy,iz)%Vertex
             END IF
             IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex%VertX)) THEN
               nrCutP=Vertices(ix,iy,iz)%Vertex%VertX%nrCutP
               VertSoilOut(nrCutP)%Vertex=>Vertices(ix,iy,iz)%Vertex%VertX
             END IF  
             IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex%VertY)) THEN
               nrCutP=Vertices(ix,iy,iz)%Vertex%VertY%nrCutP
               VertSoilOut(nrCutP)%Vertex=>Vertices(ix,iy,iz)%Vertex%VertY
             END IF  
             IF (ASSOCIATED(Vertices(ix,iy,iz)%Vertex%VertZ)) THEN
               nrCutP=Vertices(ix,iy,iz)%Vertex%VertZ%nrCutP
               VertSoilOut(nrCutP)%Vertex=>Vertices(ix,iy,iz)%Vertex%VertZ
             END IF  
           END DO
         END DO
       END DO
     END DO   ! ib
     V(1:2)%Point%x=Domain%xP(0)
     V(1:2)%Point%y=Domain%yP(0)
     V(1)%Point%z=Domain%zP(0)
     V(2)%Point%z=Domain%zP(1)
     dzi_soil=zParametricOut(V(2))-zParametricOut(V(1))
     DO i=2,nr_lsoil+1
       DO j=1,nr_cutplane
         ALLOCATE(VertSoilOut(((i-1)*nr_cutplane)+j)%Vertex)
         VertSoilOut(((i-1)*nr_cutplane)+j)%Vertex=VertSoilOut(j)%Vertex
         VertSoilOut(((i-1)*nr_cutplane)+j)%Vertex%Point%z= &
          VertSoilOut(((i-1)*nr_cutplane)+j)%Vertex%Point%z-(i-1)*dzi_soil !dz(1) !dz(1)/2 !200.0d0 !dz(1)/7
         VertSoilOut(((i-1)*nr_cutplane)+j)%Vertex%nrCutP= &
          VertSoilOut(((i-1)*nr_cutplane)+j)%Vertex%nrCutP+((i-1)*nr_cutplane)
       END DO 
     END DO 
   END IF ! if(nr_sb>0)

END SUBROUTINE SortVertexSoil
END MODULE GridNeu_Mod
