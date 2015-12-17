MODULE Neighbor_Mod


  USE Kind_Mod
  USE Geometry_Mod
  USE Polygon_Mod

  IMPLICIT NONE

  TYPE Nachbar_T
    CHARACTER(2) :: nTYPE
    INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1    ! Coordinates of border in block coordinates
    INTEGER :: iNx0,iNx1,iNy0,iNy1,iNz0,iNz1 ! Coordinates of border in neighbour coordinates
    INTEGER :: ixO,ixI,iyO,iyI,izO,izI
    INTEGER :: iNxO,iNxI,iNyO,iNyI,iNzO,iNzI
    INTEGER :: ib
    INTEGER :: Refine
    INTEGER :: RefineX
    INTEGER :: RefineY
    INTEGER :: RefineZ
    INTEGER :: IncrX
    INTEGER :: IncrY
    INTEGER :: IncrZ
    INTEGER :: CopyCase
    TYPE(VertexP_T), POINTER :: Vertices(:,:,:)=>NULL()
    TYPE (EdgeP_T), POINTER :: Edges_X(:,:,:)=>NULL()
    TYPE (EdgeP_T), POINTER :: Edges_Y(:,:,:)=>NULL()
    TYPE (EdgeP_T), POINTER :: Edges_Z(:,:,:)=>NULL()
    TYPE (FaceP_T), POINTER :: Faces_XY(:,:,:)=>NULL()
    TYPE (FaceP_T), POINTER :: Faces_YZ(:,:,:)=>NULL()
    TYPE (FaceP_T), POINTER :: Faces_ZX(:,:,:)=>NULL()
    TYPE (CellP_T), POINTER :: Cells(:,:,:)=>NULL()
  END TYPE Nachbar_T

  INTEGER :: jx0,jx1,jy0,jy1,jz0,jz1
  INTEGER :: jNx0,jNx1,jNy0,jNy1,jNz0,jNz1
  INTEGER :: jxO,jxI,jyO,jyI,jzO,jzI
  INTEGER :: jNxO,jNxI,jNyO,jNyI,jNzO,jNzI
  INTEGER :: ibn
  INTEGER :: RefineNachbar
  INTEGER :: RefineNachbarX
  INTEGER :: RefineNachbarY
  INTEGER :: RefineNachbarZ
  INTEGER :: IncrX,IncrY,IncrZ
  INTEGER :: CopyCase
  CHARACTER(2) :: nType
  TYPE(VertexP_T), POINTER :: VerticesNachbar(:,:,:)=>NULL()
  TYPE (EdgeP_T), POINTER :: Edges_XNachbar(:,:,:)=>NULL()
  TYPE (EdgeP_T), POINTER :: Edges_YNachbar(:,:,:)=>NULL()
  TYPE (EdgeP_T), POINTER :: Edges_ZNachbar(:,:,:)=>NULL()
  TYPE (FaceP_T), POINTER :: Faces_XYNachbar(:,:,:)=>NULL()
  TYPE (FaceP_T), POINTER :: Faces_YZNachbar(:,:,:)=>NULL()
  TYPE (FaceP_T), POINTER :: Faces_ZXNachbar(:,:,:)=>NULL()
  TYPE (CellP_T), POINTER :: CellsNachbar(:,:,:)=>NULL()

  INTERFACE Set
    MODULE PROCEDURE SetNachbar
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE NachbarAllocate
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE NachbarDeallocate
  END INTERFACE

CONTAINS

SUBROUTINE SetNachbar(Nachbar)

  TYPE (Nachbar_T) :: Nachbar

  ibn=Nachbar%ib
  jx0=Nachbar%ix0
  jx1=Nachbar%ix1
  jy0=Nachbar%iy0
  jy1=Nachbar%iy1
  jz0=Nachbar%iz0
  jz1=Nachbar%iz1
  jNx0=Nachbar%iNx0
  jNx1=Nachbar%iNx1
  jNy0=Nachbar%iNy0
  jNy1=Nachbar%iNy1
  jNz0=Nachbar%iNz0
  jNz1=Nachbar%iNz1
  jxO=Nachbar%ixO
  jxI=Nachbar%ixI
  jyO=Nachbar%iyO
  jyI=Nachbar%iyI
  jzO=Nachbar%izO
  jzI=Nachbar%izI
  jNxO=Nachbar%iNxO
  jNxI=Nachbar%iNxI
  jNyO=Nachbar%iNyO
  jNyI=Nachbar%iNyI
  jNzO=Nachbar%iNzO
  jNzI=Nachbar%iNzI
  RefineNachbar=Nachbar%Refine
  RefineNachbarX=Nachbar%RefineX
  RefineNachbarY=Nachbar%RefineY
  RefineNachbarZ=Nachbar%RefineZ
  IncrX=Nachbar%IncrX
  IncrY=Nachbar%IncrY
  IncrZ=Nachbar%IncrZ
  CopyCase=Nachbar%CopyCase
  nType=Nachbar%nType
  VerticesNachbar=>Nachbar%Vertices
  Edges_XNachbar=>Nachbar%Edges_X
  Edges_YNachbar=>Nachbar%Edges_Y
  Edges_ZNachbar=>Nachbar%Edges_Z
  Faces_XYNachbar=>Nachbar%Faces_XY
  Faces_YZNachbar=>Nachbar%Faces_YZ
  Faces_ZXNachbar=>Nachbar%Faces_ZX
  CellsNachbar=>Nachbar%Cells

END SUBROUTINE SetNachbar

SUBROUTINE NachbarAllocate(Nachbar)
  
  TYPE (Nachbar_T) :: Nachbar
  CALL Set(Nachbar)
  IF (nType=='iw'.OR.nType=='ie') THEN
    ALLOCATE(Nachbar%Vertices(jx0:jx1,jy0:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Edges_X(jx0+1:jx1,jy0:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Edges_Y(jx0:jx1,jy0+1:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Edges_Z(jx0:jx1,jy0:jy1,jz0+1:jz1))
    ALLOCATE(Nachbar%Faces_XY(jx0+1:jx1,jy0+1:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Faces_YZ(jx0:jx1,jy0+1:jy1,jz0+1:jz1))
    ALLOCATE(Nachbar%Faces_ZX(jx0+1:jx1,jy0:jy1,jz0+1:jz1))
    ALLOCATE(Nachbar%Cells(jx0+1:jx1,jy0+1:jy1,jz0+1:jz1))
  ELSE IF (nType=='is'.OR.nType=='in') THEN
    ALLOCATE(Nachbar%Vertices(jx0:jx1,jy0:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Edges_X(jx0+1:jx1,jy0:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Edges_Y(jx0:jx1,jy0+1:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Edges_Z(jx0:jx1,jy0:jy1,jz0+1:jz1))
    ALLOCATE(Nachbar%Faces_XY(jx0+1:jx1,jy0+1:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Faces_YZ(jx0:jx1,jy0+1:jy1,jz0+1:jz1))
    ALLOCATE(Nachbar%Faces_ZX(jx0+1:jx1,jy0:jy1,jz0+1:jz1))
    ALLOCATE(Nachbar%Cells(jx0+1:jx1,jy0+1:jy1,jz0+1:jz1))
  ELSE IF (nType=='ib'.OR.nType=='it') THEN
    ALLOCATE(Nachbar%Vertices(jx0:jx1,jy0:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Edges_X(jx0+1:jx1,jy0:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Edges_Y(jx0:jx1,jy0+1:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Edges_Z(jx0:jx1,jy0:jy1,jz0+1:jz1))
    ALLOCATE(Nachbar%Faces_XY(jx0+1:jx1,jy0+1:jy1,jz0:jz1))
    ALLOCATE(Nachbar%Faces_YZ(jx0:jx1,jy0+1:jy1,jz0+1:jz1))
    ALLOCATE(Nachbar%Faces_ZX(jx0+1:jx1,jy0:jy1,jz0+1:jz1))
    ALLOCATE(Nachbar%Cells(jx0+1:jx1,jy0+1:jy1,jz0+1:jz1))
  END IF
    
END SUBROUTINE NachbarAllocate

SUBROUTINE NachbarDeallocate(Nachbar)

  TYPE (Nachbar_T) :: Nachbar

  INTEGER :: jx,jy,jz

  CALL Set(Nachbar)
  IF (ASSOCIATED(VerticesNachbar)) THEN
    DEALLOCATE(VerticesNachbar)
  END IF
  IF (ASSOCIATED(Edges_XNachbar)) THEN
    DO jx=jx0+1,jx1
      DO jy=jy0,jy1
        DO jz=jz0,jz1
          IF (ASSOCIATED(Edges_XNachbar(jx,jy,jz)%Edge)) THEN
            DEALLOCATE(Edges_XNachbar(jx,jy,jz)%Edge)
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(Edges_XNachbar)
  END IF  
  IF (ASSOCIATED(Edges_YNachbar)) THEN
    DO jx=jx0,jx1
      DO jy=jy0+1,jy1
        DO jz=jz0,jz1
          IF (ASSOCIATED(Edges_YNachbar(jx,jy,jz)%Edge)) THEN
            DEALLOCATE(Edges_YNachbar(jx,jy,jz)%Edge)
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(Edges_YNachbar)
  END IF  
  IF (ASSOCIATED(Edges_ZNachbar)) THEN
    DO jx=jx0,jx1
      DO jy=jy0,jy1
        DO jz=jz0+1,jz1
          IF (ASSOCIATED(Edges_ZNachbar(jx,jy,jz)%Edge)) THEN
            DEALLOCATE(Edges_ZNachbar(jx,jy,jz)%Edge)
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(Edges_ZNachbar)
  END IF  
  IF (ASSOCIATED(Faces_XYNachbar)) THEN
    DO jx=jx0+1,jx1
      DO jy=jy0+1,jy1
        DO jz=jz0,jz1
          IF (ASSOCIATED(Faces_XYNachbar(jx,jy,jz)%Face)) THEN
            DEALLOCATE(Faces_XYNachbar(jx,jy,jz)%Face)
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(Faces_XYNachbar)
  END IF  
  IF (ASSOCIATED(Faces_ZXNachbar)) THEN
    DO jx=jx0+1,jx1
      DO jy=jy0,jy1
        DO jz=jz0+1,jz1
          IF (ASSOCIATED(Faces_ZXNachbar(jx,jy,jz)%Face)) THEN
            DEALLOCATE(Faces_ZXNachbar(jx,jy,jz)%Face)
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(Faces_ZXNachbar)
  END IF  
  IF (ASSOCIATED(Faces_YZNachbar)) THEN
    DO jx=jx0,jx1
      DO jy=jy0+1,jy1
        DO jz=jz0+1,jz1
          IF (ASSOCIATED(Faces_YZNachbar(jx,jy,jz)%Face)) THEN
            DEALLOCATE(Faces_YZNachbar(jx,jy,jz)%Face)
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(Faces_YZNachbar)  
  END IF  
  IF (ASSOCIATED(CellsNachbar)) THEN
    DO jx=jx0+1,jx1
      DO jy=jy0+1,jy1
        DO jz=jz0+1,jz1
          IF (ASSOCIATED(CellsNachbar(jx,jy,jz)%Cell)) THEN
            DEALLOCATE(CellsNachbar(jx,jy,jz)%Cell)
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(CellsNachbar)
  END IF  

END SUBROUTINE NachbarDeallocate

END MODULE Neighbor_Mod

