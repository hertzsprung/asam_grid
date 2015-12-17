!################################################################################
! OutputWeightBlk_Mod.f90
!  (Ausgabe Koordinaten,Weight-,MittelPoint-, Faces-Cellen zum Input ASAM-Model)
!     -> NrW_*  : alle angrenzenden am Berg (geschnittene) und unterhalb Berg
!     -> NrMP_* : alle geschnittenen am Berg
!     -> NrB_   : alle Cellen, geschnittene und komplett mit Face1 angrenzende,
!                 (NrB==NrMP) der Cellen
!     -> Soil-Struct Output
!
!------------------------------------------------------------------------------
! Description:  -Ausgabe dx,dy,dz,Faces,Cellen
!                     -> jeweils Block/Blöcke
!                     -> In-/Out Border east,-west,-north,-south,-top 
!                     -> Border 'ob', 'out bottom' nicht
!               -Ausgabe Numbers Weight:
!                     -> NrW_FacesYZ, NrW_FacesZX, NrW_FacesXY
!                     -> NrW_Cells 
!               -Ausgabe Numbers MittelPoint, Bound, Bound_Struct
!                     -> NrMP_FacesYZ, NrMP_FacesZX, NrMP_FacesXY
!                     -> NrB_Cells, struct_bound
!------------------------------------------------------------------------------
!################################################################################

!------------------------------------------------------------------------------
!################################################################################

MODULE OutputWeightBlk_Mod

  USE Floor_Mod

  USE Function_Mod ,  ONLY:   &
               conv_gk

  USE Parametric_Mod ,  ONLY:  &
               indata_type,OutGrid

  USE GridNeu_Mod, ONLY: &
               VolFace_XY,VolFace_YZ,VolFace_ZX

  USE IOControl_Mod

  IMPLICIT NONE

!================================================================================
CONTAINS
!================================================================================

SUBROUTINE WeightNachbar(WeightIn,WeightBoundary,NrW)
  REAL(RealKind) :: WeightIn(:,:)
  REAL(RealKind) :: WeightBoundary(:,:)
  INTEGER :: NrW

  INTEGER :: i,j

  DO i=1,SIZE(WeightIn,1)
    DO j=1,SIZE(WeightIn,2)
      IF (WeightIn(i,j)>=0.0d0) THEN
        WeightBoundary(i,j)=WeightIn(i,j)
        NrW=NrW+1
      END IF  
    END DO  
  END DO  
END SUBROUTINE WeightNachbar

SUBROUTINE WriteWeight(ix0,ix1,iy0,iy1,iz0,iz1,Weight)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  REAL(RealKind) :: Weight(ix0:ix1,iy0:iy1,iz0:iz1)

  INTEGER :: ix,iy,iz

  DO ix=ix0,ix1
    DO iy=iy0,iy1
      DO iz=iz0,iz1
        IF (Weight(ix,iy,iz)>=0.0d0) THEN
          WRITE(10,*) ix,iy,iz
          WRITE(10,*) Weight(ix,iy,iz)
        ELSE  
          Weight(ix,iy,iz)=0.0d0 !OSSI
        END IF
      END DO  
    END DO  
  END DO  
END SUBROUTINE WriteWeight

SUBROUTINE WriteWeightAll(ix0,ix1,iy0,iy1,iz0,iz1,Weight)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  REAL(RealKind) :: Weight(ix0:ix1,iy0:iy1,iz0:iz1)

  INTEGER :: ix,iy,iz

  WRITE(11,*) 'ix0....',ix0,ix1,iy0,iy1,iz0,iz1
  DO ix=ix0,ix1
    DO iy=iy0,iy1
      DO iz=iz0,iz1
        IF (ix<=ix0+1.OR.ix>=ix1-1) THEN
        WRITE(11,*) ix,iy,iz
        WRITE(11,*) Weight(ix,iy,iz)
        END IF
      END DO  
    END DO  
  END DO  
END SUBROUTINE WriteWeightAll
SUBROUTINE WriteMidPointFaces(ix0,ix1,iy0,iy1,iz0,iz1,Faces,Name)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  TYPE(FaceP_T) :: Faces(ix0:ix1,iy0:iy1,iz0:iz1)
  CHARACTER*15 :: Name

  INTEGER :: ix,iy,iz
  INTEGER :: NrMP

  NrMP=0
  DO iz=iz0,iz1
    DO iy=iy0,iy1
      DO ix=ix0,ix1
        IF (ASSOCIATED(Faces(ix,iy,iz)%Face).AND.Faces(ix,iy,iz)%Face%mp>0) THEN
          NrMP=NrMP+1
        END IF
      END DO
    END DO
  END DO
  WRITE(10,*) NrMP,Name
  DO iz=iz0,iz1
    DO iy=iy0,iy1
      DO ix=ix0,ix1
        IF (ASSOCIATED(Faces(ix,iy,iz)%Face).AND.Faces(ix,iy,iz)%Face%mp>0) THEN
          WRITE(10,*) ix,iy,iz
          WRITE(10,*) Faces(ix,iy,iz)%Face%MidPoint
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE WriteMidPointFaces

SUBROUTINE CountWeightFacesYZ(ix0,ix1,iy0,iy1,iz0,iz1,Vertices,Faces,Weight)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  TYPE (VertexP_T) :: Vertices(ix0:ix1,iy0:iy1,iz0:iz1)
  TYPE (FaceP_T) :: Faces(ix0:ix1,iy0+1:iy1,iz0+1:iz1)
  REAL(RealKind) :: Weight(ix0:ix1,iy0+1:iy1,iz0+1:iz1)

  INTEGER :: ix,iy,iz
  INTEGER :: in_out
  TYPE(Face_T), POINTER :: Face

  DO ix=ix0,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        IF (ASSOCIATED(Faces(ix,iy,iz)%Face)) THEN
          Face=>Faces(ix,iy,iz)%Face
          IF ((Face%in_out<4.AND.Face%ec>0).OR. &
             (Face%in_out<0.AND.Face%ec==-1).OR. &
              Face%Vol==0.0d0) THEN
            NrRW_FacesYZ=NrRW_FacesYZ+1
            Weight(ix,iy,iz)=Face%Vol
          END IF
        ELSE
          in_out=Vertices(ix,iy-1,iz-1)%Vertex%in_out &
                +Vertices(ix,iy,iz-1)%Vertex%in_out &
                +Vertices(ix,iy-1,iz)%Vertex%in_out &
                +Vertices(ix,iy,iz)%Vertex%in_out
          IF (in_out<=0) THEN
            Weight(ix,iy,iz)=0.0d0
            NrRW_FacesYZ=NrRW_FacesYZ+1
          END IF
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE CountWeightFacesYZ

SUBROUTINE CountWeightFacesZX(ix0,ix1,iy0,iy1,iz0,iz1,Vertices,Faces,Weight)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  TYPE (VertexP_T) :: Vertices(ix0:ix1,iy0:iy1,iz0:iz1)
  TYPE (FaceP_T) :: Faces(ix0+1:ix1,iy0:iy1,iz0+1:iz1)
  REAL(RealKind) :: Weight(ix0+1:ix1,iy0:iy1,iz0+1:iz1)

  INTEGER :: ix,iy,iz
  INTEGER :: in_out
  TYPE(Face_T), POINTER :: Face

  DO ix=ix0+1,ix1
    DO iy=iy0,iy1
      DO iz=iz0+1,iz1
        IF (ASSOCIATED(Faces(ix,iy,iz)%Face)) THEN
          Face=>Faces(ix,iy,iz)%Face
          IF ((Face%in_out<4.AND.Face%ec>0).OR. &
              (Face%in_out<0.AND.Face%ec==-1).OR. &
               Face%Vol==0.0d0) THEN
            Weight(ix,iy,iz)=Face%Vol
            NrRW_FacesZX=NrRW_FacesZX+1
          END IF
        ELSE
          in_out=Vertices(ix-1,iy,iz-1)%Vertex%in_out &
                +Vertices(ix,iy,iz-1)%Vertex%in_out &
                +Vertices(ix-1,iy,iz)%Vertex%in_out &
                +Vertices(ix,iy,iz)%Vertex%in_out
          IF (in_out<=0) THEN
            Weight(ix,iy,iz)=0.0d0
            NrRW_FacesZX=NrRW_FacesZX+1
          END IF
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE CountWeightFacesZX

SUBROUTINE CountWeightFacesXY(ix0,ix1,iy0,iy1,iz0,iz1,Vertices,Faces,Weight)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  TYPE (VertexP_T) :: Vertices(ix0:ix1,iy0:iy1,iz0:iz1)
  TYPE (FaceP_T) :: Faces(ix0+1:ix1,iy0+1:iy1,iz0:iz1)
  REAL(RealKind) :: Weight(ix0+1:ix1,iy0+1:iy1,iz0:iz1)

  INTEGER :: ix,iy,iz
  INTEGER :: in_out
  TYPE(Face_T), POINTER :: Face

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0,iz1
        IF (ASSOCIATED(Faces(ix,iy,iz)%Face)) THEN
          Face=>Faces(ix,iy,iz)%Face
          IF ((Face%in_out<4.AND.Face%ec>0).OR. &
              (Face%in_out<0.AND.Face%ec==-1).OR. & 
               Face%Vol==0.0d0) THEN
            Weight(ix,iy,iz)=Face%Vol
            NrRW_FacesXY=NrRW_FacesXY+1
          END IF
        ELSE
          in_out=Vertices(ix-1,iy-1,iz)%Vertex%in_out &
                +Vertices(ix,iy-1,iz)%Vertex%in_out &
                +Vertices(ix-1,iy,iz)%Vertex%in_out &
                +Vertices(ix,iy,iz)%Vertex%in_out
          IF (in_out<=0) THEN
            Weight(ix,iy,iz)=0.0d0
            NrRW_FacesXY=NrRW_FacesXY+1
          END IF
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE CountWeightFacesXY

SUBROUTINE WriteCutCells(ix0,ix1,iy0,iy1,iz0,iz1,Cells,Name)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  TYPE(CellP_T) :: Cells(ix0:ix1,iy0:iy1,iz0:iz1)
  CHARACTER*15 :: Name

  INTEGER :: ix,iy,iz
  INTEGER :: nl
  INTEGER :: NrMP
  INTEGER :: i

  NrMP=0
  DO iz=iz0,iz1
    DO iy=iy0,iy1
      DO ix=ix0,ix1
        IF (ASSOCIATED(Cells(ix,iy,iz)%Cell).AND.Cells(ix,iy,iz)%Cell%mp>0) THEN
          NrMP=NrMP+1
        END IF
      END DO
    END DO
  END DO
  WRITE(10,'(I12,A16,2(I6,A16))') NrMP,    "      NrB_Cells" &
                              &  , struct_bound, "   Struct_Bound" &
                              &  , nr_lsoil,     "  Nr_SoilLayers"
  DO iz=iz0,iz1
    DO iy=iy0,iy1
      DO ix=ix0,ix1
        IF (ASSOCIATED(Cells(ix,iy,iz)%Cell).AND.Cells(ix,iy,iz)%Cell%mp>0) THEN
          WRITE(10,*) ix,iy,iz
          WRITE(10,*) Cells(ix,iy,iz)%Cell%MidPoint
          WRITE(10,'(6d15.7)') zRauh,Albedo,Emissivity,Cells(ix,iy,iz)%Cell%CutF_MidP
          IF (nr_soildef>0 .AND. nr_landdef==0) THEN ! Bound-Struct: 4
            DO nl=1,nr_soildef
              IF (ix>(s_ixa(nl)*2.e0**RefineX).AND.ix<=(s_ixe(nl)*2.e0**RefineX) .AND. &
                 iy>(s_iya(nl)*2.e0**RefineY).AND.iy<=(s_iye(nl)*2.e0**RefineY)) THEN
                WRITE(10,*) soil_type(nl,1:nr_lsoil)
              END IF
            END DO
          ELSE IF (nr_soildef==0 .AND. nr_landdef>0) THEN  ! Bound-Struct: 5
            WRITE(10,'(I3)') Cells(ix,iy,iz)%Cell%LandClass
          ELSE IF (nr_soildef>0 .AND. nr_landdef>0) THEN  ! Bound-Struct: 6
            DO nl=1,nr_soildef
              IF (ix>(s_ixa(nl)*2.e0**RefineX).AND.ix<=(s_ixe(nl)*2.e0**RefineX) .AND. &
                  iy>(s_iya(nl)*2.e0**RefineY).AND.iy<=(s_iye(nl)*2.e0**RefineY)) THEN
                WRITE(10,*) Cells(ix,iy,iz)%Cell%LandClass,soil_type(nl,1:nr_lsoil)
              END IF
            END DO
          END IF
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE WriteCutCells

SUBROUTINE CountWeightCells(ix0,ix1,iy0,iy1,iz0,iz1,Vertices,Cells,VolC)
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  TYPE (VertexP_T) :: Vertices(ix0:ix1,iy0:iy1,iz0:iz1)
  TYPE (CellP_T) :: Cells(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)
  REAL(RealKind) :: VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)
  
  INTEGER :: ix,iy,iz
  INTEGER :: in_out
  TYPE(Cell_T), POINTER :: Cell

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        IF (ASSOCIATED(Cells(ix,iy,iz)%Cell)) THEN
          Cell=>Cells(ix,iy,iz)%Cell
          IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
          ELSE IF (Cell%in_out<6.AND.Cell%Vol>0.0d0.OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
            VolC(ix,iy,iz)=Cell%Vol
            NrRW_Cells=NrRW_Cells+1
          ELSE
            VolC(ix,iy,iz)=Cell%Vol
            NrRW_Cells=NrRW_Cells+1
          END IF
        ELSE
          in_out=Vertices(ix-1,iy-1,iz-1)%Vertex%in_out &
                +Vertices(ix,iy-1,iz-1)%Vertex%in_out &
                +Vertices(ix-1,iy-1,iz)%Vertex%in_out &
                +Vertices(ix-1,iy,iz-1)%Vertex%in_out &
                +Vertices(ix-1,iy,iz)%Vertex%in_out   &
                +Vertices(ix,iy-1,iz)%Vertex%in_out   &
                +Vertices(ix,iy,iz-1)%Vertex%in_out   &
                +Vertices(ix,iy,iz)%Vertex%in_out
          IF (in_out<=0) THEN
            VolC(ix,iy,iz)=0.0d0
            NrRW_Cells=NrRW_Cells+1
          END IF
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE CountWeightCells

SUBROUTINE CountWeightCell(Cell,i,j,k,ib)
  TYPE(Cell_T) :: Cell
  INTEGER :: i,j,k,ib

  INTEGER :: in_out
  INTEGER :: yes_v,no_v

   yes_v=1    ! yes_Volume
   no_v=0     ! no_volume
   ! RandCellen, alternative da Rand-Faces noch nicht bearbeitet
   IF (ASSOCIATED(Floor(ib)%Cells(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
         !Cellen mit Edge od. Point an Berg
         !Cellen mit MaxVol, nicht unter NrW_ ausgeben
         !maxVol je Face für RandCellen zu checken noch nicht möglich
     ELSE IF (Cell%in_out<6.AND.Cell%Vol>0.0d0.OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       NrRW_Cells=NrRW_Cells+1
     ELSE
        !Bsp./s1agi/reutgen/WRK_GITTER/TEST_GRID_GMV_8.1.1.2/Agi_2Blk_H5000_test.grid
        !    (2.Block Celle 6,1,1 unterhalb am Berg wird allokiert)
        ! Cellen mit Kanten unterhalb am Berg noch filtern in AnalyzeCell
       NrRW_Cells=NrRW_Cells+1
     END IF
   ELSE
        in_out=Floor(ib)%Vertices(i-1,j-1,k-1)%Vertex%in_out &
              +Floor(ib)%Vertices(i,j-1,k-1)%Vertex%in_out &
              +Floor(ib)%Vertices(i-1,j-1,k)%Vertex%in_out &
              +Floor(ib)%Vertices(i-1,j,k-1)%Vertex%in_out &
              +Floor(ib)%Vertices(i-1,j,k)%Vertex%in_out   &
              +Floor(ib)%Vertices(i,j-1,k)%Vertex%in_out   &
              +Floor(ib)%Vertices(i,j,k-1)%Vertex%in_out   &
              +Floor(ib)%Vertices(i,j,k)%Vertex%in_out
        IF (in_out<=0) THEN
          NrRW_Cells=NrRW_Cells+1
        END IF
   END IF

END SUBROUTINE CountWeightCell

SUBROUTINE CountWeightCellBotton(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out
INTEGER :: yes_v,no_v

   yes_v=1    ! yes_Volume
   no_v=0     ! no_volume
   ! RandCellen, alternative da Rand-Faces noch nicht bearbeitet
   IF (ASSOCIATED(Floor(ib)%Cells(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
         !Cellen mit Edge od. Point an Berg
         !Cellen mit MaxVol, nicht unter NrW_ ausgeben
         !maxVol je Face für RandCellen zu checken noch nicht möglich
     ELSE IF ((Cell%in_out<6.AND.Cell%Vol>0.0d0.AND.Cell%vc>0).OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       NrRW_Cells=NrRW_Cells+1
     ELSE
      !wenn celle Vol=0.d0 ist, ist für Botton nicht von Interesse im ASAM 
        !Bsp./s1agi/reutgen/WRK_GITTER/TEST_GRID_GMV_8.1.1.2/Agi_2Blk_H5000_test.grid
        !    (2.Block Celle 6,1,1 unterhalb am Berg wird allokiert)
        ! Cellen mit Kanten unterhalb am Berg noch filtern in AnalyzeCell
      ! NrRW_Cells=NrRW_Cells+1
     END IF
   ELSE
      !wenn celle Vol=0.d0 ist, ist für Botton nicht von Interesse im ASAM 
      !ELSE auch weg lassen
      !  in_out=Floor(ib)%Vertices(i-1,j-1,k-1)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i,j-1,k-1)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j-1,k)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k-1)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j-1,k)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j,k-1)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j,k)%Vertex%in_out
      !  IF (in_out<=0) THEN
      !    !checke Cell z0+1, da Grenzfall 'ob'
      !    in_out=Floor(ib)%Vertices(i-1,j-1,k)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i,j-1,k)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j-1,k+1)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k+1)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j-1,k+1)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j,k)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j,k+1)%Vertex%in_out
      !    IF (in_out<=0) THEN
      !      NrRW_Cells=NrRW_Cells+1
      !    END IF
      !  END IF
   END IF

END SUBROUTINE CountWeightCellBotton

SUBROUTINE OutWeightCell(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out
INTEGER :: yes_v,no_v

   yes_v=1    ! yes_Volume
   no_v=0     ! no_volume
   ! RandCellen, alternative da Rand-Faces noch nicht bearbeitet
   IF (ASSOCIATED(Floor(ib)%Cells(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
         !Cellen mit Edge od. Point an Berg
         !Cellen mit MaxVol, nicht unter NrW_ ausgeben
         !maxVol je Face für RandCellen zu checken noch nicht möglich
     ELSE IF (Cell%in_out<6.AND.Cell%Vol>0.0d0.OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) Cell%Vol
       !WRITE(10,*) yes_v,Cell%Vol
       !WRITE(10,*) Cell(i,j,k)%Cell%Vol,Cell(i,j,k)%Cell%MidPoint
     ELSE
        !Bsp./s1agi/reutgen/WRK_GITTER/TEST_GRID_GMV_8.1.1.2/Agi_2Blk_H5000_test.grid
        !    (2.Block Celle 6,1,1 unterhalb am Berg wird allokiert)
        ! Cellen mit Kanten unterhalb am Berg noch filtern in AnalyzeCell
       WRITE(10,*) i,j,k
       WRITE(10,*) 0.0d0
       !WRITE(10,*) no_v,0.0d0
     END IF
   ELSE
        in_out=Floor(ib)%Vertices(i-1,j-1,k-1)%Vertex%in_out &
              +Floor(ib)%Vertices(i,j-1,k-1)%Vertex%in_out &
              +Floor(ib)%Vertices(i-1,j-1,k)%Vertex%in_out &
              +Floor(ib)%Vertices(i-1,j,k-1)%Vertex%in_out &
              +Floor(ib)%Vertices(i-1,j,k)%Vertex%in_out   &
              +Floor(ib)%Vertices(i,j-1,k)%Vertex%in_out   &
              +Floor(ib)%Vertices(i,j,k-1)%Vertex%in_out   &
              +Floor(ib)%Vertices(i,j,k)%Vertex%in_out
        IF (in_out<=0) THEN
          WRITE(10,*) i,j,k
          !WRITE(10,*) no_v,0.0d0
          WRITE(10,*) 0.0d0
        END IF
   END IF

END SUBROUTINE OutWeightCell

SUBROUTINE OutWeightCellBotton(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out
INTEGER :: yes_v,no_v

   yes_v=1    ! yes_Volume
   no_v=0     ! no_volume
   ! RandCellen, alternative da Rand-Faces noch nicht bearbeitet
   IF (ASSOCIATED(Floor(ib)%Cells(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
         !Cellen mit Edge od. Point an Berg
         !Cellen mit MaxVol, nicht unter NrW_ ausgeben
         !maxVol je Face für RandCellen zu checken noch nicht möglich
     ELSE IF ((Cell%in_out<6.AND.Cell%Vol>0.0d0.AND.Cell%vc>0).OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       WRITE(10,*) i,j,k-1
       WRITE(10,*) Cell%Vol
       !WRITE(10,*) yes_v,Cell%Vol
       !WRITE(10,*) Cell(i,j,k)%Cell%Vol,Cell(i,j,k)%Cell%MidPoint
     ELSE
      !wenn celle Vol=0.d0 ist, ist für Botton nicht von Interesse im ASAM 
        !Bsp./s1agi/reutgen/WRK_GITTER/TEST_GRID_GMV_8.1.1.2/Agi_2Blk_H5000_test.grid
        !    (2.Block Celle 6,1,1 unterhalb am Berg wird allokiert)
        ! Cellen mit Kanten unterhalb am Berg noch filtern in AnalyzeCell
      ! WRITE(10,*) i,j,k-1
      ! WRITE(10,*) 0.0d0
     END IF
   ELSE
      !wenn celle Vol=0.d0 ist, ist für Botton nicht von Interesse im ASAM 
      !ELSE auch weg lassen
      !  in_out=Floor(ib)%Vertices(i-1,j-1,k-1)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i,j-1,k-1)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j-1,k)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k-1)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j-1,k)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j,k-1)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j,k)%Vertex%in_out
      !  IF (in_out<=0) THEN
      !    !WRITE(10,*) i,j,k-1
      !    !WRITE(10,*) 0.0d0
      !    !checke Cell z0+1, da Grenzfall 'ob'
      !    in_out=Floor(ib)%Vertices(i-1,j-1,k)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i,j-1,k)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j-1,k+1)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k)%Vertex%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k+1)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j-1,k+1)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j,k)%Vertex%in_out   &
      !        +Floor(ib)%Vertices(i,j,k+1)%Vertex%in_out
      !    IF (in_out<=0) THEN
      !       WRITE(10,*) i,j,k-1
      !       WRITE(10,*) Cell%Vol
      !    END IF
      !  END IF
   END IF

END SUBROUTINE OutWeightCellBotton

SUBROUTINE OutMPCellBorder(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

   IF (ASSOCIATED(Floor(ib)%Cells(i,j,k)%Cell)) THEN
     IF (Cell%in_out<7.AND.Cell%Vol>0.0d0.AND.Cell%vc>0 ) THEN
       !ausschalten Celle mit VolMax,Edge-Grenze,1Point Grenze
       WRITE(10,*) i,j,k
       WRITE(10,*) Cell%MidPoint
     ELSE IF (Cell%in_out==6.AND.Cell%vc>0) THEN
       !Cell%in_out==6 mit Schnittfläche
       WRITE(10,*) i,j,k
       WRITE(10,*) Cell%MidPoint
     END IF
   END IF

END SUBROUTINE OutMPCellBorder


!SUBROUTINE WriteVerticesNachbar(Vertices,ib,Type)
!  TYPE(VertexP_T) :: Vertices(:,:,:) 
!  INTEGER :: ib
!  CHARACTER :: Type
!  
!  INTEGER :: i,j,k
!  TYPE(Vertex_T), POINTER :: Vertex
!  WRITE(300,*) '  Block ib Type',ib,Type
!  DO i=LBOUND(Vertices,1),UBOUND(Vertices,1)
!    DO j=LBOUND(Vertices,2),UBOUND(Vertices,2)
!     DO k=LBOUND(Vertices,3),UBOUND(Vertices,3)
!       IF (ASSOCIATED(Vertices(i,j,k)%Vertex)) THEN
!         Vertex=>Vertices(i,j,k)%Vertex
!         IF (Vertex%in_out<=0) THEN
!         WRITE(300,*) '   i,j,k',i,j,k
!         WRITE(300,*) 'ix,iy,iz',Vertex%ix,Vertex%iy,Vertex%iz
!         WRITE(300,*) '  in_out',Vertex%in_out
!         END IF
!       END IF  
!     END DO
!   END DO
! END DO
!END SUBROUTINE WriteVerticesNachbar

SUBROUTINE WriteWeightBlk(FileName)
  CHARACTER*50 :: FileName

! -----------------------------
  INTEGER :: nl,k,j,i,izv,widthx,widthy,widthz
  INTEGER :: ib,in,posb
  INTEGER :: nx_ges,ny_ges,nz_ges
  INTEGER :: in_out
  INTEGER :: yes_v,no_v
  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  REAL(RealKind) :: zH(domain%igx0+1:domain%igx1 &
                      ,domain%igy0+1:domain%igy1 &
                      ,domain%igz0+1:domain%igz1)
!------------------------------------------------------------------------------
  zRauh=0.1d0
  !Albedo=0.5
  Albedo=0.08
  Emissivity=0.8
  yes_v=1
  no_v=0

  OPEN(UNIT=10,FILE=TRIM(FileName)//'.Weight2',STATUS='unknown')
  CALL WriteAuswWeightToProt 

  IF((indata_type=='gk'.AND.conv_gk=='s').OR. &
     (indata_type=='gk'.AND.OutGrid=="Globe") .OR. &
      indata_type=='geo' .OR. indata_type=='rad' .OR. &
      indata_type=='wgs84') THEN
    WRITE(10,*) 'Globe'
  ELSE IF (indata_type=='cyl') THEN
    WRITE(10,*) 'Cyl' 
    !WRITE(10,*) OutGrid(1:3)
  ELSE
    If(OutGrid(1:4)/='Cart') THEN
       WRITE(*,*) "Output Grid stimmt nicht mit definition Cart überein!"
    END IF
    WRITE(10,*) 'Cart'
  END IF

  WRITE(10,*) nb,"  ! Number of blocks"  
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO in=1,AnzahlNachbar
      IF (Nachbars(in)%nType(2:2)=="b") THEN
        posb=in
      END IF
    END DO
    WRITE(10,*)  "#.... Block",ib,"....................................................."
    WRITE(10,*)  xP(ix0),yP(iy0),zP(iz0)
    CALL WriteBlkNrToProt(ib)

    WRITE(10,*) nx
    DO i=ix0+1,ix1
       WRITE(10,*) dx(i)
    END DO

    WRITE(10,*) ny
    DO j=iy0+1,iy1
       WRITE(10,*) dy(j)
    END DO

    WRITE(10,*) nz
    DO k=iz0+1,iz1
       WRITE(10,*) dz(k)
    END DO

    FU=-1.d40
    FV=-1.d40
    FW=-1.d40
    VolC=-1.d40
    NrRW_FacesYZ=0
    NrRW_FacesZX=0
    NrRW_FacesXY=0
    CALL CountWeightFacesYZ(ix0,ix1,iy0,iy1,iz0,iz1,Vertices,Faces_YZ &
                           ,FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1))
    CALL CountWeightFacesZX(ix0,ix1,iy0,iy1,iz0,iz1,Vertices,Faces_ZX &
                           ,FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1))
    CALL CountWeightFacesXY(ix0,ix1,iy0,iy1,iz0,iz1,Vertices,Faces_XY &
                           ,FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1))
    CALL CountWeightCells(ix0,ix1,iy0,iy1,iz0,iz1,Vertices,Cells &
                         ,VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
    DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      IF (nType=='iw') THEN
        CALL CountWeightFacesYZ(jx0,jx0,jy0,jy1,jz0,jz1 &
                               ,VerticesNachbar(jx0:jx0,jy0:jy1,jz0:jz1) &
                               ,Faces_YZNachbar(jx0:jx0,jy0+1:jy1,jz0+1:jz1) &
                               ,FU(ix0-1:ix0-1,jy0+1:jy1,jz0+1:jz1))
        CALL CountWeightCells(jx0,jx1,jy0,jy1,jz0,jz1,VerticesNachbar,CellsNachbar &
                             ,VolC(ix0:ix0,jy0+1:jy1,jz0+1:jz1))
      ELSE IF (nType=='ie') THEN
        CALL CountWeightFacesYZ(jx1,jx1,jy0,jy1,jz0,jz1 &
                               ,VerticesNachbar(jx1:jx1,jy0:jy1,jz0:jz1) &
                               ,Faces_YZNachbar(jx1:jx1,jy0+1:jy1,jz0+1:jz1) &
                               ,FU(ix1+1:ix1+1,jy0+1:jy1,jz0+1:jz1))
        CALL CountWeightCells(jx0,jx1,jy0,jy1,jz0,jz1,VerticesNachbar,CellsNachbar &
                             ,VolC(ix1+1:ix1+1,jy0+1:jy1,jz0+1:jz1))
      ELSE IF (nType=='is') THEN
        CALL CountWeightFacesZX(jx0,jx1,jy0,jy0,jz0,jz1 &
                               ,VerticesNachbar(jx0:jx1,jy0:jy0,jz0:jz1) &
                               ,Faces_ZXNachbar(jx0+1:jx1,jy0:jy0,jz0+1:jz1) &
                               ,FV(jx0+1:jx1,iy0-1:iy0-1,jz0+1:jz1))
        CALL CountWeightCells(jx0,jx1,jy0,jy1,jz0,jz1,VerticesNachbar,CellsNachbar &
                             ,VolC(jx0+1:jx1,iy0:iy0,jz0+1:jz1))
      ELSE IF (nType=='in') THEN
        CALL CountWeightFacesZX(jx0,jx1,jy1,jy1,jz0,jz1 &
                               ,VerticesNachbar(jx0:jx1,jy1:jy1,jz0:jz1) &
                               ,Faces_ZXNachbar(jx0+1:jx1,jy1:jy1,jz0+1:jz1) &
                               ,FV(jx0+1:jx1,iy1+1:iy1+1,jz0+1:jz1))
        CALL CountWeightCells(jx0,jx1,jy0,jy1,jz0,jz1,VerticesNachbar,CellsNachbar &
                             ,VolC(jx0+1:jx1,iy1+1:iy1+1,jz0+1:jz1))
      ELSE IF (nType=='ib') THEN
        CALL CountWeightFacesXY(jx0,jx1,jy0,jy1,jz0,jz0 &
                               ,VerticesNachbar(jx0:jx1,jy0:jy1,jz0:jz0) &
                               ,Faces_XYNachbar(jx0+1:jx1,jy0+1:jy1,jz0:jz0) &
                               ,FW(jx0+1:jx1,jy0+1:jy1,iz0-1:iz0-1))
        CALL CountWeightCells(jx0,jx1,jy0,jy1,jz0,jz1-1,VerticesNachbar,CellsNachbar &
                             ,VolC(jx0+1:jx1,jy0+1:jy1,iz0:iz0))
      ELSE IF (nType=='it') THEN
        CALL CountWeightFacesXY(jx0,jx1,jy0,jy0,jz1,jz1 &
                               ,VerticesNachbar(jx0:jx1,jy0:jy1,jz1:jz1) &
                               ,Faces_XYNachbar(jx0+1:jx1,jy0+1:jy1,jz1:jz1) &
                               ,FW(jx0+1:jx1,jy0+1:jy1,iz1+1:iz1+1))
        CALL CountWeightCells(jx0,jx1,jy0,jy1,jz0+1,jz1,VerticesNachbar,CellsNachbar &
                             ,VolC(jx0+1:jx1,jy0+1:jy1,iz1+1:iz1+1))
      ELSE IF (nType=='ow') THEN
        CALL WeightNachbar(FU(ix0,iy0+1:iy1,iz0+1:iz1),FU(jx0-1,jy0+1:jy1,jz0+1:jz1),NrRW_FacesYZ)       
        CALL WeightNachbar(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1),VolC(ix0,iy0+1:iy1,iz0+1:iz1),NrRW_Cells)       
      ELSE IF (nType=='oe') THEN
        CALL WeightNachbar(FU(ix1,iy0+1:iy1,iz0+1:iz1),FU(jx1+1,jy0+1:jy1,jz0+1:jz1),NrRW_FacesYZ)       
        CALL WeightNachbar(VolC(ix1,iy0+1:iy1,iz0+1:iz1),VolC(ix1+1,iy0+1:iy1,iz0+1:iz1),NrRW_Cells)       
      ELSE IF (nType=='os') THEN
        CALL WeightNachbar(FV(ix0+1:ix1,iy0,iz0+1:iz1),FV(jx0+1:jx1,jy0-1,jz0+1:jz1),NrRW_FacesZX)       
        CALL WeightNachbar(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1),VolC(ix0+1:ix1,iy0,iz0+1:iz1),NrRW_Cells)       
      ELSE IF (nType=='on') THEN
        CALL WeightNachbar(FV(ix0+1:ix1,iy1,iz0+1:iz1),FV(jx0+1:jx1,jy1+1,jz0+1:jz1),NrRW_FacesZX)       
        CALL WeightNachbar(VolC(ix0+1:ix1,iy1,iz0+1:iz1),VolC(ix0+1:ix1,iy1+1,iz0+1:iz1),NrRW_Cells)       
      ELSE IF (nType=='ob') THEN
        CALL WeightNachbar(FW(ix0+1:ix1,iy0+1:iy1,iz0),FW(jx0+1:jx1,jy0+1:jy1,jz0-1),NrRW_FacesXY)       
        CALL WeightNachbar(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1),VolC(ix0+1:ix1,iy0+1:iy1,iz0),NrRW_Cells)       
      ELSE IF (nType=='ot') THEN
        CALL WeightNachbar(FW(ix0+1:ix1,iy0+1:iy1,iz1),FW(jx0+1:jx1,jy0+1:jy1,jz1+1),NrRW_FacesXY)       
        CALL WeightNachbar(VolC(ix0+1:ix1,iy0+1:iy1,iz1),VolC(ix0+1:ix1,iy0+1:iy1,iz1+1),NrRW_Cells)       
      END IF  
    END DO
    WRITE(10,*) NrRW_FacesYZ,'   NrW_FacesYZ'
    CALL WriteWeight(ix0-1,ix1+1,iy0+1,iy1,iz0+1,iz1,FU)
    CALL WriteMidPointFaces(ix0,ix1,iy0+1,iy1,iz0+1,iz1 &
                           ,Faces_YZ,'   NrMP_FacesYZ')
    WRITE(10,*) NrRW_FacesZX,'   NrW_FacesZX'
    CALL WriteWeight(ix0+1,ix1,iy0-1,iy1+1,iz0+1,iz1,FV)
    CALL WriteMidPointFaces(ix0+1,ix1,iy0,iy1,iz0+1,iz1 &
                           ,Faces_ZX,'   NrMP_FacesZX')
    WRITE(10,*) NrRW_FacesXY,'   NrW_FacesXY'
    CALL WriteWeight(ix0+1,ix1,iy0+1,iy1,iz0-1,iz1+1,FW)
    CALL WriteMidPointFaces(ix0+1,ix1,iy0+1,iy1,iz0,iz1 &
                           ,Faces_XY,'   NrMP_FacesXY')
    WRITE(10,*) NrRW_Cells, '    NrW_Cells'
    CALL WriteWeight(ix0,ix1+1,iy0,iy1+1,iz0,iz1+1,VolC)
!   WRITE(11,*) NrRW_Cells, '    NrW_Cells'
!   CALL WriteWeightAll(ix0,ix1+1,iy0,iy1+1,iz0,iz1+1,VolC)
    CALL WriteCutCells(ix0+1,ix1,iy0+1,iy1,iz0+1,iz1 &
                      ,Cells,'     NrMP_Cells')
  END DO
  CLOSE(10)

END SUBROUTINE WriteWeightBlk


SUBROUTINE OutWeightFaceYZ(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Floor(ib)%Faces_YZ(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. &
          Face%Vol==0.0d0) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) Face%Vol
     END IF
   ELSE
        in_out=Vertices(i,j-1,k-1)%Vertex%in_out &
              +Vertices(i,j,k-1)%Vertex%in_out &
              +Vertices(i,j-1,k)%Vertex%in_out &
              +Vertices(i,j,k)%Vertex%in_out
        IF (in_out<=0) THEN
          WRITE(10,*) i,j,k
          WRITE(10,*) 0.0d0
        END IF
   END IF
END SUBROUTINE OutWeightFaceYZ
SUBROUTINE OutWeightFaceZX(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out

   IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. &
          Face%Vol==0.0d0) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) Face%Vol
     END IF
   ELSE
        in_out=Vertices(i-1,j,k-1)%Vertex%in_out &
              +Vertices(i,j,k-1)%Vertex%in_out &
              +Vertices(i-1,j,k)%Vertex%in_out &
              +Vertices(i,j,k)%Vertex%in_out
        IF (in_out<=0) THEN
          WRITE(10,*) i,j,k
          WRITE(10,*) 0.0d0
        END IF
   END IF
END SUBROUTINE OutWeightFaceZX
SUBROUTINE OutWeightFaceXY(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. & 
          Face%Vol==0.0d0) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) Face%Vol
     END IF
   ELSE
     in_out=Vertices(i-1,j-1,k)%Vertex%in_out &
           +Vertices(i,j-1,k)%Vertex%in_out &
           +Vertices(i-1,j,k)%Vertex%in_out &
           +Vertices(i,j,k)%Vertex%in_out
     IF (in_out<=0) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) 0.0d0
     END IF
   END IF
END SUBROUTINE OutWeightFaceXY

SUBROUTINE CountWeightFaceXYBotton(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. & 
          Face%Vol==0.0d0) THEN
       NrRW_FacesXY=NrRW_FacesXY+1
     END IF
   !ELSE
   !   in_out=Vertices(i-1,j-1,k)%Vertex%in_out &
   !         +Vertices(i,j-1,k)%Vertex%in_out &
   !         +Vertices(i-1,j,k)%Vertex%in_out &
   !         +Vertices(i,j,k)%Vertex%in_out
   !   IF (in_out<=0) THEN
   !     NrRW_FacesXY=NrRW_FacesXY+1
   !   END IF
   END IF
END SUBROUTINE CountWeightFaceXYBotton

SUBROUTINE OutWeightFaceXYBotton(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. & 
          Face%Vol==0.0d0) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) Face%Vol
     END IF
   !ELSE
   !  in_out=Vertices(i-1,j-1,k)%Vertex%in_out &
   !        +Vertices(i,j-1,k)%Vertex%in_out &
   !        +Vertices(i-1,j,k)%Vertex%in_out &
   !        +Vertices(i,j,k)%Vertex%in_out
   !  IF (in_out<=0) THEN
   !    WRITE(10,*) i,j,k
   !    WRITE(10,*) 0.0d0
   !  END IF
   END IF
END SUBROUTINE OutWeightFaceXYBotton
END MODULE OutputWeightBlk_Mod
