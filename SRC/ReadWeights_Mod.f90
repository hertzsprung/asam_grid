MODULE ReadWeights_Mod
  USE Kind_Mod
  USE Floor_Mod
  USE Control_Mod
  USE Physics_Mod
  USE Parameter_Mod
  USE Parallel_Mod
  USE SoilData_Mod
  USE Canopy_Mod

  IMPLICIT NONE

  TYPE(PressureVelocity), POINTER :: VolFace(:)

CONTAINS 

SUBROUTINE CompareEW(F1,F2)

  TYPE(Domain_T) :: F1,F2

  INTEGER :: iy,iz

  DO iy=F1%iy0+1,F1%iy1
    DO iz=F1%iz0+1,F1%iz1
      IF (F1%WeiFU(F1%ix1-1,iy,iz)/=F2%WeiFU(F2%ix0-1,iy,iz))THEN
        WRITE(*,*) 'I',iy
        WRITE(*,*) F1%WeiFU(F1%ix1-1,iy,iz),F2%WeiFU(F2%ix0-1,iy,iz)
      END IF  
      IF (F1%WeiFU(F1%ix1,iy,iz)/=F2%WeiFU(F2%ix0,iy,iz))THEN
        WRITE(*,*) 'R',iy
        WRITE(*,*) F1%WeiFU(F1%ix1,iy,iz),F2%WeiFU(F2%ix0,iy,iz)
      END IF  
      IF (F1%WeiFU(F1%ix1+1,iy,iz)/=F2%WeiFU(F2%ix0+1,iy,iz))THEN
        WRITE(*,*) 'O',iy
        WRITE(*,*) F1%WeiFU(F1%ix1+1,iy,iz),F2%WeiFU(F2%ix0+1,iy,iz)
      END IF  
    END DO
  END DO
  WRITE(*,*) 'CompareEW Vol',F1%ib,F2%ib
  DO iy=F1%iy0+1,F1%iy1
    IF (F1%VolC(F1%ix1,iy,F1%iz0+1)/=F2%VolC(F2%ix0,iy,F1%iz0+1))THEN
      WRITE(*,*) 'I',iy
      WRITE(*,*) F1%VolC(F1%ix1,iy,F1%iz0+1),F2%VolC(F2%ix0,iy,F1%iz0+1)
    END IF  
    IF (F1%VolC(F1%ix1+1,iy,F1%iz0+1)/=F2%VolC(F2%ix0+1,iy,F1%iz0+1))THEN
      WRITE(*,*) 'O',iy
      WRITE(*,*) F1%VolC(F1%ix1+1,iy,F1%iz0+1),F2%VolC(F2%ix0+1,iy,F1%iz0+1)
    END IF  
  END DO
END SUBROUTINE CompareEW

SUBROUTINE TestWeight

  REAL(RealKind) :: VL(nb,nb)
  REAL(RealKind) :: VR(nb,nb)
  REAL(RealKind) :: FUL(nb,nb)
  REAL(RealKind) :: FUR(nb,nb)
  REAL(RealKind) :: FVL(nb,nb)
  REAL(RealKind) :: FVR(nb,nb)
  REAL(RealKind) :: FWL(nb,nb)
  REAL(RealKind) :: FWR(nb,nb)

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: ib,ibLoc,in

  VL=0.0d0
  VR=0.0d0
  FUL=0.0d0
  FUR=0.0d0
  FVL=0.0d0
  FVR=0.0d0
  FWL=0.0d0
  FWR=0.0d0
  DO ib=1,nb
    CALL Set(Floor(ib))
    IF (blMPI(ib)%Proc==MyId) THEN
      CALL Set(Floor(ib)) 
      DO in=1,AnzahlNachbar
        CALL Set(Nachbars(in))
        IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              VR(ib,ibn)=VR(ib,ibn)+VolC(ix0+1,jy,jz)
              VL(ib,ibn)=VL(ib,ibn)+VolC(ix0,jy,jz)
              FUR(ib,ibn)=FUR(ib,ibn)+FU(ix0+1,jy,jz)
              FUL(ib,ibn)=FUL(ib,ibn)+FU(ix0,jy,jz)
            END DO  
          END DO
        END IF
        IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              VL(ib,ibn)=VL(ib,ibn)+VolC(ix1,jy,jz)
              VR(ib,ibn)=VR(ib,ibn)+VolC(ix1+1,jy,jz)
              FUL(ib,ibn)=FUL(ib,ibn)+FU(ix1,jy,jz)
              FUR(ib,ibn)=FUR(ib,ibn)+FU(ix1+1,jy,jz)
            END DO  
          END DO
        END IF
        IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              VR(ib,ibn)=VR(ib,ibn)+VolC(jx,iy0+1,jz)
              VL(ib,ibn)=VL(ib,ibn)+VolC(jx,iy0,jz)
              FVR(ib,ibn)=FVR(ib,ibn)+FV(jx,iy0+1,jz)
              FVL(ib,ibn)=FVL(ib,ibn)+FV(jx,iy0,jz)
            END DO  
          END DO
        END IF
        IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              VL(ib,ibn)=VL(ib,ibn)+VolC(jx,iy1,jz)
              VR(ib,ibn)=VR(ib,ibn)+VolC(jx,iy1+1,jz)
              FVL(ib,ibn)=FVL(ib,ibn)+FV(jx,iy1,jz)
              FVR(ib,ibn)=FVR(ib,ibn)+FV(jx,iy1+1,jz)
            END DO  
          END DO
        END IF
        IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              VR(ib,ibn)=VR(ib,ibn)+VolC(jx,jy,iz0+1)
              VL(ib,ibn)=VL(ib,ibn)+VolC(jx,jy,iz0)
              FWR(ib,ibn)=FWR(ib,ibn)+FW(jx,jy,iz0+1)
              FWL(ib,ibn)=FWL(ib,ibn)+FW(jx,jy,iz0)
            END DO  
          END DO
        END IF
        IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              VL(ib,ibn)=VL(ib,ibn)+VolC(jx,jy,iz1)
              VR(ib,ibn)=VR(ib,ibn)+VolC(jx,jy,iz1+1)
              FWL(ib,ibn)=FWL(ib,ibn)+FW(jx,jy,iz1)
              FWR(ib,ibn)=FWR(ib,ibn)+FW(jx,jy,iz1+1)
            END DO  
          END DO
        END IF
      END DO
    END IF
    CALL MPI_Bcast(VL(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(VR(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FUL(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FUR(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FVL(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FVR(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FWL(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FWR(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
  END DO
  IF (MyId==0) THEN
    WRITE(*,*) 'VL'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (VL(ib,ibn)/=VL(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,VL(ib,ibn),VL(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'VR'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (VR(ib,ibn)/=VR(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,VR(ib,ibn),VR(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FUL'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FUL(ib,ibn)/=FUL(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FUL(ib,ibn),FUL(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FUR'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FUR(ib,ibn)/=FUR(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FUR(ib,ibn),FUR(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FVL'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FVL(ib,ibn)/=FVL(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FVL(ib,ibn),FVL(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FVR'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FVR(ib,ibn)/=FVR(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FVR(ib,ibn),FVR(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FWL'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FWL(ib,ibn)/=FWL(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FWL(ib,ibn),FWL(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FWR'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FWR(ib,ibn)/=FWR(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FWR(ib,ibn),FWR(ibn,ib)
        END IF  
      END DO  
    END DO  
  END IF  

END SUBROUTINE TestWeight

SUBROUTINE TestWeight2

  REAL(RealKind) :: VL(nb,nb)
  REAL(RealKind) :: VR(nb,nb)
  REAL(RealKind) :: FUL(nb,nb)
  REAL(RealKind) :: FUR(nb,nb)
  REAL(RealKind) :: FVL(nb,nb)
  REAL(RealKind) :: FVR(nb,nb)
  REAL(RealKind) :: FWL(nb,nb)
  REAL(RealKind) :: FWR(nb,nb)

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: ib,ibLoc,in

  VL=0.0d0
  VR=0.0d0
  FUL=0.0d0
  FUR=0.0d0
  FVL=0.0d0
  FVR=0.0d0
  FWL=0.0d0
  FWR=0.0d0
  DO ib=1,nb
    CALL Set(Floor(ib))
    IF (blMPI(ib)%Proc==MyId) THEN
      CALL Set(Floor(ib)) 
      DO in=1,AnzahlNachbar
        CALL Set(Nachbars(in))
        IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              VR(ib,ibn)=VR(ib,ibn)+VolC(ix0,jy,jz)
              VL(ib,ibn)=VL(ib,ibn)+VolC(ix0-1,jy,jz)
              FUR(ib,ibn)=FUR(ib,ibn)+FU(ix0,jy,jz)
              FUL(ib,ibn)=FUL(ib,ibn)+FU(ix0-1,jy,jz)
            END DO  
          END DO
        END IF
        IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              VL(ib,ibn)=VL(ib,ibn)+VolC(ix1-1,jy,jz)
              VR(ib,ibn)=VR(ib,ibn)+VolC(ix1,jy,jz)
              FUL(ib,ibn)=FUL(ib,ibn)+FU(ix1-1,jy,jz)
              FUR(ib,ibn)=FUR(ib,ibn)+FU(ix1,jy,jz)
            END DO  
          END DO
        END IF
        IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
!             VR(ib,ibn)=VR(ib,ibn)+VolC(jx,iy0,jz)
!             VL(ib,ibn)=VL(ib,ibn)+VolC(jx,iy0-1,jz)
              FVR(ib,ibn)=FVR(ib,ibn)+FV(jx,iy0,jz)
              FVL(ib,ibn)=FVL(ib,ibn)+FV(jx,iy0-1,jz)
            END DO  
          END DO
        END IF
        IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
!             VL(ib,ibn)=VL(ib,ibn)+VolC(jx,iy1-1,jz)
!             VR(ib,ibn)=VR(ib,ibn)+VolC(jx,iy1,jz)
              FVL(ib,ibn)=FVL(ib,ibn)+FV(jx,iy1-1,jz)
              FVR(ib,ibn)=FVR(ib,ibn)+FV(jx,iy1,jz)
            END DO  
          END DO
        END IF
        IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
!             VR(ib,ibn)=VR(ib,ibn)+VolC(jx,jy,iz0+1)
!             VL(ib,ibn)=VL(ib,ibn)+VolC(jx,jy,iz0)
              FWR(ib,ibn)=FWR(ib,ibn)+FW(jx,jy,iz0)
              FWL(ib,ibn)=FWL(ib,ibn)+FW(jx,jy,iz0-1)
            END DO  
          END DO
        END IF
        IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
!             VL(ib,ibn)=VL(ib,ibn)+VolC(jx,jy,iz1)
!             VR(ib,ibn)=VR(ib,ibn)+VolC(jx,jy,iz1+1)
              FWL(ib,ibn)=FWL(ib,ibn)+FW(jx,jy,iz1-1)
              FWR(ib,ibn)=FWR(ib,ibn)+FW(jx,jy,iz1)
            END DO  
          END DO
        END IF
      END DO
    END IF
    CALL MPI_Bcast(VL(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(VR(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FUL(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FUR(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FVL(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FVR(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FWL(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
    CALL MPI_Bcast(FWR(ib,1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
&                 MPI_COMM_WORLD, MPIErr) 
  END DO
  IF (MyId==0) THEN
    WRITE(*,*) 'VL'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (VL(ib,ibn)/=VL(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,VL(ib,ibn),VL(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'VR'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (VR(ib,ibn)/=VR(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,VR(ib,ibn),VR(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FUL'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FUL(ib,ibn)/=FUL(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FUL(ib,ibn),FUL(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FUR'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FUR(ib,ibn)/=FUR(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FUR(ib,ibn),FUR(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FVL'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FVL(ib,ibn)/=FVL(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FVL(ib,ibn),FVL(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FVR'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FVR(ib,ibn)/=FVR(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FVR(ib,ibn),FVR(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FWL'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FWL(ib,ibn)/=FWL(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FWL(ib,ibn),FWL(ibn,ib)
        END IF  
      END DO  
    END DO  
    WRITE(*,*) 'FWR'
    DO ib=1,nb
      DO ibn=ib,nb
        IF (FWR(ib,ibn)/=FWR(ibn,ib)) THEN
          WRITE(*,*) ib,ibn,FWR(ib,ibn),FWR(ibn,ib)
        END IF  
      END DO  
    END DO  
  END IF  

END SUBROUTINE TestWeight2


SUBROUTINE ReadWeights(FileName)

  CHARACTER(300) :: Line
  CHARACTER(*) :: FileName
  CHARACTER(80) :: mi_wgts2
  CHARACTER(16) :: Str_Read
  LOGICAL :: Cartesian,Globe,CylinderT

  INTEGER :: i,ib,ibLoc,ibInput,in,j,k
  INTEGER :: nLines,nxx,nyy,nzz
  INTEGER :: nxWeight,nyWeight,nzWeight 
  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: ixref, iyref, izref
  INTEGER :: jxref, jyref, jzref
  INTEGER :: levX,levY,levZ,shX,shY,shZ
  INTEGER :: levNX,levNY,levNZ,shNX,shNY,shNZ
  INTEGER :: icut 
  INTEGER :: nr_structbound
  INTEGER :: nLocBoundary
  INTEGER, ALLOCATABLE :: BoundaryBlk(:,:)
  

  REAL(RealKind) :: n1,n2,n3,FL
  REAL(RealKind) :: n1G,n2G,n3G
  REAL(RealKind) :: n1G1,n2G1,n3G1
  REAL(RealKind) :: n1G2,n2G2,n3G2
  REAL(RealKind) :: phiL,phiR
  REAL(RealKind) :: lamL,lamR
  REAL(RealKind) :: CosLam,SinPhi,CosPhi
  REAL(RealKind) :: MaxHeightLoc
  REAL(RealKind) :: xFLSG,yFLSG,zFLSG
  REAL(RealKind) :: xFLSGmin
  REAL(RealKind) :: dxLoc,dyLoc,dzLoc
  REAL(RealKind) :: xw0,xw1
  REAL(RealKind) :: yw0,yw1
  REAL(RealKind) :: zw0,zw1

  ! Domain-Part-------------------------------------------------------------------------
  ALLOCATE(domain%dx(domain%ix0+1:domain%ix1))
  ALLOCATE(domain%dy(domain%iy0+1:domain%iy1))
  ALLOCATE(domain%dz(domain%iz0+1:domain%iz1))
  ALLOCATE(domain%xP(domain%ix0:domain%ix1))
  ALLOCATE(domain%yP(domain%iy0:domain%iy1))
  ALLOCATE(domain%zP(domain%iz0:domain%iz1))

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'#OutputDomain')>0) THEN
      READ(InputUnit,*) Line
      IF (Line=='XYZ_Number') THEN 
        READ(InputUnit,*) ix0Out,ix1Out,iy0Out,iy1Out,iz0Out,iz1Out
      ELSE IF (Line=='XYZ_Coord') THEN 
      END IF  
      EXIT
    END IF
  END DO
1 CONTINUE
  REWIND(InputUnit)
  DO 
    READ(InputUnit,*,END=2) Line
    !.......................................................................
    IF (INDEX(Line,'#xGrid')>0) THEN   ! set different distance x-direction
      READ(InputUnit,*) nlines
      nxx=0
      domain%xP(domain%ix0)=domain%x0
      DO i=1,nlines
        READ(InputUnit,*) xw0,xw1,domain%nx
        DO ix=nxx+1,nxx+domain%nx
          domain%dx(ix)=(xw1-xw0)/domain%nx
          domain%xP(ix)=domain%xP(ix-1)+domain%dx(ix)
        END DO
        nxx=nxx+domain%nx
      END DO
      domain%nx=nxx
      EXIT
    ELSE
      domain%xP(domain%ix0)=domain%x0
      DO ix=domain%ix0+1,domain%ix1
        domain%dx(ix)=(domain%x1-domain%x0)/domain%nx
        domain%xP(ix)=domain%xP(ix-1)+domain%dx(ix)
      END DO
    END IF
  END DO
2 CONTINUE
  REWIND(InputUnit)
  DO 
    READ(InputUnit,*,END=3) Line
    !.......................................................................
    IF (INDEX(Line,'#yGrid')>0) THEN   ! set different distance y-direction
      READ(InputUnit,*) nlines
      nyy=0
      domain%yP(domain%iy0)=domain%y0
      DO j=1,nlines
        READ(InputUnit,*) yw0,yw1,domain%ny
        DO iy=nyy+1,nyy+domain%ny
          domain%dy(iy)=(yw1-yw0)/domain%ny
          domain%yP(iy)=domain%yP(iy-1)+domain%dy(iy)
        END DO
        nyy=nyy+domain%ny
      END DO
      domain%ny=nyy
      EXIT
    ELSE
      domain%yP(domain%iy0)=domain%y0
      DO iy=domain%iy0+1,domain%iy1
        domain%dy(iy)=(domain%y1-domain%y0)/domain%ny
        domain%yP(iy)=domain%yP(iy-1)+domain%dy(iy)
      END DO
    END IF
  END DO
3 CONTINUE
  REWIND(InputUnit)
  DO 
    READ(InputUnit,*,END=4) Line
    !........................................................................
    IF (INDEX(Line,'#zGrid')>0) THEN   ! set different distance z-direction
      READ(InputUnit,*) nlines
      nzz=0
      domain%zP(domain%iz0)=domain%z0
      DO k=1,nlines
        READ(InputUnit,*) zw0,zw1,domain%nz
        DO iz=nzz+1,nzz+domain%nz
          domain%dz(iz)=(zw1-zw0)/domain%nz
          domain%zP(iz)=domain%zP(iz-1)+domain%dz(iz)
        END DO
        nzz=nzz+domain%nz
      END DO
      domain%nz=nzz
      EXIT
    ELSE IF (INDEX(Line,'#DCMIPzGrid')>0) THEN
      domain%zP(domain%iz0)=domain%z0
      DO iz=domain%iz0+1,domain%iz1
        zw0=Domain%z1*(SQRT(15.d0*((iz-1)/30.0d0)**2.0d0+1.d0)-1)/(SQRT(15.d0+1)-1)
        zw1=Domain%z1*(SQRT(15.d0*((iz)/30.0d0)**2.0d0+1.d0)-1)/(SQRT(15.d0+1)-1)
        domain%dz(iz)=(zw1-zw0)
        domain%zP(iz)=domain%zP(iz-1)+domain%dz(iz)
      END DO
      EXIT
    ELSE
      domain%zP(domain%iz0)=domain%z0
      DO iz=domain%iz0+1,domain%iz1
        domain%dz(iz)=(domain%z1-domain%z0)/domain%nz
        domain%zP(iz)=domain%zP(iz-1)+domain%dz(iz)
      END DO
    END IF
    !........................................................................
  END DO
4 CONTINUE
  CLOSE(UNIT=InputUnit)

  mi_wgts2 = TRIM(FileName(1:INDEX(FileName,'.grid')-1))//'.Weight2'
  OPEN(UNIT=InputUnit,FILE=mi_wgts2,STATUS='old')
  READ(InputUnit,*) GridType
  GridType=TRIM(GridType)
  Cartesian=GridType=='Cart'
  Globe=GridType=='Globe'
  CylinderT=GridType=='Cyl'
  READ(InputUnit,*) nb ! Number of blocks
  MaxHeightLoc=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    S:DO
      READ(InputUnit,'(a300)') Line
      IF (INDEX(Line,'#.... Block')>0) THEN
        READ(Line(INDEX(Line,'#.... Block')+LEN('#.... Block'):),*) ibInput ! the indices of blocks
        IF(ib==ibInput) THEN 
          READ(InputUnit,*) Floor(ib)%x0,Floor(ib)%y0,Floor(ib)%z0 ! the indices of the westest, southest and lowest boundaries of domain
          !.......................................
          ! x-direction
          READ(InputUnit,*) nxWeight ! number of cells in x-direction, including ghost cells
          IF (nxWeight/=Floor(ib)%ix1-Floor(ib)%ix0) STOP 'Error WeightFile (nx)'
          Floor(ib)%xP(Floor(ib)%ix0)=Floor(ib)%x0
          DO ix=ix0+1,ix1
            READ(InputUnit,*) Floor(ib)%dx(ix)
            Floor(ib)%xP(ix)=Floor(ib)%xP(ix-1)+Floor(ib)%dx(ix)
          END DO
          Floor(ib)%x1=Floor(ib)%xP(ix1)
          !.......................................
          ! y-direction
          READ(InputUnit,*) nyWeight ! number of cells in y-direction, including ghost cells
          IF (nyWeight/=Floor(ib)%iy1-Floor(ib)%iy0) STOP 'Error WeightFile (ny)'
          Floor(ib)%yP(Floor(ib)%iy0)=Floor(ib)%y0
          DO iy=Floor(ib)%iy0+1,Floor(ib)%iy1
            READ(InputUnit,*) Floor(ib)%dy(iy)
            Floor(ib)%yP(iy)=Floor(ib)%yP(iy-1)+Floor(ib)%dy(iy)
          END DO
          Floor(ib)%y1=Floor(ib)%yP(iy1)
          !.......................................
          ! z-direction
          READ(InputUnit,*) nzWeight ! number of cells in z-direction, including ghost cells
          IF (nzWeight/=Floor(ib)%iz1-Floor(ib)%iz0) STOP 'Error WeightFile (nz)'
          Floor(ib)%zP(Floor(ib)%iz0)=Floor(ib)%z0
          DO iz=Floor(ib)%iz0+1,Floor(ib)%iz1
            READ(InputUnit,*) Floor(ib)%dz(iz)
            Floor(ib)%zP(iz)=Floor(ib)%zP(iz-1)+Floor(ib)%dz(iz)
          END DO
          Floor(ib)%z1=Floor(ib)%zP(iz1)
          DO ix=Floor(ib)%ix0,Floor(ib)%ix1
            DO iy=Floor(ib)%iy0,Floor(ib)%iy1
              DO iz=Floor(ib)%iz0,Floor(ib)%iz1
                IF (Sphere) THEN
                  Floor(ib)%xPG(ix,iy,iz)=GlobCartXP(Floor(ib)%zP(iz)+RadEarth,Floor(ib)%xP(ix),Floor(ib)%yP(iy))
                  Floor(ib)%yPG(ix,iy,iz)=GlobCartYP(Floor(ib)%zP(iz)+RadEarth,Floor(ib)%xP(ix),Floor(ib)%yP(iy))
                  Floor(ib)%zPG(ix,iy,iz)=GlobCartZP(Floor(ib)%zP(iz)+RadEarth,Floor(ib)%yP(iy))
                ELSE
                  Floor(ib)%xPG(ix,iy,iz)=Floor(ib)%xP(ix)
                  Floor(ib)%yPG(ix,iy,iz)=Floor(ib)%yP(iy)
                  Floor(ib)%zPG(ix,iy,iz)=Floor(ib)%zP(iz)
                END IF
              END DO
            END DO
          END DO
          !
          !.........................................................................
          !--- Weights at x-boundaries
          IF (Cartesian) THEN
            DO iz=iz0+1,iz1
              DO iy=iy0+1,iy1
                DO ix=ix0-1,ix1+1          ! with boundary coordinate
                  FU(ix,iy,iz)=dy(iy)*dz(iz)
                END DO
              END DO
            END DO
          ELSE IF (Globe) THEN
            DO iz=iz0+1,iz1
              DO iy=iy0+1,iy1
                DO ix=ix0-1,ix1+1          ! boundary coordinate
                  FU(ix,iy,iz)=FaceGlobeYZ(zP(iz-1)+RadEarth &
                                          ,zP(iz)+RadEarth   &
                                          ,yP(iy-1) &
                                          ,yP(iy))
                 END DO
              END DO
            END DO
          ELSE IF (CylinderT) THEN
            DO iz=iz0+1,iz1
              DO iy=iy0+1,iy1
                DO ix=ix0-1,ix1+1          ! boundary coordinate
                  FU(ix,iy,iz)=FaceCylinderYZ(zP(iz-1) &
                                             ,zP(iz)   &
                                             ,yP(iy-1) &
                                             ,yP(iy))
                 END DO
              END DO
            END DO
          END IF
          ! Weight-FacesYZ, Schnitte-/unterhalb Berg
!         FU=0.0d0 !OSSI
          READ(InputUnit,*) icut    !NrW_FacesYZ
          DO i=1,icut
            READ(InputUnit,*) ix,iy,iz
            READ(InputUnit,*) FU(ix,iy,iz)
          END DO

          ! MP_FacesYZ
          READ(InputUnit,*) Floor(ib)%NumBoundFaceU ! NrMP_FacesYZ
          ALLOCATE(Floor(ib)%BoundFaceU(Floor(ib)%NumBoundFaceU))
          DO i=1,Floor(ib)%NumBoundFaceU
             READ(InputUnit,*) Floor(ib)%BoundFaceU(i)%ix &
                              ,Floor(ib)%BoundFaceU(i)%iy &
                              ,Floor(ib)%BoundFaceU(i)%iz 
             READ(InputUnit,*) Floor(ib)%BoundFaceU(i)%xS &
                              ,Floor(ib)%BoundFaceU(i)%yS &
                              ,Floor(ib)%BoundFaceU(i)%zS 
          END DO
          !.........................................................................
          !--- Weights at y-boundaries
          IF (Cartesian) THEN
            DO iz=iz0+1,iz1
              DO iy=iy0-1,iy1+1               ! with boundary coordinate
                DO ix=ix0+1,ix1
                   FV(ix,iy,iz)=dx(ix)*dz(iz)
                END DO
              END DO
            END DO
          ELSE IF (Globe) THEN
            DO iz=iz0+1,iz1
              DO iy=iy0,iy1               ! boundary coordinate
                DO ix=ix0+1,ix1
                   FV(ix,iy,iz)=FaceGlobeZX(zP(iz-1)+RadEarth &
                                           ,zP(iz)+RadEarth &
                                           ,yP(iy) &
                                           ,xP(ix-1) &
                                           ,xP(ix))
                END DO
              END DO
            END DO
            DO in=1,AnzahlNachbar
              CALL Set(Nachbars(in))
              IF (Nachbars(in)%nType(2:2)=='s') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dyLoc=dy(iy0+1)
                ELSE  
                  IF (igy0==Domain%igy0) THEN
                    igy0=Domain%igy1
                  END IF
                  dyLoc=Domain%yP(igy0)-Domain%yP(igy0-2**(-RefineNachbarY))
                END IF
                DO jx=jx0+1,jx1
                  DO jz=jz0+1,jz1
                    FV(jx,iy0-1,jz)=FaceGlobeZX(zP(jz-1)+RadEarth &
                                           ,zP(jz)+RadEarth &
                                           ,yP(iy0)-dyLoc &
                                           ,xP(jx-1) &
                                           ,xP(jx))
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='n') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dyLoc=dy(iy1)
                ELSE  
                  IF (igy1==Domain%igy1) THEN
                    igy1=Domain%igy0
                  END IF  
                  dyLoc=Domain%yP(igy1+2**(-RefineNachbarY))-Domain%yP(igy1)
                END IF
                DO jx=jx0+1,jx1
                  DO jz=jz0+1,jz1
                    FV(jx,iy1+1,jz)=FaceGlobeZX(zP(jz-1)+RadEarth &
                                           ,zP(jz)+RadEarth &
                                           ,yP(iy1)+dyLoc &
                                           ,xP(jx-1) &
                                           ,xP(jx))
                  END DO
                END DO
              END IF  
            END DO  
          ELSE IF (CylinderT) THEN
            DO iz=iz0+1,iz1
              DO iy=iy0-1,iy1+1               ! boundary coordinate
                DO ix=ix0+1,ix1
                   FV(ix,iy,iz)=FaceCylinderZX(zP(iz-1) &
                                              ,zP(iz)   &
                                              ,yP(iy) &
                                              ,xP(ix-1) &
                                              ,xP(ix))
                END DO
              END DO
            END DO
          END IF 
          READ(InputUnit,*) icut     !NrW_FacesZX
          DO i=1,icut
            READ(InputUnit,*) ix,iy,iz
            READ(InputUnit,*) FV(ix,iy,iz)
          END DO
          IF (Globe) THEN
            IF (yP(iy0)<=-Half*Pi+1.d-10) THEN
              DO iz=iz0+1,iz1
                DO ix=ix0+1,ix1
                 FV(ix,iy0,iz)=Zero
                END DO
              END DO
              FV(:,iy0-1,:)=FV(:,iy0,:)
            END IF
            IF (yP(iy1)>=Half*Pi-1.d-10) THEN
              DO iz=iz0+1,iz1
                DO ix=ix0+1,ix1
                  FV(ix,iy1,iz)=Zero
                END DO
              END DO
              FV(:,iy1+1,:)=FV(:,iy1,:)
            END IF
          END IF

          ! MP_FacesZX
          READ(InputUnit,*) Floor(ib)%NumBoundFaceV ! NrMP_FacesZX
          ALLOCATE(Floor(ib)%BoundFaceV(Floor(ib)%NumBoundFaceV))
          DO i=1,Floor(ib)%NumBoundFaceV
            READ(InputUnit,*) Floor(ib)%BoundFaceV(i)%ix &
                             ,Floor(ib)%BoundFaceV(i)%iy &
                             ,Floor(ib)%BoundFaceV(i)%iz
            READ(InputUnit,*) Floor(ib)%BoundFaceV(i)%xS &
                             ,Floor(ib)%BoundFaceV(i)%yS &
                             ,Floor(ib)%BoundFaceV(i)%zS
          END DO

          !.........................................................................
          !--- Weights at z-boundaries
          IF (Cartesian) THEN
            DO iz=iz0-1,iz1+1                  ! with boundary coordinate
              DO iy=iy0+1,iy1
                DO ix=ix0+1,ix1
                  FW(ix,iy,iz)=dx(ix)*dy(iy)
                END DO
              END DO
            END DO
          ELSE IF (Globe) THEN
            DO iz=iz0,iz1                  ! boundary coordinate
              DO iy=iy0+1,iy1
                DO ix=ix0+1,ix1
                  FW(ix,iy,iz)=FaceGlobeXY(zP(iz)+RadEarth &
                                          ,yP(iy-1) &
                                          ,yP(iy) &
                                          ,xP(ix-1) &
                                          ,xP(ix))
                END DO
              END DO
            END DO
            DO in=1,AnzahlNachbar
              CALL Set(Nachbars(in))
              IF (Nachbars(in)%nType(2:2)=='b') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dzLoc=dz(iz0+1)
                ELSE  
                  IF (igz0==Domain%igz0) THEN
                    igz0=Domain%igz1
                  END IF
                  dzLoc=Domain%zP(igz0)-Domain%zP(igz0-2**(-RefineNachbarZ))
                END IF
                DO jx=jx0+1,jx1
                  DO jy=jy0+1,jy1
                    FW(jx,jy,iz0-1)=FaceGlobeXY(zP(iz0)+RadEarth-dzLoc &
                                          ,yP(jy-1) &
                                          ,yP(jy) &
                                          ,xP(jx-1) &
                                          ,xP(jx))
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='t') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dzLoc=dz(iz1)
                ELSE  
                  IF (igz1==Domain%igz1) THEN
                    igz1=Domain%igz0
                  END IF  
                  dzLoc=Domain%zP(igz1+2**(-RefineNachbarZ))-Domain%zP(igz1)
                END IF
                DO jx=jx0+1,jx1
                  DO jy=jy0+1,jy1
                    FW(jx,jy,iz1+1)=FaceGlobeXY(zP(iz1)+RadEarth+dzLoc &
                                          ,yP(jy-1) &
                                          ,yP(jy) &
                                          ,xP(jx-1) &
                                          ,xP(jx))
                  END DO
                END DO
              END IF  
            END DO
          ELSE IF (CylinderT) THEN
            DO iz=iz0-1,iz1+1                  ! boundary coordinate
              DO iy=iy0+1,iy1
                DO ix=ix0+1,ix1
                  FW(ix,iy,iz)=FaceCylinderXY(yP(iy-1) &
                                             ,yP(iy) &
                                             ,xP(ix-1) &
                                             ,xP(ix))
                END DO
              END DO
            END DO
          END IF
          READ(InputUnit,*) icut
          DO i=1,icut       !NrW_FacesXY
            READ(InputUnit,*) ix,iy,iz
            READ(InputUnit,*) FW(ix,iy,iz)
          END DO

          ! MP_FacesXY, Schnitte
          READ(InputUnit,*) Floor(ib)%NumBoundFaceW ! NrMP_FacesXY
          ALLOCATE(Floor(ib)%BoundFaceW(Floor(ib)%NumBoundFaceW))
          DO i=1,Floor(ib)%NumBoundFaceW
             READ(InputUnit,*) Floor(ib)%BoundFaceW(i)%ix &
                              ,Floor(ib)%BoundFaceW(i)%iy &
                              ,Floor(ib)%BoundFaceW(i)%iz
             READ(InputUnit,*) Floor(ib)%BoundFaceW(i)%xS &
                              ,Floor(ib)%BoundFaceW(i)%yS &
                              ,Floor(ib)%BoundFaceW(i)%zS
          END DO
          Floor(ib)%zH(:,:)=0.d0
          DO iy=iy0+1,iy1
             DO ix=ix0+1,ix1
                Floor(ib)%zH(ix,iy)=zP(iz0)
                DO iz=iz0+1,iz1
                   IF (FW(ix,iy,iz-1)<FW(ix,iy,iz)) THEN
                     Floor(ib)%zH(ix,iy)=zP(iz-1)
                     MaxHeightLoc=MAX(MaxHeightLoc,zP(iz-1))
                   END IF
                END DO
             END DO
          END DO

          !.........................................................................
          !--- Weights for cell volume

          IF (Cartesian) THEN
            VolC=1.0d0
            DO iz=iz0+1,iz1
              DO iy=iy0+1,iy1
                DO ix=ix0+1,ix1
                  VolC(ix,iy,iz)=dx(ix)*dy(iy)*dz(iz)
                END DO
              END DO
            END DO
            DO in=1,AnzahlNachbar
              CALL Set(Nachbars(in))
              IF (Nachbars(in)%nType(2:2)=='w') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dxLoc=dx(ix0+1)
                ELSE  
                  IF (igx0==Domain%igx0) THEN
                    igx0=Domain%igx1
                  END IF  
                  dxLoc=Domain%xP(igx0)-Domain%xP(igx0-2**(-RefineNachbarX))
                END IF
                Nachbars(in)%dLoc=dxLoc
                DO jy=jy0+1,jy1
                  DO jz=jz0+1,jz1
                    VolC(ix0,jy,jz)=dxLoc*dy(jy)*dz(jz)
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='e') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dxLoc=dx(ix1)
                ELSE  
                  IF (igx1==Domain%igx1) THEN
                    igx1=Domain%igx0
                  END IF
                  dxLoc=Domain%xP(igx1+2**(-RefineNachbarX))-Domain%xP(igx1)
                END IF
                Nachbars(in)%dLoc=dxLoc
                DO jy=jy0+1,jy1
                  DO jz=jz0+1,jz1
                    VolC(ix1+1,jy,jz)=dxLoc*dy(jy)*dz(jz)
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='s') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dyLoc=dy(iy0+1)
                ELSE  
                  IF (igy0==Domain%igy0) THEN
                    igy0=Domain%igy1
                  END IF
                  dyLoc=Domain%yP(igy0)-Domain%yP(igy0-2**(-RefineNachbarY))
                END IF
                Nachbars(in)%dLoc=dyLoc
                DO jx=jx0+1,jx1
                  DO jz=jz0+1,jz1
                    VolC(jx,iy0,jz)=dx(jx)*dyLoc*dz(jz)
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='n') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dyLoc=dy(iy1)
                ELSE  
                  IF (igy1==Domain%igy1) THEN
                    igy1=Domain%igy0
                  END IF  
                  dyLoc=Domain%yP(igy1+2**(-RefineNachbarY))-Domain%yP(igy1)
                END IF
                Nachbars(in)%dLoc=dyLoc
                DO jx=jx0+1,jx1
                  DO jz=jz0+1,jz1
                    VolC(jx,iy1+1,jz)=dx(jx)*dyLoc*dz(jz)
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='b') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dzLoc=dz(iz0+1)
                ELSE  
                  IF (igz0==Domain%igz0) THEN
                    igz0=Domain%igz1
                  END IF
                  dzLoc=Domain%zP(igz0)-Domain%zP(igz0-2**(-RefineNachbarZ))
                END IF
                Nachbars(in)%dLoc=dzLoc
                DO jx=jx0+1,jx1
                  DO jy=jy0+1,jy1
                    VolC(jx,jy,iz0)=dx(jx)*dy(jy)*dzLoc
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='t') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dzLoc=dz(iz1)
                ELSE  
                  IF (igz1==Domain%igz1) THEN
                    igz1=Domain%igz0
                  END IF  
                  dzLoc=Domain%zP(igz1+2**(-RefineNachbarZ))-Domain%zP(igz1)
                END IF
                Nachbars(in)%dLoc=dzLoc
                DO jx=jx0+1,jx1
                  DO jy=jy0+1,jy1
                    VolC(jx,jy,iz1+1)=dx(jx)*dy(jy)*dzLoc
                  END DO
                END DO
              END IF  
            END DO
          ELSE IF (Globe) THEN
            DO iz=iz0+1,iz1
              DO iy=iy0+1,iy1
                DO ix=ix0+1,ix1
                  VolC(ix,iy,iz)=VolCellGlobe(zP(iz-1)+RadEarth &
                                             ,zP(iz)+RadEarth &
                                             ,yP(iy-1) &
                                             ,yP(iy) &
                                             ,xP(ix-1) &
                                             ,xP(ix))
                END DO
              END DO
            END DO
            DO in=1,AnzahlNachbar
              CALL Set(Nachbars(in))
              IF (Nachbars(in)%nType(2:2)=='w') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dxLoc=dx(ix0+1)
                ELSE  
                  IF (igx0==Domain%igx0) THEN
                    igx0=Domain%igx1
                  END IF  
                  dxLoc=Domain%xP(igx0)-Domain%xP(igx0-2**(-RefineNachbarX))
                END IF
                Nachbars(in)%dLoc=dxLoc
                DO jy=jy0+1,jy1
                  DO jz=jz0+1,jz1
                    VolC(ix0,jy,jz)=VolCellGlobe(zP(jz-1)+RadEarth &
                                               ,zP(jz)+RadEarth &
                                               ,yP(jy-1) &
                                               ,yP(jy) &
                                               ,xP(ix0)-dxLoc &
                                               ,xP(ix0))
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='e') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dxLoc=dx(ix1)
                ELSE  
                  IF (igx1==Domain%igx1) THEN
                    igx1=Domain%igx0
                  END IF
                  dxLoc=Domain%xP(igx1+2**(-RefineNachbarX))-Domain%xP(igx1)
                END IF
                Nachbars(in)%dLoc=dxLoc
                DO jy=jy0+1,jy1
                  DO jz=jz0+1,jz1
                    VolC(ix1+1,jy,jz)=VolCellGlobe(zP(jz-1)+RadEarth &
                                                  ,zP(jz)+RadEarth &
                                                  ,yP(jy-1) &
                                                  ,yP(jy) &
                                               ,xP(ix1) &
                                               ,xP(ix1)+dxLoc)
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='s') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dyLoc=dy(iy0+1)
                ELSE  
                  IF (igy0==Domain%igy0) THEN
                    igy0=Domain%igy1
                  END IF
                  dyLoc=Domain%yP(igy0)-Domain%yP(igy0-2**(-RefineNachbarY))
                END IF
                Nachbars(in)%dLoc=dyLoc
                DO jx=jx0+1,jx1
                  DO jz=jz0+1,jz1
                    VolC(jx,iy0,jz)=VolCellGlobe(zP(jz-1)+RadEarth &
                                                ,zP(jz)+RadEarth &
                                                ,yP(iy0)-dyLoc &
                                                ,yP(iy0) &
                                                ,xP(jx-1) &
                                                ,xP(jx))
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='n') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dyLoc=dy(iy1)
                ELSE  
                  IF (igy1==Domain%igy1) THEN
                    igy1=Domain%igy0
                  END IF  
                  dyLoc=Domain%yP(igy1+2**(-RefineNachbarY))-Domain%yP(igy1)
                END IF
                Nachbars(in)%dLoc=dyLoc
                DO jx=jx0+1,jx1
                  DO jz=jz0+1,jz1
                    VolC(jx,iy1+1,jz)=VolCellGlobe(zP(jz-1)+RadEarth &
                                                  ,zP(jz)+RadEarth &
                                                  ,yP(iy1) &
                                                  ,yP(iy1)+dyLoc &
                                                  ,xP(jx-1) &
                                                  ,xP(jx))
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='b') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dzLoc=dz(iz0+1)
                ELSE  
                  IF (igz0==Domain%igz0) THEN
                    igz0=Domain%igz1
                  END IF
                  dzLoc=Domain%zP(igz0)-Domain%zP(igz0-2**(-RefineNachbarZ))
                END IF
                Nachbars(in)%dLoc=dzLoc
                DO jx=jx0+1,jx1
                  DO jy=jy0+1,jy1
                    VolC(jx,jy,iz0)=VolCellGlobe(zP(iz0)+RadEarth-dzLoc &
                                                ,zP(iz0)+RadEarth &
                                                ,yP(jy-1) &
                                                ,yP(jy) &
                                                ,xP(jx-1) &
                                                ,xP(jx))
                  END DO
                END DO
              ELSE IF (Nachbars(in)%nType(2:2)=='t') THEN
                IF (Nachbars(in)%nType(1:1)=='o') THEN
                  dzLoc=dz(iz1)
                ELSE  
                  IF (igz1==Domain%igz1) THEN
                    igz1=Domain%igz0
                  END IF  
                  dzLoc=Domain%zP(igz1+2**(-RefineNachbarZ))-Domain%zP(igz1)
                END IF
                Nachbars(in)%dLoc=dzLoc
                DO jx=jx0+1,jx1
                  DO jy=jy0+1,jy1
                    VolC(jx,jy,iz1+1)=VolCellGlobe(zP(iz1)+RadEarth &
                                                  ,zP(iz1)+RadEarth+dzLoc &
                                                  ,yP(jy-1) &
                                                  ,yP(jy) &
                                                  ,xP(jx-1) &
                                                  ,xP(jx))
                  END DO
                END DO
              END IF  
            END DO
          ELSE IF (CylinderT) THEN
            DO iz=iz0+1,iz1
              DO iy=iy0+1,iy1
                DO ix=ix0+1,ix1
                  VolC(ix,iy,iz)=VolCellCylinder(zP(iz-1) &
                                                ,zP(iz) &
                                                ,yP(iy-1) &
                                                ,yP(iy) &
                                                ,xP(ix-1) &
                                                ,xP(ix))
                END DO
              END DO
            END DO
          END IF
          READ(InputUnit,*) icut     !NrW_Cells
          DO i=1,icut
              READ(InputUnit,*) ix,iy,iz
              READ(InputUnit,*) VolC(ix,iy,iz)
          END DO

          IF (Globe) THEN
            IF (igy0==domain%igy0) THEN
              VolC(:,iy0,:)=Zero
            END IF  
            IF (iy1==domain%igy1) THEN
              VolC(:,iy1+1,:)=Zero
            END If  
          END IF
          IF (Coriolis) THEN
            IF (CoriolisFree.AND.Sphere) THEN
              DO ix=ix0+1,ix1
                lamL=xP(ix-1) 
                lamR=xP(ix) 
                CosLam=(SIN(lamR)-SIN(lamL))/(lamR-lamL)
                CosLam=COS(Half*(lamL+lamR))
                DO iy=iy0+1,iy1
                  phiL=yP(iy-1)
                  phiR=yP(iy)
                  SinPhi=-(COS(phiR)-COS(phiL))/(phiR-phiL)
                  SinPhi=SIN(Half*(phiL+phiR))
                  CosPhi=(SIN(phiR)-SIN(phiL))/(phiR-phiL)
                  CosPhi=COS(Half*(phiL+phiR))
                  fCor(ix,iy)=Two*Omega*(-CosLam*CosPhi*SIN(RotAngle)+SinPhi*COS(RotAngle))
                END DO
              END DO
            ELSE IF (CoriolisFree.AND.Cylinder) THEN
              fCor=Two*Omega
            ELSE IF (CoriolisFree.AND.BetaPlane) THEN
              DO ix=ix0+1,ix1
                DO iy=iy0+1,iy1
                  fCor(ix,iy)=Half*OmegaPlane*(yP(iy-1)+yP(iy))
                END DO
              END DO
            ELSE
              fCor=Two*Omega*SIN(Pi*PhiCor/180.0d0)
            END IF
          END IF
    
          
          MetrXY=One
          MetrXZ=One
          MetrYX=One
          MetrYZ=One
          MetrZX=One
          MetrZY=One
          IF (Sphere) THEN
            DO iy=iy0+1,iy1
              phiL=yP(iy-1)
              phiR=yP(iy)
              MetrXY(iy)=RadEarth*((SIN(phiR)-SIN(phiL))/(phiR-phiL))
            END DO
            DO ix=ix0+1,ix1
              MetrYX(ix)=RadEarth
            END DO
          ELSE IF (Cylinder) THEN
            DO iy=iy0+1,iy1
              MetrXY(iy)=Half*(yP(iy)+yP(iy-1))
            END DO
          END IF
    
          !.........................................................................
          !--- Determination of boundary cells

          READ(InputUnit,*)  NumBoundCell,Str_Read,nr_structbound                    !>>>>>><<<<<
          Floor(ib)%NumBoundCell=NumBoundCell
          !..........................................
          ALLOCATE(Floor(ib)%BoundCell(1:NumBoundCell))
          ALLOCATE(Floor(ib)%BoundCell3D(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1))
          Floor(ib)%BoundCell3d(:,:,:)=0
          DO i=1,NumBoundCell
            Floor(ib)%BoundCell3d(ix,iy,iz)=i
            ALLOCATE(Floor(ib)%BoundCell(i)%SoilType(Domain%nrsoillayers))  !!! 7 ersetzen
            READ(InputUnit,*) Floor(ib)%BoundCell(i)%ix &
                             ,Floor(ib)%BoundCell(i)%iy &
                             ,Floor(ib)%BoundCell(i)%iz
            READ(InputUnit,*) Floor(ib)%BoundCell(i)%xS &
                             ,Floor(ib)%BoundCell(i)%yS &
                             ,Floor(ib)%BoundCell(i)%zS
            ix=Floor(ib)%BoundCell(i)%ix
            iy=Floor(ib)%BoundCell(i)%iy
            iz=Floor(ib)%BoundCell(i)%iz
            SELECT CASE (nr_structbound)  
              CASE (4)  !OnlySoilType
               !Floor(ib)%BoundCell(i)%LandClass=0
                Floor(ib)%BoundCell(i)%LandClass=LandClassDef
                READ(InputUnit,*) Floor(ib)%BoundCell(i)%zRauh &
                                 ,Floor(ib)%BoundCell(i)%alb &
                                 ,Floor(ib)%BoundCell(i)%ee &
                                 ,Floor(ib)%BoundCell(i)%xFLS &
                                 ,Floor(ib)%BoundCell(i)%yFLS &
                                 ,Floor(ib)%BoundCell(i)%zFLS &
                                 ,Floor(ib)%BoundCell(i)%SoilType
                Floor(ib)%BoundCell(i)%zRauh=0.05
              CASE (5)  !OnlyLandType (only for LandClass=9 ("sea"))
                READ(InputUnit,*) Floor(ib)%BoundCell(i)%zRauh &
                                 ,Floor(ib)%BoundCell(i)%alb &
                                 ,Floor(ib)%BoundCell(i)%ee &
                                 ,Floor(ib)%BoundCell(i)%xFLS &
                                 ,Floor(ib)%BoundCell(i)%yFLS &
                                 ,Floor(ib)%BoundCell(i)%zFLS &
                                 ,Floor(ib)%BoundCell(i)%LandClass 
                IF (Floor(ib)%BoundCell(i)%LandClass/=9) Floor(ib)%BoundCell(i)%LandClass=9
                Floor(ib)%BoundCell(i)%zRauh=zrough(Floor(ib)%BoundCell(i)%LandClass)
              CASE (6)  !LandAndSoilType
                READ(InputUnit,*) Floor(ib)%BoundCell(i)%zRauh &
                                 ,Floor(ib)%BoundCell(i)%alb &
                                 ,Floor(ib)%BoundCell(i)%ee &
                                 ,Floor(ib)%BoundCell(i)%xFLS &
                                 ,Floor(ib)%BoundCell(i)%yFLS &
                                 ,Floor(ib)%BoundCell(i)%zFLS &
                                 ,Floor(ib)%BoundCell(i)%LandClass &
                                 ,Floor(ib)%BoundCell(i)%SoilType
                Floor(ib)%BoundCell(i)%zRauh=zrough(Floor(ib)%BoundCell(i)%LandClass)
              CASE DEFAULT  !(CASE (3))
               !Floor(ib)%BoundCell(i)%LandClass=0
                Floor(ib)%BoundCell(i)%LandClass=LandClassDef
                Floor(ib)%BoundCell(i)%SoilType(Domain%nrsoillayers)=5
                Floor(ib)%BoundCell(i)%zRauh=zrough(Floor(ib)%BoundCell(i)%LandClass)
                READ(InputUnit,*) Floor(ib)%BoundCell(i)%zRauh &
                                 ,Floor(ib)%BoundCell(i)%alb &
                                 ,Floor(ib)%BoundCell(i)%ee &
                                 ,Floor(ib)%BoundCell(i)%xFLS &
                                 ,Floor(ib)%BoundCell(i)%yFLS &
                                 ,Floor(ib)%BoundCell(i)%zFLS 
            END SELECT
            ! MJ
            IF (VolC(ix,iy,iz)>=(dx(ix)*dy(iy)*dz(iz)*0.9999d0).AND.RealIsland) Floor(ib)%BoundCell(i)%LandClass=9
         !  WRITE (*,*) 'ReadWeights_Mod'
         !  WRITE (*,*) 'i,ix,iy,iz,LandClass,SoilType',i,ix,iy,iz,Floor(ib)%BoundCell(i)%LandClass,Floor(ib)%BoundCell(i)%SoilType
          END DO

          DO i=1,Floor(ib)%NumBoundCell
            ix=Floor(ib)%BoundCell(i)%ix
            iy=Floor(ib)%BoundCell(i)%iy
            iz=Floor(ib)%BoundCell(i)%iz
            Floor(ib)%BoundCell3d(ix,iy,iz)=i
            !..............................
            n1=FU(ix,iy,iz)-FU(ix-1,iy,iz)
            n2=FV(ix,iy,iz)-FV(ix,iy-1,iz)
            n3=FW(ix,iy,iz)-FW(ix,iy,iz-1)
            ! Hinneburg: falls gegenueberliegende Waende geschlossen
            !            --> FL entsprechend vergroessert
            IF (FU(ix,iy,iz)==0.d0.AND.FU(ix-1,iy,iz)==0.d0) n1=dy(iy)*dz(iz)*2.d0 ! Hinneburg
            IF (FV(ix,iy,iz)==0.d0.AND.FV(ix,iy-1,iz)==0.d0) n2=dx(ix)*dz(iz)*2.d0 ! Hinneburg
            IF (FW(ix,iy,iz)==0.d0.AND.FW(ix,iy,iz-1)==0.d0) n3=dx(ix)*dy(iy)*2.d0 ! Hinneburg
            FL=SQRT(n1*n1+n2*n2+n3*n3)
            IF (FL>0.d0) THEN
              Floor(ib)%BoundCell(i)%n1=n1/FL ! x component of the normal vector
              Floor(ib)%BoundCell(i)%n2=n2/FL ! y component of the normal vector
              Floor(ib)%BoundCell(i)%n3=n3/FL ! z component of the normal vector
              Floor(ib)%BoundCell(i)%FL=FL ! surface area of the cut cross section in a boundary cell
            ELSE
              ! Hinneburg: falls alle gegenueberliegenden Seiten gleich gross (offen bis geschlossen)
              !            --> keine Reibung
              Floor(ib)%BoundCell(i)%n1=0.d0 ! Hinneburg
              Floor(ib)%BoundCell(i)%n2=0.d0 ! Hinneburg
              Floor(ib)%BoundCell(i)%n3=1.d0 ! Hinneburg
              Floor(ib)%BoundCell(i)%FL=1.d-10 ! Hinneburg
            END IF
            IF (Sphere) THEN
              n1G1=VecCartGlobX(Floor(ib)%BoundCell(i)%xFLS & 
                               ,Floor(ib)%BoundCell(i)%yFLS & 
                               ,Floor(ib)%BoundCell(i)%n1   &
                               ,Floor(ib)%BoundCell(i)%n2   &
                               ,Floor(ib)%BoundCell(i)%n3)
              n2G1=VecCartGlobY(Floor(ib)%BoundCell(i)%xFLS & 
                               ,Floor(ib)%BoundCell(i)%yFLS & 
                               ,Floor(ib)%BoundCell(i)%n1   &
                               ,Floor(ib)%BoundCell(i)%n2   &
                               ,Floor(ib)%BoundCell(i)%n3)
              n3G1=VecCartGlobZ(Floor(ib)%BoundCell(i)%yFLS & 
                               ,Floor(ib)%BoundCell(i)%n1   &
                               ,Floor(ib)%BoundCell(i)%n2   &
                               ,Floor(ib)%BoundCell(i)%n3)
              Floor(ib)%BoundCell(i)%n1G=n1G1
              Floor(ib)%BoundCell(i)%n2G=n2G1
              Floor(ib)%BoundCell(i)%n3G=n3G1
              xFLSG=GlobCartXP(RadEarth &
                              +Floor(ib)%BoundCell(i)%zFLS &
                              ,Floor(ib)%BoundCell(i)%xFLS &
                              ,Floor(ib)%BoundCell(i)%yFLS)
              yFLSG=GlobCartYP(RadEarth &
                              +Floor(ib)%BoundCell(i)%zFLS &
                              ,Floor(ib)%BoundCell(i)%xFLS &
                              ,Floor(ib)%BoundCell(i)%yFLS)
              zFLSG=GlobCartZP(RadEarth &
                              +Floor(ib)%BoundCell(i)%zFLS &
                              ,Floor(ib)%BoundCell(i)%yFLS)
              xFLSGmin=ABS(GlobCartXP(RadEarth &
                              ,Floor(ib)%BoundCell(i)%xFLS &
                              ,Floor(ib)%BoundCell(i)%yFLS))
              Floor(ib)%BoundCell(i)%xFLS=xFLSG
              Floor(ib)%BoundCell(i)%yFLS=yFLSG
              Floor(ib)%BoundCell(i)%zFLS=zFLSG
              MaxHeightLoc=MAX(MaxHeightLoc,1.5d0*(ABS(Floor(ib)%BoundCell(i)%xFLS)-xFLSGmin))
            ELSE
              Floor(ib)%BoundCell(i)%n1G=Floor(ib)%BoundCell(i)%n1
              Floor(ib)%BoundCell(i)%n2G=Floor(ib)%BoundCell(i)%n2
              Floor(ib)%BoundCell(i)%n3G=Floor(ib)%BoundCell(i)%n3
              MaxHeightLoc=MAX(MaxHeightLoc,1.5d0*Floor(ib)%BoundCell(i)%zFLS)
            END IF
            !............................................
            Floor(ib)%BoundCell(i)%zRauhT=Floor(ib)%BoundCell(i)%zRauh
            !............................................
            Floor(ib)%BoundCell(i)%Tes=290
          END DO
         
          EXIT S
        END IF   ! ibInput
      END IF   ! INDEX...Block.. 
    END DO S
  END DO  ! ibLoc



  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    IF (GradFull) THEN
      FUG=0.0d0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1-1
            FUG(ix,iy,iz)=(VolC(ix,iy,iz)+VolC(ix+1,iy,iz))/(MetrXY(iy)*(dx(ix+1)+dx(ix)+Eps)) &
                          *FU(ix,iy,iz)/(FU(ix,iy,iz)+Eps)          
          END DO
        END DO
      END DO
      DO in=1,AnzahlNachbar
        CALL Set(Nachbars(in))
        IF (Nachbars(in)%nType(2:2)=='w') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              FUG(ix0,jy,jz)=(VolC(ix0,jy,jz)+VolC(ix0+1,jy,jz))/(MetrXY(jy)*(dx(ix0+1)+dLoc+Eps)) &
                            *FU(ix0,jy,jz)/(FU(ix0,jy,jz)+Eps)          
            END DO
          END DO
        ELSE IF (Nachbars(in)%nType(2:2)=='e') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              FUG(ix1,jy,jz)=(VolC(ix1,jy,jz)+VolC(ix1+1,jy,jz))/(MetrXY(jy)*(dLoc+dx(ix1)+Eps)) &
                            *FU(ix1,jy,jz)/(FU(ix1,jy,jz)+Eps)          
            END DO
          END DO
        END IF  
      END DO
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1-1
          DO ix=ix0+1,ix1
            FVG(ix,iy,iz)=(VolC(ix,iy,iz)+VolC(ix,iy+1,iz))/(MetrYX(ix)*(dy(iy+1)+dy(iy)+Eps)) &
                          *FV(ix,iy,iz)/(FV(ix,iy,iz)+Eps)          
          END DO
        END DO
      END DO
      DO in=1,AnzahlNachbar
        CALL Set(Nachbars(in))
        IF (Nachbars(in)%nType(2:2)=='s') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              FVG(jx,iy0,jz)=(VolC(jx,iy0,jz)+VolC(jx,iy0+1,jz))/(MetrYX(jx)*(dy(iy0+1)+dLoc+Eps)) &
                            *FV(jx,iy0,jz)/(FV(jx,iy0,jz)+Eps)          
            END DO
          END DO
        ELSE IF (Nachbars(in)%nType(2:2)=='n') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              FVG(jx,iy1,jz)=(VolC(jx,iy1,jz)+VolC(jx,iy1+1,jz))/(MetrYX(jx)*(dy(iy1)+dLoc+Eps)) &
                            *FV(jx,iy1,jz)/(FV(jx,iy1,jz)+Eps)          
            END DO
          END DO
        END IF
      END DO
      DO iz=iz0+1,iz1-1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            FWG(ix,iy,iz)=(VolC(ix,iy,iz)+VolC(ix,iy,iz+1))/(dz(iz+1)+dz(iz)+Eps) &
                          *FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps)          
          END DO
        END DO
      END DO
      DO in=1,AnzahlNachbar
        CALL Set(Nachbars(in))
        IF (Nachbars(in)%nType(2:2)=='b') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              FWG(jx,jy,iz0)=(VolC(jx,jy,iz0)+VolC(jx,jy,iz0+1))/(dz(iz0+1)+dLoc+Eps) &
                            *FW(jx,jy,iz0)/(FW(jx,jy,iz0)+Eps)          
            END DO
          END DO
        ELSE IF (Nachbars(in)%nType(2:2)=='t') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              FWG(jx,jy,iz1)=(VolC(jx,jy,iz1)+VolC(jx,jy,iz1+1))/(dz(iz1)+dLoc+Eps) &
                            *FW(jx,jy,iz1)/(FW(jx,jy,iz1)+Eps)          
            END DO
          END DO
        END IF
      END DO
    ELSE
      FUG=FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)
      FVG=FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)
      FWG=FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)
    END IF  
  END DO
  DO ib=1,nb
    Floor(ib)%FreeCells=0
    CALL Set(Floor(ib))
    IF (blMPI(ib)%Proc==MyId) THEN
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              IF (ix0Out+1<=ix.AND.ix<=ix1Out.AND. &
                  iy0Out+1<=iy.AND.iy<=iy1Out.AND. &
                  iz0Out+1<=iz.AND.iz<=iz1Out) THEN
                 Floor(ib)%FreeCells=Floor(ib)%FreeCells+1
              END IF
            END IF
          END DO
        END DO
      END DO
    END IF
    CALL MPI_Bcast(Floor(ib)%FreeCells, 1, MPI_INTEGER, blMPI(ib)%Proc,  &
&                  MPI_COMM_WORLD, MPIErr)
    CALL MPI_Bcast(Floor(ib)%NumBoundCell, 1, MPI_INTEGER, blMPI(ib)%Proc,  &
&                  MPI_COMM_WORLD, MPIErr)
  END DO
  CLOSE(InputUnit)
  CALL MPI_Allreduce(MaxHeightLoc,MaxHeight,1,MPI_RealKind, &
&                    MPI_MAX,MPI_Comm_World,MPIErr)

! Volume and Face check  
  IF (VolumeCheck) THEN
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      WRITE(*,*) 'Block ',ib
      WRITE(*,*) 'WestFace  ',SUM(FU(ix0-1,:,:)),SUM(FU(ix0,:,:)),SUM(FU(ix0+1,:,:))
      WRITE(*,*) 'EastFace  ',SUM(FU(ix1-1,:,:)),SUM(FU(ix1,:,:)),SUM(FU(ix1+1,:,:))
      WRITE(*,*) 'SouthFace ',SUM(FV(:,iy0-1,:)),SUM(FV(:,iy0,:)),SUM(FV(:,iy0+1,:))
      WRITE(*,*) 'NorthFace ',SUM(FV(:,iy1-1,:)),SUM(FV(:,iy1,:)),SUM(FV(:,iy1+1,:))
      WRITE(*,*) 'BottomFace',SUM(FW(:,:,iz0))
      WRITE(*,*) 'TopFace   ',SUM(FW(:,:,iz1))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)==Zero) THEN
              IF (FU(ix-1,iy,iz)+FU(ix,iy,iz) &
                 +FV(ix,iy-1,iz)+FV(ix,iy,iz) &
                 +FW(ix,iy,iz-1)+FV(ix,iy,iz)>Zero) THEN
                WRITE(*,*) 'False Cell',ix,iy,iz 
              END IF  
            END IF  
            IF (FU(ix-1,iy,iz)+FU(ix,iy,iz) &
               +FV(ix,iy-1,iz)+FV(ix,iy,iz) &
               +FW(ix,iy,iz-1)+FV(ix,iy,iz)==Zero) THEN
              IF (VolC(ix,iy,iz)>Zero) THEN
                WRITE(*,*) 'False Face in Cell',ix,iy,iz 
              END IF  
            END IF  
          END DO  
        END DO  
      END DO  
    END DO  
  END IF  
! WRITE(*,*) 'CompareEW' 
! CALL CompareEW(Floor(1),Floor(2))
! CALL CompareEW(Floor(2),Floor(3))
! CALL CompareEW(Floor(1),Floor(1))

  CALL Allocate_Velocity(VolFace)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
!       VolFace(ib)%u_w(iy,iz)=Half*SUM(VolC(ix0+1:ix1,iy,iz))
!       VolFace(ib)%u_e(iy,iz)=Half*SUM(VolC(ix0+1:ix1,iy,iz))
!       VolFace(ibLoc)%u_w(iy,iz)=Half*VolC(ix0+1,iy,iz)
!       VolFace(ibLoc)%u_e(iy,iz)=Half*VolC(ix1,iy,iz)
        VolFace(ibLoc)%u_w(iy,iz)=Half*dx(ix0+1)*dy(iy)*dz(iz)
        VolFace(ibLoc)%u_e(iy,iz)=Half*dx(ix1)*dy(iy)*dz(iz)
      END DO  
    END DO  
    DO ix=ix0+1,ix1
      DO iz=iz0+1,iz1
!       VolFace(ibLoc)%v_s(ix,iz)=Half*SUM(VolC(ix,iy0+1:iy1,iz))
!       VolFace(ibLoc)%v_n(ix,iz)=Half*SUM(VolC(ix,iy0+1:iy1,iz))
        VolFace(ibLoc)%v_s(ix,iz)=Half*VolC(ix,iy0+1,iz)
        VolFace(ibLoc)%v_n(ix,iz)=Half*VolC(ix,iy1,iz)
        VolFace(ibLoc)%v_s(ix,iz)=Half*dx(ix)*dy(iy0+1)*dz(iz)
        VolFace(ibLoc)%v_n(ix,iz)=Half*dx(ix)*dy(iy1)*dz(iz)
      END DO  
    END DO  
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
!       VolFace(ibLoc)%w_b(ix,iy)=Half*SUM(VolC(ix,iy,iz0+1:iz1))
!       VolFace(ibLoc)%w_t(ix,iy)=Half*SUM(VolC(ix,iy,iz0+1:iz1))
        VolFace(ibLoc)%w_b(ix,iy)=Half*VolC(ix,iy,iz0+1)
        VolFace(ibLoc)%w_t(ix,iy)=Half*VolC(ix,iy,iz1)
      END DO  
    END DO  
  END DO
  CALL Exchange(VolFace)

  IF (Canopy) THEN
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      CALL ReadCanopy(FileName)
    END DO
  END IF

  WRITE(*,*) 'TestWeight'
  CALL TestWeight
! CALL TestWeight2

CONTAINS

FUNCTION FaceCartXY(r,phi1,phi2,lam1,lam2)
  REAL(RealKind) :: FaceCartXY
  REAL(RealKind) :: r
  REAL(RealKind) :: phi1,phi2,lam1,lam2
  !Wertebereich: phi [-Pi/2;Pi/2]
  !Wertebereich: lambda [0;2*Pi]
  FaceCartXY=MAX(r**2*(SIN(phi2)-SIN(phi1))*(lam2-lam1),Zero)
END FUNCTION FaceCartXY

FUNCTION FaceCartYZ(r1,r2,phi1,phi2)
  REAL(RealKind) :: FaceCartYZ
  REAL(RealKind) :: r1,r2,phi1,phi2
  FaceCartYZ=MAX(0.5*(phi2-phi1)*(r2**2-r1**2),Zero)
END FUNCTION FaceCartYZ

FUNCTION FaceCartZX(r1,r2,phi,lam1,lam2)
  REAL(RealKind) :: FaceCartZX
  REAL(RealKind) :: r1,r2,phi,lam1,lam2
  FaceCartZX=MAX(0.5*COS(phi)*(lam2-lam1)*(r2**2-r1**2),Zero)
END FUNCTION FaceCartZX

FUNCTION VolCellCart(r1,r2,phi1,phi2,lam1,lam2)
  REAL(RealKind) :: VolCellCart
  REAL(RealKind) :: r1,r2,phi1,phi2,lam1,lam2
  VolCellCart=MAX(1.0d0/3.0d0*((lam2-lam1)*(r2**3-r1**3)*(SIN(phi2)-SIN(phi1))),Zero)
END FUNCTION VolCellCart


FUNCTION FaceGlobeXY(r,phi1,phi2,lam1,lam2)
  REAL(RealKind) :: FaceGlobeXY
  REAL(RealKind) :: r
  REAL(RealKind) :: phi1,phi2,lam1,lam2
  !Wertebereich: phi [-Pi/2;Pi/2]
  !Wertebereich: lambda [0;2*Pi]
  FaceGlobeXY=MAX(r**2*(SIN(phi2)-SIN(phi1))*(lam2-lam1),Zero)
END FUNCTION FaceGlobeXY

FUNCTION FaceGlobeYZ(r1,r2,phi1,phi2)
  REAL(RealKind) :: FaceGlobeYZ
  REAL(RealKind) :: r1,r2,phi1,phi2
  FaceGlobeYZ=MAX(0.5*(phi2-phi1)*(r2**2-r1**2),Zero)
END FUNCTION FaceGlobeYZ

FUNCTION FaceGlobeZX(r1,r2,phi,lam1,lam2)
  REAL(RealKind) :: FaceGlobeZX
  REAL(RealKind) :: r1,r2,phi,lam1,lam2
  FaceGlobeZX=MAX(0.5*COS(phi)*(lam2-lam1)*(r2**2-r1**2),Zero)
END FUNCTION FaceGlobeZX

FUNCTION VolCellGlobe(r1,r2,phi1,phi2,lam1,lam2)
  REAL(RealKind) :: VolCellGlobe
  REAL(RealKind) :: r1,r2,phi1,phi2,lam1,lam2
  VolCellGlobe=MAX(1.0d0/3.0d0*((lam2-lam1)*(r2**3-r1**3)*(SIN(phi2)-SIN(phi1))),Zero)
END FUNCTION VolCellGlobe

FUNCTION FaceCylinderXY(r1,r2,lam1,lam2)
  REAL(RealKind) :: FaceCylinderXY
  REAL(RealKind) :: r1,r2,lam1,lam2
  !Wertebereich: lambda [0;2*Pi]
  FaceCylinderXY=ABS((r2+r1)*(r2-r1)*(lam2-lam1)/2.0d0)
END FUNCTION FaceCylinderXY

FUNCTION FaceCylinderYZ(z1,z2,r1,r2)
  REAL(RealKind) :: FaceCylinderYZ
  REAL(RealKind) :: z1,z2,r1,r2
  FaceCylinderYZ=(z2-z1)*(r2-r1)
END FUNCTION FaceCylinderYZ

FUNCTION FaceCylinderZX(z1,z2,r,lam1,lam2)
  REAL(RealKind) :: FaceCylinderZX
  REAL(RealKind) :: z1,z2,r,lam1,lam2
  FaceCylinderZX=ABS(r*(lam2-lam1)*(z2-z1))
END FUNCTION FaceCylinderZX

FUNCTION VolCellCylinder(z1,z2,r1,r2,lam1,lam2)
  REAL(RealKind) :: VolCellCylinder
  REAL(RealKind) :: z1,z2,r1,r2,lam1,lam2

  VolCellCylinder=0.5d0*ABS((lam2-lam1) &
                 *(r2**2-r1**2) &
                 *(z2-z1))
END FUNCTION VolCellCylinder

FUNCTION GlobCartXP(r1,lam,phi)
  REAL(RealKind) :: GlobCartXP
  REAL(RealKind) :: r1,phi,lam
  GlobCartXP=r1*COS(lam)*COS(phi) 
!  GlobCartXP=r1*SIN(lam)*COS(phi) 
END FUNCTION GlobCartXP

FUNCTION GlobCartYP(r1,lam,phi)
  REAL(RealKind) :: GlobCartYP
  REAL(RealKind) :: r1,phi,lam
  GlobCartYP=r1*SIN(lam)*COS(phi) 
!  GlobCartYP=r1*COS(lam)*COS(phi) 
END FUNCTION GlobCartYP

FUNCTION GlobCartZP(r1,phi)
  REAL(RealKind) :: GlobCartZP
  REAL(RealKind) :: r1,phi
  GlobCartZP=r1*SIN(phi) 
END FUNCTION GlobCartZP

FUNCTION VecCartGlobX(lam,phi,e1,e2,e3)
  REAL(RealKind) :: VecCartGlobX
  REAL(RealKind) :: lam,phi
  REAL(RealKind) :: e1,e2,e3
  VecCartGlobX=-SIN(lam)*e1           &
               -SIN(phi)*COS(lam)*e2  &
               +COS(phi)*COS(lam)*e3
END FUNCTION VecCartGlobX

FUNCTION VecCartGlobY(lam,phi,e1,e2,e3)
  REAL(RealKind) :: VecCartGlobY
  REAL(RealKind) :: lam,phi
  REAL(RealKind) :: e1,e2,e3
  VecCartGlobY= COS(lam)*e1           &
               -SIN(phi)*SIN(lam)*e2  &
               +COS(phi)*SIN(lam)*e3
END FUNCTION VecCartGlobY

FUNCTION VecCartGlobZ(phi,e1,e2,e3)
  REAL(RealKind) :: VecCartGlobZ
  REAL(RealKind) :: phi
  REAL(RealKind) :: e1,e2,e3
  VecCartGlobZ= COS(phi)*e2  &
               +SIN(phi)*e3
END FUNCTION VecCartGlobZ

END SUBROUTINE ReadWeights
END MODULE ReadWeights_Mod
