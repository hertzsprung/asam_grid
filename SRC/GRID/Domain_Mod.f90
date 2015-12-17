MODULE Domain_Mod


  USE Kind_Mod
  USE Geometry_Mod
  USE Polygon_Mod
  USE Neighbor_Mod

  IMPLICIT NONE

  TYPE Topo_T
    TYPE (Vertex_T), POINTER :: Vert=>NULL()
    INTEGER :: ecnrxy,ecnryz,ecnrzx,exnr,eynr,eznr ! Edges Nummern
    INTEGER :: iecnrxy,iecnryz,iecnrzx,iexnr,ieynr,ieznr ! Edges Nummern inside
    INTEGER :: FCutNr              ! FaceCut-Nummer
    INTEGER :: fxynr,fzxnr,fyznr   ! Face(XY,YZ,ZX)-Nummern
    INTEGER :: ifxynr,ifzxnr,ifyznr   ! Face(XY,YZ,ZX)-Nummern inside
    INTEGER :: cnr                 ! Cell-Nummer
    INTEGER :: icnr                 ! Cell-Nummer inside
    INTEGER :: ctp                 ! Celle-Topo-Index
    ! später zur Cell-Struktur hinzu 
    INTEGER :: nfn                 ! Celle-nFaces
    INTEGER :: cfn(1:7)            ! Celle-Face-Liste mit FCut
    INTEGER :: fceg_nr             ! anz. EdgeCut des FaceCut
    INTEGER :: fceg_nrcut(1:8)     ! Liste Nummer EdgeCut des FaceCut
    INTEGER :: infn                ! Celle-nFaces
    INTEGER :: icfn(1:7)           ! Celle-Face-Liste mit FCut
    ! in Bearbeitung
  END TYPE

  INTEGER :: topo_cnr,topo_cvcnr                    ! all ib
  INTEGER :: topo_enr,topo_ecnr,topo_fnr,topo_fcnr  ! all ib
  !----
  INTEGER :: stopo_cnr, stopo_cv                    ! summe je ib
  INTEGER :: stopo_enr,stopo_ecnr                   ! ""
  INTEGER :: topo_pos_ec1,topo_pos_ecn              ! temp je block
  INTEGER :: stopo_exnr,stopo_eynr,stopo_eznr       ! summe je block 
  INTEGER :: stopo_fnr,stopo_fcnr                   ! "" 
  INTEGER :: stopo_fxynr,stopo_fzxnr,stopo_fyznr    ! ""
  !----
  INTEGER :: para_oro_cnr,para_grid_cnr, & !pview_oro_cnr...
      &      para_oro_enr,para_grid_enr, &
      &      para_oro_fnr,para_grid_fnr
  !----
  INTEGER, PARAMETER :: topo_border=1,  &
               &        topo_inside=0
  INTEGER, PARAMETER :: tropos_inside=1, &
               &        tropos_border=0,  &
               &        soil_inside=-1

  TYPE TopoV_T
    TYPE(Topo_T), POINTER :: Topo=>NULL()
  END TYPE

  ! Upper Bounds Layer Celle
  TYPE UpBounds_T
    INTEGER :: i,j,k     !Index Topo(:,:,:)
    INTEGER :: nrP       !nrP-Index_bounds and border
    INTEGER :: nrEcx,nrEcy,nrEcz  ! nrP-Verts-Edges
    !..............................
    INTEGER :: bexnr,beynr,beznr  ! Edges Nummern
    INTEGER :: becnrxy,becnryz,becnrzx ! Edges_ec Nummern   
    !..............................
    INTEGER :: bfxynr,bfzxnr,bfyznr   ! Face(XY,YZ,ZX)-Nummern
    INTEGER :: flisteg_nr(1:8)        ! liste-edge-zur Ebene
    INTEGER :: feg_nr                 ! anzahl edge-nummer zur Ebene
    !..............................
    INTEGER :: FCutNr              ! FaceCut-Nummer
    INTEGER :: fceg_nr             ! anz. EdgeCut des FaceCut
    INTEGER :: fceg_nrcut(1:8)     ! Liste Nummer EdgeCut des FaceCut
    !..............................
    INTEGER :: ub_cnr        ! NummerBoundCelle    
    CHARACTER :: ubc
    INTEGER :: nfn           ! Celle-nFaces
    INTEGER :: cfn(1:7)      ! Celle-Face-Liste mit FCut
    !da auf Ebene ausgerichtet auch Value mit "0" oder nicht belegt
  END TYPE

  TYPE UpBoundsV_T
    TYPE(UpBounds_T), POINTER :: BCell=>NULL()
  END TYPE

  TYPE DimUpBound_T
    INTEGER :: uZ
    INTEGER :: bZ
  END TYPE

  TYPE DimUpBoundV_T
    TYPE(DimUpBound_T), POINTER :: DimG
  END TYPE
  INTEGER :: minZ_UBL,maxZ_UBL
  INTEGER :: upbound_p,upbound_ps,upbound_enr,upbound_ecnr     ! summe all
  INTEGER :: gl_benr       ! sum identical number benr in i-w,e,s,n,b,t 
  INTEGER :: upbound_fnr,upbound_fcnr,upbound_cnr,upbound_cvc  ! ""
  INTEGER :: upb_pos_ec1,upb_pos_ecn   ! temp je block
  INTEGER :: sbexnr,sbeynr,sbeznr,sbecnr,sbenr      ! counter je block
  INTEGER :: sbfcnr,sbfxynr,sbfyznr,sbfzxnr,sbfnr    ! ""
  INTEGER :: sbcvcnr,sbchnr,sbcnr                    ! ""
  INTEGER :: sbenr_bdr(1:6),sbfnr_bdr(1:6),sbcnr_bdr(1:6) 
  !          number edges-faces-cells border: w,e,s,n,b,t je block

  !only grid
  INTEGER ::  upboundgi_enr,upboundgi_fnr,upboundgi_cnr
!.....................................................  
  CHARACTER(8), PARAMETER :: gmvinput ='gmvinput'
  CHARACTER(8), PARAMETER :: ascii    ='   ascii'
  CHARACTER(8), PARAMETER :: ieeei4r4 ='ieeei4r4'
  CHARACTER(8), PARAMETER :: ieeei4r8 ='ieeei4r8'
  CHARACTER(8), PARAMETER :: iecxi4r4 ='iecxi4r4'
  CHARACTER(8), PARAMETER :: iecxi4r8 ='iecxi4r8'
  CHARACTER(8), PARAMETER :: nodes    ='nodes   '
  CHARACTER(8), PARAMETER :: nodev    ='nodev   '
  CHARACTER(8), PARAMETER :: cellsGMV ='cells   '
  CHARACTER(8), PARAMETER :: material ='material'
  CHARACTER(8), PARAMETER :: velocity ='velocity'
  CHARACTER(8), PARAMETER :: hex      ='hex     '
  CHARACTER(8), PARAMETER :: general  ='general '
  CHARACTER(8), PARAMETER :: polygons ='polygons'
  CHARACTER(8), PARAMETER :: endpoly  ='endpoly '
  CHARACTER(8), PARAMETER :: endflag  ='endflag '
  CHARACTER(8), PARAMETER :: fromfile ='fromfile'
  CHARACTER(8), PARAMETER :: endgmv   ='endgmv  '
  !..................................................
  CHARACTER(10), PARAMETER :: def_out_index='XYZ_Number' !'Gr_Index'
  CHARACTER( 9), PARAMETER :: def_out_coord='XYZ_Coord'  !'Gr_Coord'
  !..................................................
  TYPE Boundary_T
    CHARACTER(10) :: West=''
    CHARACTER(10) :: East=''
    CHARACTER(10) :: South=''
    CHARACTER(10) :: North=''
    CHARACTER(10) :: Bottom=''
    CHARACTER(10) :: Top=''
  END TYPE Boundary_T

  TYPE(Boundary_T), SAVE :: BC

  TYPE Domain_T
    CHARACTER(2) :: TypeW,TypeE,TypeS,TypeN,TypeB,TypeT
    REAL(RealKind) :: Boundary
    INTEGER :: nx,ny,nz            ! numbers def. x-,y-,z-direction the domain-grid
    INTEGER :: nc                  ! number of cells the domain-grid
    INTEGER :: nr_out              ! number points output for gmv 
    INTEGER :: nr_view_out         ! number points view-domain output grid for gmv
    INTEGER :: WriteOffsetC        ! Cell Offset for parallel I/O
    INTEGER :: WriteOffsetN        ! Node Offset for parallel I/O
    INTEGER :: igx0,igx1,igy0,igy1,igz0,igz1
    INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1,ncz
    !......................................................................
    !Output-Area-ViewGMV
    CHARACTER(10)   :: def_out_domain
    INTEGER        :: view_ixa,view_ixe,view_iya,view_iye,view_iza,view_ize
    REAL(RealKind) :: view_xa,view_xe,view_ya,view_ye,view_za,view_ze
    !......................................................................
    INTEGER :: ib
    INTEGER :: RefLevel
    INTEGER :: Refine
    INTEGER :: RefineX
    INTEGER :: RefineY
    INTEGER :: RefineZ
    INTEGER :: Xshift
    INTEGER :: Yshift
    INTEGER :: Zshift
    REAL(RealKind) :: x0,x1,y0,y1,z0,z1,calc_z
    REAL(RealKind) :: x0View,x1View,y0View,y1View,z0View,z1View !Domain
    REAL(RealKind), POINTER :: dx(:),dy(:),dz(:)
    REAL(RealKind), POINTER :: xP(:),yP(:),zP(:)
    TYPE(VertexP_T), POINTER :: Vertices(:,:,:)
    TYPE (EdgeP_T), POINTER :: Edges_X(:,:,:)
    TYPE (EdgeP_T), POINTER :: Edges_Y(:,:,:)
    TYPE (EdgeP_T), POINTER :: Edges_Z(:,:,:)
    TYPE (FaceP_T), POINTER :: Faces_XY(:,:,:)
    TYPE (FaceP_T), POINTER :: Faces_YZ(:,:,:)
    TYPE (FaceP_T), POINTER :: Faces_ZX(:,:,:)
    TYPE (CellP_T), POINTER :: Cells(:,:,:)
    INTEGER          :: nr_soildef
    INTEGER, POINTER :: soil_type(:,:)
    INTEGER, POINTER :: s_ixa(:),s_ixe(:),s_iya(:),s_iye(:),s_iza(:),s_ize(:)
    REAL(RealKind), POINTER :: WeiFU(:,:,:),WeiFV(:,:,:),WeiFW(:,:,:)
    REAL(RealKind), POINTER :: VolC(:,:,:)
    INTEGER :: NrW_Cells,NrR_Cells,NrRN_Cells,NrRW_Cells
    !INTEGER :: NrW_Cells1,NrW_Cells2
    INTEGER :: NrMP_Cells,NrB_Cells
    INTEGER :: NrW_FacesXY,NrW_FacesYZ,NrW_FacesZX
    INTEGER :: NrR_FacesXY,NrR_FacesYZ,NrR_FacesZX
    INTEGER :: NrRN_FacesXY,NrRN_FacesYZ,NrRN_FacesZX
    INTEGER :: NrRW_FacesXY,NrRW_FacesYZ,NrRW_FacesZX
    INTEGER :: NrMP_FacesXY,NrMP_FacesYZ,NrMP_FacesZX,NrMP_F
    !INTEGER :: sub_face_gr(1:6)  !Index Grenze Face subMountain
    INTEGER :: AnzahlNachbar
    TYPE (Nachbar_T), POINTER :: Nachbars(:)
    !......................................................................
    ! Topo
    TYPE (TopoV_T), POINTER :: Topo(:,:,:)
    INTEGER :: stopo_cnr,stopo_cv
    INTEGER :: stopo_enr,stopo_ecnr
    INTEGER :: topo_pos_ec1,topo_pos_ecn
    INTEGER :: stopo_exnr,stopo_eynr,stopo_eznr
    INTEGER :: stopo_fxynr,stopo_fzxnr,stopo_fyznr
    INTEGER :: stopo_fnr,stopo_fcnr
    !......................................................................
    ! Upper Bounds Layer
    TYPE (UpBoundsV_T), POINTER :: UpBoundsLayer(:,:,:)
    TYPE (DimUpBoundV_T),  POINTER :: DimUpBounds(:,:)
    INTEGER, POINTER :: DimUpBoundsArray(:,:,:)
    INTEGER, POINTER :: DimUpBoundsN(:,:,:)
    INTEGER :: sbcvcnr,sbchnr,sbcnr,sbcnr_bdr(1:6)
    INTEGER :: sbexnr,sbeynr,sbeznr,sbecnr,sbenr,sbenr_bdr(1:6)
    INTEGER :: upb_pos_ec1,upb_pos_ecn
    INTEGER :: sbfxynr,sbfzxnr,sbfyznr,sbfcnr,sbfnr,sbfnr_bdr(1:6)
  END TYPE Domain_T
  
  TYPE (Domain_T) :: Domain
  INTEGER :: nx,ny,nz 
  INTEGER :: WriteOffsetC
  INTEGER :: WriteOffsetN
  INTEGER :: igx0,igx1,igy0,igy1,igz0,igz1
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1,ncz
  REAL(RealKind) :: calc_z
  REAL(RealKind), POINTER :: dx(:),dy(:),dz(:)
  REAL(RealKind), POINTER :: xP(:),yP(:),zP(:)
  TYPE(VertexP_T), POINTER :: Vertices(:,:,:)
  TYPE (EdgeP_T), POINTER :: Edges_X(:,:,:)
  TYPE (EdgeP_T), POINTER :: Edges_Y(:,:,:)
  TYPE (EdgeP_T), POINTER :: Edges_Z(:,:,:)

  TYPE (FaceP_T), POINTER :: Faces_XY(:,:,:)
  TYPE (FaceP_T), POINTER :: Faces_YZ(:,:,:)
  TYPE (FaceP_T), POINTER :: Faces_ZX(:,:,:)
 
  TYPE (CellP_T), POINTER :: Cells(:,:,:)

  REAL(RealKind), POINTER :: FU(:,:,:),FV(:,:,:),FW(:,:,:)
  REAL(RealKind), POINTER :: VolC(:,:,:)
  
  TYPE (TopoV_T), POINTER :: Topo(:,:,:)

  TYPE (UpBoundsV_T), POINTER :: UpBoundsLayer(:,:,:)
  TYPE (DimUpBoundV_T), POINTER :: DimUpBounds(:,:)
  INTEGER, POINTER :: DimUpBoundsArray(:,:,:)
  INTEGER, POINTER :: DimUpBoundsN(:,:,:)

  INTEGER :: nr_sb, struct_bound, std_bound=3
  INTEGER :: nr_lsoil !=7   !akt. def. Bodenschichten Cell-Soil-Darstellung
  INTEGER :: nr_soildef
  REAL(RealKind) :: dzi_soil
  INTEGER, POINTER :: soil_type(:,:)
  INTEGER, POINTER :: s_ixa(:),s_ixe(:),s_iya(:),s_iye(:),s_iza(:),s_ize(:)
  
  INTEGER :: nr_landdef
  INTEGER, POINTER :: LandDef(:,:,:)

  INTEGER :: AnzahlNachbar
  TYPE (Nachbar_T), POINTER :: Nachbars(:)
  TYPE(Nachbar_T), POINTER :: Nachbar
  INTEGER :: RefLevel
  INTEGER :: Refine
  INTEGER :: RefineX
  INTEGER :: RefineY
  INTEGER :: RefineZ
  REAL(RealKind) :: Boundary
  CHARACTER(2) :: TypeE,TypeW,TypeS,TypeN,TypeB,TypeT
  
  INTEGER :: ibc

  INTEGER :: nr_out=0,nr_view_out,nr_EgVerts,nr_viewEgVerts
! TYPE (Vertex_T), POINTER :: VertOutView(:)=>NULL()
  INTEGER, POINTER         :: VertViewScale(:)=>NULL()
  INTEGER :: nr_cutplane=0
  TYPE (Vertex_T), POINTER :: VertCutPOut(:)=>NULL()
  INTEGER , POINTER        :: VertNrCutPOut(:)=>NULL()
  INTEGER :: nr_out_soil
! TYPE (Vertex_T), POINTER :: VertSoilOut(:)=>NULL()


  INTEGER :: nr_out2
  TYPE (Vertex_T), POINTER :: VertOut2(:)=>NULL() !  zu ParaView-PartOut2-Oro

  INTEGER :: nr_upbounds,znr_upbounds ! Summe Points im Bound-area
  INTEGER :: nrp_upb_nc,nrp_upb_gi    ! nr-Points-Schnitt, nr-Points 
                                      ! auf Gitterpoint
  INTEGER :: out_bounds               ! Summe Points Bound-Out_W_E z.Zt.
  INTEGER :: per_bounds               ! Summe Points Bound-Periode z.Zt.
  INTEGER :: b_upbounds=0             ! Summe Points Bound-Periode z.Zt.
  INTEGER :: i_upbounds=0             ! Summe Points inside BorderBounds indentical
  INTEGER, POINTER  :: IndexOutBounds(:,:)=>NULL()  ! Index NrP auf VertOut(:),
                                                    ! and i,j,k
  INTEGER, POINTER  :: NrPOutBounds(:)=>NULL()      ! NrP-VertOut-Position-
                                                    ! NrP-Bounds
  TYPE (Vertex_T), POINTER  :: VertOutBounds(:)=>NULL()
                                                    ! VertBound auf VertOut(:)
  TYPE (Vertex_T), POINTER  :: VertOutBoundsBord(:)=>NULL() ! Border
  INTEGER, POINTER  :: IndexBordBounds(:,:)=>NULL()  ! Index NrP auf VertOutiBounds(:),
                                                    ! and i,j,k


  !Cells counter for GMV
  INTEGER :: nr_cutinside
  INTEGER :: nr_cells
  INTEGER :: nr_viewcells                              ! outside Topogr.
  INTEGER :: nr_cutplanecells,nr_viewcutplanecells     ! cut Topogr.
  INTEGER :: nr_cutplanecells1,nr_viewcutplanecells1   ! only 'cut' Topogr.
  INTEGER :: nr_cutplanecells2,nr_viewcutplanecells2   ! only Celle k==1 'cut' Topogr.
  INTEGER :: nr_soilplanecells,nr_viewsoilplanecells   ! cut=soilplane Topogr.
  INTEGER :: nr_soilplanecells1,nr_viewsoilplanecells1 !  
  INTEGER :: nr_soilplanecells2,nr_viewsoilplanecells2 ! only Celle k==1 'soilplane'
  INTEGER :: nr_cutcells,nr_cutf1,nr_incells,nr_orocells        ! inside Topogr.
  INTEGER :: nr_viewcutcells,nr_viewcutf1,nr_viewincells,nr_vieworocells
  INTEGER :: nr_CellOroVolNull
  ! into  Oro.out
  INTEGER :: nr_gen=0
  INTEGER :: nr_grenzeF2=0 
  INTEGER :: nr_grenzeF1=0
  INTEGER :: nr_grenzeSp3=0
  INTEGER :: nr_insidehex=0
  INTEGER :: nr_cell_vc_oro=0
  INTEGER :: nr_cell_novc_oro=0
  INTEGER :: nr_ges_oro_out=0

  !Cells counter for Weight
  INTEGER :: NrAll_Cells,NrAll_FacesXY,NrAll_FacesYZ,NrAll_FacesZX
  INTEGER :: NrW_Cells,NrW_FacesXY,NrW_FacesYZ,NrW_FacesZX
  INTEGER :: NrRW_Cells,NrRW_FacesXY,NrRW_FacesYZ,NrRW_FacesZX
  INTEGER :: NrW_All_Cells,NrW_All_FXY,NrW_All_FYZ,NrW_All_FZX
  INTEGER :: NrW_Cells1=0
  INTEGER :: NrW_Cells2=0
  INTEGER :: NrMP_Cells,NrMP_FacesXY,NrMP_FacesYZ,NrMP_FacesZX,NrMP_F
  INTEGER :: NrMP_All_Cells,NrMP_All_FXY,NrMP_All_FYZ,NrMP_All_FZX
  INTEGER :: NrB_Cells,NrB_All_Cells
  INTEGER :: NrR_Cells,NrR_FacesXY,NrR_FacesYZ,NrR_FacesZX
  INTEGER :: NrRN_Cells,NrRN_FacesXY,NrRN_FacesYZ,NrRN_FacesZX
  INTEGER :: sub_cell_gr,sub_face_gr
  !INTEGER :: sub_face_gr(1:6)  !Index Grenze Face subMountain
  !Face Output Weight
  CHARACTER(6), DIMENSION(6):: FName_MP=(/"MP_FB","MP_FT",&
                                          "MP_FS","MP_FN",&
                                          "MP_FW","MP_FE"/)
  CHARACTER(6) :: CName_MP="C_MP"
  TYPE CF_MP_T
    CHARACTER(6)   :: FN_MP
    TYPE (Point_T) :: F_MP
  END TYPE CF_MP_T 
  TYPE(CF_MP_T) :: CF_MP(1:6)
  INTEGER :: nfmp

  !
  CHARACTER(8) :: OutputType
  INTEGER :: nr_wahlfkt
  INTEGER :: v_x0,v_x1,v_y0,v_y1,v_z0,v_z1,v_zw
  CHARACTER :: out_area,view_cell,view_block

  REAL(4) :: tempw
  INTEGER :: nRec
  INTEGER, PARAMETER :: SizeOfReal=4
  
  REAL(8) :: zRauh=0.01d0
  REAL(8) :: Albedo=0.0d0
  REAL(8) :: Emissivity=0.0d0

  !Input GridDistsCtrl
  REAL(8) :: dist_fscv=1.0d-12     ! for fine scaling dist to in_out-def, CheckVertex 
  REAL(8) :: distx_coeff=0.01      ! Distance Point(x) coefficient
  REAL(8) :: disty_coeff=0.01      ! Distance Point(y) coefficient
  REAL(8) :: distz_coeff=0.01      ! Distance Point(z) coefficient
  REAL(8) :: dxViewLoc=1.0d-8      ! Distance Point(x) coefficient of model border
  REAL(8) :: dyViewLoc=1.0d-8      ! Distance Point(y) coefficient of model border
  REAL(8) :: dzViewLoc=1.0d-8      ! Distance Point(z) coefficient of model border
  INTEGER :: IncrVol=10            ! counter for cuts of volume-analysis (1 für Einheitswürfel 
  REAL(8) :: dist_scMaxCell=1.0d-12 ! adjustment value to the filter(screen) of cells 
                                    ! with roughly max. Vol
  REAL(8) :: ShiftPoint=1.0d0


  INTERFACE Set
    MODULE PROCEDURE SetDomain
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE DomainAllocate
  END INTERFACE
  INTERFACE DeAllocate
    MODULE PROCEDURE DomainDeAllocate
  END INTERFACE

CONTAINS

SUBROUTINE SetDomain(Domain)

  TYPE (Domain_T) :: Domain

  nx=Domain%nx
  ny=Domain%ny
  nz=Domain%nz
  WriteOffsetC=Domain%WriteOffsetC
  WriteOffsetN=Domain%WriteOffsetN
  ix0=Domain%ix0
  ix1=Domain%ix1
  iy0=Domain%iy0
  iy1=Domain%iy1
  iz0=Domain%iz0
  iz1=Domain%iz1
  igx0=Domain%igx0
  igx1=Domain%igx1
  igy0=Domain%igy0
  igy1=Domain%igy1
  igz0=Domain%igz0
  igz1=Domain%igz1
  dx=>Domain%dx
  dy=>Domain%dy
  dz=>Domain%dz
  xP=>Domain%xP
  yP=>Domain%yP
  zP=>Domain%zP
  !.......................
  Vertices=>Domain%Vertices
  Edges_X=>Domain%Edges_X
  Edges_Y=>Domain%Edges_Y
  Edges_Z=>Domain%Edges_Z
  Faces_XY=>Domain%Faces_XY
  Faces_YZ=>Domain%Faces_YZ
  Faces_ZX=>Domain%Faces_ZX
  Cells=>Domain%Cells
  !.............................
  Topo=>Domain%Topo
  !cells
  stopo_cnr=Domain%stopo_cnr
  stopo_cv =Domain%stopo_cv
  !edges
  stopo_exnr=Domain%stopo_exnr
  stopo_eynr=Domain%stopo_eynr
  stopo_eznr=Domain%stopo_eznr
  stopo_enr =Domain%stopo_enr
  stopo_ecnr=Domain%stopo_ecnr
  topo_pos_ec1=Domain%topo_pos_ec1
  topo_pos_ecn=Domain%topo_pos_ecn
  !faces
  stopo_fxynr=Domain%stopo_fxynr
  stopo_fzxnr=Domain%stopo_fzxnr
  stopo_fyznr=Domain%stopo_fyznr
  stopo_fnr  =Domain%stopo_fnr
  stopo_fcnr =Domain%stopo_fcnr
  !...............................
  UpBoundsLayer=>Domain%UpBoundsLayer
  DimUpBounds=>Domain%DimUpBounds
  DimUpBoundsArray=>Domain%DimUpBoundsArray
  DimUpBoundsN=>Domain%DimUpBoundsN
  !cells
  sbcvcnr  =Domain%sbcvcnr
  sbchnr   =Domain%sbchnr
  sbcnr_bdr=Domain%sbcnr_bdr
  sbcnr    =Domain%sbcnr
  !edges
  sbexnr   =Domain%sbexnr
  sbeynr   =Domain%sbeynr
  sbeznr   =Domain%sbeznr
  sbecnr   =Domain%sbecnr
  upb_pos_ec1=Domain%upb_pos_ec1
  upb_pos_ecn=Domain%upb_pos_ecn
  sbenr_bdr=Domain%sbenr_bdr
  sbenr    =Domain%sbenr
  !faces
  sbfcnr   =Domain%sbfcnr
  sbfxynr  =Domain%sbfxynr
  sbfzxnr  =Domain%sbfzxnr
  sbfyznr  =Domain%sbfyznr
  sbfnr_bdr=Domain%sbfnr_bdr
  sbfnr    =Domain%sbfnr
  !........................
  nr_soildef=Domain%nr_soildef
  soil_type=>Domain%soil_type
  s_ixa=>Domain%s_ixa
  s_ixe=>Domain%s_ixe
  s_iya=>Domain%s_iya
  s_iye=>Domain%s_iye
  s_iza=>Domain%s_iza
  s_ize=>Domain%s_ize
  !........................
  NrW_FacesXY=Domain%NrW_FacesXY
  NrW_FacesYZ=Domain%NrW_FacesYZ
  NrW_FacesZX=Domain%NrW_FacesZX
  NrMP_FacesXY=Domain%NrMP_FacesXY
  NrMP_FacesYZ=Domain%NrMP_FacesYZ
  NrMP_FacesZX=Domain%NrMP_FacesZX
  NrMP_F=Domain%NrMP_F ! Allg. für Rand->NrMP_F* Splittung später 
  NrR_FacesXY=Domain%NrR_FacesXY
  NrR_FacesYZ=Domain%NrR_FacesYZ
  NrR_FacesZX=Domain%NrR_FacesZX
  NrRN_FacesXY=Domain%NrRN_FacesXY
  NrRN_FacesYZ=Domain%NrRN_FacesYZ
  NrRN_FacesZX=Domain%NrRN_FacesZX
  NrRW_FacesXY=Domain%NrRW_FacesXY
  NrRW_FacesYZ=Domain%NrRW_FacesYZ
  NrRW_FacesZX=Domain%NrRW_FacesZX
  !sub_face_gr=Domain%sub_face_gr
  !...............................
  NrW_Cells=Domain%NrW_Cells
  NrB_Cells=Domain%NrB_Cells
  NrMP_Cells=Domain%NrMP_Cells
  NrR_Cells=Domain%NrR_Cells
  NrRN_Cells=Domain%NrRN_Cells
  NrRW_Cells=Domain%NrRW_Cells
  !............................
  FU=>Domain%WeiFU
  FV=>Domain%WeiFV
  FW=>Domain%WeiFW
  VolC=>Domain%VolC
  !................................
  AnzahlNachbar=Domain%AnzahlNachbar
  Nachbars=>Domain%Nachbars
  RefLevel=Domain%RefLevel
  Refine=Domain%Refine
  RefineX=Domain%RefineX
  RefineY=Domain%RefineY
  RefineZ=Domain%RefineZ
  Boundary=Domain%Boundary
  TypeW=Domain%TypeW
  TypeE=Domain%TypeE
  TypeS=Domain%TypeS
  TypeN=Domain%TypeN
  TypeT=Domain%TypeT
  TypeB=Domain%TypeB
END SUBROUTINE SetDomain

SUBROUTINE Set_EnvBlkFace(Domain,ib)

  TYPE (Domain_T) :: Domain 
  INTEGER :: ib

  ix0=Domain%ix0
  ix1=Domain%ix1
  iy0=Domain%iy0
  iy1=Domain%iy1
  iz0=Domain%iz0
  iz1=Domain%iz1
  dx=>Domain%dx
  dy=>Domain%dy
  dz=>Domain%dz
  xP=>Domain%xP
  yP=>Domain%yP
  zP=>Domain%zP
  !.......................
  Vertices=>Domain%Vertices
  Faces_XY=>Domain%Faces_XY
  Faces_YZ=>Domain%Faces_YZ
  Faces_ZX=>Domain%Faces_ZX
  !........................
  NrW_FacesXY=Domain%NrW_FacesXY
  NrW_FacesYZ=Domain%NrW_FacesYZ
  NrW_FacesZX=Domain%NrW_FacesZX
  NrMP_FacesXY=Domain%NrMP_FacesXY
  NrMP_FacesYZ=Domain%NrMP_FacesYZ
  NrMP_FacesZX=Domain%NrMP_FacesZX
  NrMP_F=Domain%NrMP_F ! Allg. für Rand->NrMP_F* Splittung später 
  NrR_FacesXY=Domain%NrR_FacesXY
  NrR_FacesYZ=Domain%NrR_FacesYZ
  NrR_FacesZX=Domain%NrR_FacesZX
  NrRN_FacesXY=Domain%NrRN_FacesXY
  NrRN_FacesYZ=Domain%NrRN_FacesYZ
  NrRN_FacesZX=Domain%NrRN_FacesZX
  !...............................
END SUBROUTINE Set_EnvBlkFace

SUBROUTINE DomainAllocate(Domain)
   TYPE (Domain_T)           ::  &
      Domain
   ! Local variables
   INTEGER                   ::  &
      ix,iy,iz
   INTEGER :: in   
   INTEGER                   ::  & 
      stat_P, stat_d, stat_V,    & ! error status variables Points,d[x,y,z],Verts
      stat_E, stat_F, stat_C,    & ! error status variables Edges,Faces,Cells
      stat_W, stat_VolC,         & ! error status variables Wei*,VolC
      stat_Tp,                   & ! error status variables Topo 
      stat_dimB                    ! error status varaibles DimUpBoundsArray
   CHARACTER (LEN=65)        ::  &
      nerrmsg                      ! error message
   CHARACTER (LEN=15)        ::  &
      nroutine                     ! name of this subroutine

   CALL Set(Domain)

!  Points mit Rand
   ALLOCATE(Domain%xP(ix0:ix1),STAT=stat_P)
   ALLOCATE(Domain%yP(iy0:iy1),STAT=stat_P)
   ALLOCATE(Domain%zP(iz0:iz1),STAT=stat_P)
!  d() mit Rand
   ALLOCATE(Domain%dx(ix0+1:ix1),STAT=stat_d)
   ALLOCATE(Domain%dy(iy0+1:iy1),STAT=stat_d)
   ALLOCATE(Domain%dz(iz0+1:iz1),STAT=stat_d)

!  Vertices mit Rand
   ALLOCATE(Domain%Vertices(ix0:ix1,iy0:iy1,iz0:iz1),STAT=stat_V)

!  Edges ohne Rand
   ALLOCATE(Domain%Edges_X(ix0+1:ix1,iy0:iy1,iz0:iz1),STAT=stat_E)
   ALLOCATE(Domain%Edges_Y(ix0:ix1,iy0+1:iy1,iz0:iz1),STAT=stat_E)
   ALLOCATE(Domain%Edges_Z(ix0:ix1,iy0:iy1,iz0+1:iz1),STAT=stat_E)

!  Faces mit Rand
   ALLOCATE(Domain%Faces_XY(ix0+1:ix1,iy0+1:iy1,iz0:iz1),STAT=stat_F)
   ALLOCATE(Domain%Faces_ZX(ix0+1:ix1,iy0:iy1,iz0+1:iz1),STAT=stat_F)
   ALLOCATE(Domain%Faces_YZ(ix0:ix1,iy0+1:iy1,iz0+1:iz1),STAT=stat_F)

!  Cells mit Rand
   ALLOCATE(Domain%Cells(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1),STAT=stat_C)

!  Weights mit Rand
   ALLOCATE(Domain%WeiFU(ix0-1:ix1+1,iy0+1:iy1,iz0+1:iz1),STAT=stat_W) ! [x]-+
   ALLOCATE(Domain%WeiFV(ix0+1:ix1,iy0-1:iy1+1,iz0+1:iz1),STAT=stat_W) ! [y]-+
   ALLOCATE(Domain%WeiFW(ix0+1:ix1,iy0+1:iy1,iz0-1:iz1+1),STAT=stat_W) ! [z]-+

!  VolC mit Rand
   ALLOCATE(Domain%VolC(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1),STAT=stat_VolC)

!  Topography All
   ALLOCATE(Domain%Topo(ix0:ix1,iy0:iy1,iz0:iz1),STAT=stat_Tp)
   DO ix=ix0,ix1
     DO iy=iy0,iy1
       DO iz=iz0,iz1
          ALLOCATE(Domain%Topo(ix,iy,iz)%Topo)
       END DO
     END DO
   END DO
!..Allocate-status
   IF(stat_P/=0 .OR. stat_d/=0 .OR. stat_V/=0) THEN
     nroutine= 'DomainAllocate'
     nerrmsg = 'Allocate-Error: "stat_P, stat_d, stat_V"  aus DomainAllocate'
     Write(*,*) nerrmsg
     !CALL tool_break (my_cart_id, 1001, nerrmsg, nroutine)
   END IF
   IF(stat_E/=0 .OR. stat_F/=0 .OR. stat_C/=0) THEN
     nroutine= 'DomainAllocate'
     nerrmsg = 'Allocate-Error:  "stat_E, stat_F, stat_C" aus DomainAllocate'
     Write(*,*) nerrmsg
     !CALL tool_break (my_cart_id, 1002, nerrmsg, nroutine)
   END IF
   IF(stat_W/=0 .OR. stat_VolC/=0 ) THEN
     nroutine= 'DomainAllocate'
     nerrmsg = 'Allocate-Error:  "stat_W, stat_VolC"  aus DomainAllocate'
     Write(*,*) nerrmsg
     !CALL tool_break (my_cart_id, 1003, nerrmsg, nroutine)
   END IF
   IF(stat_Tp/=0 ) THEN
     nroutine= 'DomainAllocate'
     nerrmsg = 'Allocate-Error:  "stat_Tp for Array Topo"  aus DomainAllocate'
     !CALL tool_break (my_cart_id, 1004, nerrmsg, nroutine)
   END IF 
!  Upper Bounds Layer
   ALLOCATE(Domain%DimUpBoundsArray(ix0:ix1+1,iy0:iy1+1,1:2),STAT=stat_dimB)
   ALLOCATE(Domain%DimUpBoundsN(ix0:ix1,iy0:iy1,1:2))
   ALLOCATE(Domain%UpBoundsLayer(ix0-1:ix1+1,iy0-1:iy1+1,iz0:iz1))
   DO ix=ix0-1,ix1+1
     DO iy=iy0-1,iy1+1
       DO iz=iz0,iz1
          ALLOCATE(Domain%UpBoundsLayer(ix,iy,iz)%BCell)
       END DO
     END DO
   END DO

   DO in=1,AnzahlNachbar
     CALL Allocate(Nachbars(in))
   END DO
   
END SUBROUTINE DomainAllocate

SUBROUTINE DomainDeAllocate(Domain)

   TYPE (Domain_T) :: Domain

   INTEGER :: i,j,k
   CALL Set(Domain)

   DEALLOCATE(Domain%xP)
   DEALLOCATE(Domain%yP)
   DEALLOCATE(Domain%zP)

   DEALLOCATE(Domain%dx)
   DEALLOCATE(Domain%dy)
   DEALLOCATE(Domain%dz)

   DEALLOCATE(Domain%Vertices)

   DO i=ix0+1,ix1
     DO j=iy0,iy1
       DO k=iz0,iz1
         IF (ASSOCIATED(Domain%Edges_X(i,j,k)%Edge)) THEN
           DEALLOCATE(Domain%Edges_X(i,j,k)%Edge)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Edges_X)
   DO i=ix0,ix1
     DO j=iy0+1,iy1
       DO k=iz0,iz1
         IF (ASSOCIATED(Domain%Edges_Y(i,j,k)%Edge)) THEN
           DEALLOCATE(Domain%Edges_Y(i,j,k)%Edge)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Edges_Y)
   DO i=ix0,ix1
     DO j=iy0,iy1
       DO k=iz0+1,iz1
         IF (ASSOCIATED(Domain%Edges_Z(i,j,k)%Edge)) THEN
           DEALLOCATE(Domain%Edges_Z(i,j,k)%Edge)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Edges_Z)

   DO i=ix0,ix1+1
     DO j=iy0,iy1+1
       DO k=iz0-1,iz1+1
         IF (ASSOCIATED(Domain%Faces_XY(i,j,k)%Face)) THEN
           DEALLOCATE(Domain%Faces_XY(i,j,k)%Face)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Faces_XY)
   DO i=ix0,ix1+1
     DO j=iy0-1,iy1+1
       DO k=iz0,iz1+1
         IF (ASSOCIATED(Domain%Faces_ZX(i,j,k)%Face)) THEN
           DEALLOCATE(Domain%Faces_ZX(i,j,k)%Face)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Faces_ZX)
   DO i=ix0-1,ix1+1
     DO j=iy0,iy1+1
       DO k=iz0,iz1+1
         IF (ASSOCIATED(Domain%Faces_YZ(i,j,k)%Face)) THEN
           DEALLOCATE(Domain%Faces_YZ(i,j,k)%Face)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Faces_YZ)  

   DO i=ix0,ix1+1
     DO j=iy0,iy1+1
       DO k=iz0,iz1+1
         IF (ASSOCIATED(Domain%Cells(i,j,k)%Cell)) THEN
          ! CALL DEALLOC_Cell_Parts(Domain%Cell(i,j,k)%Cell)
           !If (i /= ix0 .AND. i /= ix1+1 ) THEN
           !CALL DEALLOC_Cell_CutsListP(Domain%Cell(i,j,k)%Cell)
           !CALL DEALLOC_Cell_FaceP(Domain%Cell(i,j,k)%Cell)
           !END IF
           DEALLOCATE(Domain%Cells(i,j,k)%Cell)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Cells)

   DEALLOCATE(Domain%WeiFU)
   DEALLOCATE(Domain%WeiFV)
   DEALLOCATE(Domain%WeiFW)
   DEALLOCATE(Domain%VolC)
  
   DEALLOCATE(Domain%Nachbars)

   !Topography and UpBoundsLayer
   DO i=ix0,ix1
     DO j=iy0,iy1
       DO k=iz0,iz1
         IF (ASSOCIATED(Domain%Topo(i,j,k)%Topo)) THEN
           DEALLOCATE(Domain%Topo(i,j,k)%Topo)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Topo)

   DO i=ix0-1,ix1+1
     DO j=iy0-1,iy1+1
       DO k=iz0,iz1
         IF (ASSOCIATED(Domain%UpBoundsLayer(i,j,k)%BCell)) THEN
           DEALLOCATE(Domain%UpBoundsLayer(i,j,k)%BCell)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%UpBoundsLayer)

END SUBROUTINE DomainDeAllocate

SUBROUTINE CheckWriteCellFace(Cell)
  TYPE(Cell_T), POINTER :: Cell
  nfmp=0
  CName_MP="MP_C"
  IF(Cell%Face1%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(1)
     CF_MP(nfmp)%F_MP=Cell%Face1%MidPoint
  END IF  
  IF(Cell%Face2%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(2)
     CF_MP(nfmp)%F_MP=Cell%Face2%MidPoint
  END IF  
  IF(Cell%Face3%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(3)
     CF_MP(nfmp)%F_MP=Cell%Face3%MidPoint
  END IF  
  IF(Cell%Face4%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(4)
     CF_MP(nfmp)%F_MP=Cell%Face4%MidPoint
  END IF  
  IF(Cell%Face5%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(5)
     CF_MP(nfmp)%F_MP=Cell%Face5%MidPoint
  END IF  
  IF(Cell%Face6%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(6)
     CF_MP(nfmp)%F_MP=Cell%Face6%MidPoint
  END IF  
END SUBROUTINE CheckWriteCellFace

SUBROUTINE WriteVerticesNachbar(Vertices,ib,Type,Unit)
  TYPE(VertexP_T) :: Vertices(:,:,:) 
  INTEGER :: ib
  CHARACTER*2 :: Type
  INTEGER :: Unit
  
  INTEGER :: i,j,k
  TYPE(Vertex_T), POINTER :: Vertex
  REAL(RealKind) :: Sum
  Sum=0.0d0
  DO i=LBOUND(Vertices,1),UBOUND(Vertices,1)
    DO j=LBOUND(Vertices,2),UBOUND(Vertices,2)
      DO k=LBOUND(Vertices,3),UBOUND(Vertices,3)
        Sum=Sum+Vertices(i,j,k)%Vertex%in_out
      END DO
   END DO
 END DO
 WRITE(*,*) 'Sum  Block ib Type',ib,Type,Sum
END SUBROUTINE WriteVerticesNachbar
END MODULE Domain_Mod

