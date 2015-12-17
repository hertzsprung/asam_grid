  MODULE Parametric_Mod

! Declarations:
!
! Modules used:
!------------------------------------------------------------------------------
  USE Function_Mod  , ONLY :    &
      diffgx, diffgy,    &       ! Differenz: Gauss-grid-min(x,y) and *.dat-min(x,y)
      out_surf,conv_gk,  &
      gaussx_min,gaussy_min, &
      findata_type, &             ! for SRTM_Oro-Funktion
      fconv_ingrid,fout_wahl      ! for SRTM_Oro-Funktion
      
!------------------------------------------------------------------------------
  USE Domain_Mod 
  USE Geometry_Mod
!------------------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: OutUnitCheck1=21
  INTEGER :: OutUnitCheck2=22
  INTEGER :: OutUnitCheck3=23
  !..............................................................................
  CHARACTER(5) :: indata_type
  CHARACTER(2) :: gausskrueger  ="gk"
  CHARACTER(3) :: geographical  ="geo"
  CHARACTER(3) :: radian_measure="rad"
  CHARACTER(4) :: cartesian     ="cart"
  CHARACTER(5) :: wgs84         ="wgs84"  !World Geodetic System,WGS84
  !.............................................................................. 
  CHARACTER :: out_type              ! a/b -> ascii,binary
  CHARACTER :: out_wahlgrid          ! C/G -> Cartesien,Globe
  CHARACTER :: invalue_to_out        ! convert 'rad' to input value (z.Zt.'geo') 
                                     ! for gmv-Output  
  CHARACTER(9) :: conv_ingrid        ! convert input grid -> spherical,
                                     !                    -> nonspherical (cartesian)
  CHARACTER(5) :: OutGrid='Cart'     ! Cart,Globe
  !..............................................................................

  REAL(8) :: RadOutput     !Erdradius
                           !RadOutput, value to numerator for parametrization
  REAL(8) :: ScaleRad      !ScaleRad, value to denominator for parametrization
  REAL(8) :: ScaleSoil     !ScaleSoil, value scaling dzi_soil /part, std. 1%
  INTEGER :: MoveXYGrid    !MoveXYGrid, move grid to scale-x-y, 0 .or. 1  as factor

  LOGICAL :: WNull         ! -> *.WNull
  LOGICAL :: GCut          ! -> *.Cut.out.gmvG
  LOGICAL :: GCut2         ! -> *.Cut2.out.gmvG
  LOGICAL :: GSoil         ! -> *.Soil.out.gmvG
  LOGICAL :: GONull        ! -> *.ONull.out.gmvG
  LOGICAL :: GOro          ! -> *.Oro.out.gmvG
  LOGICAL :: Bound         ! -> *.bound
  LOGICAL :: Pbound        ! -> *.pva.bound
  LOGICAL :: Pgall         ! -> *.pva.gall
  LOGICAL :: Ptropo        ! -> *.pva.tropo   
 
  NAMELIST /GridDistsCtrl/ dist_fscv, distx_coeff, disty_coeff, distz_coeff, &
                           dxViewLoc, dyViewLoc, dzViewLoc, IncrVol, &
                           dist_scMaxCell,ShiftPoint

  NAMELIST /OutGMVControl/ out_wahlgrid, out_type, conv_ingrid, invalue_to_out, &
                           RadOutput, ScaleRad, ScaleSoil, MoveXYGrid

  NAMELIST /GridFileOut/   WNull,GCut,GCut2,GSoil,GONull,GOro, &
                           Bound,Pbound,Pgall,Ptropo

CONTAINS


!..................................................................................
! Input Routine
!..................................................................................

SUBROUTINE input_GMVControl(nuin, ierrstat)
 INTEGER, INTENT (IN)      :: nuin         ! Unit number for Namelist INPUT file
 INTEGER (KIND=IntKind), INTENT (INOUT)   :: ierrstat     ! error status variable 

 INTEGER (KIND=IntKind) :: iz_err

  !- Section 1: Initialize
  iz_err=0    ! Wichtig
  ierrstat=0
 
  !- Section 2: Initialize variables with defaults
  out_wahlgrid='C'
  out_type='b'
  conv_ingrid='n'
  invalue_to_out='n'
  RadOutput=0.0d0
  ScaleRad=1.0d0 
  ScaleSoil=1.0d0
  MoveXYGrid=0
 
  READ (nuin,NML=OutGMVControl,IOSTAT=iz_err)
 ! READ (nuin,NML=OutGMVControl)
  !.........................................................................
  SELECT CASE (out_wahlgrid)
    CASE ('C')
        OutGrid="Cart"
    CASE ('G')
        OutGrid="Globe"
    CASE ('Z')
        OutGrid="Cyl"
    CASE DEFAULT
        OutGrid="Cart"
  END SELECT
  !...........................................................................
 ! IF((indata_type=='gk'.AND.conv_gk=='s').OR. &      !PointParametricEarth
 !    (indata_type=='gk'.AND.OutGrid=="Globe").OR. &
 !    indata_type=='geo'.OR.indata_type=='rad') THEN
  !...........................................................................
  
  IF(out_wahlgrid=='G')THEN          ! for Surfit-approximation
     out_surf='G'
  ELSE
     out_surf='C'
  END IF
 
  IF(conv_ingrid=='spherical') THEN  ! for Surfit-approximation
   conv_gk='s'
  ELSE
   conv_gk='n'
  END IF
  fconv_ingrid=conv_ingrid           ! interface for SRTM_Oro-Function (Wgs84)
  fout_wahl=out_wahlgrid             ! interface for " "
  findata_type=indata_type           ! interface for " "
  !...........................................................................
  IF (iz_err /= 0) THEN
     ierrstat = -1
     RETURN
  END IF

END SUBROUTINE input_GMVControl

SUBROUTINE input_dists_ctrl(nuin, ierrstat)
 INTEGER, INTENT (IN)      :: nuin         ! Unit number for Namelist INPUT file
 INTEGER (KIND=IntKind), INTENT (INOUT)   :: ierrstat     ! error status variable 

 INTEGER (KIND=IntKind) :: iz_err

  !- Section 1: Initialize
  iz_err=0    ! Wichtig
  ierrstat=0
 
  !- Section 2: Initialize variables with defaults
  !  > in Domain_Mod.f90
  !  dist_fscv=1.0d-12     ! for fine scaling dist to in_out-def, CheckVertex 
  !  distx_coeff=0.01      ! Distance Point(x) coefficient
  !  disty_coeff=0.01      ! Distance Point(y) coefficient
  !  distz_coeff=0.01      ! Distance Point(z) coefficient
  !  dxViewLoc=1.0d-8      ! Distance Point(x) coefficient of model border
  !  dyViewLoc=1.0d-8      ! Distance Point(y) coefficient of model border
  !  dzViewLoc=1.0d-8      ! Distance Point(z) coefficient of model border
  !  IncrVol=10            ! counter for cuts of volume-analysis (1 für Einheitswürfel 
  !  dist_scMaxCell=1.0d-12  ! adjustment value to the filter(screen) of cells
                             ! with roughly max. Vol

  !- Section 3: Read
  READ (nuin,NML=GridDistsCtrl,IOSTAT=iz_err)
  !.........................................................................
  IF (iz_err /= 0) THEN
     ierrstat = -1
     RETURN
  END IF

END SUBROUTINE input_dists_ctrl


SUBROUTINE input_GridFileOut(nuin, ierrstat)
 INTEGER, INTENT (IN)      :: nuin         ! Unit number for Namelist INPUT file
 INTEGER (KIND=IntKind), INTENT (INOUT)   :: ierrstat     ! error status variable 

 INTEGER (KIND=IntKind) :: iz_err
  !- Section 1: Initialize
  iz_err=0    ! Wichtig
  ierrstat=0
 
  !- Section 2: Initialize variables with defaults
   WNull  =.TRUE.      ! -> *.WNull
   GCut   =.TRUE.      ! -> *.Cut.out.gmvG
   GCut2  =.TRUE.      ! -> *.Cut2.out.gmvG
   GSoil  =.TRUE.      ! -> *.Soil.out.gmvG
   GONull =.TRUE.      ! -> *.ONull.out.gmvG
   GOro   =.TRUE.      ! -> *.Oro.out.gmvG
   Bound  =.TRUE.      ! -> *.bound
   Pbound =.TRUE.      ! -> *.pva.bound
   Pgall  =.TRUE.      ! -> *.pva.gall
   Ptropo =.TRUE.      ! -> *.pva.tropo   

  READ (nuin,NML=GridFileOut,IOSTAT=iz_err)

  IF (iz_err /= 0) THEN
     ierrstat = -1
     RETURN
  END IF

END SUBROUTINE input_GridFileOut


!.........................................................................
SUBROUTINE GKtoRad(y,x,yr,xr)
  REAL(8) :: y,x,yr,xr
  REAL(8) :: breite,laenge

  Call OGKGEO(y/1000, x/1000, breite, laenge)
  yr=(breite*4.0d0*ATAN(1.0d0))/180
  xr=(laenge*4.0d0*ATAN(1.0d0))/180
END SUBROUTINE GKtoRad

SUBROUTINE GeotoRad(breite,laenge,yr,xr)
  REAL(8) :: breite,laenge,yr,xr
  yr=(breite*4.0d0*ATAN(1.0d0))/180
  xr=(laenge*4.0d0*ATAN(1.0d0))/180
END SUBROUTINE GeotoRad

!...................................................................................
!Parametric auf Vertex_T-Struktur
!...................................................................................

FUNCTION xParametricOutGlobe(Vertex)
  REAL(8) :: xParametricOutGlobe
  TYPE(Vertex_T) :: Vertex
  xParametricOutGlobe=(Vertex%Point%z+RadOutput)*SIN(Vertex%Point%x)*COS(Vertex%Point%y)/ScaleRad
END FUNCTION xParametricOutGlobe

FUNCTION yParametricOutGlobe(Vertex)
  REAL(8) :: yParametricOutGlobe
  TYPE(Vertex_T) :: Vertex
  yParametricOutGlobe=(Vertex%Point%z+RadOutput)*COS(Vertex%Point%x)*COS(Vertex%Point%y)/ScaleRad
END FUNCTION yParametricOutGlobe

FUNCTION zParametricOutGlobe(Vertex)
  REAL(8) :: zParametricOutGlobe
  TYPE(Vertex_T) :: Vertex
  zParametricOutGlobe=(Vertex%Point%z+RadOutput)*SIN(Vertex%Point%y)/ScaleRad
END FUNCTION zParametricOutGlobe

FUNCTION xParametricOutCyl(Vertex)
  REAL(8) :: xParametricOutCyl
  TYPE(Vertex_T) :: Vertex
  xParametricOutCyl=SIN(Vertex%Point%x)*Vertex%Point%y
END FUNCTION xParametricOutCyl

FUNCTION yParametricOutCyl(Vertex)
  REAL(8) :: yParametricOutCyl
  TYPE(Vertex_T) :: Vertex
  yParametricOutCyl=COS(Vertex%Point%x)*Vertex%Point%y
END FUNCTION yParametricOutCyl
  
FUNCTION zParametricOutCyl(Vertex)
  REAL(8) :: zParametricOutCyl
  TYPE(Vertex_T) :: Vertex
  zParametricOutCyl=Vertex%Point%z
END FUNCTION zParametricOutCyl


FUNCTION xParametricOutCart(Vertex)
  REAL(8) :: xParametricOutCart
  TYPE(Vertex_T) :: Vertex
  xParametricOutCart=Vertex%Point%x
END FUNCTION xParametricOutCart

FUNCTION yParametricOutCart(Vertex)
  REAL(8) :: yParametricOutCart
  TYPE(Vertex_T) :: Vertex
  yParametricOutCart=Vertex%Point%y
END FUNCTION yParametricOutCart

FUNCTION zParametricOutCart(Vertex)
  REAL(8) :: zParametricOutCart
  TYPE(Vertex_T) :: Vertex
  zParametricOutCart=Vertex%Point%z
END FUNCTION zParametricOutCart

FUNCTION xParametricRadToGeo(Vertex)
  REAL(8) :: xParametricRadToGeo
  TYPE(Vertex_T) :: Vertex
  xParametricRadToGeo=(Vertex%Point%x*180)/(ATAN(1.0d0)*4.0d0)
END FUNCTION xParametricRadToGeo

FUNCTION yParametricRadToGeo(Vertex)
  REAL(8) :: yParametricRadToGeo
  TYPE(Vertex_T) :: Vertex
  yParametricRadToGeo=(Vertex%Point%y*180)/(ATAN(1.0d0)*4.0d0)
END FUNCTION yParametricRadToGeo

FUNCTION x_GeoToRad(x)
REAL(8) :: x_GeoToRad
REAL(8) :: x
  x_GeoToRad=(x*4.0d0*ATAN(1.0d0))/180.0d0
END FUNCTION x_GeoToRad

FUNCTION y_GeoToRad(y)
REAL(8) :: y_GeoToRad
REAL(8) :: y
  y_GeoToRad=(y*4.0d0*ATAN(1.0d0))/180.0d0
END FUNCTION y_GeoToRad

FUNCTION val_RadToGeo(value)
  REAL(8) :: val_RadToGeo
  REAL(8) :: value
  val_RadToGeo=(value*180)/(ATAN(1.0d0)*4.0d0)
END FUNCTION val_RadToGeo

FUNCTION xParametricOut(Vertex)
  REAL(8) :: xParametricOut
  TYPE(Vertex_T) :: Vertex
  REAL(8) :: re_x0
  SELECT CASE (OutGrid)
    CASE ('Globe')
      xParametricOut=xParametricOutGlobe(Vertex)
    CASE ('Cyl')
      xParametricOut=xParametricOutCyl(Vertex)
    CASE ('Cart')
      IF(conv_gk=='s') THEN
         xParametricOut=xParametricOutCart(Vertex)
      ELSE IF (invalue_to_out=='y'.AND. &
             &  (TRIM(indata_type)==geographical) )THEN
           re_x0=val_RadToGeo(Domain%x0)
           xParametricOut=xParametricRadToGeo(Vertex)-MoveXYGrid*re_x0
      ELSE
        !xParametricOut=xParametricOutCart(Vertex)-diffgx
         !xParametricOut=xParametricOutCart(Vertex)-gaussx_min
          xParametricOut=xParametricOutCart(Vertex)-MoveXYGrid*Domain%x0
      END IF
    CASE DEFAULT
      xParametricOut=Vertex%Point%x
   END SELECT
END FUNCTION xParametricOut

FUNCTION yParametricOut(Vertex)
  REAL(8) :: yParametricOut
  TYPE(Vertex_T) :: Vertex
  REAL(8) :: re_y0
  SELECT CASE (OutGrid)
    CASE ('Globe')
      yParametricOut=yParametricOutGlobe(Vertex)
    CASE ('Cyl')
      yParametricOut=yParametricOutCyl(Vertex)
    CASE ('Cart')
      IF(conv_gk=='s') THEN
        yParametricOut=yParametricOutCart(Vertex)
      ELSE IF (invalue_to_out=='y'.AND. &
             &  (TRIM(indata_type)==geographical)) THEN
           re_y0=val_RadToGeo(Domain%y0)
           yParametricOut=yParametricRadToGeo(Vertex)-MoveXYGrid*re_y0
      ELSE
      !yParametricOut=yParametricOutCart(Vertex)-diffgy
      !yParametricOut=yParametricOutCart(Vertex)-gaussy_min
      yParametricOut=yParametricOutCart(Vertex)-MoveXYGrid*Domain%y0
      END IF
    CASE DEFAULT
      yParametricOut=Vertex%Point%y
   END SELECT
END FUNCTION yParametricOut

FUNCTION zParametricOut(Vertex)
  REAL(8) :: zParametricOut
  TYPE(Vertex_T) :: Vertex
  SELECT CASE (OutGrid)
    CASE ('Globe')
      zParametricOut=zParametricOutGlobe(Vertex)
    CASE ('Cyl')
      zParametricOut=zParametricOutCyl(Vertex)
    CASE ('Cart')
      zParametricOut=zParametricOutCart(Vertex)
    CASE DEFAULT
      zParametricOut=Vertex%Point%z
   END SELECT
END FUNCTION zParametricOut

!...................................................................................
!Parametric auf Point_T-Struktur
!...................................................................................
FUNCTION xPointParametricOutGlobe(Point)
  REAL(8) :: xPointParametricOutGlobe
  TYPE(Point_T) :: Point
  xPointParametricOutGlobe=(Point%z+RadOutput)*SIN(Point%x)*COS(Point%y)/RadOutput
END FUNCTION xPointParametricOutGlobe

FUNCTION yPointParametricOutGlobe(Point)
  REAL(8) :: yPointParametricOutGlobe
  TYPE(Point_T) :: Point
  yPointParametricOutGlobe=(Point%z+RadOutput)*COS(Point%x)*COS(Point%y)/RadOutput
END FUNCTION yPointParametricOutGlobe

FUNCTION zPointParametricOutGlobe(Point)
  REAL(8) :: zPointParametricOutGlobe
  TYPE(Point_T) :: Point
  zPointParametricOutGlobe=(Point%z+RadOutput)*SIN(Point%y)/RadOutput
END FUNCTION zPointParametricOutGlobe

!...................................................................................
FUNCTION xPointParametricOutCyl(Point)
  REAL(8) :: xPointParametricOutCyl
  TYPE(Point_T) :: Point
  xPointParametricOutCyl=SIN(Point%x)*Point%y
END FUNCTION xPointParametricOutCyl

FUNCTION yPointParametricOutCyl(Point)
  REAL(8) :: yPointParametricOutCyl
  TYPE(Point_T) :: Point
  yPointParametricOutCyl=COS(Point%x)*Point%y
END FUNCTION yPointParametricOutCyl

FUNCTION zPointParametricOutCyl(Point)
  REAL(8) :: zPointParametricOutCyl
  TYPE(Point_T) :: Point
  zPointParametricOutCyl=Point%z
END FUNCTION zPointParametricOutCyl

FUNCTION xPointParametricOutCart(Point)
  REAL(8) :: xPointParametricOutCart
  TYPE(Point_T) :: Point
  xPointParametricOutCart=Point%x
END FUNCTION xPointParametricOutCart

FUNCTION yPointParametricOutCart(Point)
  REAL(8) :: yPointParametricOutCart
  TYPE(Point_T) :: Point
  yPointParametricOutCart=Point%y
END FUNCTION yPointParametricOutCart

FUNCTION zPointParametricOutCart(Point)
  REAL(8) :: zPointParametricOutCart
  TYPE(Point_T) :: Point
  zPointParametricOutCart=Point%z
END FUNCTION zPointParametricOutCart

FUNCTION xPointParametricOut(Point)
  REAL(8) :: xPointParametricOut
  TYPE(Point_T) :: Point
  SELECT CASE (OutGrid)
    CASE ('Globe')
      xPointParametricOut=xPointParametricOutGlobe(Point)
    CASE ('Cyl')
      xPointParametricOut=xPointParametricOutCyl(Point)
    CASE ('Cart')
      xPointParametricOut=xPointParametricOutCart(Point)
    CASE DEFAULT
      xPointParametricOut=Point%x
   END SELECT
END FUNCTION xPointParametricOut

FUNCTION yPointParametricOut(Point)
  REAL(8) :: yPointParametricOut
  TYPE(Point_T) :: Point
  SELECT CASE (OutGrid)
    CASE ('Globe')
      yPointParametricOut=yPointParametricOutGlobe(Point)
    CASE ('Cyl')
      yPointParametricOut=yPointParametricOutCyl(Point)
    CASE ('Cart')
      yPointParametricOut=yPointParametricOutCart(Point)
    CASE DEFAULT
      yPointParametricOut=Point%y
   END SELECT
END FUNCTION yPointParametricOut

FUNCTION zPointParametricOut(Point)
  REAL(8) :: zPointParametricOut
  TYPE(Point_T) :: Point
  SELECT CASE (OutGrid)
    CASE ('Globe')
      zPointParametricOut=zPointParametricOutGlobe(Point)
    CASE ('Cyl')
      zPointParametricOut=zPointParametricOutGlobe(Point)
    CASE ('Cart')
      zPointParametricOut=zPointParametricOutCart(Point)
    CASE DEFAULT
      zPointParametricOut=Point%z
   END SELECT
END FUNCTION zPointParametricOut

FUNCTION PointParametricEarth(Point)
   TYPE(Point_T) :: PointParametricEarth
   TYPE(Point_T) :: Point

      IF((indata_type=='gk'.AND.conv_gk=='s').OR. &
         (indata_type=='gk'.AND.OutGrid=="Globe").OR. &
         indata_type=='geo'.OR.indata_type=='rad') THEN
        PointParametricEarth%x=xParametricS(Point%x,Point%y,Point%z)
        PointParametricEarth%y=yParametricS(Point%x,Point%y,Point%z)
        PointParametricEarth%z=zParametricS(Point%x,Point%y,Point%z)
      ELSE IF(indata_type=='cyl') THEN
        PointParametricEarth%x=xParametricC(Point%x,Point%y,Point%z)
        PointParametricEarth%y=yParametricC(Point%x,Point%y,Point%z)
        PointParametricEarth%z=zParametricC(Point%x,Point%y,Point%z)
      ELSE
        PointParametricEarth%x=Point%x
        PointParametricEarth%y=Point%y
        PointParametricEarth%z=Point%z
      END IF
END FUNCTION PointParametricEarth

!..............................................................
! Spherical coordinates
!..............................................................
FUNCTION xParametricS(x,y,z)

  REAL(8) :: xParametricS
  REAL(8) :: x,y,z
  xParametricS=(z+RadEarth)*SIN(x)*COS(y)

END FUNCTION xParametricS

FUNCTION yParametricS(x,y,z)

  REAL(8) :: yParametricS
  REAL(8) :: x,y,z
  yParametricS=(z+RadEarth)*COS(x)*COS(y)

END FUNCTION yParametricS

FUNCTION zParametricS(x,y,z)

  REAL(8) :: zParametricS
  REAL(8) :: x,y,z
  zParametricS=(z+RadEarth)*SIN(y)

END FUNCTION zParametricS

!..............................................................
! Cylindrical coordinates
!..............................................................
FUNCTION xParametricC(x,y,z)
   
  REAL(8) :: xParametricC
  REAL(8) :: x,y,z
  xParametricC=SIN(x)*y

END FUNCTION xParametricC

FUNCTION yParametricC(x,y,z)

  REAL(8) :: yParametricC
  REAL(8) :: x,y,z
  yParametricC=COS(x)*y

END FUNCTION yParametricC

FUNCTION zParametricC(x,y,z)

  REAL(8) :: zParametricC
  REAL(8) :: x,y,z
  zParametricC=z

END FUNCTION zParametricC



END MODULE Parametric_Mod

