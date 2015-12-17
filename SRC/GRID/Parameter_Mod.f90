MODULE Parameter_Mod

  IMPLICIT NONE

  REAL(8) :: RadEarth  !=6371.229d3 ! Erdradius
  REAL(8) :: Pi
  REAL(8) :: Omega
  REAL(8),PARAMETER :: ScaleFactor=1.0d0 !500.0d0
CONTAINS  

SUBROUTINE ComputeParameter
  
  Pi=ATAN(1.0d0)*4.0d0
  Omega=2.0d0*Pi/(3600.0d0*24.0d0)
  RadEarth=6371.229d3/ScaleFactor
END SUBROUTINE ComputeParameter

SUBROUTINE regrot(pxreg,pyreg,pxrot,pyrot,pxcen,pycen,kcall)

!----------------------------------------------------------------------
!
!*    conversion between regular and rotated spherical coordinates.
!*
!*    pxreg     longitudes of the regular coordinates
!*    pyreg     latitudes of the regular coordinates
!*    pxrot     longitudes of the rotated coordinates
!*    pyrot     latitudes of the rotated coordinates
!*              all coordinates given in degrees n (negative for s)
!*              and degrees e (negative values for w)
!*    pxcen     regular longitude of the south pole of the rotated grid
!*    pycen     regular latitude of the south pole of the rotated grid
!*
!*    kcall=-1: find regular as functions of rotated coordinates.
!*    kcall= 1: find rotated as functions of regular coordinates.
!
!-----------------------------------------------------------------------
!
  INTEGER :: kxdim,kydim,kx,ky,kcall
  REAL(8) :: pxreg,pyreg,&
                    pxrot,pyrot,&
                    pxcen,pycen
  REAL(8) :: zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg,zcyreg, &
                    zsyrot,zcyrot,zcxrot,zsxrot,zpi,zpih
  INTEGER  :: jy,jx

  zpih = Pi*0.5d0
  zsycen = SIN((pycen+zpih))
  zcycen = COS((pycen+zpih))
!
  IF (kcall.EQ.1) THEN
     zxmxc  = pxreg - pxcen
     zsxmxc = SIN(zxmxc)
     zcxmxc = COS(zxmxc)
     zsyreg = SIN(pyreg)
     zcyreg = COS(pyreg)
     zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc
     zsyrot = max(zsyrot,-1.d0)
     zsyrot = min(zsyrot,+1.d0)
     
     pyrot = ASIN(zsyrot)
    
     zcyrot = COS(pyrot)
     zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot
     zcxrot = max(zcxrot,-1.d0)
     zcxrot = min(zcxrot,+1.d0)
     zsxrot = zcyreg*zsxmxc/zcyrot
   
     pxrot = ACOS(zcxrot)
  
     IF (zsxrot<0.0) pxrot = -pxrot
  ELSE IF (kcall.EQ.-1) THEN
     zsxrot = SIN(pxrot)
     zcxrot = COS(pxrot)
     zsyrot = SIN(pyrot)
     zcyrot = COS(pyrot)
     zsyreg = zcycen*zsyrot + zsycen*zcyrot*zcxrot
     zsyreg = max(zsyreg,-1.d0)
     zsyreg = min(zsyreg,+1.d0)
       
     pyreg = ASIN(zsyreg)
         
     zcyreg = COS(pyreg)
     zcxmxc = (zcycen*zcyrot*zcxrot -&
              zsycen*zsyrot)/zcyreg
     zcxmxc = max(zcxmxc,-1.d0)
     zcxmxc = min(zcxmxc,+1.d0)
     zsxmxc = zcyrot*zsxrot/zcyreg
     zxmxc  = ACOS(zcxmxc)
     IF (zsxmxc<0.0) zxmxc = -zxmxc
         
     pxreg = zxmxc + pxcen
          
  ELSE
     WRITE(6,'(1x,''invalid kcall in regrot'')')
     STOP
  ENDIF
END SUBROUTINE regrot

SUBROUTINE Rotate(lam,phi,rot_lam,rot_phi,rotation_angle)

  REAL(8) :: lam,phi,rot_lam,rot_phi,rotation_angle

  IF (ABS(rotation_angle)<1.0E-8) THEN
    rot_lam = lam
    rot_phi = phi
  ELSE
     CALL regrot(lam,phi,rot_lam,rot_phi,0.0d0,-0.5d0*pi+rotation_angle,1)
  ENDIF
END SUBROUTINE Rotate

END MODULE Parameter_Mod

