MODULE Function_Mod

  USE Geometry_Mod
  USE Haus_Mod
  USE Tree_Mod
  USE Domain_Mod

  IMPLICIT NONE

  INTEGER, PRIVATE :: InputUnit=20
  INTEGER, PRIVATE :: InputUnitData=21
  CHARACTER*40 :: NameFunction
  CHARACTER*40 :: FileData

  REAL(8), PRIVATE :: H
  REAL(8), PRIVATE :: a
  REAL(8), PRIVATE :: aSchaer
  REAL(8), PRIVATE :: aBell
  REAL(8), PRIVATE :: HBell
  REAL(8), PRIVATE :: aUP
  REAL(8), PRIVATE :: aDOWN
  REAL(8), PRIVATE :: xM
  REAL(8), PRIVATE :: yM
  REAL(8), PRIVATE :: Lambda
  REAL(8), PRIVATE :: xStart
  REAL(8), PRIVATE :: xEnd
  REAL(8), PRIVATE :: xL,xC,xR,hMin
  REAL(8), PRIVATE :: r1,r2
! Parameter HalfCircle
  REAL(8), PRIVATE :: RadHalfCircle
  REAL(8), PRIVATE :: xMHalfCircle
  REAL(8), PRIVATE :: zMHalfCircle
  REAL(8), PRIVATE :: xMetHalfCircle
  REAL(8), PRIVATE :: zBottomHalfCircle
! Parameter Valley
  REAL(8), PRIVATE :: hP,Vx,Sx,Px,Sy,Py
! Parameter Inlet
  REAL(8), PRIVATE :: h1Inlet,h2Inlet,h3Inlet,h4Inlet,h5Inlet
  REAL(8), PRIVATE :: x1Inlet,x2Inlet,x3Inlet,x4Inlet,x5Inlet
  REAL(8), PRIVATE :: xOut1,xOut2,xOut3
  REAL(8), PRIVATE :: hOut1,hOut2
! Parameter BaroIn
  REAL(8), PRIVATE :: RotAngle


  REAL(8), PRIVATE :: b
  REAL(8), PRIVATE :: c
  REAL(8), PRIVATE :: d
  REAL(8), PRIVATE :: alpha
  REAL(8), PRIVATE :: beta
  REAL(8), PRIVATE :: nn
  REAL(8), PRIVATE :: mu
  REAL(8), PRIVATE :: phi

  INTEGER :: NumberHaus
  TYPE(Haus_T),POINTER :: Haus(:)
  TYPE(Box_T), POINTER :: Root
  INTEGER :: Depth=8
  
  CHARACTER(1)  :: leerz =" "
  CHARACTER(11) :: leer11="           "
  CHARACTER(20) :: leer20="                    "
  CHARACTER(20) :: name_fkt
  CHARACTER(50) :: read_namefkt
  CHARACTER(50) :: file_namefkt
  INTEGER, PARAMETER :: anz_fkt=19
  INTEGER, PARAMETER :: nr_fkt_haus=10
  INTEGER, PARAMETER :: nr_fkt_gml=11
  INTEGER, PARAMETER :: nr_fkt_oro=12
  INTEGER, PARAMETER :: nr_fkt_kugel=13
  INTEGER, PARAMETER :: nr_fkt_sierra=14
  INTEGER, PARAMETER :: nr_fkt_surfit=15
  INTEGER, PARAMETER :: nr_fkt_wgs_oro=17
  INTEGER, PARAMETER :: nr_fkt_wgs_surf=18
  INTEGER, PARAMETER :: nr_fkt_utm=19
  INTEGER :: OutUnitCheck24=24
  !.......................................................................
  ! def 'Oro', 'UTM'
  REAL(8), ALLOCATABLE :: xPH(:),yPH(:),Height(:,:)
  REAL(8) :: x0H,y0H,x1H,y1H,dxH,dyH,maxHOro
  INTEGER :: nxH,nyH,nrows
  !........................................................................
  ! def 'WGS84_Oro', 'WGS84_Surf'
  REAL(8) :: grdy1,miny1,seky1
  REAL(8) :: grdy0,miny0,seky0
  REAL(8) :: grdx1,minx1,sekx1
  REAL(8) :: grdx0,minx0,sekx0
  INTEGER :: f_nlat,f_slat,f_elon,f_wlon
  INTEGER :: nxWgs,nyWgs
  CHARACTER    :: fout_wahl      ! Schnittstelle: out_wahlgrid -> C/G -> Cartesien,Globe
  CHARACTER(9) :: fconv_ingrid   ! Schnittstelle: conv_ingrid
  CHARACTER(5) :: findata_type   ! Schnittstelle: indata_type
  ! def 'Corine'
  CHARACTER(8)  :: str_corine
  CHARACTER(8)  :: logical_corine
  CHARACTER(50) :: file_ncorine
  CHARACTER(25)  :: file_clc_list="CLC_IfT.csv"    !"CLC.csv"  ! Def.: CLC1990
  LOGICAL :: ifcorine=.FALSE.
  INTEGER, PARAMETER :: Nr_Land_Cl=44 ! Anz.def.LandNutzungsKlassen, CLC1990
  INTEGER, ALLOCATABLE :: LandCover(:,:)
  TYPE rgb_T
    INTEGER :: nr_r
    INTEGER :: nr_g
    INTEGER :: nr_b
    INTEGER :: mat_rgb
  END TYPE
  TYPE Corine_T
    INTEGER :: grid_c
    INTEGER :: clc_c
    INTEGER :: mat_c          ! Erweitert für File=CLC_IfT.csv
    CHARACTER(32) :: lab1
    CHARACTER(50) :: lab2
    CHARACTER(90) :: lab3
    CHARACTER(11) :: lrgb
  END TYPE
  TYPE (Corine_T), POINTER :: v_clc(:)    ! (1:Nr_Land_Cl)
  TYPE (rgb_T),    POINTER :: v_rgb(:)    ! (1:Nr_Land_Cl)
  !........................................................................
  !def 'Kugel'
  REAL(8) :: phi0,lam0
  CHARACTER*20 :: NameFun
  REAL(8) :: RIn,ROut
  !........................................................................
  ! for sierra (approximation with curfit,splev)
  REAL(8), DIMENSION(:), ALLOCATABLE :: x, y, w      ! with 'm' allocate
  REAL(8), DIMENSION(:), ALLOCATABLE :: wrk          ! with 'lwrk' allocate
  !REAL(8), DIMENSION(:), ALLOCATABLE :: sp           ! with 'm' allocate
  INTEGER, DIMENSION(:), ALLOCATABLE :: iwrk         ! with 'nest' allocate
  REAL(8), DIMENSION (:), ALLOCATABLE :: t, b_coef   ! with 'nest' allocate
  INTEGER :: n   ! total number of knots
  INTEGER :: k   ! degree of s(x) smoothing
  REAL(8) :: s   ! specify the smoothing factor
  REAL(8) :: lo  ! lowering orographie >0m, when oscillation around z=0
  INTEGER :: iopt  ! specify whether a weighted least-square spline variation (-1)
                   !             or  a smoothing spline must be determined (0 or 1)
  INTEGER :: fst   ! factor for struct of knots
  !.........................................................................
  ! for        (approximation with surfit,bispev)
  REAL(8), DIMENSION(:), ALLOCATABLE :: z    !,x, y, w
  REAL(8), DIMENSION(:), ALLOCATABLE :: tx,ty,c_spl,wrk1,wrk2  !,iwrk
  INTEGER :: kx,ky,tx_nx,ty_ny
  INTEGER :: kwrk,lwrk2
  CHARACTER :: out_surf,conv_gk
  REAL(8) :: gaussx_min=0.0,gaussy_min=0.0
  REAL(8) :: diffgx=0.0,diffgy=0.0  ! Differenz: Gauss-grid-min(x,y) and *.dat-min(x,y)
             ! !!diffgx,diffgy also summand at cart-parametric-calculations !!!
  !.........................................................................

CONTAINS

FUNCTION Level(x,y,z)
  REAL(8) :: Level
  REAL(8) :: x,y,z

  REAL(8) :: r
  REAL(8) :: hx,hy
  REAL(8) :: fOut
  CHARACTER*20 :: Casex,Casey
  REAL(8) :: Grav
  REAL(8) :: eta,eta0,U0,eta_vertical,SurfGeoPot,SurfGeoPot0
  REAL(8) :: lam,phi,rot_lam,rot_phi
  REAL(8) :: t0,r0,r1,t,Xphi,Yphi,x0,y0,omega
  REAL(8) :: l,hh,zr
  REAL(8) :: qO,gO,RO,dO
  REAL(8) :: RS
  REAL(8) :: Pi
  REAL(8) :: d,xi,h0,as,lambdac,phic,zs,Rc,zetaC
  REAL(8) :: sin_tmp,cos_tmp
  REAL(8) :: x1,y1,xc
  REAL(8) :: HillHeight,xLoc,zb

  Pi=ATAN(1.0d0)*4.0d0
  SELECT CASE(NameFunction)
    CASE ('Flat')
      Level=z-(Domain%z0 + 1.0d-10)
      Level=z-(Domain%z0 + 100.0d0)
    CASE ('UTMFlat')
      Level=BiQuadratic(x,y,z)
    CASE ('UTM')
      Level=BiQuadratic(x,y,z)
    CASE ('UTM1')
      Level=f_bispev(x,y,z)
    CASE ('Haus') 
     Level=DistanceHaus(x,y,z)
    CASE ('HalfCircle')
      Level=(z-zMHalfCircle)*(z-zMHalfCircle) &
           +(x-xMHalfCircle)*(x-xMHalfCircle) &
                -RadHalfCircle*RadHalfCircle
      Level=MIN(Level,MAX(x-xMetHalfCircle,-x-xMetHalfCircle,z,-zBottomHalfCircle+z))
    CASE ('Agnesi2D')
      IF (x<xM) THEN
        Level=z-H/(1.0d0+((x-xM)/aDown)**2)
      ELSE
        Level=z-H/(1.0d0+((x-xM)/aUp)**2)
      END IF
    CASE ('Agnesi2DP')
      xLoc=x
      IF (xLoc-xM>Domain%x1) THEN
        xLoc=xLoc-(Domain%x1-Domain%x0)
      END IF  
      IF (xLoc<xM) THEN
        Level=z-H/(1.0d0+((xLoc-xM)/aDown)**2)
      ELSE
        Level=z-H/(1.0d0+((xLoc-xM)/aUp)**2)
      END IF
    CASE ('Agnesi3D')
      Level=z-(H/(1.0d0+((x-xM)/a)**2+ ((y-yM)/b)**2))
    CASE ('Bell')
      Level=z-(HBell/(((x**2.0d0+y**2.0d0)/aBell**2.d0+1.0d0)**1.5d0))
    CASE ('AgnesiDouble')
      Level=MAX(z-1000.0/(1.0+(x/10000.0)**2.0+(y/1.0d15)**2.0) &
             -1000.0/(1+(y/10000.0)**2.0),0.0d0)
!     Level=z-(H/(1.0+ ((x-d/2.0)/a)**2.0+ (y/b)**2.0) &
!            +H/(1.0+ ((x+d/2.0)/a)**2.0+ (y/b)**2.0))
    CASE ('SchaerHill')
      Level=z-H*EXP(-x**2/aSchaer**2)*(COS(Pi*x/Lambda))**2-1.d-4
    Case ('SchaerCos')
      IF (ABS(x) < aSchaer) THEN
        Level=z-H*COS(x*Pi/(2.0*aSchaer))**2 * COS(x*Pi/Lambda)**2
      ELSE
        Level=z
      END IF
    CASE ('Concave')
!     Aus Zängl MWR    
      IF (y>=d) THEN
        x1=x*COS(alpha)+(y-b)*SIN(alpha)
        y1=(y-b)*COS(alpha)-x*SIN(alpha)
      ELSE IF (y<=-d) THEN
        x1=x*COS(alpha)-(y+b)*SIN(alpha)
        y1=(y+b)*COS(alpha)+x*SIN(alpha)
      ELSE
        xc=(b-d/2.0d0-y**2/(2.0d0*d))*TAN(alpha)
        x1=(x-xc)*COS(alpha)
        y1=0
      END IF
      IF (ABS(y)>c) THEN
        IF (y*y1>=0) THEN
          Level=H/(1.0d0+x1**2/a**2+y1**2/a**2)
        ELSE
          Level=H/(1.0d0+x1**2/a**2)
        END IF
      ELSE
        Level=H*(1.0d0+beta*COS(Pi*y/c)**2)/(1.0d0+x1**2/a**2)
      END IF
    CASE ('Bannon2D')
!     Level=z-H*(1.0d0+mu*COS(n*Pi*x/a+phi))/(1.0d0+(x/a)**2.0d0)  !!!ORIGINAL
!     Level=z-H*(1.0d0+mu*COS(n*Pi*x/a+(phi*Pi)))/(1.0d0+(x/a)**Pot)
      Level=z-H*(1.0d0+mu*COS(nn*Pi*x/a+phi))/(1.0d0+(x/a)**2.0d0)
    CASE ('Bannon3D')
!     Level=z-H*(1.0d0+mu*COS(n*Pi*x/a+phi))/(1.0d0+(x/a)**2.0d0+(y/b)**2.0d0) !!!ORIGINAL
!     Level=z-H*(1.0d0+mu*COS(n*Pi*x/a+(phi*Pi)))/(1.0d0+(x/a)**2.0d0+(y/b)**Pot)
      Level=z-H*(1.0d0+mu*COS(nn*Pi*x/a+phi))/(1.0d0+(x/a)**2.0d0+(y/b)**2.0d0)
    CASE ('HillFroehlich')
      HillHeight=28.0d0
      xLoc=x*1000.0d0
      IF (xLoc>=4.5d0*HillHeight) THEN
        xLoc=MAX(9.0d0*HillHeight-xLoc,0.0d0)
      END IF

      IF (xLoc>=0.0d0.AND.xLoc<9.0d0) THEN
        zb=2.800000000000d+01 &
          +0.000000000000d+00*xLoc &
          +6.775070969851d-03*(xLoc**2.0) &
          -2.124527775800d-03*(xLoc**3.0)
        zb=MIN(zb,28.0d0)
      ELSE IF (xLoc>=9.0d0.AND.xLoc<14.0d0) THEN
        zb=2.507355893131d+01 &
          +9.754803562315d-01*xLoc &
          -1.016116352781d-01*(xLoc**2.0d0) &
          +1.889794677828d-03*(xLoc**3.0d0)
      ELSE IF (xLoc>=14.0d0.AND.xLoc<20.0d0) THEN
        zb=2.579601052357d+01 &
          +8.206693007457d-01*xLoc &
          -9.055370274339d-02*(xLoc**2.0d0) &
          +1.626510569859d-03*(xLoc**3.0d0)
      ELSE IF (xLoc>=20.0d0.AND.xLoc<30.0d0) THEN
        zb=4.046435022819d+01 &
          -1.379581654948d+00*xLoc &
          +1.945884504128d-02*(xLoc**2.0d0) &
          -2.070318932190d-04*(xLoc**3.0d0)
      ELSE IF (xLoc>=30.0d0.AND.xLoc<40.0d0) THEN
        zb=1.792461334664d+01  &
          +8.743920332081d-01*xLoc &
          -5.567361123058d-02*(xLoc**2.0d0) &
          +6.277731764683d-04*(xLoc**3.0d0)
      ELSE IF (xLoc>=40.0d0.AND.xLoc<54.0d0) THEN
        zb=5.639011190988d+01 &
          -2.010520359035d+00*xLoc &
          +1.644919857549d-02*(xLoc**2.0d0) &
          +2.674976141766d-05*(xLoc**3.0d0)
        zb=MAX(0.0d0,zb)
      ELSE
        zb=0.0
      END IF
      zb=zb/1000.0d0
      Level=z-zb
    CASE ('Zeppelin')
      l=0.5d0*4660.0d0
      hh=0.5d0*3330.0d0
      zr=5000.0d0
      Level=SQRT((x/l)**2.0d0+((z-zr)/hh)**2.0d0)-1.0d0
    CASE ('Cardioid')  
      Level=-(3.0d0*(x*x+y*y)-y)**2.0d0+x*x+y*y
    CASE ('Sphere')  
      RS=0.5d0
      Level=-x*x-y*y-z*z+RS*RS
    CASE ('OakAcorn')  
      qO=-6.0d0/7.0d0
      gO=0.5d0
      RO=15.0d0/7.0d0
      dO=SQRT((RO*RO-gO*gO)/(qO*qO))
      IF (z>0) THEN
        Level=-(x/dO)**2.0d0-(y/dO)**2.0d0-(z+qO)*ABS(z+qO)
      ELSE
        Level=-x*x-y*y-(z-gO)**2.0d0+RO*RO
      END IF
    CASE ('StarLi')  
      r0=0.5d0
      r1=0.2d0
      omega=5.0d0
      x0=0.2d0/SQRT(20.0d0)
      y0=0.2d0/SQRT(20.0d0)
      phi=ATAN2(y-y0,x-x0)
      IF (phi<0.0d0) phi=phi+2.0d0*Pi
      r=r0+r1*SIN(omega*phi)
      Xphi=r*COS(phi)! +x0
      Yphi=r*SIN(phi)! +y0
      Level=-SQRT((x-x0)*(x-x0)+(y-y0)*(y-y0))+SQRT(Xphi*Xphi+Yphi*Yphi)
    CASE ('Star1')
      t0=0.00132d0
      r0=0.02d0*SQRT(5.0d0)
      phi=ATAN2(y,x)
      IF (phi<0.0d0) phi=phi+2.0d0*Pi
      t=phi+t0
      Xphi=r0+(0.5d0+.2d0*SIN(5.0d0*t))*COS(t)
      Yphi=r0+(0.5d0+.2d0*SIN(5.0d0*t))*SIN(t)
      Level=-SQRT(x*x+y*y)+SQRT(Xphi*Xphi+Yphi*Yphi)
    CASE ('Star2')
      t0=0.45234d0
      t=ATAN2(y,x)
      IF (t<0.0d0) t=t+2.0d0*Pi
      phi=t+t0+SIN(4.0d0*(t+t0))
      r=.60125d0+0.24012d0*COS(4.0d0*(t+t0)+0.5d0*Pi)
      Xphi=r*COS(phi)
      Yphi=r*SIN(phi)
      Level=-SQRT(x*x+y*y)+SQRT(Xphi*Xphi+Yphi*Yphi)
    CASE ('Rhombus')
      Level=MIN(y+0.5d0,0.5d0-y)+MIN(x+0.5d0,0.5d0-x)
      Level=MIN(MIN(y+0.5d0,0.5d0-y),MIN(x+0.5d0,0.5d0-x))
    CASE ('Annulus')
      r=SQRT(x*x+y*y)
      Level=MIN(r-RIn,Rout-r)
    CASE ('Leer')
      Level=MIN(y-r1,r2-y)
    CASE ('ValleyTwo')
      IF (ABS(x)<=Vx) THEN
        hx=0.0d0
        Casex='1x'
      ELSE IF (ABS(x)<=Vx+Sx) THEN
        hx=0.5d0-0.5d0*COS(Pi*(ABS(x)-Vx)/Sx)
        Casex='2x'
      ELSE IF (ABS(x)<=Vx+Sx+Px) THEN
        hx=1.0d0
        Casex='3x'
      ELSE IF (ABS(x)<=Vx+2.0d0*Sx+Px) THEN
        hx=0.5d0+0.5d0*COS(Pi*(ABS(x)-(Vx+Sx+Px))/Sx)
        Casex='4x'
      ELSE
        hx=0.0d0
        Casex='5x'
      END IF
      IF (ABS(y)<=Py) THEN
        hy=1.0d0
        Casey='1y'
      ELSE IF (ABS(y)<=Py+Sy) THEN
        hy=0.5d0+0.5d0*COS(Pi*(ABS(y)-Py)/Sy)
        Casey='2y'
      ELSE
        hy=0.0d0
        Casey='3y'
      END IF
      Level=z-hP*hx*hy-1.d-10
    CASE ('Valley3D')
      IF (ABS(x)<=Vx) THEN
        hx=0.0d0
      ELSE IF (ABS(x)<=Vx+Sx) THEN
        hx=0.5d0-0.5d0*COS(Pi*(ABS(x)-Vx)/Sx)
      ELSE
        hx=1.0d0
      END IF
      hy=0.5d0+0.5d0*TANH(y/Sy)
      Level=z-hP*hx*hy
    CASE ('Valley2D')
      hy=0.5d0+0.5d0*TANH(y/Sy)
      Level=z-hP*hy
    CASE ('Kegel')
      IF (x<=xL) THEN
        Level=z
      ELSE IF (x<=xC.AND.xC>xL) THEN
        Level=z-MIN((x-xL)/(xC-xL)*h,hMin)
      ELSE IF (x<=xR.AND.xR>xC) THEN
        Level=z-MIN((xR-x)/(xR-xC)*h,hMin)
      ELSE
        Level=z
      END IF
    CASE ('Inlet')
      IF (x<=x1Inlet) THEN
        Level=0.0d0
        Level=z-Level
      ELSE IF (x<=x2Inlet) THEN
        Level=h2Inlet/(x2Inlet-x1Inlet)**2.0d0*(x-x1Inlet)**2.0d0
        Level=z-Level
      ELSE IF (x<=x4Inlet) THEN
        Level=h3Inlet-(h3Inlet-h2Inlet)/(x2Inlet-x3Inlet)**2.0d0*(x-x3Inlet)**2.0d0
        Level=z-Level
      ELSE
        Level=h3Inlet-(h3Inlet-h2Inlet)/(x2Inlet-x3Inlet)**2.0d0*(x4Inlet-x3Inlet)**2.0d0 &
           -h4Inlet*(x-x4Inlet)/(x5Inlet-x4Inlet)
        Level=z-Level
      END IF
      IF (x>=xOut1.AND.x<=xOut3) THEN
        fOut=hOut1
        Level=MIN(Level,fOut-z)
      END IF
      IF (x>=xOut2.AND.x<=xOut3) THEN
        fOut=hOut2
        Level=MAX(Level,z-fOut)
      END IF  
      Level=MIN(Level,h5Inlet-z)
    CASE ('BaroIn')
      Grav=9.81d0
      eta=1.0d0
      eta0=0.252d0
      U0=35.0d0
      phi=0.5d0*Pi
      SurfGeoPot0   = U0*1.5d0*(COS(eta_vertical))**1.5d0*                                          &
                     ((-2.d0*(SIN(phi))**6.0d0*((COS(phi))**2.0d0+1.d0/3.d0)+10.d0/63.d0)* &
                      U0*(COS(eta_vertical))**1.5d0  +                                             &
                     (8.d0/5.d0*(COS(phi))**3.0d0*((SIN(phi))**2.0d0+2.d0/3.d0)            &
                     - pi/4.d0)*Omega*RadEarth)
      lam=x
      phi=y
      CALL Rotate(lam,phi,rot_lam,rot_phi,RotAngle)
      eta_vertical = (eta - eta0) * 0.5d0*pi
      SurfGeoPot   = U0*1.5d0*(COS(eta_vertical))**1.5d0*                                          &
                     ((-2.d0*(SIN(rot_phi))**6.0d0*((COS(rot_phi))**2.0d0+1.d0/3.d0)+10.d0/63.d0)* &
                      U0*(COS(eta_vertical))**1.5d0  +                                             &
                     (8.d0/5.d0*(COS(rot_phi))**3.0d0*((SIN(rot_phi))**2.0d0+2.d0/3.d0)            &
                     - pi/4.d0)*Omega*RadEarth)-SurfGeoPot0
      Level=z-SurfGeoPot/Grav
    CASE ('Schar')
      lambdac=pi/4.0d0
      phic=0.0d0
      d=5000.0d0
      h0=250.0d0
      xi=4000.0d0
      lam=x
      phi=y
      sin_tmp=SIN(phi)*SIN(phic)
      cos_tmp=COS(phi)*COS(phic)
      as=RadEarth
    ! great circle distance with 'a/X'  
      r=as*ACOS(sin_tmp+cos_tmp*COS(lam-lambdac))
      zs =h0*exp(-(r**2)/(d**2))*(cos(pi*r/xi)**2)
      Level=z-zs

   CASE ('ScharCloud')
      lambdac=3.0d0*pi/2.0d0
      phic=0.0d0
      Rc=3.0d0*Pi/4.0d0
      zetaC=Pi/16.0d0
      d=5000.0d0
      h0=2000.0d0
      xi=4000.0d0
      lam=x
      phi=y
      sin_tmp=SIN(phi)*SIN(phic)
      cos_tmp=COS(phi)*COS(phic)
      as=RadEarth
    ! great circle distance with 'a/X'  
      r=ACOS(sin_tmp+cos_tmp*COS(lam-lambdac))
      IF (r<Rc) THEN
        zs=0.5d0*h0*(1.0d0+COS(Pi*r/Rc))*COS(Pi*r/zetaC)**2
      ELSE
        zs=0.0d0
      END IF
      Level=z-zs
    CASE Default
      Level=z+1.d20
  END SELECT
END FUNCTION Level

SUBROUTINE InputFunction(Filename)
  CHARACTER(*) :: FileName

  REAL(8) :: A,B,C,D
  REAL(8) :: x4Inlet1
  REAL(8) :: x4Inlet2
  CHARACTER(300) :: Line

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'#Funktion')>0) THEN
      READ(InputUnit,*) NameFunction
      SELECT CASE(NameFunction)
        CASE ('HalfCircle')
          READ(InputUnit,*) RadHalfCircle
          READ(InputUnit,*) xMHalfCircle
          READ(InputUnit,*) zMHalfCircle
          READ(InputUnit,*) xMetHalfCircle
          READ(InputUnit,*) zBottomHalfCircle
        CASE ('Agnesi2D')
          READ(InputUnit,*) H
          READ(InputUnit,*) aUP   ! Halbwertsbreite vor Berg
          READ(InputUnit,*) aDOWN ! Halbwertsbreite nach Berg
          READ(InputUnit,*) xM 
        CASE ('Agnesi2DP')
          READ(InputUnit,*) H
          READ(InputUnit,*) aUP   ! Halbwertsbreite vor Berg
          READ(InputUnit,*) aDOWN ! Halbwertsbreite nach Berg
          READ(InputUnit,*) xM 
          !Level=z-H/(1.0d0+((x-xM)/aUP)**2)
          !h/(1.0d0+((x0-xM)/aDown)**2)=h/(1.0d0+((x1-xM)/aUp)**2)
          !((x0-xM)/aDown)**2=((x1-xM)/aUp)**2
          !(xM-x0)/aDown=(x1-xM)/aUp
          !aDown=aUp*(xM-Domain%x0)/(Domain%x1-xM)
        CASE ('Agnesi3D')
          READ(InputUnit,*) H
          READ(InputUnit,*) a
          READ(InputUnit,*) b
          READ(InputUnit,*) xM 
          READ(InputUnit,*) yM 
        CASE ('Bell')
          READ(InputUnit,*) HBell
          READ(InputUnit,*) aBell
        CASE ('AgnesiDouble')
          READ(InputUnit,*) H
          READ(InputUnit,*) a
          READ(InputUnit,*) b
          READ(InputUnit,*) d  !!!Abstand Bergmaxima
        CASE ('SchaerHill')
          READ(InputUnit,*) H
          READ(InputUnit,*) aSchaer
          READ(InputUnit,*) Lambda
        CASE ('SchaerCos')
          READ(InputUnit,*) H
          READ(InputUnit,*) aSchaer
          READ(InputUnit,*) Lambda
        CASE ('Concave')
          READ(InputUnit,*) H
          READ(InputUnit,*) a
          READ(InputUnit,*) b
          READ(InputUnit,*) c
          READ(InputUnit,*) d
          READ(InputUnit,*) alpha
          READ(InputUnit,*) beta
        CASE ('Bannon2D')
          READ(InputUnit,*) H
          READ(InputUnit,*) a
          READ(InputUnit,*) nn
          READ(InputUnit,*) mu
          READ(InputUnit,*) phi
        CASE ('Bannon3D')
          READ(InputUnit,*) H
          READ(InputUnit,*) a
          READ(InputUnit,*) b
          READ(InputUnit,*) nn
          READ(InputUnit,*) mu
          READ(InputUnit,*) phi
        CASE ('Star')
        CASE ('Annulus')
          READ(InputUnit,*) ROut
          READ(InputUnit,*) RIn
        CASE ('Kegel')
          READ(InputUnit,*) H
          READ(InputUnit,*) xL
          READ(InputUnit,*) xC
          READ(InputUnit,*) xR
          READ(InputUnit,*) hMin
        CASE ('Leer')
          READ(InputUnit,*) r1
          READ(InputUnit,*) r2
        CASE ('ValleyTwo')
          WRITE(*,*) 'ValleyTwo'
          READ(InputUnit,*) hP
          READ(InputUnit,*) Vx
          READ(InputUnit,*) Sx
          READ(InputUnit,*) Px
          READ(InputUnit,*) Sy
          READ(InputUnit,*) Py
        CASE ('Valley3D')
          WRITE(*,*) 'Valley3D'
          READ(InputUnit,*) hP
          READ(InputUnit,*) Vx
          READ(InputUnit,*) Sx
          READ(InputUnit,*) Sy
        CASE ('Valley2D')
          WRITE(*,*) 'Valley2D'
          READ(InputUnit,*) hP
          READ(InputUnit,*) Sy
        CASE ('Inlet')
          WRITE(*,*) 'Inlet'
          READ(InputUnit,*) h1Inlet
          READ(InputUnit,*) h2Inlet
          READ(InputUnit,*) h3Inlet
          READ(InputUnit,*) h4Inlet
          READ(InputUnit,*) h5Inlet
          READ(InputUnit,*) x1Inlet
          READ(InputUnit,*) x2Inlet
          READ(InputUnit,*) x3Inlet
          READ(InputUnit,*) x4Inlet
          READ(InputUnit,*) x5Inlet
          READ(InputUnit,*) xOut1
          READ(InputUnit,*) xOut2
          READ(InputUnit,*) xOut3
          READ(InputUnit,*) hOut1
          READ(InputUnit,*) hOut2
!           left derivative
!           f9Prime=-2.0d0*(h3Inlet-h2Inlet)/(x2Inlet-x3Inlet)**2.0d0*(x4Inlet-x3Inlet)
!           right derivative
!           f9Prime=-h4Inlet/(x5Inlet-x4Inlet)
!           A*(x4-B)=C/(D-x4)
!           (x4-B)*(D-x4)=C/A
!           x4^2-(B+D)*x4+BD+C/A=0
!           x4=(B+D)/2+-SQRT(1/4*(B*D)^2-BD-C/A)
          A=-2.0d0*(h3Inlet-h2Inlet)/(x2Inlet-x3Inlet)**2.0d0
          B=x3Inlet
          C=-h4Inlet
          D=x5Inlet
          x4Inlet=0.5d0*(B+D)-SQRT(0.25d0*(B+D)**2.0d0-B*D-C/A)
        CASE ('BaroIn')
          READ(InputUnit,*) RotAngle
          RotAngle=RotAngle/180.0d0*Pi
        CASE ('Kugel')
          Pi=ATAN(1.0d0)*4.0d0
          READ(InputUnit,*) H
          READ(InputUnit,*) d
          READ(InputUnit,*) lam0
          READ(InputUnit,*) phi0
          READ(InputUnit,*) NameFun
          phi0=phi0/180.0d0*Pi
          lam0=lam0/180.0d0*Pi
        CASE ('Haus')  
          READ(InputUnit,*) FileData
          WRITE(*,*) 'FileData  ',FileData
          CALL ReadHaus(FileData)
        CASE ('UTMFlat')  
          READ(InputUnit,*) FileData
          CALL ReadUTMFlat(FileData)
        CASE ('UTM')  
          READ(InputUnit,*) FileData
          CALL ReadUTM(FileData)
        CASE ('UTM1')  
          READ(InputUnit,*) FileData
          CALL ReadUTM1(FileData)
        CASE ('Oro')  
          READ(InputUnit,*) FileData
          CALL ReadOro(FileData)
        CASE ('WGS84Koord')  
          READ(InputUnit,*) FileData
          CALL ReadWGS84Koord(FileData)
        CASE ('WGS84Oro')  
          READ(InputUnit,*) FileData
          CALL ReadWGS84Oro(FileData)
        CASE ('WGS84Surf')  
          READ(InputUnit,*) FileData
          CALL ReadWGS84Surf(FileData)
      END SELECT  
    END IF  
  END DO  
1 CONTINUE
  CLOSE(UNIT=InputUnit)
END SUBROUTINE InputFunction

SUBROUTINE  ReadHaus(FileName)
  CHARACTER(*) :: FileName
  INTEGER:: h,i,j,k,l,nr_p
  TYPE(Point_T) ::Pdis
  REAL(8):: disp
  !Read-Syntax:
  !numberHaus
  !numberFlächen
  !numberEckenJeFläche ..
  !array's x- y- z-KoordinatenFläche

  OPEN(UNIT=InputUnitData,FILE=TRIM(FileName),STATUS='OLD')

  READ(InputUnitData,*) NumberHaus
  ALLOCATE(Haus(1:NumberHaus))
  DO h=1,NumberHaus
    READ(InputUnitData,*) Haus(h)%Number,Haus(h)%NumberOfFaces
    ALLOCATE(Haus(h)%Faces(Haus(h)%NumberOfFaces))
    READ(InputUnitData,*) (Haus(h)%Faces(i)%NumberOfPoints, i=1,Haus(h)%NumberOfFaces)
    READ(InputUnitData,*) (Haus(h)%Faces(i)%type, i=1,Haus(h)%NumberOfFaces)
    DO i=1,Haus(h)%NumberOfFaces
      ALLOCATE(Haus(h)%Faces(i)%Points(Haus(h)%Faces(i)%NumberOfPoints))
      READ(InputUnitData,*)  (Haus(h)%Faces(i)%Points(j)%x, j=1,Haus(h)%Faces(i)%NumberOfPoints) &
                       ,(Haus(h)%Faces(i)%Points(k)%y, k=1,Haus(h)%Faces(i)%NumberOfPoints) &
                       ,(Haus(h)%Faces(i)%Points(l)%z, l=1,Haus(h)%Faces(i)%NumberOfPoints)
    END DO
    CALL NormalForm(Haus(h))
    CALL BoundingBox(Haus(h))
  END DO
  CLOSE(UNIT=InputUnitData)
  ALLOCATE(Root)
  NULLIFY(Root%List)
  Root%P0%x=Domain%x0
  Root%P0%y=Domain%y0
  Root%P0%z=Domain%z0
  Root%P1%x=Domain%x1
  Root%P1%y=Domain%y1
  Root%P1%z=Domain%z1
  Root%TypeCut='x'
  Root%Index=0
  CALL CreateTree(Root,Depth)
  DO i=1,NumberHaus
    CALL InsertHaus(Root,Haus(i),i)
  END DO
END SUBROUTINE ReadHaus

SUBROUTINE ReadUTMFlat(FileName)
  CHARACTER(*) :: FileName

  INTEGER :: i,j
  INTEGER :: Nxnew,Nynew
  REAL(8) :: Dummy,Easting,Northing,HeightLoc
  INTEGER:: InUnitOro=8
  CHARACTER(1) :: DummyRow
  REAL(8), ALLOCATABLE :: Work(:,:)
  INTEGER :: iSmooth,NumSmooth
  REAL(8) :: SFac
  IF (ALLOCATED(Height)) THEN
    RETURN
  END IF
  OPEN(UNIT=InputUnitData,FILE=TRIM(FileName),STATUS='OLD')
  READ(InputUnitData,*) nxH
  READ(InputUnitData,*) nyH
  READ(InputUnitData,*) x0H
  READ(InputUnitData,*) y0H
  READ(InputUnitData,*) dxH
  READ(InputUnitData,*) dyH
  READ(InputUnitData,*) nrows
  READ(InputUnitData,'(A1)') DummyRow
  ALLOCATE(xPH(nxH))
  ALLOCATE(yPH(nyH))
  ALLOCATE(Height(nxH,nyH))
  DO i=1,nxH
    DO j=1,nyH
      Height(i,j)=0.0
    END DO
  END DO
  DO i=1,nrows
    READ(InputUnitData,*) Easting,Northing,HeightLoc
    IF (HeightLoc>0.0) THEN
      nxNew=INT((Easting-x0H)/dxH+1.0) 
      nyNew=INT((Northing-y0H)/dyH+1.0) 
      Height(nxNew,nyNew)=1
    END IF
  END DO
  SFac=0.5d0
  NumSmooth=25
  WRITE(*,*) 'Smoothing data (smoothing factor =',SFac,', iterations: ',NumSmooth,')'
  ALLOCATE(Work(nxH,nyH))
  WRITE(*,*) 'MaxHeight',MAXVAL(Height)
  DO iSmooth=1,NumSmooth
    Work=Height
    DO j=2,nyH-1
      DO i=2,nxH-1
        Height(i,j)=Work(i,j)+SFac*(Work(i+1,j)+Work(i-1,j)+Work(i,j+1)+Work(i,j-1))
        Height(i,j)=Height(i,j)/(1.0d0+4.0d0*SFac)
      END DO  
    END DO  
    WRITE(*,*) 'MaxHeight',MAXVAL(Height)
  END DO  
  DEALLOCATE(Work)
  CLOSE(UNIT=InputUnitData)
  CLOSE(UNIT=InputUnitData)
END SUBROUTINE ReadUTMFlat

SUBROUTINE ReadUTM(FileName)
  CHARACTER(*) :: FileName

  INTEGER :: i,j
  INTEGER :: Nxnew,Nynew
  REAL(8) :: Dummy,Easting,Northing,HeightLoc
  INTEGER:: InUnitOro=8
  CHARACTER(1) :: DummyRow
  REAL(8), ALLOCATABLE :: Work(:,:)
  INTEGER :: iSmooth,NumSmooth
  REAL(8) :: SFac
  IF (ALLOCATED(Height)) THEN
    RETURN
  END IF
  OPEN(UNIT=InputUnitData,FILE=TRIM(FileName),STATUS='OLD')
  READ(InputUnitData,*) nxH
  READ(InputUnitData,*) nyH
  READ(InputUnitData,*) x0H
  READ(InputUnitData,*) y0H
  READ(InputUnitData,*) dxH
  READ(InputUnitData,*) dyH
  READ(InputUnitData,*) nrows
  READ(InputUnitData,'(A1)') DummyRow
  ALLOCATE(xPH(nxH))
  ALLOCATE(yPH(nyH))
  ALLOCATE(Height(nxH,nyH))
  DO i=1,nxH
    DO j=1,nyH
      Height(i,j)=0.0
    END DO
  END DO
  DO i=1,nrows
    READ(InputUnitData,*) Easting,Northing,HeightLoc
    IF (HeightLoc>0.0) THEN
      nxNew=INT((Easting-x0H)/dxH+1.0) 
      nyNew=INT((Northing-y0H)/dyH+1.0) 
      Height(nxNew,nyNew)=HeightLoc
    END IF
  END DO
  xPH(1)=x0H
! xPH(1)=0.0
  DO i=2,nxH
    xPH(i)=xPH(i-1)+dxH
  END DO
  yPH(1)=y0H
! yPH(1)=0.0
  DO j=2,nyH
    yPH(j)=yPH(j-1)+dyH
  END DO
  SFac=0.1d0
  NumSmooth=4
  WRITE(*,*) 'Smoothing data (smoothing factor =',SFac,', iterations: ',NumSmooth,')'
  ALLOCATE(Work(nxH,nyH))
  WRITE(*,*) 'MaxHeight',MAXVAL(Height)
  DO iSmooth=1,NumSmooth
    Work=Height
    DO j=2,nyH-1
      DO i=2,nxH-1
        Height(i,j)=Work(i,j)+SFac*(Work(i+1,j)+Work(i-1,j)+Work(i,j+1)+Work(i,j-1))
        Height(i,j)=Height(i,j)/(1.0d0+4.0d0*SFac)
      END DO  
    END DO  
    WRITE(*,*) 'MaxHeight',MAXVAL(Height)
  END DO  
  DEALLOCATE(Work)
  CLOSE(UNIT=InputUnitData)
END SUBROUTINE ReadUTM

SUBROUTINE ReadOro(FileName)
  CHARACTER(*) :: FileName

  INTEGER :: i,j
  REAL(8) :: Dummy
  INTEGER :: InUnitOro=8
  IF (ALLOCATED(Height)) THEN
    RETURN
  END IF
  OPEN(UNIT=InUnitOro,FILE=TRIM(file_namefkt),STATUS='OLD')
  READ(InUnitOro,*) nxH
  READ(InUnitOro,*) nyH
  READ(InUnitOro,*) x0H
  READ(InUnitOro,*) y0H
  READ(InUnitOro,*) dxH
  READ(InUnitOro,*) dyH
  ALLOCATE(xPH(nxH))
  ALLOCATE(yPH(nyH))
  ALLOCATE(Height(nxH,nyH))
  DO i=1,nxH
    DO j=1,nyH
      READ(InUnitOro,*) Dummy,Dummy,Height(i,j)
    END DO
  END DO
  xPH(1)=x0H
  DO i=2,nxH
    xPH(i)=xPH(i-1)+dxH
  END DO
  yPH(1)=y0H
  DO j=2,nyH
    yPH(j)=yPH(j-1)+dyH
  END DO
  CLOSE(UNIT=InUnitOro)
END SUBROUTINE ReadOro

SUBROUTINE ReadWGS84Koord(FileName)
  CHARACTER(*) :: FileName
  CHARACTER(20) :: def_str

  OPEN(UNIT=InputUnitData,FILE=TRIM(file_namefkt),STATUS='OLD')
  READ(InputUnitData,*) def_str,grdy1,miny1,seky1
  READ(InputUnitData,*) def_str,grdy0,miny0,seky0
  READ(InputUnitData,*) def_str,grdx1,minx1,sekx1
  READ(InputUnitData,*) def_str,grdx0,minx0,sekx0
  READ(InputUnitData,*) def_str,nyWgs
  READ(InputUnitData,*) def_str,nxWgs

  Domain%x0=grdx0+minx0/60+sekx0/3600
  Domain%x1=grdx1+minx1/60+sekx1/3600
  Domain%y0=grdy0+miny0/60+seky0/3600
  Domain%y1=grdy1+miny1/60+seky1/3600
  Domain%nx=nxWgs-1  !nxWGS = number value direction
  Domain%ny=nyWgs-1  !nyWGS = number value direction
  CLOSE(UNIT=InputUnitData)
END SUBROUTINE ReadWGS84Koord

SUBROUTINE ReadWGS84Oro(FileName)
  CHARACTER(*) :: FileName

  INTEGER   :: iy,ix,i,j
  REAL(8)   :: Dummy
  REAL(8)   :: grdy1,miny1,seky1,gy1H,gk_y1
  REAL(8)   :: grdy0,miny0,seky0,gy0H,gk_y0
  REAL(8)   :: grdx1,minx1,sekx1,gx1H,gk_x1,gxzw,xzw
  REAL(8)   :: grdx0,minx0,sekx0,gx0H,gk_x0
  REAL(8)   :: breite,laenge
  INTEGER   :: rows,cols
  CHARACTER (LEN=20)  :: def_str,l_str
  CHARACTER (LEN=512) :: read_line
  INTEGER   :: InputUnitData=8

  IF (ALLOCATED(Height)) THEN
    RETURN
  END IF
  OPEN(UNIT=InputUnitData,FILE=TRIM(file_namefkt),STATUS='OLD')
  IF(TRIM(findata_type)=="gk") THEN
     READ(InputUnitData,*) def_str,gk_y1
     READ(InputUnitData,*) def_str,gk_y0
     READ(InputUnitData,*) def_str,gk_x1
     READ(InputUnitData,*) def_str,gk_x0
     READ(InputUnitData,*) def_str,nyH
     READ(InputUnitData,*) def_str,nxH
     IF(fconv_ingrid=='spherical'.OR.fout_wahl=='G') THEN
        Call OGKGEO(gk_y0/1000, gk_x0/1000, breite, laenge)
        y0H=(breite*4.0d0*ATAN(1.0d0))/180
        x0H=(laenge*4.0d0*ATAN(1.0d0))/180
        Call OGKGEO(gk_y1/1000, gk_x1/1000, breite, laenge)
        y1H=(breite*4.0d0*ATAN(1.0d0))/180
        x1H=(laenge*4.0d0*ATAN(1.0d0))/180
        IF(x1H<x0H) THEN
          xzw=x1H
          x1H=x0H
          x0H=xzw
        END IF
     ELSE
        y1H=gk_y1; y0H=gk_y0
        x1H=gk_x1; x0H=gk_x0
     END IF
  ELSE
     !IF(TRIM(indata_type)=="wgs84") THEN
     IF (TRIM(findata_type)=="geo") THEN
       READ(InputUnitData,*) def_str,grdy1,miny1,seky1
       READ(InputUnitData,*) def_str,grdy0,miny0,seky0
       READ(InputUnitData,*) def_str,grdx1,minx1,sekx1
       READ(InputUnitData,*) def_str,grdx0,minx0,sekx0
       READ(InputUnitData,*) def_str,nyH
       READ(InputUnitData,*) def_str,nxH
       !wgs84-->geo
       gy1H=grdy1+miny1/60+seky1/3600
       gy0H=grdy0+miny0/60+seky0/3600
       gx1H=grdx1+minx1/60+sekx1/3600
       gx0H=grdx0+minx0/60+sekx0/3600
       IF(gx1H<gx0H) THEN
         gxzw=gx1H
         gx1H=gx0H
         gx0H=gxzw
       END IF
       !geo-->rad
       y1H=(gy1H*4.0d0*ATAN(1.0d0))/180.0d0
       y0H=(gy0H*4.0d0*ATAN(1.0d0))/180.0d0
       x1H=(gx1H*4.0d0*ATAN(1.0d0))/180.0d0
       x0H=(gx0H*4.0d0*ATAN(1.0d0))/180.0d0
     ELSE
       Write(*,*) "Funktion Oro: Spezifiziert zu 'gk', 'geo', 'wgs84'"
       STOP "Input Koordinate nicht spezifiziert zur Funktion Oro!"
     END IF
  END IF

  IF(Domain%x0<x0H.OR.Domain%x1>x1H.OR. &
     Domain%y0<y0H.OR.Domain%y1>y1H )THEN
     Write(*,*) "---------------------------------------------------------"
     Write(*,*) "Wertebereich x-y"
     Write(*,*) "****************"
     Write(*,*) "Input-WGS84 :","  x0H  = ",x0H      ,"  x1H = ",x1H
     Write(*,*) "Def.-Domain :","  x0   = ",Domain%x0,"  x1  = ",Domain%x1
     Write(*,*) "---"
     Write(*,*) "Input-WGS84 :","  y0H  = ",y0H      ,"  y1H = ",y1H
     Write(*,*) "Def.-Domain :","  y0   = ",Domain%y0,"  y1  = ",Domain%y1
     Write(*,*) "---------------------------------------------------------"
     IF(x0H>Domain%x0) THEN
         Write(*,*) "x0H-Domain%x0 = ",x0H-Domain%x0, &
             &      "!Domain%x0 nicht im Input-WGS84-Area!"
     END IF
     IF(x1H<Domain%x1) THEN
         Write(*,*) "x1H-Domain%x1 = ",x1H-Domain%x1, &
             &      "!Domain%x1 nicht im Input-WGS84-Area!"
     END IF
     IF(y0H>Domain%y0) THEN
         Write(*,*) "y0H-Domain%y0 = ",y0H-Domain%y0, &
             &      "!Domain%y0 nicht im Input-WGS84-Area!"
     END IF
     IF(y1H<Domain%y1) THEN
         Write(*,*) "Domain%y1-y1H = ",Domain%y1-y1H, &
             &      "!Domain%y1 nicht im Input-WGS84-Area!"
     END IF
     STOP "STOP Read_WGS84_Oro"
  ELSE
    dyH=(y1H-y0H)/(nyH-1)
    dxH=(x1H-x0H)/(nxH-1)
  
    ALLOCATE(xPH(nxH))
    ALLOCATE(yPH(nyH))
    ALLOCATE(Height(nxH,nyH))
    DO iy=nyH,1,-1
        READ(InputUnitData,*) (Height(ix,iy),ix=1,nxH)
    END DO
    maxHOro=MAXVAL(Height)
  
    xPH(1)=x0H
    DO i=2,nxH
      xPH(i)=xPH(i-1)+dxH
    END DO
    yPH(1)=y0H
    DO j=2,nyH
      yPH(j)=yPH(j-1)+dyH
    END DO
  END IF

  CLOSE(UNIT=InputUnitData)
END SUBROUTINE ReadWGS84Oro


SUBROUTINE Read_Corine_List
  INTEGER             :: InUnitCorine=8
  CHARACTER(100)      :: read_line
  CHARACTER (LEN=20)  :: nroutine     ! name of this subroutine
  CHARACTER (LEN=40)  :: nerrmsg      ! error message
  INTEGER             :: nstat,cstat  ! error status variables
  INTEGER             :: i

  ! Datei-Orig: CLC.csv geändert!!!, Label 1-3 teilweise Hochkomma Nachtrag !!!
  ! Datei     : CLC_IfT.csv , erweitert um Spalte mat_IfT (1-9mat)
  OPEN(UNIT=InUnitCorine,FILE=TRIM(file_clc_list),STATUS='OLD',IOSTAT=nstat)
  IF(nstat/=0) THEN
     nroutine='READ_Corine_List'
     nerrmsg = "OPENING OF FILE \'file_clc_list\' FAILED"
     Write(*,*) nerrmsg, nroutine
     !CALL tool_break (my_cart_id, 1001, nerrmsg, nroutine)
  END IF

  ALLOCATE(v_clc(1:Nr_Land_Cl),STAT=cstat)
  IF(cstat/=0) THEN
     nroutine='READ_Corine_List'
     nerrmsg = 'Allocate-Error:  "v_clc" '
     Write(*,*) nerrmsg, nroutine
     !CALL tool_break (my_cart_id, 1003, nerrmsg, nroutine)
  END IF
  READ(InUnitCorine,*) read_line
  DO i=1,Nr_Land_Cl
      READ(InUnitCorine,*) v_clc(i)
  END DO
  CLOSE(UNIT=InUnitCorine)
END SUBROUTINE Read_Corine_List

SUBROUTINE Read_Corine_Data

  REAL(8)    :: grdy1,miny1,seky1  ! read
  REAL(8)    :: grdy0,miny0,seky0  ! read
  REAL(8)    :: grdx1,minx1,sekx1  ! read
  REAL(8)    :: grdx0,minx0,sekx0  ! read
  INTEGER    :: nxLC,nyLC          ! read
  INTEGER    :: InUnitCorine=8
  CHARACTER (LEN=20)  :: def_str,l_str
  CHARACTER (LEN=512) :: read_line
  CHARACTER (LEN=20)  :: nroutine     ! name of this subroutine
  CHARACTER (LEN=40)  :: nerrmsg      ! error message
  INTEGER             :: nstat,cstat  ! error status variables
  INTEGER             :: ix,iy

  IF (ALLOCATED(LandCover)) THEN
    RETURN
  END IF
  OPEN(UNIT=InUnitCorine,FILE=TRIM(file_ncorine),STATUS='OLD',IOSTAT=nstat)
  IF(nstat/=0) THEN
     nroutine='READ_Corine_Data'
     nerrmsg = "OPENING OF FILE \'file_ncorine\' FAILED"
     Write(*,*) nerrmsg
     !CALL tool_break (my_cart_id, 1001, nerrmsg, nroutine)
  ELSE
    READ(InUnitCorine,*) read_line !to skip sth.: def_str,grdy1,miny1,seky1,l_str
    READ(InUnitCorine,*) read_line !to skip sth.: def_str,grdy0,miny0,seky0,l_str
    READ(InUnitCorine,*) read_line !to skip sth.: def_str,grdx1,minx1,sekx1,l_str
    READ(InUnitCorine,*) read_line !to skip sth.: def_str,grdx0,minx0,sekx0,l_str
    READ(InUnitCorine,*) def_str,nyLC
    READ(InUnitCorine,*) def_str,nxLC

    IF(nxH/=nxLC.OR.nyH/=nyLC)THEN
      Write(*,*) "Keine Übereinstimmung Anzahl: 'nxH'  'nyH'   aus Read_WGS84_Oro   "
      Write(*,*) "                  mit Anzahl: 'nxLC' 'nyLC'  aus Read_Corine_Data ! "
      !CALL tool_message(my_cart_id,2001,nwarning,nroutine))
    END IF

    ALLOCATE(LandCover(nxLC,nyLC),STAT=cstat)
    DO iy=nyLC,1,-1
      READ(InUnitCorine,*) (LandCover(ix,iy),ix=1,nxLC)
    END DO
    CLOSE(UNIT=InUnitCorine)
  END IF
END SUBROUTINE Read_Corine_Data


SUBROUTINE View_Def_Surf(viopt,vkx,vky,vs,vfst,vlo)
  INTEGER :: viopt,vkx,vky
  REAL(8) :: vs
  INTEGER :: vfst
  REAL(8) :: vlo
  !...
  CHARACTER (LEN=20) :: stxt_i = "(A12,A7,I12,A3,A18)"
  CHARACTER (LEN=14) :: stxt_2 = "(A12,A26,A46)"
  CHARACTER (LEN=14) :: stxt_3 = "(A12,A26,A49)"
  CHARACTER (LEN=20) :: stxt_x = "(A12,A7,I12,A3,A35)"
  CHARACTER (LEN=20) :: stxt_y = "(A12,A7,I12,A3,A35)"
  CHARACTER (LEN=22) :: stxt_s = "(A12,A7,F12.2,A3,A19)"
  CHARACTER (LEN=20) :: stxt_f = "(A12,A7,I12,A3,A29)"
  CHARACTER (LEN=22) :: stxt_l = "(A12,A7,F12.2,A3,A30)"
  CHARACTER (LEN=14) :: stxt_9 = "(A12,A22,A30)"
  !...
  Write(*,FMT=stxt_i)  leerz, "iopt = ", viopt," ","(specify whether,"
  Write(*,FMT=stxt_2)  leerz, leerz, "- weighted least-square spline variation (-1)"
  Write(*,FMT=stxt_3)  leerz, leerz, "- smoothing spline must be determined (0 or 1) )"
  Write(*,FMT=stxt_x)  leerz, "kx   = ", vkx, " ","(degree of the spline x-direction)"
  Write(*,FMT=stxt_y)  leerz, "ky   = ", vky, " ","(degree of the spline y-direction)"
  Write(*,FMT=stxt_s)  leerz, "s    = ", vs,  " ","(smoothing factor)"
  Write(*,FMT=stxt_f)  leerz, "fst  = ", vfst," ","(factor for struct of knots)"
  Write(*,FMT=stxt_l)  leerz, "lo   = ", vlo, " ","(lowering orographie in m >0,"
  Write(*,FMT=stxt_9)  leerz, leerz,              " when oscillation around z=0)"

END SUBROUTINE View_Def_Surf

SUBROUTINE View_Para_Surf(m,xb,xe,yb,ye,nxest,nyest,nmax,l1,l2,tx_nx,ty_ny)
  INTEGER :: m
  REAL(8) :: xb,xe,yb,ye
  INTEGER :: nxest,nyest,nmax
  INTEGER :: l1,l2,tx_nx,ty_ny
  !...
  CHARACTER (LEN=17) :: form_m  = "(A12,A8,I12,A47)"
  CHARACTER (LEN=19) :: form_xy = "(A12,A8,F12.2,A42)"
  CHARACTER (LEN=17) :: form_s  = "(A12,A8,I12,A35)"
  CHARACTER (LEN=17) :: form_n  = "(A12,A8,I12,A39)"
  CHARACTER (LEN=17) :: form_l  = "(A12,A8,I12,A36)"
  CHARACTER (LEN=24) :: form_t  = "(A12,A8,I12,A37)"
  CHARACTER (LEN=14) :: form_f  = "(A12,A20,A40)"
  !...
  Write(*,*)
  Write(*,*) leer11,"Computed Parameter :"
  Write(*,FMT=form_m )leerz,"m     = ", m, " (denotes the number of data points (nyH*nxH))"

  ! set up the boundaries of the approximation domain.
  Write(*,FMT=form_xy)leerz,"xb    = ",xb, " (boundaries of the approximation domain)"
  Write(*,FMT=form_xy)leerz,"xe    = ",xe, " (                ::                    )"
  Write(*,FMT=form_xy)leerz,"yb    = ",yb, " (                ::                    )"
  Write(*,FMT=form_xy)leerz,"ye    = ",ye, " (                ::                    )"
  ! set up the dimension information.
  !    -must specify an upper bound for the number of knots required
  Write(*,FMT=form_s )leerz,"nxest = ",nxest, " (upper bound for number of knots)"
  Write(*,FMT=form_s )leerz,"nyest = ",nyest, " (              ::               )"
  Write(*,FMT=form_n )leerz,"nmax  = ",nmax," (actual dimension of the array tx,ty)"
  !    -wrk-arrays
  Write(*,FMT=form_l )leerz,"lwrk1 = ",l1," (actual dimensions the wrk-arrays)"
  Write(*,FMT=form_l )leerz,"lwrk2 = ",l2," (               ::               )"
  !    -iopt >= 0 will contain the number of knots with respect to the x-/y-variable,
  !               of the spline approximation returned
  !    -iopt =-1 should be specified on entry
  Write(*,FMT=form_t )leerz,"tx_nx = ",tx_nx, " (will contain the number of knots, "
  Write(*,FMT=form_t )leerz,"ty_ny = ",ty_ny, "  with respect to the x-/y-variable,"
  Write(*,FMT=form_f )leerz,leerz,            "  of the spline approximation returned)"
 
END SUBROUTINE View_Para_Surf

SUBROUTINE View_Value_Area
     Write(*,*) "---------------------------------------------------------"
     Write(*,*) "Wertebereich x-y"
     Write(*,*) "****************"
     !Write(*,*) "WGS84->geo  :"," gx0H  = ",gx0H     ," gx1H = ",gx1H
     !Write(*,*) "WGS84->rad  :","  x0H  = ",x0H      ,"  x1H = ",x1H
     !Write(*,*) "Grid-Domain :","  x0   = ",Domain%x0,"  x1  = ",Domain%x1
     !Write(*,*) "---"
     !Write(*,*) "WGS84->rad  :","  y0H  = ",y0H      ,"  y1H = ",y1H
     !Write(*,*) "Grid-Domain :","  y0   = ",Domain%y0,"  y1  = ",Domain%y1
     !Write(*,*) "---------------------------------------------------------"

     !Write(*,*) leerz,"Wertebereich x-y:" 
     !Write(*,'(a12,a51)') leerz,"  ---> 1) WGS84->geo  2)WGS84->rad  3) Grid-Domain"
     !Write(*,'(a12,a5,a10)') leerz,"  1)"," gx0H  = ",gx0H     ," gx1H = ",gx1H 
     !Write(*,'(a12,a5,a10)') leerz,"  2)","  x0H  = ",x0H      ,"  x1H = ",x1H
     !Write(*,'(a12,a5,a10)') leerz,"  3)","  x0   = ",Domain%x0,"  x1  = ",Domain%x1
     !Write(*,'(a12,a5)') leerz," ---"
     !Write(*,'(a12,a5,a10)') leerz,"  1)"," gy0H  = ",gy0H     ," gy1H = ",gy1H
     !Write(*,'(a12,a5,a10)') leerz,"  2)","  y0H  = ",y0H      ,"  y1H = ",y1H
     !Write(*,'(a12,a5,a10)') leerz,"  3)","  y0   = ",Domain%y0,"  y1  = ",Domain%y1
     !Write(*,*) ""
END SUBROUTINE View_Value_Area

SUBROUTINE  ReadWGS84Surf(FileName)
  CHARACTER(*) :: FileName
  INTEGER :: i,j,ix,iy
  INTEGER :: iopt,m,nxest,nyest,nmax,is,lwrk1,ier
  INTEGER :: mx,my
  REAL(8) :: xb,xe,yb,ye,s,eps,fp
  REAL(8) :: xx(11),yy(11),zz(121) !Bsp: Surfit
  !
  INTEGER :: u,v,km,ne,bx,by,b1,b2
  REAL(8) :: ai,delta,ww,yp
  REAL(8) :: m_half,xmin,xmax,ymin,ymax
  REAL(8) :: diffmx,diffmy
  INTEGER :: InUnitV2=7
  REAL(8) :: grdy1,miny1,seky1,gy1H,gk_y1
  REAL(8) :: grdy0,miny0,seky0,gy0H,gk_y0
  REAL(8) :: grdx1,minx1,sekx1,gx1H,gk_x1,gxzw,xzw
  REAL(8) :: grdx0,minx0,sekx0,gx0H,gk_x0
  REAL(8) :: breite,laenge
  REAL(8) :: x_datmin,y_datmin,min_px,min_py
  CHARACTER (LEN=20)  :: def_str
  CHARACTER           :: nlat_id,slat_id,elon_id,wlon_id
  REAL(8) :: zw(1:1024,1:1024)
  INTEGER :: ii
  !....................................................
  !Initialisierung:
  lo=0.0d0
  f_nlat=1;f_slat=1;f_elon=1;f_wlon=1
  IF(ALLOCATED(c_spl)) THEN
    RETURN
  END IF
  !....................................................
  Write(*,*) ""
  Write(*,*) leer11, "Read :"
  READ(10,*) iopt   ! specify whether a weighted least-square spline variation (-1)
                    !             or  a smoothing spline must be determined (0 or 1)
  READ(10,*) kx     ! degree of the spline x-direction
  READ(10,*) ky     ! degree of the spline y-direction
  READ(10,*) s      ! smoothing factor
  READ(10,*) fst    ! factor for strut of knots
  READ(10,*) lo     ! lowering orographie > null, when oscillation (z=0)

  OPEN(UNIT=InputUnitData,FILE=TRIM(FileName),STATUS='OLD')
  IF(TRIM(findata_type)=="gk") THEN
     READ(InputUnitData,*) def_str,gk_y1
     READ(InputUnitData,*) def_str,gk_y0
     READ(InputUnitData,*) def_str,gk_x1
     READ(InputUnitData,*) def_str,gk_x0
     READ(InputUnitData,*) def_str,nyH
     READ(InputUnitData,*) def_str,nxH
     IF(fconv_ingrid=='spherical'.OR.fout_wahl=='G') THEN
        Call OGKGEO(gk_y0/1000, gk_x0/1000, breite, laenge)
        y0H=(breite*4.0d0*ATAN(1.0d0))/180
        x0H=(laenge*4.0d0*ATAN(1.0d0))/180
        Call OGKGEO(gk_y1/1000, gk_x1/1000, breite, laenge)
        y1H=(breite*4.0d0*ATAN(1.0d0))/180
        x1H=(laenge*4.0d0*ATAN(1.0d0))/180
        IF(gx1H<gx0H) THEN
          xzw=x1H
          x1H=x0H
          x0H=xzw
        END IF
     ELSE
        y1H=gk_y1; y0H=gk_y0
        x1H=gk_x1; x0H=gk_x0
     END IF
  ELSE
     !IF(TRIM(indata_type)=="wgs84") THEN
     !im grid-file geo-Angabe umgestellt
     IF (TRIM(findata_type)=="geo") THEN
       READ(InputUnitData,*) def_str,grdy1,miny1,seky1,nlat_id
       READ(InputUnitData,*) def_str,grdy0,miny0,seky0,slat_id
       READ(InputUnitData,*) def_str,grdx1,minx1,sekx1,elon_id
       READ(InputUnitData,*) def_str,grdx0,minx0,sekx0,wlon_id
       READ(InputUnitData,*) def_str,nyH
       READ(InputUnitData,*) def_str,nxH
       IF(nlat_id=='S') f_nlat=-1
       IF(slat_id=='S') f_slat=-1
       IF(elon_id=='W') f_elon=-1
       IF(wlon_id=='W') f_wlon=-1
       !wgs84-->geo
       gy1H=grdy1+miny1/60+seky1/3600
       gy0H=grdy0+miny0/60+seky0/3600
       gx1H=grdx1+minx1/60+sekx1/3600
       gx0H=grdx0+minx0/60+sekx0/3600
       !geo-->rad
       y1H=((gy1H*4.0d0*ATAN(1.0d0))/180.0d0)*f_nlat
       y0H=((gy0H*4.0d0*ATAN(1.0d0))/180.0d0)*f_slat
       x1H=((gx1H*4.0d0*ATAN(1.0d0))/180.0d0)*f_elon
       x0H=((gx0H*4.0d0*ATAN(1.0d0))/180.0d0)*f_wlon
     ELSE
       Write(*,*) "Read_WGS84_Surf() : Spezifiziert 'gk', 'geo', 'wgs84'"
       STOP "Input Koordinaten nicht spezifiziert !"
     END IF
  END IF

     Write(*,*) ""
     Write(*,*) "            ","Wertebereich x-y:"
     Write(*,*) "            ","  ---> 1) WGS84->geo  2)WGS84->rad  3) Grid-Domain"
     Write(*,*) "            ","  1)"," gx0H  = ",gx0H     ," gx1H = ",gx1H
     Write(*,*) "            ","  2)","  x0H  = ",x0H      ,"  x1H = ",x1H
     Write(*,*) "            ","  3)","  x0   = ",Domain%x0,"  x1  = ",Domain%x1
     Write(*,*) "            "," ---"
     Write(*,*) "            ","  1)"," gy0H  = ",gy0H     ," gy1H = ",gy1H
     Write(*,*) "            ","  2)","  y0H  = ",y0H      ,"  y1H = ",y1H
     Write(*,*) "            ","  3)","  y0   = ",Domain%y0,"  y1  = ",Domain%y1
     Write(*,*) ""
     
  IF(Domain%x0<x0H.OR.Domain%x1>x1H.OR. &
     Domain%y0<y0H.OR.Domain%y1>y1H )THEN
     IF(x0H>Domain%x0) THEN
         Write(*,*) "x0H-Domain%x0 = ",x0H-Domain%x0, &
             &      "!Domain%x0 nicht im Input-WGS84-Area!"
     END IF
     IF(x1H<Domain%x1) THEN
         Write(*,*) "x1H-Domain%x1 = ",x1H-Domain%x1, &
             &      "!Domain%x1 nicht im Input-WGS84-Area!"
     END IF
     IF(y0H>Domain%y0) THEN
         Write(*,*) "y0H-Domain%y0 = ",y0H-Domain%y0, &
             &      "!Domain%y0 nicht im Input-WGS84-Area!"
    END IF
     IF(y1H<Domain%y1) THEN
         Write(*,*) "Domain%y1-y1H = ",Domain%y1-y1H, &
             &      "!Domain%y1 nicht im Input-WGS84-Area!"
     END IF
    STOP "STOP Read_WGS84_Surf"
 ELSE

  ! m denotes the number of data points
  ! m>=(kx+1)*(ky+1) !!!Wichtig!!!
      !read(InputUnitData,*) m
      m=nyH*nxH
      ALLOCATE(x(1:m))
      ALLOCATE(y(1:m))
      ALLOCATE(z(1:m))
      ALLOCATE(w(1:m))

  ! fetch the co-ordinate and function values of each data point.
      !In Datei WGS84-Daten entspricht 1.Zeile/1.Zahl der Höhendaten
      !die linke obere Ecke (north/west) des Koordinatensystems
      !for INPUT east (W) . west (W)  Bsp kapverden.dat-Datensatz
      ! V1 Zuordnung: Z(i) --> iy=(nyH-1)-0, ix=1-nxH, linke_obere Ecke
      DO  iy=nyH,1,-1
         i=(iy-1)*nxH
         read(InputUnitData,*) (z(i+ix),ix=1,nxH)
      END DO  !--> z1

      maxHOro=MAXVAL(z)
      Write(*,*) leer11, "z-Max-Orographie =",maxHOro
      Write(*,*) ""
      dyH=(y1H-y0H)/(nyH-1)
      dxH=(x1H-x0H)/(nxH-1)
      ! nur Zwischenergebnis
      ! ......................
      ALLOCATE(xPH(nxH))
      ALLOCATE(yPH(nyH))
      xPH(1)=x0H
      DO i=2,nxH
        xPH(i)=xPH(i-1)+dxH
      END DO
      yPH(1)=y0H
      DO j=2,nyH
       yPH(j)=yPH(j-1)+dyH
      END DO
      ! ......................
      DO iy=1,nyH
        DO  ix=1,nxH
         i=(iy-1)*nxH+ix
         x(i)=x0H+dxH*(ix-1)
        END DO
      END DO
      DO iy=1,nyH
        yp=y0H+(dyH*(iy-1))
        DO  ix=1,nxH
         i=(iy-1)*nxH+ix
         y(i)=yp
        END DO
      END DO
      x_datmin=x(1)
      y_datmin=y(1)
 
  ! fetch an estimate of the standard deviation of the data values.
      !read(InputUnitData,*) delta
      delta=1.0d-12
  ! the weights are set equal to delta**(-1)
      ww= 97.49  !ww = 1./delta
      DO  i=1,m
        w(i) = ww
      END DO

  ! set up the boundaries of the approximation domain.
      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
      DO i=2,m
        xmin=MIN(xmin,x(i))
        xmax=MAX(xmax,x(i))
        ymin=MIN(ymin,y(i))
        ymax=MAX(ymax,y(i))
      END DO
      diffmx=xmax-xmin
      diffmy=ymax-ymin
      IF(conv_gk=='s'.OR.out_surf=='G') THEN
        xb=xmin-0.01*diffmx
        xe=xmax+0.01*diffmx
        yb=ymin-0.01*diffmy
        ye=ymax+0.01*diffmy
      ELSE IF (TRIM(findata_type)=="geo") THEN
        xb=xmin-1.0d-12*diffmx
        xe=xmax+1.0d-12*diffmx
        yb=ymin-1.0d-12*diffmy
        ye=ymax+1.0d-12*diffmy
      ELSE
        xb=REAL(FLOOR(xmin))   ! nicht fuer Bogenmass
        xe=REAL(CEILING(xmax)) ! nicht fuer Bogenmass
        yb=REAL(FLOOR(ymin))   ! nicht fuer Bogenmass
        ye=REAL(CEILING(ymax)) ! nicht fuer Bogenmass
      END IF

  ! integer flag. on entry iopt must specify whether
  !     - a weighted least-squares spline (iopt=-1)
  !     - a smoothing spline (iopt= 0 or 1) must be determined.
     ! iopt = 0    ! extern read
     ! kx = 3      ! extern read
     ! ky = 3      ! extern read
     ! s = 500.0   ! extern read ,s= 30.0-900000.0
     CALL View_Def_Surf(iopt,kx,ky,s,fst,lo)

  ! set up the dimension information.
  !    -must specify an upper bound for the number of knots required
  !     in the x- and y-directions respect.
  !    -nxest >= 2*(kx+1), nyest >= 2*(ky+1)
  !    -in most practical situation
  !     nxest = kx+1+sqrt(m/2), nyest = ky+1+sqrt(m/2) will be sufficient.
      m_half=MIN(m/2,5000)
      m_half=fst
      nxest = kx+1+CEILING(SQRT(m_half))
      nyest = ky+1+CEILING(SQRT(m_half))
      IF (nxest<(2*kx+2)) THEN
        nxest=2*kx+2
      END IF
      IF (nyest<(2*ky+2)) THEN
        nyest=2*ky+2
      END IF
      nmax = MAX(nxest,nyest)
      ALLOCATE(tx(1:nmax))
      ALLOCATE(ty(1:nmax))
      ALLOCATE(c_spl((nxest-kx-1)*(nyest-ky-1)))
      kwrk = m+(nxest-2*kx-1)*(nyest-2*ky-1)
      ALLOCATE(iwrk(1:kwrk))
      u = nxest-kx-1
      v = nyest-ky-1
      km = max(kx,ky)+1
      ne = nmax
      bx = kx*v+ky+1
      by = ky*u+kx+1
      if(bx.le.by) THEN
           b1 = bx
           b2 = b1+v-ky
      END IF
      if(bx.gt.by) THEN
           b1 = by
           b2 = b1+u-kx
      END IF
      lwrk1=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
      ALLOCATE(wrk1(1:lwrk1))
      lwrk2=u*v*(b2+1)+b2
      ALLOCATE(wrk2(1:lwrk2))
 
  ! choose a value for eps
      eps=0.1e-05

  ! integer flag. on entry iopt must specify whether a weighted
  ! least-squares spline (iopt=-1)
 
      IF(iopt==-1) THEN
        tx_nx=11
        ty_ny=11
        j=kx+2
        DO i=1,3
          ai = i-2
          tx(j)=ai
          ty(j)=ai
          j=j+1
        END DO
      ELSE
        tx_nx=nmax
        ty_ny=nmax
      END IF
  !   CALL View_Para_Surf(m,xb,xe,yb,ye,nxest,nyest,nmax,lwrk1,lwrk2,tx_nx,ty_ny)

  ! spline approximations of degree k
      Write(*,*) leer11,"running spline approximation...."
      call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
        nmax,eps,tx_nx,tx,ty_ny,ty,c_spl,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)

      Write(*,*) leer11,"ende spline approximation  ier = ",ier
  END IF

  CLOSE(UNIT=InputUnitData)
END SUBROUTINE ReadWGS84Surf

FUNCTION DistanceHaus(x,y,z)
  REAL(8) :: DistanceHaus
  REAL(8) :: x,y,z
  INTEGER :: i,h
  TYPE(Point_T) ::Pdis
  TYPE(Box_T), POINTER :: Box
  REAL(8):: disp

  Pdis%x=x
  Pdis%y=y
  Pdis%z=z
  DistanceHaus=+1.0d99
  CALL FindBoxPoint(Root,Box,PDis)
  !Distance
  IF (z==0.0d0) THEN
    DistanceHaus=-1.d-2
  END IF
  DO
    IF (ASSOCIATED(Box)) THEN
      IF (ASSOCIATED(Box%List)) THEN
        DO i=1,SIZE(Box%List)
          disp=Dist(Pdis,Haus(Box%List(i)))
          DistanceHaus=MIN(DistanceHaus,disp)
        END DO
      END IF
      Box=>Box%Parent
    ELSE
      EXIT
    END IF
  END DO
  DistanceHaus=DistanceHaus-1.d-2
END FUNCTION DistanceHaus

FUNCTION BiQuadratic(x,y,z)
  REAL(8) :: BiQuadratic
  REAL(8) :: x,y,z

  INTEGER :: iPos,jPos
  REAL(8) :: x1,x2,y1,y2
  REAL(8) :: h11,h12,h21,h22
  REAL(8) :: Dist

  iPos=(x-x0H)/dxH+1
  jPos=(y-y0H)/dyH+1

  x1=xPH(iPos)
  x2=xPH(iPos+1)
  y1=yPH(jPos)
  y2=yPH(jPos+1)
  h11=Height(iPos,jPos)
  h21=Height(iPos+1,jPos)
  h12=Height(iPos,jPos+1)
  h22=Height(iPos+1,jPos+1)
  IF (iPos==nxH) THEN
    x2=xPH(iPos)+dxH/2
    h21=Height(iPos,jPos)
    IF (jPos==nyH) THEN
      h22=Height(iPos,jPos)
    ELSE
      h22=Height(iPos,jPos+1)
    END IF
  END IF
  IF (jPos==nyH) THEN
    y2=yPH(jPos)+dyH/2
    h12=Height(iPos,jPos)
    IF (iPos/=nxH)  h22=Height(iPos+1,jPos)
  END IF

  Dist=(x2-x1)*(y2-y1)
  BiQuadratic = &
        (x2-x)*(y2-y)*h11 &
       -(x1-x)*(y2-y)*h21 &
       -(x2-x)*(y1-y)*h12 &
       +(x1-x)*(y1-y)*h22
  BiQuadratic=z-(BiQuadratic/Dist+1.d-10)
END FUNCTION BiQuadratic


FUNCTION fOroWgs(x,y,z)
  REAL(8) :: fOroWgs
  REAL(8) :: x,y,z

  INTEGER :: iPos,jPos
  REAL(8) :: x1,x2,y1,y2
  REAL(8) :: h11,h12,h21,h22
  REAL(8) :: Dist
  !search pos into Height
  !Write(*,*)"x0H=", x0H, " x=",x  !8942-WGSOro-Test Kapverde 75x75x100
  !Write(*,*)"y0H=", y0H, " y=",y
  
  IF(INT(x-x0H)==0) THEN
    iPos=1
    !Write(*,*)"1", " iPos=",iPos
  ELSE
    iPos=INT((x-x0H)/dxH)
    !Write(*,*)"2", " iPos=",iPos
  END IF
  IF(INT(y-y0H)==0)THEN
    jPos=1
    !Write(*,*)"1", " jPos=",jPos
  ELSE
    jPos=INT((y-y0H)/dyH)
    !Write(*,*)"2", " jPos=",jPos
  END IF
  !iPos=(x-x0H)/dxH+1             !8942 kommentiert, Kapverde 75x75x100
  !  Write(*,*)"3", " iPos=",iPos
  !jPos=(y-y0H)/dyH+1
  !  Write(*,*)"3", " jPos=",jPos

  x1=xPH(iPos)
  x2=xPH(iPos+1)
  y1=yPH(jPos)
  y2=yPH(jPos+1)
  h11=Height(iPos,jPos)
  h21=Height(iPos+1,jPos)
  h12=Height(iPos,jPos+1)
  h22=Height(iPos+1,jPos+1)
  !wenn iPos,jPos Grenze: +dxH/2 od. +dyH/2
  ! ---> für Skalierung damit Dist != 0 werden kann
  !wenn im grid-File nx,ny Grenzen '0 bis n(x,y)WGS-1' angegeben wird
  !-->  'IF' nicht benutzt
  IF(iPos==nxH) THEN
     x2=xPH(iPos)+dxH/2
     h21=Height(iPos,jPos)
     IF(jPos==nyH) THEN
        h22=Height(iPos,jPos)
     ELSE
        h22=Height(iPos,jPos+1)
     END IF
  END IF
  IF(jPos==nyH) THEN
     y2=yPH(jPos)+dyH/2
     h12=Height(iPos,jPos)
     IF(iPos/=nxH)  h22=Height(iPos+1,jPos)
  END IF

  ! 2-dimensionale Interpolation
  Dist=(x2-x1)*(y2-y1)
  fOroWgs = &
        (x2-x)*(y2-y)*h11 &
       -(x1-x)*(y2-y)*h21 &
       -(x2-x)*(y1-y)*h12 &
       +(x1-x)*(y1-y)*h22
  fOroWgs = z-fOroWgs/Dist

END FUNCTION fOroWgs


FUNCTION f_r(lambda0,lambda1,phi0,phi1)
    REAL(8) :: f_r
    REAL(8) :: lambda0,lambda1,phi0,phi1

    f_r = RadEarth*ACOS(SIN(phi0)*SIN(phi1)+COS(phi0)*COS(phi1)*COS(lambda1-lambda0))
END FUNCTION f_r


FUNCTION f_h(lam,phi,z)  !shape of mountain
     REAL(8) :: f_h
     REAL(8) :: lam,phi,z
     REAL(8) :: r,PiNinth
     !h0, height at the center of the mountain
     !d, half-width of the mountain
     !r, distance from the center

     IF (TRIM(NameFun)=='Agnesi') THEN
       r=f_r(lam0,lam,phi0,phi)
       f_h = z-H/(1.0d0+(r/d)**2) !h-h0
     ELSE IF (TRIM(NameFun)=='Cone') THEN
       !PiNinth=ATAN(1.0d0)*4.0d0/9.0d0
       PiNinth=Pi/9.0d0
       r=MIN(SQRT((lam-lam0)**2+(phi-phi0)**2),PiNinth)
       f_h=z-H*(PiNinth-r)/PiNinth
     ELSE
       STOP 'False NameFun'
     END IF
END FUNCTION f_h


FUNCTION f_splev(x,y,z)
   REAL(8)  :: f_splev
   REAL(8)  :: x,y,z
   INTEGER, PARAMETER :: m=1
   INTEGER  :: ier
   REAL(8)  :: xs(m),sp(m)
   xs(1)=x
   call splev(t,n,b_coef,k,xs,sp,m,ier)
   f_splev=z-sp(m)
END FUNCTION f_splev

FUNCTION f_bispev(x,y,z)
   REAL(8)  :: f_bispev
   REAL(8)  :: x,y,z
   INTEGER, PARAMETER :: m=1
   INTEGER  :: ier
   REAL(8)  :: xg(1),yg(1),b_sp(1)
   xg(1)=x
   yg(1)=y
   call bispev(tx,tx_nx,ty,ty_ny,c_spl,kx,ky,xg,m,yg,m,b_sp, &
        wrk2,lwrk2,iwrk,kwrk,ier)
   b_sp(m)=b_sp(m)-10.0d0
!  f_bispev=z-b_sp(m)
   f_bispev=z-MAX(b_sp(m),0.0d0)
END FUNCTION f_bispev

SUBROUTINE  ReadUTM1(FileName)
  CHARACTER(*) :: FileName

  REAL(8) :: xx(11),yy(11),zz(121) !Bsp:
  REAL(8) :: x_datmin,y_datmin,min_px,min_py
  INTEGER :: i,ier,iopt,is,j,lwrk1,m,mx,my,nc
  INTEGER :: nmax,nxest,nyest
  INTEGER :: u,v,km,ne,bx,by,b1,b2
  REAL(8) :: ai,delta,eps,fp,s,ww,xb,xe,yb,ye
  REAL(8) :: m_half,xmin,xmax,ymin,ymax,breite,laenge
  REAL(8) :: diffmx,diffmy
  REAL(8) :: x0H,y0H,dxH,dyH
  CHARACTER*1 :: DummyRow
  INTEGER(16) :: lwrk1_long
  REAL(8) :: zMin
  INTEGER :: iFlat

  READ(InputUnit,*) iopt   ! specify whether a weighted least-square spline variation (-1)
                    !             or  a smoothing spline must be determined (0 or 1)
  READ(InputUnit,*) kx     ! degree of the spline x-direction
  READ(InputUnit,*) ky     ! degree of the spline y-direction
  READ(InputUnit,*) s      ! smoothing factor
  READ(InputUnit,*) fst    ! factor for strut of knots
  READ(InputUnit,*) iFlat  ! 0 = flat island, 1 = original topography

  OPEN(UNIT=InputUnitData,FILE=TRIM(FileName),STATUS='OLD')
  READ(InputUnitData,*) nxH
  READ(InputUnitData,*) nyH
  READ(InputUnitData,*) x0H
  READ(InputUnitData,*) y0H
  READ(InputUnitData,*) dxH
  READ(InputUnitData,*) dyH
  READ(InputUnitData,*) m
  READ(InputUnitData,'(A1)') DummyRow
  ALLOCATE(x(1:m))
  ALLOCATE(y(1:m))
  ALLOCATE(z(1:m))
  ALLOCATE(w(1:m))
  zMin=0.0d0
  IF (iFlat==0) WRITE(*,*) 'Flat island!'
  DO i=1,m
    READ(InputUnitData,*) x(i),y(i),z(i)
    IF (iFlat==0) THEN 
      z(i)=1.0d0
    END IF
    IF (z(i)>0.0d0) THEN
      z(i)=z(i)+10.0d0
    END IF
    zMin=MIN(zMin,z(i))
  END DO
  WRITE(*,*) 'zMin',zMin
  delta=1.0d0

! the weights are set equal to delta**(-1)
   ww=1./delta
   DO i=1,m
     w(i)=ww
   END DO

! set up the boundaries of the approximation domain.
  xmin=x(1)
  xmax=x(1)
  ymin=y(1)
  ymax=y(1)
  DO i=2,m
    xmin=MIN(xmin,x(i))
    xmax=MAX(xmax,x(i))
    ymin=MIN(ymin,y(i))
    ymax=MAX(ymax,y(i))
  END DO
  diffmx=xmax-xmin
  diffmy=ymax-ymin
  xb=xmin-0.01*diffmx
  xe=xmax+0.01*diffmx
  yb=ymin-0.01*diffmy
  ye=ymax+0.01*diffmy
 
! integer flag. on entry iopt must specify whether a weighted
! least-squares spline (iopt=-1) or a smoothing spline (iopt=
! 0 or 1) must be determined.
     ! iopt = 0    ! read
     ! kx = 3      ! read
     ! ky = 3      ! read
     ! s = 900000.   !s= 30 !read


! set up the dimension information.
!    -must specify an upper bound for the number of knots required
!     in the x- and y-directions respect.
!    -nxest >= 2*(kx+1), nyest >= 2*(ky+1)
!    -in most practical situation
!     nxest = kx+1+sqrt(m/2), nyest = ky+1+sqrt(m/2) will be sufficient.
  m_half=MIN(m/2,fst)
  WRITE(*,*) 'm_half,m,fst',m_half,m,fst
  nxest=kx+1+CEILING(SQRT(m_half))
  nyest=ky+1+CEILING(SQRT(m_half))
  WRITE(*,*) 'nxest',nxest
  WRITE(*,*) 'nyest',nyest
  IF (nxest<(2*kx+2)) THEN
    nxest=2*kx+2
  END IF
  IF (nyest<(2*ky+2)) THEN
    nyest=2*ky+2
  END IF
  nmax = MAX(nxest,nyest)
  ALLOCATE(tx(1:nmax))
  ALLOCATE(ty(1:nmax))
  ALLOCATE(c_spl((nxest-kx-1)*(nyest-ky-1)))
  kwrk=m+(nxest-2*kx-1)*(nyest-2*ky-1)
  ALLOCATE(iwrk(1:kwrk))
! computation for 'lwrk1'-> must specify the actual dimension of wrk1
  u=nxest-kx-1
  v=nyest-ky-1
  km=MAX(kx,ky)+1
  ne=nmax
  bx=kx*v+ky+1
  by=ky*u+kx+1
  IF (bx<=by) THEN
    b1=bx
    b2=b1+v-ky
  END IF
  IF (bx>by) THEN
    b1=by
    b2=b1+u-kx
  END IF
  lwrk1_long=16222049172_16
  lwrk1=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
  lwrk1_long=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
  ALLOCATE(wrk1(1:lwrk1))
  lwrk2=u*v*(b2+1)+b2
  ALLOCATE(wrk2(1:lwrk2))
! im Bsp surfit def. : nxest = 15,nyest = 15,nmax = 15
!                     !kwrk = 300, lwrk1 = 12000, lwrk2 = 6000
 
! choose a value for eps
  eps=0.1d-01

! integer flag. on entry iopt must specify whether a weighted
! least-squares spline (iopt=-1)
 
  IF (iopt==-1) THEN
!if the computation mode iopt=-1 is used, the values tx(kx+2),
!    c          ...tx(nx-kx-1) must be supplied by the user, before entry.
!    c          see also the restrictions (ier=10).
!if iopt=-1: 2*kx+2<=nx<=nxest
!    c                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
!    c                        2*ky+2<=ny<=nyest
!    c                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
!kx = 3
!ky = 3
    tx_nx=nmax
    ty_ny=nmax
    ai=(xe-xb)/(tx_nx-kx-1-kx-2+2)
    tx(kx+2)=xb+ai
    DO i=kx+3,tx_nx-kx-1
      tx(i)=tx(i-1)+ai
    END DO
    ai=(ye-yb)/(ty_ny-ky-1-ky-2+2)
    ty(ky+2)=yb+ai
    DO i=ky+3,ty_ny-ky-1
      ty(i)=ty(i-1)+ai
    END DO
  ELSE
    tx_nx=nmax
    ty_ny=nmax
  END IF

! spline approximations of degree k
  WRITE(*,*) 'vor Surfit'
  CALL surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
              nmax,eps,tx_nx,tx,ty_ny,ty,c_spl,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)

  WRITE(*,*) 'nach Surfit'
  DO i=1,m
    WRITE(*,*) f_bispev(x(i),y(i),0.0d0),z(i)
  END DO  
! evaluation of the spline approximation
!      call bispev(tx,tx_nx,ty,ty_ny,c_spl,kx,ky,xx,mx,yy,my,zz,   &
!         wrk2,lwrk2,iwrk,kwrk,ier)

  CLOSE(UNIT=InputUnitData)
END SUBROUTINE ReadUTM1

END MODULE Function_Mod
