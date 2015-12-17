
! File utilities.f90 from Library lm_f90: Nonhydrostatic local model - preoperational testversion
! Version 1.35 from 1/4/0 extracted: 1/4/0
!+ Source module for utility routines
!-------------------------------------------------------------------------------

! MODULE  utilities
!-------------------------------------------------------------------------------
!
! Description:
!   This module provides service utilities for the model. All routines are 
!   written in a manner that also other models can use it. That means:
!     - no routine uses other modules, except the declarations for the 
!       KIND-type parameter; the data access is by parameter list only
!     - no routine allocates dynamic memory; work space needed is
!       provided via the parameter list
!     - no derived data types are used
!
!   Routines (module procedures) currently contained:
!
!     - dfilt4:
!       digital filter of length 4
!
!     - dfilt8:
!       digital filter of length 8
!
!     - dolph:
!       calculates the Dolph-Chebyshev window for the initialization
!
!     - elapsed_time:
!       Returns the elapsed wall-clock time in seconds since the last call.
!       On the first call the variables are only initialized. If no system
!       clock is present, an error-value will be returned
!
!     - get_utc_date:
!       Calculates the actual date using the date of the forecast-start and 
!       the number of timesteps performed.
!
!     - istringlen:
!       Determines the length of a textstring without tailing blanks.
!
!     - phirot2phi:
!       Converts phi from the rotated system to phi in the real
!       geographical system.
!
!     - phi2phirot:
!       Converts phi from the real geographical system to phi
!       in the rotated system.
!
!     - rlarot2rla:
!       Converts lambda from the rotated system to lambda in the real
!       geographical system.
!
!     - rla2rlarot:
!       Converts lambda from the real geographical system to lambda 
!       in the rotated system.
!
!     - smoother:
!       smoothes a 2-D field by applying digital filters
!
!     - tautsp:
!       computes tension splines
!
!     - uvrot2uv:
!       Converts the wind components u and v from the rotated system
!       to the real geographical system.
!
!     - uv2uvrot:
!       Converts the wind components u and v from the real geographical
!       system to the rotated system.
!
!     - uv2df:
!       Converts the wind components u and v to wind direction and speed.
!
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8236 1493
!  email:  uschaettler@dwd.d400.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.2        1998/03/30 Ulrich Schaettler
!  Introduction of subroutine dolph used during the initialization
! 1.9        1998/09/16 Guenther Doms
!  Introduction of a smoothing routine 'smoother' which uses digital
!  filters 'dfilt4' and 'dfilt8'.
! 1.10       1998/09/29 Ulrich Schaettler
!  Routine remark eliminated and put to parallel_utilities.
!  Routines uv2uvrot and uv2df introduced
! 1.16       1998/11/02 Guenther Doms
!  Correction of filter processing in routine 'smoother'.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adaptations to use this module also in GME2LM
! 1.32       1999/08/24 Guenther Doms
!  some _ireals declarations added.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!
! Declarations:
!

!=======================================================================

! CONTAINS

!*******************************************************************************

!-------------------------------------------------------------------------------

SUBROUTINE dfilt4 (fin, idim, fhelp, fout, nfilt)

!-------------------------------------------------------------------------------
!
! Description:
!   This routine smoothes an arbitrary field (fin) of length idim by applying
!   a digital filters of length nlength 4 nfilt times. The filterd field
!   is written on fout.
!
! Method:
!   Digital filter according to Shapiro
!
!-------------------------------------------------------------------------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)          ::    &
  idim,           & ! Dimension of the field
  nfilt             ! Number of iterative filerings
REAL (KIND=ireals), INTENT (IN)          ::    &
  fin (idim)        ! input field (unfilterd)
REAL (KIND=ireals), INTENT (OUT)         ::    &
  fout (idim)       ! smoothed output field (filtered)
REAL (KIND=ireals), INTENT (INOUT)       ::    &
  fhelp(idim)       ! additional storage supplied by the calling routine

! Local variables
INTEGER (KIND=iintegers) ::    &
  i,m,            & ! loop indicees
  nf_o2             ! nfilt/2

REAL (KIND=ireals)  ::    & 
  fw(5)             ! filter weights

!-------------------------------------------------------------------------------
  DATA fw / -0.00390625, 0.03125, -0.109375, 0.21875, 0.7265625 /


! begin subroutine dfilt4

  nf_o2 = (nfilt+1)/2

  fout (:) = fin(:)
  fhelp(:) = fin(:)

  DO i = 2, idim-1
    fhelp(i) = 0.15*fout (i-1) + 0.7*fout (i) + 0.15*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) = 0.15*fhelp(i-1) + 0.7*fhelp(i) + 0.15*fhelp(i+1)
  ENDDO

  DO m = 1, nf_o2
    DO i = 5, idim-4
      fhelp(i) =  fw(5)*fout(i) &
                + fw(4)*(fout(i-1)+fout(i+1)) + fw(3)*(fout(i-2)+fout(i+2)) &
                + fw(2)*(fout(i-3)+fout(i+3)) + fw(1)*(fout(i-4)+fout(i+4))  
    ENDDO
    DO i = 5, idim-4
      fout(i) = fw(5)*fhelp(i) &
              + fw(4)*(fhelp(i-1)+fhelp(i+1)) + fw(3)*(fhelp(i-2)+fhelp(i+2)) &
              + fw(2)*(fhelp(i-3)+fhelp(i+3)) + fw(1)*(fhelp(i-4)+fhelp(i+4))  
    ENDDO
  ENDDO

  DO i = 2, idim-1
    fhelp(i) = 0.15*fout (i-1) + 0.7*fout (i) + 0.15*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) = 0.15*fhelp(i-1) + 0.7*fhelp(i) + 0.15*fhelp(i+1)
  ENDDO

END SUBROUTINE dfilt4

!-------------------------------------------------------------------------------

!*******************************************************************************
!*******************************************************************************
!-------------------------------------------------------------------------------

SUBROUTINE dfilt8 (fin, idim, fhelp, fout, nfilt)

!-------------------------------------------------------------------------------
!
! Description:
!   This routine smoothes an arbitrary field (fin) of length idim by applying
!   a digital filters of length nlength 8 nfilt times. The filterd field
!   is written on fout.
!
! Method:
!   Digital filter according to Shapiro
!
!-------------------------------------------------------------------------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)          ::    &
  idim,           & ! Dimension of the field
  nfilt             ! Number of iterative filerings
REAL (KIND=ireals), INTENT (IN)          ::    &
  fin (idim)        ! input field (unfilterd)
REAL (KIND=ireals), INTENT (OUT)         ::    &
  fout (idim)       ! smoothed output field (filtered)
REAL (KIND=ireals), INTENT (INOUT)       ::    &
  fhelp(idim)       ! additional storage supplied by the calling routine

! Local variables
INTEGER (KIND=iintegers) ::    &
  i,m,            & ! loop indicees
  nf_o2             ! nfilt/2

REAL (KIND=ireals)       ::  & 
  fw(9)  ! filter weights

!-------------------------------------------------------------------------------
  DATA fw / -0.000015259, 0.0002441406, -0.0018310546, 0.0085449218, &
            -0.027770996, 0.0666503906, -0.1221923828, 0.1745605469, &
             0.8036193848 /

! begin subroutine dfilt8

  nf_o2 = (nfilt+1)/2

  fout (:) = fin(:)
  fhelp(:) = fin(:)

  DO i = 2, idim-1
    fhelp(i) = 0.25*fout (i-1) + 0.5*fout (i) + 0.25*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) = 0.25*fhelp(i-1) + 0.5*fhelp(i) + 0.25*fhelp(i+1)
  ENDDO

  DO m = 1, nf_o2
    DO i = 9, idim-8
      fhelp(i) = fw(9)*fout(i) &
               + fw(8)*(fout(i-1)+fout(i+1)) + fw(7)*(fout(i-2)+fout(i+2)) &
               + fw(6)*(fout(i-3)+fout(i+3)) + fw(5)*(fout(i-4)+fout(i+4)) &
               + fw(4)*(fout(i-5)+fout(i+5)) + fw(3)*(fout(i-6)+fout(i+6)) &
               + fw(2)*(fout(i-7)+fout(i+7)) + fw(1)*(fout(i-8)+fout(i+8))  
    ENDDO
    DO i = 9, idim-8
      fout(i) = fw(9)*fhelp(i) &
              + fw(8)*(fhelp(i-1)+fhelp(i+1)) + fw(7)*(fhelp(i-2)+fhelp(i+2)) &
              + fw(6)*(fhelp(i-3)+fhelp(i+3)) + fw(5)*(fhelp(i-4)+fhelp(i+4)) &
              + fw(4)*(fhelp(i-5)+fhelp(i+5)) + fw(3)*(fhelp(i-6)+fhelp(i+6)) &
              + fw(2)*(fhelp(i-7)+fhelp(i+7)) + fw(1)*(fhelp(i-8)+fhelp(i+8))
    ENDDO
  ENDDO

  DO i = 2, idim-1
    fhelp(i) = 0.25*fout (i-1) + 0.5*fout (i) + 0.25*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) = 0.25*fhelp(i-1) + 0.5*fhelp(i) + 0.25*fhelp(i+1)
  ENDDO

END SUBROUTINE dfilt8
!-------------------------------------------------------------------------------

!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

SUBROUTINE dolph (deltat, taus, m, window, t, time, time2, w, w2)

!-------------------------------------------------------------------------------
!
! Description:
!  Calculation of Dolph-Chebyshev window or, for short, Dolph Window, using
!  the expression in the reference:
!    Antoniou, Andreas, 1993: Digital Filters: Analysis,
!    Design and Applications. McGraw-Hill, Inc., 689pp.
!
!  The Dolph window is optimal in the following sense:
!  For a given main-lobe width, the stop-band attenuation is minimal;
!  for a given stop-band level, the main-lobe width is minimal.
!
! Method:
!
! Modules used:    NONE
!
!-------------------------------------------------------------------------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

! Parameter List:
! ---------------

INTEGER (KIND=iintegers), INTENT (IN)             ::  &
  m                   ! for dimensioning the work arrays

REAL  (KIND=ireals), INTENT (IN)                  ::  &
  deltat, taus        ! time step and cutoff period for filtering

REAL  (KIND=ireals), INTENT (OUT)                 ::  &
  window(0:2*m)       ! result

! The following variables are only used for work space
REAL  (KIND=ireals), INTENT (OUT)                 ::  &
  t(0:2*m), time(0:2*m), time2(0:2*m), w(0:2*m), w2(0:2*m)

! Local Variables:
! ----------------

INTEGER (KIND=iintegers)  :: nt, i, n, nm1, nn
REAL    (KIND=ireals)     :: zpi, zthetas, zx0, zarg, zterm1, zterm2, zrr,   &
                             zr, zdb, zsum, zsumw

!------------ End of header ----------------------------------------------------

! Begin subroutine dolph

  zpi = 4.0_ireals * ATAN(1.0_ireals)

  n = 2*m+1
  nm1 = n-1
  zthetas = 2.0_ireals*zpi*deltat/taus
  zx0 = 1.0_ireals / COS(zthetas/2.0_ireals)
  zterm1 = (zx0 + SQRT(zx0**2-1))**(FLOAT(N-1))
  zterm2 = (zx0 - SQRT(zx0**2-1))**(FLOAT(N-1))
  zrr = 0.5*(zterm1 + zterm2)
  zr = 1/zrr
  zdb = 20.0_ireals * LOG10(zr)

!------------------------------------------------------------

  DO nt = 0, M
    zsum = 1
    DO i = 1, M
      zarg = zx0 * cos(i*zpi/N)
      ! Calculate the Chebyshev polynomials
      ! Reference: Numerical Recipes, Page 184, recurrence
      !   T_n(x) = 2xT_{n-1}(x) - T_{n-2}(x) ,  n>=2.
      T(0) = 1
      T(1) = zarg
      DO nn=2,nm1
        T(nn) = 2*zarg*T(nn-1) - T(nn-2)
      ENDDO
      zterm1 = T(nm1)
      zterm2 = cos(2*nt*zpi*i/n)
      zsum   = zsum + zr*2 * zterm1 * zterm2
    ENDDO
    w(nt) = zsum / n
    TIME(nt) = nt
  ENDDO

  ! Fill in the negative-time values by symmetry.
  DO nt = 0, m
    w2(m+nt) = w(nt)
    w2(m-nt) = w(nt)
    time2(m+nt) =  time(nT)
    time2(m-nt) = -time(nT)
  ENDDO

  ! Fill up the array for return
  zsumw = 0.0_ireals
  DO nt = 0, 2*m
    zsumw = zsumw + w2(nt)
  ENDDO

  DO nt=0,2*m
    WINDOW(nt) = w2(nt)
  ENDDO
!
!
!----------------------------------------------------------
!       PRINT *, (w2(nT),    nT=0,2*M)
!----------------------------------------------------------
!

END SUBROUTINE dolph

!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

SUBROUTINE elapsed_time    (realtimedif, istat)

!-------------------------------------------------------------------------------
!
! Description:
!   Returns the elapsed wall-clock time in seconds since the last call. On
!   the first call the variables are only initialized. If no system clock is
!   present, an error value of istat=1 will be returned, if the optional
!   argument istat was passed from the calling routine. 
!   realtimedif is set to 0 then.
!
! Method:
!   The intrinsic function SYSTEM_CLOCK is used, that returns the number of
!   clock counts since some system dependent event in the past (e.g. midnight
!   for a 24-hour system clock). The difference of clock counts since the last
!   call is determined and converted into seconds. The variables "lfirst"
!   and "icountsold" (see below) have to be SAVEd for the next call.
!
! Modules used:    NONE
!
!-------------------------------------------------------------------------------
!
! Parameter List:
! ---------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

REAL  (KIND=ireals), INTENT (OUT)                 ::  &
      realtimedif     ! wall-clock time since the last call in seconds
                      ! (0 if no system-clock is available)

INTEGER (KIND=iintegers), INTENT (OUT), OPTIONAL  ::  &
      istat           ! optional argument for error value


! Local Variables:
! ----------------

LOGICAL, SAVE      :: lfirst = .TRUE.   ! determine whether first call or not

INTEGER, SAVE      :: icountsold        ! number of counts in the last call

INTEGER            :: icountsnew,     & ! number of counts in this call
                      ir, im            ! other arguments to SYSTEM_CLOCK

LOGICAL            :: lpres             ! if optional argument is present

!------------ End of header ----------------------------------------------------

! Begin subroutine elapsed_time

  lpres = PRESENT (istat)

  CALL SYSTEM_CLOCK ( COUNT=icountsnew, COUNT_RATE=ir, COUNT_MAX=im )

  IF ( ir /= 0 ) THEN
    ! system clock is present
    IF (lpres) THEN
      istat = 0
    ENDIF

    IF (lfirst) THEN
      ! first call: store value for the number of clock counts
      icountsold = icountsnew
      lfirst     = .FALSE.
    ELSE
      ! convert the clock counts to seconds
      IF ( icountsnew >= icountsold ) THEN
        realtimedif = ( REAL (icountsnew - icountsold, ireals) )      &
                      / REAL (ir,ireals)
      ELSE
        realtimedif = REAL (im- (icountsold-icountsnew ), ireals)     &
                      / REAL (ir, ireals)
      ENDIF
      icountsold = icountsnew
    ENDIF
  ELSE
    ! no system clock present: set error value
    realtimedif = 0.0
    IF ( lpres ) THEN
      istat = 1
    ENDIF
  ENDIF

END SUBROUTINE elapsed_time

!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

SUBROUTINE get_utc_date (ntsteps, ystartdate, dt,                          &
                         yactdate1, yactdate2, nactday, acthour)

!-------------------------------------------------------------------------------
!
! Description:
!   This routine determines the actual date of this forecast step.
!
! Method:
!   Using the date of the forecast-start, the number of time steps 
!   already performed and the length of the time steps, the actual
!   date is calculated taking leap-years into consideration.
!   The date is given in three different formats.
!
! Modules used:    NONE
!
!-------------------------------------------------------------------------------
!
! Input Parameter list:
! ---------------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

INTEGER   (KIND=iintegers), INTENT(IN)   ::                           &
  ntsteps    ! number of actual performed time-steps

REAL      (KIND=ireals), INTENT(IN)      ::                           &
  dt         ! time step in seconds

CHARACTER (LEN=10), INTENT(IN)           ::                           &
  ystartdate ! start date of the forecast


! Output Parameter list:
! ----------------------

CHARACTER (LEN=10), INTENT(OUT)          ::                           &
  yactdate1  ! actual date in the form   yyyymmddhh

CHARACTER (LEN=22), INTENT(OUT)          ::                           &
  yactdate2  ! actual date in the form   wd   dd.mm.yy  hh UTC


INTEGER   (KIND=iintegers), INTENT(OUT)  ::                           &
  nactday    ! day of the year

REAL      (KIND=ireals), INTENT(OUT)     ::                           &
  acthour    ! actual hour of the day

! Local variables:
INTEGER   (KIND=iintegers)   ::                                       &
  month(12), monthsum(13), ileap, iweek, iy, m,                       &
  idd, imm, iyy, ihh, iday, imonth, iyear, ihour, immhours, iyyhours
CHARACTER (LEN=3)            :: yweek(7)

!------------ End of header ----------------------------------------------------

! Begin subroutine get_utc_date

DATA         month  / 31 ,  28 ,  31 ,  30 ,  31 ,  30 ,       &
                      31 ,  31 ,  30 ,  31 ,  30 ,  31 /
DATA         yweek  /'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /


! Statementfunction: ileap(yy) = 0:  no leap year, 
!                    ileap(yy) = 1:  leap year
  ileap (iy) = IABS( MOD(iy,4) - 4) / 4

! Divide ystartdate in day, month, year and hour
! and calculate the sums of days from the beginning of the year to the 
! end of the months
  READ ( ystartdate, '(I4,3I2)' ) iyy, imm, idd, ihh

  month (2)    = 28 + ileap (iyy)
  monthsum(1) =  0
  DO m =  2 , 13
    monthsum(m) =  monthsum(m-1) + month(m-1)
  enddo

! Determine how many hours have passed in this year
  iyyhours = (idd*24) + monthsum(imm)*24 + (ihh-24)
  iyyhours = iyyhours + NINT (ntsteps*dt)/3600

! Take turning of the year into account
  IF (iyyhours < 0) THEN
    iyear    = iyy-1
    iyyhours = 8760 + ileap(iyear)*24 + iyyhours
  ELSE IF (iyyhours >= (8760+ileap(iyy)*24)) THEN
    iyear    = iyy+1
    iyyhours = iyyhours - (8760+ileap(iyy)*24)
  ELSE
    iyear    =   iyy
  ENDIF

! Determine the actual date from iyyhours
  m        = 1
  immhours = iyyhours
  DO WHILE (immhours >= 0)
    m        = m+1
    immhours = iyyhours - monthsum(m) * 24
  ENDDO
  imonth   = m-1

  immhours = iyyhours - monthsum(imonth)*24
  iday     = immhours/24 + 1
  ihour    = MOD(immhours,24)
  acthour  = FLOAT(ihour) + dt/3600.* MOD(ntsteps,INT(3600./dt+0.01)) + 0.0001
  ihour    = INT(acthour)
  nactday  = monthsum(imonth) + iday + INT(acthour/24. + 0.0001)
  iweek    = MOD(monthsum(imonth) + iday + (iyear-1981)+(iyear-1981)/4+2 , 7) + 1

  WRITE ( yactdate1(1:4) , '(I4.4)' ) iyear
  WRITE ( yactdate1(5:6) , '(I2.2)' ) imonth
  WRITE ( yactdate1(7:8) , '(I2.2)' ) iday
  WRITE ( yactdate1(9:10), '(I2.2)' ) ihour

  yactdate2 = yweek(iweek)//' '//yactdate1(7:8)//'.'// yactdate1(5:6)//'.'// &
                         yactdate1(1:4)//'  '//yactdate1(9:10)//' UTC'

END SUBROUTINE get_utc_date

!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

FUNCTION istringlen (ytext)

!-------------------------------------------------------------------------------
!
! Description:
!   This function determines the length of a textstring without tailing 
!   blanks.
!
! Method:
!   The intrinsic function LEN is used that determines the length of
!   the textstring with tailing blanks. If blanks are present at the
!   end of the string, they are not counted here.
!
!-------------------------------------------------------------------------------
!
! Parameter list:
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

CHARACTER, INTENT(IN)      :: ytext*(*)
INTEGER (KIND=iintegers)   :: istringlen

! Local variables:
INTEGER (KIND=iintegers)   :: itest             ! 

!-------------------------------------------------------------------------------

! Begin subroutine istringlen

! Determine the length of the string with tailing blanks 
  istringlen = LEN(ytext)
  itest      = 0

! Cut tailing blanks off for the counting
  DO WHILE (itest == 0) 
    IF (ytext(istringlen:istringlen) /= ' ') THEN
      itest = 1
    ELSE
      istringlen = istringlen - 1
      IF (istringlen == 0) THEN
        itest = 1
      ENDIF
    ENDIF
  ENDDO
 
END FUNCTION istringlen

!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

FUNCTION  phirot2phi ( phirot, rlarot, polphi, pollam, polgam )

!-------------------------------------------------------------------------------
!
! Description:
!   This function converts phi from one rotated system to phi in another
!   system. If the optional argument polgam is present, the other system
!   can also be a rotated one, where polgam is the angle between the two
!   north poles.
!   If polgam is not present, the other system is the real geographical
!   system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

! Parameter list:
REAL (KIND=ireals), INTENT (IN)      ::        &
  polphi,   & ! latitude of the rotated north pole
  pollam,   & ! longitude of the rotated north pole
  phirot,   & ! latitude in the rotated system
  rlarot      ! longitude in the rotated system

REAL (KIND=ireals), INTENT (IN), OPTIONAL      ::        &
  polgam      ! angle between the north poles of the systems

REAL (KIND=ireals)                   ::        &
  phirot2phi  ! latitude in the geographical system

! Local variables
REAL (KIND=ireals)                   ::        &
  zsinpol, zcospol, zphis, zrlas, zarg, zgam

REAL (KIND=ireals)                   ::        &
  zrpi18 = 57.2957795_ireals,                  &
  zpir18 = 0.0174532925_ireals

!-------------------------------------------------------------------------------

! Begin function phirot2phi

  zsinpol     = SIN (zpir18 * polphi)
  zcospol     = COS (zpir18 * polphi)
 
  zphis       = zpir18 * phirot
  IF (rlarot > 180.0_ireals) THEN
    zrlas = rlarot - 360.0_ireals
  ELSE
    zrlas = rlarot
  ENDIF
  zrlas       = zpir18 * zrlas

  IF ( PRESENT (polgam) ) THEN
    zgam  = zpir18 * polgam
    zarg  = zsinpol*SIN (zphis) +                                           &
        zcospol*COS(zphis) * ( COS(zrlas)*COS(zgam) - SIN(zgam)*SIN(zrlas) )
  ELSE
    zarg  = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
  ENDIF
 
  phirot2phi  = zrpi18 * ASIN (zarg)

END FUNCTION phirot2phi

!******************************************************************************
!******************************************************************************

!------------------------------------------------------------------------------

FUNCTION  phi2phirot ( phi, rla, polphi, pollam )

!------------------------------------------------------------------------------
! Description:
!   This routine converts phi from the real geographical system to phi
!   in the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!-------------------------------------------------------------------------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

! Parameter list:
REAL (KIND=ireals), INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in the rotated system
  rla        ! longitude in the rotated system

REAL (KIND=ireals)                   ::        &
  phi2phirot ! longitude in the rotated system

! Local variables
REAL (KIND=ireals)                       ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & !
  zpir18 = 0.0174532925_ireals

!-------------------------------------------------------------------------------

! Begin function phi2phirot

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0_ireals) THEN
    zrla1  = rla - 360.0_ireals
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = SIN (zphi) * zsinpol
  zarg2    = COS (zphi) * zcospol * COS (zrla - zlampol)

  phi2phirot = zrpi18 * ASIN (zarg1 + zarg2)

END FUNCTION phi2phirot

!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

FUNCTION  rlarot2rla (phirot, rlarot, polphi, pollam, polgam)

!-------------------------------------------------------------------------------
!
! Description:
!   This function converts lambda from one rotated system to lambda in another
!   system. If the optional argument polgam is present, the other system
!   can also be a rotated one, where polgam is the angle between the two
!   north poles.
!   If polgam is not present, the other system is the real geographical
!   system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
! Modules used:    NONE
!
!-------------------------------------------------------------------------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------

! Parameter list:
REAL (KIND=ireals), INTENT (IN)      ::        &
  polphi,   & ! latitude of the rotated north pole
  pollam,   & ! longitude of the rotated north pole
  phirot,   & ! latitude in the rotated system
  rlarot      ! longitude in the rotated system

REAL (KIND=ireals), INTENT (IN), OPTIONAL      ::        &
  polgam      ! latitude of the rotated north pole

REAL (KIND=ireals)                   ::        &
  rlarot2rla  ! latitude in the geographical system

! Local variables
REAL (KIND=ireals)                   ::        &
  zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam

REAL (KIND=ireals)                   ::        &
  zrpi18 = 57.2957795_ireals,   & !
  zpir18 = 0.0174532925_ireals

!-------------------------------------------------------------------------------

! Begin function rlarot2rla

  zsinpol = SIN (zpir18 * polphi)
  zcospol = COS (zpir18 * polphi)

  zlampol = zpir18 * pollam
  zphis   = zpir18 * phirot
  IF (rlarot > 180.0_ireals) THEN
    zrlas = rlarot - 360.0_ireals
  ELSE
    zrlas = rlarot
  ENDIF
  zrlas   = zpir18 * zrlas

  IF ( PRESENT(polgam) ) THEN
    zgam    = zpir18 * polgam
    zarg1   = SIN (zlampol) *                                                &
      (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
       + zcospol * SIN(zphis))                                               &
    - COS (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))

    zarg2   = COS (zlampol) *                                                &
      (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
       + zcospol * SIN(zphis))                                               &
    + SIN (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))
  ELSE
    zarg1   = SIN (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
                                zcospol *              SIN(zphis)) -    &
              COS (zlampol) *             SIN(zrlas) * COS(zphis)
    zarg2   = COS (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
                                zcospol *              SIN(zphis)) +   &
              SIN (zlampol) *             SIN(zrlas) * COS(zphis)
  ENDIF
 
  IF (zarg2 == 0.0) zarg2 = 1.0E-20_ireals
 
  rlarot2rla = zrpi18 * ATAN2(zarg1,zarg2)
 
END FUNCTION rlarot2rla

!**********************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

FUNCTION  rla2rlarot ( phi, rla, polphi, pollam )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine converts lambda from the real geographical system to lambda 
!   in the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!-------------------------------------------------------------------------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

!
! Parameter list:
REAL (KIND=ireals), INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in the rotated system
  rla        ! longitude in the rotated system

REAL (KIND=ireals)                   ::        &
  rla2rlarot ! latitude in the geographical system

! Local variables
REAL (KIND=ireals)                       ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & !
  zpir18 = 0.0174532925_ireals

!-------------------------------------------------------------------------------

! Begin function rla2rlarot

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0_ireals) THEN
    zrla1  = rla - 360.0_ireals
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = - SIN (zrla-zlampol) * COS(zphi)
  zarg2    = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)

  IF (zarg2 == 0.0) zarg2 = 1.0E-20_ireals

  rla2rlarot = zrpi18 * ATAN2 (zarg1,zarg2)

END FUNCTION rla2rlarot

!*******************************************************************************
!******************************************************************************

!-------------------------------------------------------------------------------

SUBROUTINE TAUTSP(TAU,GTAU,NTAU,GAMMA,S,BREAK,COEF,L,IFLAG)

!-------------------------------------------------------------------------------
!
! Description:
!   BERECHNET DEN TAUT-SPLINE FUER DIE DATEN: TAU(I),GTAU(I),I=1,.,NTAU
!
! Method:
!   WENN GAMMA GT.0  WERDEN ZUSAETZLICHE KNOTEN BERECHNET
!   GAMMA I.A.  2.5   BZW.  5.5
!
!   BREAK,COEF,L,K GEBEN DIE PP-DARSTELLUNG DER INTERPOLATIONSFUNKTION
!
!   FUER BREAK(I).LE.X.LE.BREAK(I+1) HAT DIE INTERPOLATIONSFUNKTION
!   FOLGENDE FORM
!
!   F(X)=COEF(1,I)+DX(COEF(2,I)+DX/2(COEF(3,I)+DX/3(COEF(4,I)))
!   MIT   DX=X-BREAK(I) UND I=1,...,L
!
!   IFLAG=0  KEIN FEHLER IN TAUTSP
!   IFLAG=2  FALSCHER INPUT
!
!   S(NTAU,6)  WORK-ARRAY
!
!=======================================================================
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

INTEGER (KIND=iintegers) IFLAG,L,NTAU,I,METHOD,NTAUM1
REAL (KIND=ireals)                                                 &
   BREAK(L),COEF(4,L),GAMMA,GTAU(NTAU),S(NTAU,6),TAU(NTAU),        &
   ALPHA,C,D,DEL,DENOM,DIVDIF,ENTRY,ENTRY3,FACTOR,FACTR2,GAM,      &
   ONEMG3,ONEMZT,RATIO,SIXTH,TEMP,X,Z,ZETA,ZT2,ALPH

!=======================================================================
!
      ALPH(X)= MIN (1.0_ireals,ONEMG3/X)
!
!   TEST DER INPUTPARAMETER
!
      IF (NTAU .LT. 4) THEN
        PRINT 600,NTAU
 600    FORMAT('  NTAU =',I4,'  NTAU SOLL GROESSER ALS 4 SEIN')
        GO TO 999
      ENDIF
!
!   BERECHNUNG VON DELTA TAU UND DER 1. UND 2.ABLEITUNG DER DATEN
!
      NTAUM1=NTAU-1
      DO I=1,NTAUM1
      S(I,1)=TAU(I+1)-TAU(I)
      IF (S(I,1) .LE. 0.0) THEN
        PRINT 610,I,TAU(I),TAU(I+1)
 610    FORMAT(' PUNKT ',I3,'  UND DIE NAECHSTEN',2E15.6,'SIND IN&
     &  FALSCHER REIHENFOLGE')
        GO TO 999
      ENDIF
      S(I+1,4) = (GTAU(I+1) - GTAU(I))/S(I,1)
      ENDDO
!
      DO I=2,NTAUM1
      S(I,4) = S(I+1,4) - S(I,4)
      ENDDO
!
!   2.ABLEITUNG VON GTAU AN DEN PUNKTEN TAU
!
      I=2
      S(2,2) = S(1,1)/3.0
      SIXTH = 1.0/6.0
      METHOD = 2
      GAM = GAMMA
      IF(GAM .LE. 0.0) METHOD = 1
      IF(GAM .GT. 3.0) THEN
        METHOD = 3
        GAM = GAM - 3.0
      ENDIF
      ONEMG3=1.0 - GAM/3.0
!
!   SCHLEIFE UEBER I
!
 10   CONTINUE
      Z=0.5
      IF (METHOD .EQ. 1) THEN
        GO TO 19
      ELSE IF (METHOD .EQ. 2) THEN
        GO TO 11
      ELSE IF (METHOD .EQ. 3) THEN
        GO TO 12
      ENDIF
 11   CONTINUE
      IF(S(I,4)*S(I+1,4).LT.0.0) GO TO 19
 12   CONTINUE
      TEMP = ABS(S(I+1,4))
      DENOM = ABS(S(I,4)) + TEMP
      IF(DENOM.EQ.0.0) GO TO 19
      Z = TEMP/DENOM
      IF(ABS(Z-0.5).LE.SIXTH) Z=0.5
 19   CONTINUE
      S(I,5) = Z
!
!   ERSTELLEN EINES TEILES DER I-TEN GLEICHUNG
!
      IF (Z-0.5 .LT. 0.) THEN
        ZETA = GAM*Z
        ONEMZT = 1.0 - ZETA
        ZT2 = ZETA**2
        ALPHA = ALPH(ONEMZT)
        FACTOR = ZETA/(ALPHA*(ZT2 - 1.0) + 1.0)
        S(I,6) = ZETA*FACTOR/6.0
        S(I,2) = S(I,2) + S(I,1)*((1.0-ALPHA*ONEMZT)*FACTOR*0.5-S(I,6))
        IF(S(I,2).LE.0.0) S(I,2) = 1.0
        S(I,3) = S(I,1)/6.0
!
      ELSE IF (Z-0.5 .EQ. 0.) THEN
!
        S(I,2) = S(I,2) + S(I,1)/3.0
        S(I,3) = S(I,1)/6.0
!
      ELSE
!
        ONEMZT = GAM*(1.0 - Z)
        ZETA = 1.0 - ONEMZT
        ALPHA = ALPH(ZETA)
        FACTOR = ONEMZT/(1.0 - ALPHA*ZETA*(1.0 + ONEMZT))
        S(I,6) = ONEMZT*FACTOR/6.0
        S(I,2) = S(I,2) + S(I,1)/3.0
        S(I,3) = S(I,6) * S(I,1)
      ENDIF
!
      IF (I .GT. 2) GO TO 30
      S(1,5) = 0.5
!
!   DIE ERSTEN BEIDEN GLEICHUNGEN ERZWINGEN STETIGKEIT DER 1. UND 3.AB-
!   LEITUNG IN TAU(I)
!
      S(1,2) = S(1,1)/6.0
      S(1,3) = S(2,2)
      ENTRY3 = S(2,3)

      IF (Z-0.5 .LT. 0.) THEN
        FACTR2 = ZETA*(ALPHA*(ZT2-1.0)+1.0)/(ALPHA*(ZETA*ZT2-1.0) + 1.0)
        RATIO = FACTR2*S(2,1)/S(1,2)
        S(2,2) = FACTR2*S(2,1) + S(1,1)
        S(2,3) = -FACTR2 * S(1,1)
!
      ELSE IF (Z-0.5 .EQ. 0.) THEN
!
        RATIO = S(2,1)/S(1,2)
        S(2,2) = S(2,1) + S(1,1)
        S(2,3) = -S(1,1)
!
      ELSE
!
        RATIO = S(2,1)/S(1,2)
        S(2,2) = S(2,1) + S(1,1)
        S(2,3) = -S(1,1)*6.0*ALPHA*S(2,6)
      ENDIF
!
!   ELIMINATION DER 1.UNBEKANNTEN AUS DER 2.GLEICHUNG
!
      S(2,2) = RATIO*S(1,3) + S(2,2)
      S(2,3) = RATIO*ENTRY3 + S(2,3)
      S(1,4) = S(2,4)
      S(2,4) = RATIO*S(1,4)
      GO TO 35
!
!
 30   CONTINUE
      S(I,2) = RATIO*S(I-1,3) + S(I,2)
      S(I,4) = RATIO*S(I-1,4) + S(I,4)
!
!   AUFSTELLEN DES TEILES DER NAECHSTEN GLEICHUNG,DER VOM I-TEN INTERVAL
!   ABHAENGT
!
 35   CONTINUE
      IF (Z-0.5 .LT. 0.) THEN
        RATIO = -S(I,6)*S(I,1)/S(I,2)
        S(I+1,2) = S(I,1)/3.0
!
      ELSE IF (Z-0.5 .EQ. 0.) THEN
!
        RATIO = -S(I,1)/(6.0*S(I,2))
        S(I+1,2) = S(I,1)/3.0
!
      ELSE
!
        RATIO = -(S(I,1)/6.0)/S(I,2)
        S(I+1,2) = S(I,1)*((1.0 - ZETA*ALPHA)*0.5*FACTOR-S(I,6))
      ENDIF
!
!   ENDE DER SCHLEIFE UEBER I
!
      I = I + 1
      IF(I.LT.NTAUM1) GO TO 10
      S(I,5) = 0.5
!
!   DIE BEIDEN LETZTEN GLEICHUNGEN ERZWINGEN STETIGKEIT DER
!   1. UND 3. ABLEITUNG IN TAU(NTAU-1)
!
      ENTRY = RATIO*S(I-1,3) + S(I,2) + S(I,1)/3.0
      S(I+1,2) = S(I,1)/6.0
      S(I+1,4) = RATIO*S(I-1,4) + S(I,4)
      IF (Z-0.5 .LT. 0.) THEN
        RATIO = S(I,1)*6.0*S(I-1,6)*ALPHA/S(I-1,2)
        S(I,2) = RATIO*S(I-1,3) + S(I,1) + S(I-1,1)
        S(I,3) = -S(I-1,1)
!
      ELSE IF (Z-0.5 .EQ. 0.) THEN
!
        RATIO = S(I,1)/S(I-1,2)
        S(I,2) = RATIO*S(I-1,3) + S(I,1) + S(I-1,1)
        S(I,3) = -S(I-1,1)
!
      ELSE
!
        FACTR2 = ONEMZT*(ALPHA*(ONEMZT**2-1.0)+1.0)/     &
                        (ALPHA*(ONEMZT**3-1.0)+1.0)
        RATIO = FACTR2*S(I,1)/S(I-1,2)
        S(I,2) = RATIO*S(I-1,3) + FACTR2*S(I-1,1) + S(I,1)
        S(I,3) = -FACTR2*S(I-1,1)
      ENDIF
!
!   ELIMINATION VON XI
!
      S(I,4) = RATIO*S(I-1,4)
      RATIO = -ENTRY/S(I,2)
      S(I+1,2) = RATIO*S(I,3) + S(I+1,2)
      S(I+1,4) = RATIO*S(I,4) + S(I+1,4)
!
!   RUECKSUBSTITUTION
!
      S(NTAU,4) = S(NTAU,4)/S(NTAU,2)
 50   CONTINUE
      S(I,4) = (S(I,4) - S(I,3)*S(I+1,4))/S(I,2)
      I = I - 1
      IF(I.GT.1) GO TO 50
      S(1,4) = (S(1,4) -S(1,3)*S(2,4)-ENTRY3*S(3,4))/S(1,2)
!
!   ERZEUGEN DER POLYNOM-TEILE
!
      BREAK(1) = TAU(1)
      L = 1
      DO 70 I = 1,NTAUM1
      COEF(1,L) = GTAU(I)
      COEF(3,L) = S(I,4)
      DIVDIF = (GTAU(I+1) - GTAU(I))/S(I,1)
      Z = S(I,5)
!
      IF (Z-0.5 .LT. 0.) THEN
        IF(Z.EQ.0.0) GO TO 65
        ZETA = GAM*Z
        ONEMZT = 1.0 -ZETA
        C = S(I+1,4)/6.0
        D = S(I,4)*S(I,6)
        L = L + 1
        DEL = ZETA*S(I,1)
        BREAK(L) = TAU(I) + DEL
        ZT2 = ZETA**2
        ALPHA = ALPH(ONEMZT)
        FACTOR = ONEMZT**2*ALPHA
        COEF(1,L) = GTAU(I) + DIVDIF*DEL+S(I,1)**2*(D*ONEMZT*(FACTOR- &
           1.0) + C*ZETA*(ZT2 - 1.0))
        COEF(2,L) = DIVDIF + S(I,1)*(D*(1.0-3.0*FACTOR) + C*(3.0*ZT2- &
           1.0))
        COEF(3,L) = 6.0*(D*ALPHA*ONEMZT + C*ZETA)
        COEF(4,L) = 6.0*(C - D*ALPHA)/S(I,1)
        COEF(4,L-1) = COEF(4,L) -6.0*D*(1.0-ALPHA)/(DEL*ZT2)
        COEF(2,L-1) = COEF(2,L) - DEL*(COEF(3,L)-DEL*0.5*COEF(4,L-1))
        GO TO 68
!
      ELSE IF (Z-0.5 .EQ. 0.) THEN
!
        COEF(2,L) = DIVDIF - S(I,1)*(2.0*S(I,4) + S(I+1,4))/6.0
        COEF(4,L) = (S(I+1,4) - S(I,4))/S(I,1)
        GO TO 68
!
      ELSE
!
        ONEMZT = GAM*(1.0 - Z)
        IF(ONEMZT.EQ.0.0) GO TO 65
        ZETA = 1.0 - ONEMZT
        ALPHA = ALPH(ZETA)
        C = S(I+1,4)*S(I,6)
        D = S(I,4)/6.0
        DEL = ZETA*S(I,1)
        BREAK(L+1) = TAU(I) + DEL
        COEF(2,L) = DIVDIF -S(I,1)*(2.0*D + C)
        COEF(4,L) = 6.0*(C*ALPHA - D)/S(I,1)
        L = L + 1
        COEF(4,L) = COEF(4,L-1) + 6.0*(1.0-ALPHA)*C/(S(I,1)*ONEMZT**3)
        COEF(3,L) = COEF(3,L-1) + DEL* COEF(4,L-1)
        COEF(2,L) = COEF(2,L-1) + DEL*(COEF(3,L-1)+DEL*0.5*COEF(4,L-1))
        COEF(1,L) = COEF(1,L-1) + DEL*(COEF(2,L-1)+DEL*0.5*   &
              (COEF(3,L-1) + (DEL/3.0)*COEF(4,L-1)))
        GO TO 68
      ENDIF
!
!
   65 CONTINUE
      COEF(2,L) = DIVDIF
      COEF(3,L) = 0.0
      COEF(4,L) = 0.0
   68 CONTINUE
!
      L = L + 1
      BREAK(L) = TAU(I+1)
   70 CONTINUE

      IFLAG = 0
      RETURN
!
 999  IFLAG = 2

END SUBROUTINE tautsp

!******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

SUBROUTINE smoother (fin, fout, ie, je, nlength, nfilt)

!-------------------------------------------------------------------------------
!
! Description:
!   This routine smoothes an arbitrary two-dimensional field (fin) by applying
!   digital filters of length nlength (4 or 8) nfilt times. The filterd field
!   is written on fout.
!
! Method:
!   Call of digital filters (dfilt4 or dfilt8) in each direction.    
!
!-------------------------------------------------------------------------------
!
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)      ::    &
  ie, je,         & ! Dimension of the field
  nlength,        & ! Filter lenght
  nfilt             ! Number of iterative filerings

REAL (KIND=ireals), INTENT (IN)          ::    &
  fin (ie,je)       ! 2-d input field (unfilterd)

REAL (KIND=ireals), INTENT (OUT)         ::    &
  fout (ie,je)      ! 2-d smoothed output field (filtered)

! Local variables
INTEGER (KIND=iintegers) ::    &
  i,j,            & ! loop indicees
  ix,jx,          & ! dimension of the field
  nf                ! number of filter iterations

REAL (KIND=ireals)      ::    &
  sxin(ie), sxh(ie), sxout(ie),  & ! local storage
  syin(je), syh(je), syout(je)     ! local storage

!-------------------------------------------------------------------------------
! begin subroutine smoother

  IF ( nlength /= 4 .AND. nlength /= 8 ) THEN
    PRINT*, ' CAUTION: Filterlength =',nlength,' not implemented'
    PRINT*, ' No filtering of output field done'
    RETURN
  ENDIF

  nf = nfilt
  ix = ie
  jx = je

  DO j = 1, jx
    sxin(:) = fin(:,j)
    IF(nlength==4)  CALL dfilt4 ( sxin, ix, sxh, sxout, nf )
    IF(nlength==8)  CALL dfilt8 ( sxin, ix, sxh, sxout, nf )
    fout(:,j) = sxout(:)
  ENDDO
  DO i = 1, ix
    syin(:) = fout(i,:)
    IF(nlength==4)  CALL dfilt4 ( syin, jx, syh, syout, nf )
    IF(nlength==8)  CALL dfilt8 ( syin, jx, syh, syout, nf )
    fout(i,:) = syout(:)
  ENDDO

END SUBROUTINE smoother

!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

SUBROUTINE uvrot2uv (urot, vrot, phi, rla, polphi, pollam, u, v)


!-------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the rotated system
!   to the real geographical system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!-------------------------------------------------------------------------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

! Parameter list:
REAL (KIND=ireals), INTENT (IN)          ::    &
  urot, vrot,     & ! wind components in the rotated grid
  phi, rla,       & ! latitude in the true geographical system
  polphi, pollam    ! latitude and longitude of the north pole of the
                    ! rotated grid

REAL (KIND=ireals), INTENT (OUT)         ::    &
  u, v              ! wind components in the true geographical system

! Local variables

REAL (KIND=ireals)                       ::    &
  zpolphi, zpollam, zphi, zrla, pollamd, zrlas, zarg, zbeta

REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & !
  zpir18 = 0.0174532925_ireals

REAL rla2rlarot
!-------------------------------------------------------------------------------

! Begin subroutine uvrot2uv

! Converting from degree to  BOGENMASS
  zpolphi = polphi * zpir18
  zpollam = pollam * zpir18
  zrla    = rla    * zpir18
  zphi    = phi    * zpir18
  pollamd = pollam
  IF (pollam < 0.0) pollamd = 360.0_ireals + pollam

! Longitude in the rotated system
  zrlas   = rla2rlarot (phi, rla, polphi, pollam) * zpir18

! Calculate the angle between the latitudes
  zarg   = - SIN (zpolphi) * SIN(zrla-zpollam) * SIN(zrlas)              &
                           - COS(zrla-zpollam) * COS(zrlas)
  zarg   = MIN ( 1.0_ireals, zarg)
  zarg   = MAX (-1.0_ireals, zarg)
  zbeta  = ACOS(zarg)
  zbeta  = SIGN (zbeta, - (rla - (pollamd-180.0_ireals)))

! Convert the u-component
  u      = urot * COS(zbeta) - vrot * SIN(zbeta)

!  Convert the v-component
  v      = urot * SIN(zbeta) + vrot * COS(zbeta)

END SUBROUTINE uvrot2uv

!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

SUBROUTINE uv2uvrot(u, v, phi, rla, polphi, pollam, urot, vrot)

!-------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the real
!   geographical system to the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!-------------------------------------------------------------------------------
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Parameter list:
REAL (KIND=ireals), INTENT (IN)          ::    &
  u   , v   ,     & ! wind components in the true geographical system
  phi, rla,       & ! coordinates in the true geographical system
  polphi, pollam    ! latitude and longitude of the north pole of the
                    ! rotated grid

REAL (KIND=ireals), INTENT (OUT)         ::    &
  urot, vrot        ! wind components in the rotated grid             

! Local variables

REAL (KIND=ireals)                       ::    &
  zpolphi, zpollam, zphi, zrla, pollamd, zrlas, zarg, zbeta

REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & !
  zpir18 = 0.0174532925_ireals

REAL rla2rlarot
!-------------------------------------------------------------------------------
! Begin Subroutine uv2uvrot
!-------------------------------------------------------------------------------

! Converting from degree to radians
  zpolphi = polphi * zpir18
  zpollam = pollam * zpir18
  zrla    = rla    * zpir18
  zphi    = phi    * zpir18
  pollamd = pollam
  IF (pollam < 0.0) pollamd = 360.0_ireals + pollam

! Longitude in the rotated system
  zrlas   = rla2rlarot (phi, rla, polphi, pollam) * zpir18

! Calculate the angle between the latitudes
  zarg   = - SIN (zpolphi) * SIN(zrla-zpollam) * SIN(zrlas)              &
                           - COS(zrla-zpollam) * COS(zrlas)
  zarg   = MIN ( 1.0_ireals, zarg)
  zarg   = MAX (-1.0_ireals, zarg)
  zbeta  = ACOS(zarg)
  zbeta  = SIGN (zbeta, - (rla - (pollamd-180.0_ireals)))

! Convert the u-component
  urot   =  u * COS(zbeta) + v * SIN(zbeta)

!  Convert the v-component
  vrot   = -u * SIN(zbeta) + v * COS(zbeta)

END SUBROUTINE uv2uvrot

!******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------

SUBROUTINE uv2df (u, v, d, f)

!-------------------------------------------------------------------------------
!
! Description:
!   This routine computes wind speed and wind direction from the wind
!   components.
!
! Method:
!   Straightforward.
!
!-------------------------------------------------------------------------------
!
INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers
!=======================================================================

IMPLICIT NONE

! Parameter list:
REAL (KIND=ireals), INTENT (IN)          ::    &
  u   , v           ! wind components in the true geographical system

REAL (KIND=ireals), INTENT (OUT)         ::    &
  f   , d           ! wind speed and wind direction

! Local variables

REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & ! conversion from radians to degrees
  zsmall = 0.001_ireals

!-------------------------------------------------------------------------------
! Begin Subroutine uv2df
!-------------------------------------------------------------------------------

  IF (ABS(u) > zsmall) THEN
    f  =  SQRT( u*u + v*v )
    d  =  v / u
    d  =  180.0_ireals + SIGN( 90.0_ireals , u ) - ATAN( d ) *zrpi18
  ELSEIF (ABS(v) > zsmall) THEN
    f  =  ABS( v )
    d  =  270.0_ireals - SIGN( 90.0_ireals , v )
  ELSE
    f  =    0.0_ireals
    d  =    0.0_ireals
  ENDIF

END SUBROUTINE uv2df

!******************************************************************************
! Transformationsprogramme aus METRAS
 
      SUBROUTINE OGEOGK(AB,AL,AH,AR,IREFS)
!
!---------------------------------------------------------------------
!!    *OGEOGK* : Conversion of geographic coordinates to
!!               Gauss-Krueger-coordinates
!
!     Input:  geographic coordinates AB and AL; IREFS
!             AB    = latitude  [degree]
!             AL    = longitude [degree]
!             IREFS = number of selected Gauss-Krueger-sector
!     Output: Gauss-Krueger-coordinates AH and AR
!             AH = Hochwert   [distance from equator in km]
!             AR = Rechtswert [sector-number and distance from
!                              corresponding longitude in km]
!             Example: AH = 5926. and AR = 3553. means
!                      5926.km distance from equator, corresponding
!                      longitude 9 degrees east ( = first number*3 )
!                      and 53.km east ( = 553.-500. ) from
!                      corresponding longitude
!
!     ATTENTION: Normally IREFS should have the value "-999". That means
!     ========== OGEOGK will calculate the correct Gauss-Krueger-sector
!                from the given longitude AL.
!                For special purpose the user may wish Gauss-Krueger-
!                coordinates with respect to other median longitudes 
!                (=other sectors). Then IREFS has to be the number of
!                this preferred sector.
!                Example: AL = 7. [deg] --> sector number 2 with
!                         median longitude BLH = 6. [deg].
!                         But Gauss-Krueger-coordinates are wanted with
!                         respect to sector number 3 (median longitude
!                         BLH = 9. [deg]). Then the user has to choose
!                         IREFS = 3 !
!                It should be noticed that the error of the calculated
!                Gauss-Krueger-coordinates increases with increasing
!                distance from the median longitude ! So coordinates
!                are less exact within others than the default sectors !!
!---------------------------------------------------------------------
!
!     CHANGES :
!     =========
!     NAME         DATE         COMMENT
!     ----         ----         -------
!     Nikolaus     06/12/91     new subroutine
!     Klaus        23/11/92     others than default sectors possible
!     Klaus        30/01/93     number of possible sectors expanded
!     Jerzy        09/06/94     IMPLICIT NONE. if the latitude is negative
!                               (southern from equator), the AH coordinate
!                               will be also negative. new variable KNEGAT.
!
!     METHOD :
!     ========
!
!     see: A.Schoedlbauer, Rechenformeln und Rechenbeispiele zur
!                          Landvermessung, Karlsruhe 1982, pp.33-41,53-54,
!                          93-101,113-115.
!          W.Grossmann, Geod. Rechnungen fuer die Landvermessung,
!                       Stuttgart 1964, pp.16-19.
!
!     TAPES :
!     =======
!     NAME         READ/WRITE
!     ----         ----------
!
!     --- no ---
!
!     SUBROUTINES :
!     =============
!
!     --- no ---
!
!     COMMON :
!     ========
!
!     --- no ---
!
!     ADDITIONAL VARIABLES AND CONSTANTS:
!     ===================================
!
!     Name      Type     Meaning
!     ----      ----     -------
!     AB        R        geographic latitude [degree]
!     AH        R        "Hochwert" (=north) of Gauss-Krueger-coordinates
!     AL        R        geographic longitude [degree]
!     AR        R        "Rechtswert" (=east) of Gauss-Krueger-coordinates
!     IREFS     I        number of selected Gauss-Krueger-sector
!                        = -999:  no sector selected, correct sector will 
!                                 be calculated from given longitudes
!                        = value: number of selected Gauss-Krueger-sector,
!                                 median longitude of selected sector BLH
!                                 is calculated by BLH = 3.*"value"
!     BB0       R        rounded value of latitude
!     BCPI      R        Pi
!     BDELTB    R        difference of latitude to rounded value
!     BDELTL    R        difference of longitude to rounded value
!     BDELTX    R        difference of "Hochwert" to rounded value
!     BETA0     R
!     BFXnm     R        Koeffizienten der Reihenentwicklung fuer x
!     BFYnm     R        Koeffizienten der Reihenentwicklung fuer y
!     BH0       R        rounded value of "Hochwert"
!     BLH       R        rounded value of longitude
!     BN0       R        Querkruemmungshalbmesser
!     BPCQU     R        Polkruemmungshalbmesser
!     BPE       R        Quadrat der 2. numerischen Exzentritaet
!     BPK       R        Zuschlag in Ostrichtung fuer G.-K.-Koord.
!     BPRHO     R        Umrechnungsfaktor
!     BSECT     R        index for section of longitude
!     BT0       R
!     BY        R        difference of "Rechtswert" to rounded value
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
      REAL BPCQU,BPE,BPK,BPRHO
      PARAMETER ( BPCQU = 6398786.849          &
                , BPE   = 0.006719219          &
                , BPK   = 500.                 &
                , BPRHO = 57.2957795131 )
!
      INTEGER IREFS
      REAL    AB,AH,AL,AR
!
      REAL BB0,BCPI,BDELTB,BDELTL,BDELTX,BETA0,BH0,BLH,BN0   &
          ,BSECT,BT0,BY
!
      REAL BFX02,BFX04,BFX10,BFX12,BFX14,BFX20,BFX22,BFX24   &
          ,BFX30,BFX32,BFX40                                 &
          ,BFY01,BFY03,BFY05,BFY11,BFY13,BFY15,BFY21,BFY23   &
          ,BFY31,BFY33,BFY41
!--- whether latitude is negative (southern)
      LOGICAL KNEGAT
!-----------------------------------------------------------------------
!
!---  defining PI
      BCPI   = 4.*ATAN(1.)
!
!--- check, if latitude is negative
      KNEGAT = AB .LT. 0.0
      IF (KNEGAT) AB = ABS (AB)
!
!---  rounding of latitude and longitude; median longitude of sector
      BB0    = FLOAT( NINT(AB) )
      BDELTB = AB - BB0
      IF ( IREFS.EQ.-999 ) THEN
         IF ( AL.GE.1.5  .AND.  AL.LT.4.5 ) THEN
            BLH = 3.
         ELSE IF ( AL.GE.4.5  .AND.  AL.LT.7.5 ) THEN
            BLH = 6.
         ELSE IF ( AL.GE.7.5  .AND.  AL.LT.10.5 ) THEN
            BLH = 9.
         ELSE IF ( AL.GE.10.5  .AND.  AL.LT.13.5 ) THEN
            BLH = 12.
         ELSE IF ( AL.GE.13.5  .AND.  AL.LT.16.5 ) THEN
            BLH = 15.
         ELSE IF ( AL.GE.16.5  .AND.  AL.LT.19.5 ) THEN
            BLH = 18.
         ELSE
            STOP 'OGEOGK: given longitude not implemented'
         END IF
         IREFS = INT( BLH/3. + 1.E-3 )
      ELSE
         BLH = 3.*FLOAT(IREFS)
         IF ( BLH.LT.3.  .OR.  BLH.GT.18. ) THEN
            PRINT*,'error OGEOGK: given sector not implemented'
            PRINT*,'sector number:                  IREFS = ',IREFS
            PRINT*,'corresponding median longitude: BLH   = ',BLH
            PRINT*,' '
            PRINT*,'for default sector please choose IREFS = -999'
            STOP 'error OGEOGK: given sector not implemented'
         END IF
      END IF
      BDELTL = AL - BLH
!
!---  sector of latitude between 47 and 55 degree?
      IF ( AB.GT.47.  .AND.  AB.LE.55. ) THEN
         BSECT = 0.
      ELSE
         BSECT = 1.
      ENDIF
!-----------------------------------------------------------------------
!
!---  calculation of constants
      BT0   = TAN(BB0*BCPI/180.)
      BETA0 = SQRT( BPE * (COS(BB0*BCPI/180.))**2. )
      BN0   = BPCQU / SQRT( 1. + BETA0**2. )
!
!---  calculation of coefficients
      BFX02 = ( 1./(2.*BPRHO**2.) )                   &
               * BN0 * (COS(BB0*BCPI/180.))**2. * BT0
      BFX04 = ( 1./(24.*BPRHO**4.) )                  &
               * BN0 * (COS(BB0*BCPI/180.))**4. * BT0  &
               * ( 5. - BT0**2. + 9.*BETA0**2. )
      BFX10 = ( 1./BPRHO )                                &
               * BN0                                       &
               * ( 1. - BETA0**2. + BETA0**4. - BETA0**6. )
      BFX12 = ( 1./(2.*BPRHO**3.) )                         &
               * BN0 * (COS(BB0*BCPI/180.))**2.             &
               * ( 1. - BT0**2. + BETA0**2.*BT0**2.         &
                      - BETA0**4.*BT0**2.           )
      BFX14 = ( 1./(24.*BPRHO**5.) )                        &
               * BN0 * (COS(BB0*BCPI/180.))**4.             &
               * ( 5. - 18.*BT0**2. + BT0**4.               &
                  + BSECT*BETA0**2.*(9.-40.*BT0**2.-BT0**4.) )
      BFX20 = ( 3./(2.*BPRHO**2.) )                        &
               * BN0 * BT0                                 &
               * ( BETA0**2. - 2.*BETA0**4. )
      BFX22 = ( 1./(4.*BPRHO**4. ) )                       &
               * BN0 * (COS(BB0*BCPI/180.))**2. * BT0      &
               * ( -4. + 3.*BETA0**2.*(1.-BT0**2.) )       
      BFX24 = ( 1./(6.*BPRHO**6.) )                        &
               * BN0 * (COS(BB0*BCPI/180.))**4. * BT0      &
               * ( -7. + 5.*BT0**2. )
      BFX30 = ( 1./(2.*BPRHO**3.) )                        &
               * BN0                                       &
               * ( BETA0**2.*(1.-BT0**2.) + BETA0**4.*(-2.+7.*BT0**2.) )
      BFX32 = ( 1./(12.*BPRHO**5.) )                       &
               * BN0 * (COS(BB0*BCPI/180.))**2.            &
               * ( -4. + 4.*BT0**2. )
      BFX40 = -( 1./(2.*BPRHO**4.) )                       &
               * BN0 * BETA0**2. * BT0
!
      BFY01 = ( 1./BPRHO )                                 &
               * BN0 * COS(BB0*BCPI/180.)
      BFY03 = ( 1./(6.*BPRHO**3.) )                        &
               * BN0 * (COS(BB0*BCPI/180.))**3.            &
               * ( 1. - BT0**2. + BETA0**2. )
      BFY05 = ( 1./(120.*BPRHO**5.) )                      &
               * BN0 * (COS(BB0*BCPI/180.))**5.            &
               * ( 5. - 18.*BT0**2. + BT0**4. )
      BFY11 = ( 1./BPRHO**2. )                             &
               * BN0 * COS(BB0*BCPI/180.) * BT0            &
               * ( -1. + BETA0**2. - BETA0**4. + BETA0**6. )
      BFY13 = ( 1./(6.*BPRHO**4.) )                        &
               * BN0 * (COS(BB0*BCPI/180.))**3. * BT0      &
               * ( -5. + BT0**2. - BETA0**2.*(4.+BT0**2.) )
      BFY15 = ( 1./(120.*BPRHO**6.) )                      &
               * BN0 * (COS(BB0*BCPI/180.))**5. * BT0      &
               * ( -61. + 58.*BT0**2. - BT0**4. )
      BFY21 = ( 1./(2.*BPRHO**3.) )                        &
               * BN0 * COS(BB0*BCPI/180.)                  &
               * ( -1. + BETA0**2.*( 1.-3.*BT0**2.)        &
                       + BETA0**4.*(-1.+6.*BT0**2.) )
      BFY23 = ( 1./(12.*BPRHO**5.) )                      &
               * BN0 * (COS(BB0*BCPI/180.))**3.           &
               * ( -5. + 13.*BT0**2. )
      BFY31 = ( 1./(6.*BPRHO**4.) )                       &
               * BN0 * COS(BB0*BCPI/180.) * BT0           &
               * ( 1. + BETA0**2.*(-10.+3.*BT0**2.) )
      BFY33 = ( 1./(36.*BPRHO**6.) )                      &
               * BN0 * (COS(BB0*BCPI/180.))**3. * BT0     &
               * ( 41. - 13.*BT0**2. )
      BFY41 = ( 1./(24.*BPRHO**5.) )                      &
               * BN0 * COS(BB0*BCPI/180.)
!-----------------------------------------------------------------------
!
!---  calculation of meridian arc at main longitude
      BH0 =  111120.61962 * BB0                           &
          -  15988.63853 * SIN(2.*BB0*BCPI/180.)         &
          +     16.72995 * SIN(4.*BB0*BCPI/180.)         &
          -      0.02178 * SIN(6.*BB0*BCPI/180.)         &
          +      0.00003 * SIN(8.*BB0*BCPI/180.)         
!                                                         
!---  calculation of differences of coordinates to main meridian
!---  and rounded latitude
      BDELTX =        BFX10*BDELTB                        &
             +       BFX20*BDELTB**2.                    &
             +       BFX30*BDELTB**3.                    &
             +       BFX40*BDELTB**4.                    &
             +       BFX02           *BDELTL**2.         &
             +       BFX12*BDELTB    *BDELTL**2.         &
             +       BFX22*BDELTB**2.*BDELTL**2.         &
             +       BFX32*BDELTB**3.*BDELTL**2.         &
             +       BFX04           *BDELTL**4.         &
             +       BFX14*BDELTB    *BDELTL**4.         &
             + BSECT*BFX24*BDELTB**2.*BDELTL**4.         
!                                                         
      BY =        BFY01           *BDELTL                &
        +       BFY11*BDELTB    *BDELTL                &
        +       BFY21*BDELTB**2.*BDELTL                &
        +       BFY31*BDELTB**3.*BDELTL                &
        +       BFY41*BDELTB**4.*BDELTL                &
        +       BFY03           *BDELTL**3.            &
        +       BFY13*BDELTB    *BDELTL**3.            &
        +       BFY23*BDELTB**2.*BDELTL**3.            &
        +       BFY33*BDELTB**3.*BDELTL**3.            &
        +       BFY05           *BDELTL**5.            &
        + BSECT*BFY15*BDELTB    *BDELTL**5.            
!                                                        
!---  conversion of coordinates from [m] to [km]
      BH0    = BH0    / 1000.
      BDELTX = BDELTX / 1000.
      BY     = BY     / 1000.
!
!---  calculation of Gauss-Krueger-coordinates
      AR = FLOAT(IREFS)*1.E3 + BPK + BY
      AH = BH0 + BDELTX
!
!--- if latitude is negative, correct AH coordinate
      IF (KNEGAT) AH = -AH
!
      RETURN
      END     
!----------------------------------------------------------------
      SUBROUTINE OGKGEO (AH, AR, AB, AL)
!
!---------------------------------------------------------------------
!C    *OGKGEO* : Conversion of Gauss-Krueger-coordinates to
!C               geographic coordinates
!
!     Input:  Gauss-Krueger-coordinates AH and AR
!             AH = Hochwert   [distance from equator in km]
!             AR = Rechtswert [sector-number and distance from
!                              corresponding longitude in km]
!             Example: AH = 5926. and AR = 3553. means
!                      5926.km distance from equator, corresponding
!                      longitude 9 degrees east ( = first number*3 )
!                      and 53.km east ( = 553.-500. ) from
!                      corresponding longitude
!     Output: geographic coordinates AB and AL
!             AB = latitude  [degree]
!             AL = longitude [degree]
!---------------------------------------------------------------------
!
!     CHANGES :
!     =========
!     NAME         DATE         COMMENT
!     ----         ----         -------
!     Nikolaus     06/12/91     new subroutine
!     Jerzy        29/06/93     IMPLICIT NONE
!     Jerzy        09/06/94     if AH coordinate is negative (southern
!                               from equator), latitude will be also
!                               negative. new variable KNEGAT.
!
!     METHOD :
!     ========
!
!     see: A.Schoedlbauer, Rechenformeln und Rechenbeispiele zur
!                          Landvermessung, Karlsruhe 1982, pp.33-41,53-54,
!                          93-101,113-115.
!          W.Grossmann, Geod. Rechnungen fuer die Landvermessung,
!                       Stuttgart 1964, pp.16-19.
!
!     TAPES :
!     =======
!     NAME         READ/WRITE
!     ----         ----------
!
!     --- no ---
!
!     SUBROUTINES :
!     =============
!
!     --- no ---
!
!     COMMON :
!     ========
!
!     --- no ---
!
!     ADDITIONAL VARIABLES AND CONSTANTS:
!     ===================================
!
!     Name      Type     Meaning
!     ----      ----     -------
!     AB        R        geographic latitude [degree]
!     AH        R        "Hochwert" (=north) of Gauss-Krueger-coordinates
!     AL        R        geographic longitude [degree]
!     AR        R        "Rechtswert" (=east) of Gauss-Krueger-coordinates
!     BB0       R        rounded value of latitude
!     BCPI      R        Pi
!     BDELTB    R        difference of latitude to rounded value
!     BDELTL    R        difference of longitude to rounded value
!     BDELTX    R        difference of "Hochwert" to rounded value
!     BETA0     R
!     BFBnm     R        Koeffizienten der Reihenentwicklung fuer B
!     BFLnm     R        Koeffizienten der Reihenentwicklung fuer L
!     BH0       R        rounded value of "Hochwert"
!     BLH       R        rounded value of longitude
!     BN0       R        Querkruemmungshalbmesser
!     BPCQU     R        Polkruemmungshalbmesser
!     BPE       R        Quadrat der 2. numerischen Exzentritaet
!     BPK       R        Zuschlag in Ostrichtung fuer G.-K.-Koord.
!     BR0       R        "coordinate" of main longitude G.K.-coord.
!     BPRHO     R        Umrechnungsfaktor
!     BSECT     R        index for section of longitude
!     BT0       R
!     BY        R        difference of "Rechtswert" to rounded value
!     JJ        I        loop index
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
      REAL BPCQU,BPE,BPK,BPRHO
      PARAMETER ( BPCQU = 6398786.849            &
                , BPE   = 0.006719219            &
                , BPK   = 500.                   &
                , BPRHO = 57.2957795131 )
!
      REAL AB,AH,AL,AR
!
      INTEGER JJ
      REAL BB0,BCPI,BDELTB,BDELTL,BDELTX,BETA0,BH0,BLH,BN0,BR0   &
          ,BSECT,BT0,BY
!
      REAL BFB02,BFB04,BFB06,BFB10,BFB12,BFB14,BFB20,BFB22,BFB24  &
          ,BFB30,BFB32,BFB34,BFB40,BFB42                          &
          ,BFL01,BFL03,BFL05,BFL11,BFL13,BFL15,BFL21,BFL23,BFL25  &
          ,BFL31,BFL33,BFL41
!
!--- whether AH is negative
      LOGICAL KNEGAT
!-----------------------------------------------------------------------
!
!---  defining PI
      BCPI   = 4.*ATAN(1.)
!
!--- checking the AH coordinate
      KNEGAT = AH .LT. 0.0
      IF (KNEGAT) AH = ABS (AH)
!
!---  calculation of rounded value of latitude
      BB0 = -0.5
      DO 100 JJ=0,180
         BB0 = BB0 + 0.5
         BH0 =  111120.61962 * BB0                       &
              -  15988.63853 * SIN(2.*BB0*BCPI/180.)     &
              +     16.72995 * SIN(4.*BB0*BCPI/180.)     &
              -      0.02178 * SIN(6.*BB0*BCPI/180.)     &
              +      0.00003 * SIN(8.*BB0*BCPI/180.)
         IF ( (AH*1000.-BH0).LT.0. ) GOTO 101
  100 CONTINUE
  101 CONTINUE
      DO 200 JJ=1,5
         BB0 = BB0 - 0.1
         BH0 =  111120.61962 * BB0                       &
              -  15988.63853 * SIN(2.*BB0*BCPI/180.)     &
              +     16.72995 * SIN(4.*BB0*BCPI/180.)     &
              -      0.02178 * SIN(6.*BB0*BCPI/180.)     &
              +      0.00003 * SIN(8.*BB0*BCPI/180.)
         IF ( (AH*1000.-BH0).GT.0. ) GOTO 201
  200 CONTINUE
  201 CONTINUE
      DO 300 JJ=1,10
         BB0 = BB0 + 0.01
         BH0 =  111120.61962 * BB0                      &
              -  15988.63853 * SIN(2.*BB0*BCPI/180.)    &
              +     16.72995 * SIN(4.*BB0*BCPI/180.)    &
              -      0.02178 * SIN(6.*BB0*BCPI/180.)    &
              +      0.00003 * SIN(8.*BB0*BCPI/180.)
         IF ( (AH*1000.-BH0).LT.0. ) GOTO 301
  300 CONTINUE
  301 CONTINUE
      IF ( ABS(AH*1000.-BH0).GT.50000. ) THEN
         PRINT*, 'AH  = ',AH
         PRINT*, 'BH0 = ',BH0
         PRINT*, 'BB0 = ',BB0
         STOP 'OGKGEO: no convergence for BB0'
      END IF
!
!---  differences of coordinates x and y
      BDELTX = AH - BH0/1000.
      BR0    = FLOAT( INT(AR/1.E3) )*1.E3 + BPK
      BY     = AR - BR0
!
!---  conversion of coordinates from [km] to [m]
      BDELTX = BDELTX * 1000.
      BR0    = BR0    * 1000.
      BY     = BY     * 1000.
!
!---  sector of latitude between 47 and 55 degree?
      IF ( BB0.GT.47.  .AND.  BB0.LE.55. ) THEN
         BSECT = 0.
      ELSE
         BSECT = 1.
      ENDIF
!-----------------------------------------------------------------------
!
!---  calculation of constants
      BT0   = TAN(BB0*BCPI/180.)
      BETA0 = SQRT( BPE * (COS(BB0*BCPI/180.))**2. )
      BN0   = BPCQU / SQRT( 1. + BETA0**2. )
!
!---  calculation of coefficients
      BFB02 = ( -BPRHO/(2.*BN0**2.) ) * BT0           &
               * ( 1. + BETA0**2. )
      BFB04 = ( BPRHO/(24.*BN0**4.) ) * BT0           &
               * ( 5. + 3.*BT0**2. + BETA0**2.*(6.-6.*BT0**2.) )
      BFB06 = ( -BPRHO/(720.*BN0**6.) ) * BT0         &
               * ( 61. + 90.*BT0**2. + 45.*BT0**4. )
      BFB10 = ( BPRHO/BN0 )                           &
               * ( 1. + BETA0**2. )
      BFB12 = ( BPRHO/(2.*BN0**3.) )                  &
               * ( -1. - BT0**2.                      &
                  + BETA0**2.*(-2.+2.*BT0**2.)        &
                  + BETA0**4.*(-1.+3.*BT0**2.) )
      BFB14 = ( BPRHO/(24.*BN0**5.) )                 &
               * ( 5. + 14.*BT0**2. + 9.*BT0**4.      &
                  + BSECT*BETA0**2.*(11.-30.*BT0**2.-9.*BT0**4.) )
      BFB20 = ( -3.*BPRHO/(2.*BN0**2.) ) * BT0        &
               * ( BETA0**2. + BETA0**4. )
      BFB22 = ( BPRHO/(4.*BN0**4.) ) * BT0            &
               * ( -2. - 2.*BT0**2. + BETA0**2.*(9.+BT0**2.) )

      BFB24 = ( BPRHO/(12.*BN0**6.) ) * BT0           &
               * ( 7. + 16.*BT0**2. + 9.*BT0**4. )
      BFB30 = ( BPRHO/(2.*BN0**3.) )                  &
               * ( BETA0**2.*(-1.+BT0**2.)            &
                 + BETA0**4.*(-2.+6.*BT0**2.) )
      BFB32 = ( -BPRHO/(6.*BN0**5.) )                 &
               * ( 1. + 4.*BT0**2. + 3.*BT0**4. )
      BFB34 = BSECT                                   &
              * ( BPRHO/(36.*BN0**7.) )               &
              * ( 7. + 55.*BT0**2. + 93.*BT0**4. + 45.*BT0**6. )
      BFB40 = ( BPRHO/(2.*BN0**4.) ) * BT0 * BETA0**2.
      BFB42 = ( -BPRHO/(6.*BN0**6.) ) * BT0           &
               * ( 2. + 5.*BT0 + 3.*BT0**4. )
!
      BFL01 = ( BPRHO/(BN0*COS(BB0*BCPI/180.)) )
      BFL03 = ( -BPRHO/(6.*BN0**3.*COS(BB0*BCPI/180.)) )  &
               * ( 1. + 2.*BT0**2. + BETA0**2. )
      BFL05 = ( BPRHO/(120.*BN0**5.*COS(BB0*BCPI/180.)) ) &
               * ( 5. + 28.*BT0**2. + 24.*BT0**4. )
      BFL11 = ( BPRHO/(BN0**2.*COS(BB0*BCPI/180.)) )      &
               * BT0
      BFL13 = ( -BPRHO/(6.*BN0**4.*COS(BB0*BCPI/180.)) )  &
               * BT0                                      &
               * ( 5. + 6.*BT0**2. + BETA0**2. )
      BFL15 = ( BPRHO/(120.*BN0**6.*COS(BB0*BCPI/180.)) ) &
               * BT0                                      &
               * ( 61. + 180.*BT0**2. + 120.*BT0**4. )
      BFL21 = ( BPRHO/(2.*BN0**3.*COS(BB0*BCPI/180.)) )   &
               * ( 1. + 2.*BT0**2. + BETA0**2. )
      BFL23 = ( -BPRHO/(12.*BN0**5.*COS(BB0*BCPI/180.)) ) &
               * ( 5. + 28.*BT0**2. + 24.*BT0**4. )
      BFL25 = BSECT                                        &
              * ( BPRHO/(240.*BN0**7.*COS(BB0*BCPI/180.)) ) &
              * ( 61. + 662.*BT0**2. + 1320.*BT0**4. + 720.*BT0**6. )
      BFL31 = ( BPRHO/(6.*BN0**4.*COS(BB0*BCPI/180.)) )  &
               * BT0                                     &
               * ( 5. + 6.*BT0**2. + BETA0**2. )
      BFL33 = ( -BPRHO/(36.*BN0**6.*COS(BB0*BCPI/180.)) )  &
               * BT0                                       &
               * ( 61. + 180.*BT0**2. + 120.*BT0**4. )
      BFL41 = ( BPRHO/(24.*BN0**5.*COS(BB0*BCPI/180.)) )   &
               * ( 5. + 28.*BT0**2. + 24.*BT0**4. )
!-----------------------------------------------------------------------
!
!---  calculation of differences of geographic coordinates
      BDELTB =        BFB10*BDELTX                       &
             +       BFB20*BDELTX**2.                   &
             +       BFB30*BDELTX**3.                   &
             +       BFB40*BDELTX**4.                   &
             +       BFB02           *BY**2.            &
             +       BFB12*BDELTX    *BY**2.            &
             +       BFB22*BDELTX**2.*BY**2.            &
             +       BFB32*BDELTX**3.*BY**2.            &
             +       BFB42*BDELTX**4.*BY**2.            &
             +       BFB04           *BY**4.            &
             +       BFB14*BDELTX    *BY**4.            &
             +       BFB24*BDELTX**2.*BY**4.            &
             + BSECT*BFB34*BDELTX**3.*BY**4.            &
             +       BFB06           *BY**6.            
!
      BDELTL =        BFL01           *BY               &
             +       BFL11*BDELTX    *BY               &
             +       BFL21*BDELTX**2.*BY               &
             +       BFL31*BDELTX**3.*BY               &
             +       BFL41*BDELTX**4.*BY               &
             +       BFL03           *BY**3.           &
             +       BFL13*BDELTX    *BY**3.           &
             +       BFL23*BDELTX**2.*BY**3.           &
             +       BFL33*BDELTX**3.*BY**3.           &
             +       BFL05           *BY**5.           &
             +       BFL15*BDELTX    *BY**5.           &
             + BSECT*BFL25*BDELTX**2.*BY**5.           
!
!---  calculation of geographic coordinates
      BLH = 3. * FLOAT( INT(AR/1000.) )
      AB  = BB0 + BDELTB
      AL  = BLH + BDELTL
!
!--- if AH coordinate is negative, latitude must be also negative
      IF (KNEGAT) AB = -AB
!
      RETURN
      END
