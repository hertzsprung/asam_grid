      SUBROUTINE OGKGEO (AH, AR, AB, AL)
C
C---------------------------------------------------------------------
CC    *OGKGEO* : Conversion of Gauss-Krueger-coordinates to
CC               geographic coordinates
C
C     Input:  Gauss-Krueger-coordinates AH and AR
C             AH = Hochwert   [distance from equator in km]
C             AR = Rechtswert [sector-number and distance from
C                              corresponding longitude in km]
C             Example: AH = 5926. and AR = 3553. means
C                      5926.km distance from equator, corresponding
C                      longitude 9 degrees east ( = first number*3 )
C                      and 53.km east ( = 553.-500. ) from
C                      corresponding longitude
C     Output: geographic coordinates AB and AL
C             AB = latitude  [degree]
C             AL = longitude [degree]
C---------------------------------------------------------------------
C
C     CHANGES :
C     =========
C     NAME         DATE         COMMENT
C     ----         ----         -------
C     Nikolaus     06/12/91     new subroutine
C     Jerzy        29/06/93     IMPLICIT NONE
C     Jerzy        09/06/94     if AH coordinate is negative (southern
C                               from equator), latitude will be also
C                               negative. new variable KNEGAT.
C
C     METHOD :
C     ========
C
C     see: A.Schoedlbauer, Rechenformeln und Rechenbeispiele zur
C                          Landvermessung, Karlsruhe 1982, pp.33-41,53-54,
C                          93-101,113-115.
C          W.Grossmann, Geod. Rechnungen fuer die Landvermessung,
C                       Stuttgart 1964, pp.16-19.
C
C     TAPES :
C     =======
C     NAME         READ/WRITE
C     ----         ----------
C
C     --- no ---
C
C     SUBROUTINES :
C     =============
C
C     --- no ---
C
C     COMMON :
C     ========
C
C     --- no ---
C
C     ADDITIONAL VARIABLES AND CONSTANTS:
C     ===================================
C
C     Name      Type     Meaning
C     ----      ----     -------
C     AB        R        geographic latitude [degree]
C     AH        R        "Hochwert" (=north) of Gauss-Krueger-coordinates
C     AL        R        geographic longitude [degree]
C     AR        R        "Rechtswert" (=east) of Gauss-Krueger-coordinates
C     BB0       R        rounded value of latitude
C     BCPI      R        Pi
C     BDELTB    R        difference of latitude to rounded value
C     BDELTL    R        difference of longitude to rounded value
C     BDELTX    R        difference of "Hochwert" to rounded value
C     BETA0     R
C     BFBnm     R        Koeffizienten der Reihenentwicklung fuer B
C     BFLnm     R        Koeffizienten der Reihenentwicklung fuer L
C     BH0       R        rounded value of "Hochwert"
C     BLH       R        rounded value of longitude
C     BN0       R        Querkruemmungshalbmesser
C     BPCQU     R        Polkruemmungshalbmesser
C     BPE       R        Quadrat der 2. numerischen Exzentritaet
C     BPK       R        Zuschlag in Ostrichtung fuer G.-K.-Koord.
C     BR0       R        "coordinate" of main longitude G.K.-coord.
C     BPRHO     R        Umrechnungsfaktor
C     BSECT     R        index for section of longitude
C     BT0       R
C     BY        R        difference of "Rechtswert" to rounded value
C     JJ        I        loop index
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
      REAL(8) BPCQU,BPE,BPK,BPRHO
      PARAMETER ( BPCQU = 6398786.849
     1          , BPE   = 0.006719219
     2          , BPK   = 500.
     3          , BPRHO = 57.2957795131 )
C
      REAL(8) AB,AH,AL,AR
C
      INTEGER JJ
      REAL(8) BB0,BCPI,BDELTB,BDELTL,BDELTX,BETA0,BH0,BLH,BN0,BR0
     1    ,BSECT,BT0,BY
C
      REAL(8) BFB02,BFB04,BFB06,BFB10,BFB12,BFB14,BFB20,BFB22,BFB24
     1    ,BFB30,BFB32,BFB34,BFB40,BFB42
     2    ,BFL01,BFL03,BFL05,BFL11,BFL13,BFL15,BFL21,BFL23,BFL25
     3    ,BFL31,BFL33,BFL41
C
C--- whether AH is negative
      LOGICAL KNEGAT
C-----------------------------------------------------------------------
C
C---  defining PI
      BCPI   = 4.*ATAN(1.)
C
C--- checking the AH coordinate
      KNEGAT = AH .LT. 0.0
      IF (KNEGAT) AH = ABS (AH)
C
C---  calculation of rounded value of latitude
      BB0 = -0.5
      DO 100 JJ=0,180
         BB0 = BB0 + 0.5
         BH0 =  111120.61962 * BB0
     1        -  15988.63853 * SIN(2.*BB0*BCPI/180.)
     2        +     16.72995 * SIN(4.*BB0*BCPI/180.)
     3        -      0.02178 * SIN(6.*BB0*BCPI/180.)
     4        +      0.00003 * SIN(8.*BB0*BCPI/180.)
         IF ( (AH*1000.-BH0).LT.0. ) GOTO 101
  100 CONTINUE
  101 CONTINUE
      DO 200 JJ=1,5
         BB0 = BB0 - 0.1
         BH0 =  111120.61962 * BB0
     1        -  15988.63853 * SIN(2.*BB0*BCPI/180.)
     2        +     16.72995 * SIN(4.*BB0*BCPI/180.)
     3        -      0.02178 * SIN(6.*BB0*BCPI/180.)
     4        +      0.00003 * SIN(8.*BB0*BCPI/180.)
         IF ( (AH*1000.-BH0).GT.0. ) GOTO 201
  200 CONTINUE
  201 CONTINUE
      DO 300 JJ=1,10
         BB0 = BB0 + 0.01
         BH0 =  111120.61962 * BB0
     1        -  15988.63853 * SIN(2.*BB0*BCPI/180.)
     2        +     16.72995 * SIN(4.*BB0*BCPI/180.)
     3        -      0.02178 * SIN(6.*BB0*BCPI/180.)
     4        +      0.00003 * SIN(8.*BB0*BCPI/180.)
         IF ( (AH*1000.-BH0).LT.0. ) GOTO 301
  300 CONTINUE
  301 CONTINUE
      IF ( ABS(AH*1000.-BH0).GT.50000. ) THEN
         PRINT*, 'AH  = ',AH
         PRINT*, 'BH0 = ',BH0
         PRINT*, 'BB0 = ',BB0
         STOP 'OGKGEO: no convergence for BB0'
      END IF
C
C---  differences of coordinates x and y
      BDELTX = AH - BH0/1000.
      BR0    = FLOAT( INT(AR/1.E3) )*1.E3 + BPK
      BY     = AR - BR0
C
C---  conversion of coordinates from [km] to [m]
      BDELTX = BDELTX * 1000.
      BR0    = BR0    * 1000.
      BY     = BY     * 1000.
C
C---  sector of latitude between 47 and 55 degree?
      IF ( BB0.GT.47.  .AND.  BB0.LE.55. ) THEN
         BSECT = 0.
      ELSE
         BSECT = 1.
      ENDIF
C-----------------------------------------------------------------------
C
C---  calculation of constants
      BT0   = TAN(BB0*BCPI/180.)
      BETA0 = SQRT( BPE * (COS(BB0*BCPI/180.))**2. )
      BN0   = BPCQU / SQRT( 1. + BETA0**2. )
C
C---  calculation of coefficients
      BFB02 = ( -BPRHO/(2.*BN0**2.) ) * BT0
     1         * ( 1. + BETA0**2. )
      BFB04 = ( BPRHO/(24.*BN0**4.) ) * BT0
     1         * ( 5. + 3.*BT0**2. + BETA0**2.*(6.-6.*BT0**2.) )
      BFB06 = ( -BPRHO/(720.*BN0**6.) ) * BT0
     1         * ( 61. + 90.*BT0**2. + 45.*BT0**4. )
      BFB10 = ( BPRHO/BN0 )
     1         * ( 1. + BETA0**2. )
      BFB12 = ( BPRHO/(2.*BN0**3.) )
     1         * ( -1. - BT0**2.
     2            + BETA0**2.*(-2.+2.*BT0**2.)
     3            + BETA0**4.*(-1.+3.*BT0**2.) )
      BFB14 = ( BPRHO/(24.*BN0**5.) )
     1         * ( 5. + 14.*BT0**2. + 9.*BT0**4.
     2            + BSECT*BETA0**2.*(11.-30.*BT0**2.-9.*BT0**4.) )
      BFB20 = ( -3.*BPRHO/(2.*BN0**2.) ) * BT0
     1         * ( BETA0**2. + BETA0**4. )
      BFB22 = ( BPRHO/(4.*BN0**4.) ) * BT0
     1         * ( -2. - 2.*BT0**2. + BETA0**2.*(9.+BT0**2.) )

      BFB24 = ( BPRHO/(12.*BN0**6.) ) * BT0
     1         * ( 7. + 16.*BT0**2. + 9.*BT0**4. )
      BFB30 = ( BPRHO/(2.*BN0**3.) )
     1         * ( BETA0**2.*(-1.+BT0**2.)
     2           + BETA0**4.*(-2.+6.*BT0**2.) )
      BFB32 = ( -BPRHO/(6.*BN0**5.) )
     1         * ( 1. + 4.*BT0**2. + 3.*BT0**4. )
      BFB34 = BSECT
     1        * ( BPRHO/(36.*BN0**7.) )
     2        * ( 7. + 55.*BT0**2. + 93.*BT0**4. + 45.*BT0**6. )
      BFB40 = ( BPRHO/(2.*BN0**4.) ) * BT0 * BETA0**2.
      BFB42 = ( -BPRHO/(6.*BN0**6.) ) * BT0
     1         * ( 2. + 5.*BT0 + 3.*BT0**4. )
C
      BFL01 = ( BPRHO/(BN0*COS(BB0*BCPI/180.)) )
      BFL03 = ( -BPRHO/(6.*BN0**3.*COS(BB0*BCPI/180.)) )
     1         * ( 1. + 2.*BT0**2. + BETA0**2. )
      BFL05 = ( BPRHO/(120.*BN0**5.*COS(BB0*BCPI/180.)) )
     1         * ( 5. + 28.*BT0**2. + 24.*BT0**4. )
      BFL11 = ( BPRHO/(BN0**2.*COS(BB0*BCPI/180.)) )
     1         * BT0
      BFL13 = ( -BPRHO/(6.*BN0**4.*COS(BB0*BCPI/180.)) )
     1         * BT0
     2         * ( 5. + 6.*BT0**2. + BETA0**2. )
      BFL15 = ( BPRHO/(120.*BN0**6.*COS(BB0*BCPI/180.)) )
     1         * BT0
     2         * ( 61. + 180.*BT0**2. + 120.*BT0**4. )
      BFL21 = ( BPRHO/(2.*BN0**3.*COS(BB0*BCPI/180.)) )
     1         * ( 1. + 2.*BT0**2. + BETA0**2. )
      BFL23 = ( -BPRHO/(12.*BN0**5.*COS(BB0*BCPI/180.)) )
     1         * ( 5. + 28.*BT0**2. + 24.*BT0**4. )
      BFL25 = BSECT
     1        * ( BPRHO/(240.*BN0**7.*COS(BB0*BCPI/180.)) )
     2        * ( 61. + 662.*BT0**2. + 1320.*BT0**4. + 720.*BT0**6. )
      BFL31 = ( BPRHO/(6.*BN0**4.*COS(BB0*BCPI/180.)) )
     1         * BT0
     2         * ( 5. + 6.*BT0**2. + BETA0**2. )
      BFL33 = ( -BPRHO/(36.*BN0**6.*COS(BB0*BCPI/180.)) )
     1         * BT0
     2         * ( 61. + 180.*BT0**2. + 120.*BT0**4. )
      BFL41 = ( BPRHO/(24.*BN0**5.*COS(BB0*BCPI/180.)) )
     1         * ( 5. + 28.*BT0**2. + 24.*BT0**4. )
C-----------------------------------------------------------------------
C
C---  calculation of differences of geographic coordinates
      BDELTB =        BFB10*BDELTX
     1        +       BFB20*BDELTX**2.
     2        +       BFB30*BDELTX**3.
     3        +       BFB40*BDELTX**4.
     4        +       BFB02           *BY**2.
     5        +       BFB12*BDELTX    *BY**2.
     6        +       BFB22*BDELTX**2.*BY**2.
     7        +       BFB32*BDELTX**3.*BY**2.
     8        +       BFB42*BDELTX**4.*BY**2.
     9        +       BFB04           *BY**4.
     *        +       BFB14*BDELTX    *BY**4.
     1        +       BFB24*BDELTX**2.*BY**4.
     2        + BSECT*BFB34*BDELTX**3.*BY**4.
     3        +       BFB06           *BY**6.
C
      BDELTL =        BFL01           *BY
     1        +       BFL11*BDELTX    *BY
     2        +       BFL21*BDELTX**2.*BY
     3        +       BFL31*BDELTX**3.*BY
     4        +       BFL41*BDELTX**4.*BY
     5        +       BFL03           *BY**3.
     6        +       BFL13*BDELTX    *BY**3.
     7        +       BFL23*BDELTX**2.*BY**3.
     8        +       BFL33*BDELTX**3.*BY**3.
     9        +       BFL05           *BY**5.
     9        +       BFL15*BDELTX    *BY**5.
     *        + BSECT*BFL25*BDELTX**2.*BY**5.
C
C---  calculation of geographic coordinates
      BLH = 3. * FLOAT( INT(AR/1000.) )
      AB  = BB0 + BDELTB
      AL  = BLH + BDELTL
C
C--- if AH coordinate is negative, latitude must be also negative
      IF (KNEGAT) AB = -AB
C
      RETURN
      END
