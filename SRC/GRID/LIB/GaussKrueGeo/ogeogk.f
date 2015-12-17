      SUBROUTINE ogeogk(AB,AL,AH,AR,IREFS)
C
C---------------------------------------------------------------------
CC    *OGEOGK* : Conversion of geographic coordinates to
CC               Gauss-Krueger-coordinates
C
C     Input:  geographic coordinates AB and AL; IREFS
C             AB    = latitude  [degree]
C             AL    = longitude [degree]
C             IREFS = number of selected Gauss-Krueger-sector
C     Output: Gauss-Krueger-coordinates AH and AR
C             AH = Hochwert   [distance from equator in km]
C             AR = Rechtswert [sector-number and distance from
C                              corresponding longitude in km]
C             Example: AH = 5926. and AR = 3553. means
C                      5926.km distance from equator, corresponding
C                      longitude 9 degrees east ( = first number*3 )
C                      and 53.km east ( = 553.-500. ) from
C                      corresponding longitude
C
C     ATTENTION: Normally IREFS should have the value "-999". That means
C     ========== OGEOGK will calculate the correct Gauss-Krueger-sector
C                from the given longitude AL.
C                For special purpose the user may wish Gauss-Krueger-
C                coordinates with respect to other median longitudes 
C                (=other sectors). Then IREFS has to be the number of
C                this preferred sector.
C                Example: AL = 7. [deg] --> sector number 2 with
C                         median longitude BLH = 6. [deg].
C                         But Gauss-Krueger-coordinates are wanted with
C                         respect to sector number 3 (median longitude
C                         BLH = 9. [deg]). Then the user has to choose
C                         IREFS = 3 !
C                It should be noticed that the error of the calculated
C                Gauss-Krueger-coordinates increases with increasing
C                distance from the median longitude ! So coordinates
C                are less exact within others than the default sectors !!
C---------------------------------------------------------------------
C
C     CHANGES :
C     =========
C     NAME         DATE         COMMENT
C     ----         ----         -------
C     Nikolaus     06/12/91     new subroutine
C     Klaus        23/11/92     others than default sectors possible
C     Klaus        30/01/93     number of possible sectors expanded
C     Jerzy        09/06/94     IMPLICIT NONE. if the latitude is negative
C                               (southern from equator), the AH coordinate
C                               will be also negative. new variable KNEGAT.
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
C     IREFS     I        number of selected Gauss-Krueger-sector
C                        = -999:  no sector selected, correct sector will 
C                                 be calculated from given longitudes
C                        = value: number of selected Gauss-Krueger-sector,
C                                 median longitude of selected sector BLH
C                                 is calculated by BLH = 3.*"value"
C     BB0       R        rounded value of latitude
C     BCPI      R        Pi
C     BDELTB    R        difference of latitude to rounded value
C     BDELTL    R        difference of longitude to rounded value
C     BDELTX    R        difference of "Hochwert" to rounded value
C     BETA0     R
C     BFXnm     R        Koeffizienten der Reihenentwicklung fuer x
C     BFYnm     R        Koeffizienten der Reihenentwicklung fuer y
C     BH0       R        rounded value of "Hochwert"
C     BLH       R        rounded value of longitude
C     BN0       R        Querkruemmungshalbmesser
C     BPCQU     R        Polkruemmungshalbmesser
C     BPE       R        Quadrat der 2. numerischen Exzentritaet
C     BPK       R        Zuschlag in Ostrichtung fuer G.-K.-Koord.
C     BPRHO     R        Umrechnungsfaktor
C     BSECT     R        index for section of longitude
C     BT0       R
C     BY        R        difference of "Rechtswert" to rounded value
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
      REAL BPCQU,BPE,BPK,BPRHO
      PARAMETER ( BPCQU = 6398786.849
     1          , BPE   = 0.006719219
     2          , BPK   = 500.
     3          , BPRHO = 57.2957795131 )
C
      INTEGER IREFS
      REAL    AB,AH,AL,AR
C
      REAL BB0,BCPI,BDELTB,BDELTL,BDELTX,BETA0,BH0,BLH,BN0
     1    ,BSECT,BT0,BY
C
      REAL BFX02,BFX04,BFX10,BFX12,BFX14,BFX20,BFX22,BFX24
     1    ,BFX30,BFX32,BFX40
     2    ,BFY01,BFY03,BFY05,BFY11,BFY13,BFY15,BFY21,BFY23
     3    ,BFY31,BFY33,BFY41
C--- whether latitude is negative (southern)
      LOGICAL KNEGAT
C-----------------------------------------------------------------------
C
C---  defining PI
      BCPI   = 4.*ATAN(1.)
C
C--- check, if latitude is negative
      KNEGAT = AB .LT. 0.0
      IF (KNEGAT) AB = ABS (AB)
C
C---  rounding of latitude and longitude; median longitude of sector
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
C
C---  sector of latitude between 47 and 55 degree?
      IF ( AB.GT.47.  .AND.  AB.LE.55. ) THEN
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
      BFX02 = ( 1./(2.*BPRHO**2.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**2. * BT0
      BFX04 = ( 1./(24.*BPRHO**4.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**4. * BT0
     2         * ( 5. - BT0**2. + 9.*BETA0**2. )
      BFX10 = ( 1./BPRHO )
     1         * BN0
     2         * ( 1. - BETA0**2. + BETA0**4. - BETA0**6. )
      BFX12 = ( 1./(2.*BPRHO**3.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**2.
     2         * ( 1. - BT0**2. + BETA0**2.*BT0**2.
     3                - BETA0**4.*BT0**2.           )
      BFX14 = ( 1./(24.*BPRHO**5.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**4.
     2         * ( 5. - 18.*BT0**2. + BT0**4.
     3            + BSECT*BETA0**2.*(9.-40.*BT0**2.-BT0**4.) )
      BFX20 = ( 3./(2.*BPRHO**2.) )
     1         * BN0 * BT0
     2         * ( BETA0**2. - 2.*BETA0**4. )
      BFX22 = ( 1./(4.*BPRHO**4. ) )
     1         * BN0 * (COS(BB0*BCPI/180.))**2. * BT0
     2         * ( -4. + 3.*BETA0**2.*(1.-BT0**2.) )
      BFX24 = ( 1./(6.*BPRHO**6.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**4. * BT0
     2         * ( -7. + 5.*BT0**2. )
      BFX30 = ( 1./(2.*BPRHO**3.) )
     1         * BN0
     2         * ( BETA0**2.*(1.-BT0**2.) + BETA0**4.*(-2.+7.*BT0**2.) )
      BFX32 = ( 1./(12.*BPRHO**5.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**2.
     2         * ( -4. + 4.*BT0**2. )
      BFX40 = -( 1./(2.*BPRHO**4.) )
     1         * BN0 * BETA0**2. * BT0
C
      BFY01 = ( 1./BPRHO )
     1         * BN0 * COS(BB0*BCPI/180.)
      BFY03 = ( 1./(6.*BPRHO**3.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**3.
     2         * ( 1. - BT0**2. + BETA0**2. )
      BFY05 = ( 1./(120.*BPRHO**5.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**5.
     2         * ( 5. - 18.*BT0**2. + BT0**4. )
      BFY11 = ( 1./BPRHO**2. )
     1         * BN0 * COS(BB0*BCPI/180.) * BT0
     2         * ( -1. + BETA0**2. - BETA0**4. + BETA0**6. )
      BFY13 = ( 1./(6.*BPRHO**4.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**3. * BT0
     2         * ( -5. + BT0**2. - BETA0**2.*(4.+BT0**2.) )
      BFY15 = ( 1./(120.*BPRHO**6.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**5. * BT0
     2         * ( -61. + 58.*BT0**2. - BT0**4. )
      BFY21 = ( 1./(2.*BPRHO**3.) )
     1         * BN0 * COS(BB0*BCPI/180.)
     2         * ( -1. + BETA0**2.*( 1.-3.*BT0**2.)
     3                 + BETA0**4.*(-1.+6.*BT0**2.) )
      BFY23 = ( 1./(12.*BPRHO**5.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**3.
     2         * ( -5. + 13.*BT0**2. )
      BFY31 = ( 1./(6.*BPRHO**4.) )
     1         * BN0 * COS(BB0*BCPI/180.) * BT0
     2         * ( 1. + BETA0**2.*(-10.+3.*BT0**2.) )
      BFY33 = ( 1./(36.*BPRHO**6.) )
     1         * BN0 * (COS(BB0*BCPI/180.))**3. * BT0
     2         * ( 41. - 13.*BT0**2. )
      BFY41 = ( 1./(24.*BPRHO**5.) )
     1         * BN0 * COS(BB0*BCPI/180.)
C-----------------------------------------------------------------------
C
C---  calculation of meridian arc at main longitude
      BH0 =  111120.61962 * BB0
     1     -  15988.63853 * SIN(2.*BB0*BCPI/180.)
     2     +     16.72995 * SIN(4.*BB0*BCPI/180.)
     3     -      0.02178 * SIN(6.*BB0*BCPI/180.)
     4     +      0.00003 * SIN(8.*BB0*BCPI/180.)
C
C---  calculation of differences of coordinates to main meridian
C---  and rounded latitude
      BDELTX =        BFX10*BDELTB
     1        +       BFX20*BDELTB**2.
     2        +       BFX30*BDELTB**3.
     3        +       BFX40*BDELTB**4.
     4        +       BFX02           *BDELTL**2.
     5        +       BFX12*BDELTB    *BDELTL**2.
     6        +       BFX22*BDELTB**2.*BDELTL**2.
     7        +       BFX32*BDELTB**3.*BDELTL**2.
     8        +       BFX04           *BDELTL**4.
     9        +       BFX14*BDELTB    *BDELTL**4.
     *        + BSECT*BFX24*BDELTB**2.*BDELTL**4.
C
      BY =        BFY01           *BDELTL
     1    +       BFY11*BDELTB    *BDELTL
     2    +       BFY21*BDELTB**2.*BDELTL
     3    +       BFY31*BDELTB**3.*BDELTL
     4    +       BFY41*BDELTB**4.*BDELTL
     5    +       BFY03           *BDELTL**3.
     6    +       BFY13*BDELTB    *BDELTL**3.
     7    +       BFY23*BDELTB**2.*BDELTL**3.
     8    +       BFY33*BDELTB**3.*BDELTL**3.
     9    +       BFY05           *BDELTL**5.
     *    + BSECT*BFY15*BDELTB    *BDELTL**5.
C
C---  conversion of coordinates from [m] to [km]
      BH0    = BH0    / 1000.
      BDELTX = BDELTX / 1000.
      BY     = BY     / 1000.
C
C---  calculation of Gauss-Krueger-coordinates
      AR = FLOAT(IREFS)*1.E3 + BPK + BY
      AH = BH0 + BDELTX
C
C--- if latitude is negative, correct AH coordinate
      IF (KNEGAT) AH = -AH
C
      RETURN
      END
