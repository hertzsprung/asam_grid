MODULE GridInput_Mod

  USE Function_Mod
  USE Floor_Mod
  USE Parametric_Mod
  USE IOControl_Mod
  USE BoundaryCond_Mod
  IMPLICIT NONE

  INTEGER, PRIVATE :: InputUnit=20

CONTAINS


SUBROUTINE GetNeighbor

!---------------------------------------------------------------
!---  Determine Neighbors w.r.t. Multiblock Structure
!---------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ib,ib1,in
  TYPE(Nachbar_T), ALLOCATABLE :: curr_nachb(:)

  ALLOCATE(curr_nachb(nb+5))
  DO ib=1,nb

    ix0=Floor(ib)%ix0
    ix1=Floor(ib)%ix1
    iy0=Floor(ib)%iy0
    iy1=Floor(ib)%iy1
    iz0=Floor(ib)%iz0
    iz1=Floor(ib)%iz1

    AnzahlNachbar=0

!--------------------------------------------------------
!   West Boundary

    IF (Floor(ib)%igx0==domain%ix0) THEN
      IF (BCVel%West/='Period') THEN
        AnzahlNachbar=AnzahlNachbar+1
        curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
        curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
        curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
        curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
        curr_nachb(AnzahlNachbar)%ib=ib
        curr_nachb(AnzahlNachbar)%nTYPE='ow'

      ELSE
        DO ib1=1,nb
          IF (domain%ix1==Floor(ib1)%igx1) THEN
            IF (Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                Floor(ib)%igy1>Floor(ib1)%igy0 .and. &
                Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
                Floor(ib)%igz1>Floor(ib1)%igz0) THEN
              AnzahlNachbar=AnzahlNachbar+1
              curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
              curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
              curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
              curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
              curr_nachb(AnzahlNachbar)%ib=ib1
              curr_nachb(AnzahlNachbar)%nTYPE='iw'
            END IF
          END IF
        END DO
      END IF
    ELSE
      DO ib1=1,nb
        IF (Floor(ib)%igx0==Floor(ib1)%igx1) THEN
          IF (Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
            Floor(ib)%igy1>Floor(ib1)%igy0 .and. &
            Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
            Floor(ib)%igz1>Floor(ib1)%igz0) THEN
            AnzahlNachbar=AnzahlNachbar+1
            curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
            curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
            curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
            curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
            curr_nachb(AnzahlNachbar)%ib=ib1
            curr_nachb(AnzahlNachbar)%nTYPE='iw'
          END IF
        END IF
      END DO
    END IF         ! Floor(ib)%igx0==domain%ix0

!--------------------------------------------------------
!   East Boundary

    IF (Floor(ib)%igx1==domain%ix1) THEN
      IF (BCVel%East/='Period') THEN
        AnzahlNachbar=AnzahlNachbar+1
        curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
        curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
        curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
        curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
        curr_nachb(AnzahlNachbar)%ib=ib
        curr_nachb(AnzahlNachbar)%nTYPE='oe'

      ELSE
        DO ib1=1,nb
          IF (domain%ix0==Floor(ib1)%igx0) THEN
            IF (Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                Floor(ib)%igy1>Floor(ib1)%igy0 .and. &
                Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
                Floor(ib)%igz1>Floor(ib1)%igz0) THEN
               AnzahlNachbar=AnzahlNachbar+1
               curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
               curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
               curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
               curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
               curr_nachb(AnzahlNachbar)%ib=ib1
               curr_nachb(AnzahlNachbar)%nTYPE='ie'
             END IF
           END IF
         END DO
       END IF
     ELSE
       DO ib1=1,nb
         IF (Floor(ib)%igx1==Floor(ib1)%igx0) THEN
           IF (Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
               Floor(ib)%igy1>Floor(ib1)%igy0 .and. &
               Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
               Floor(ib)%igz1>Floor(ib1)%igz0) THEN
              AnzahlNachbar=AnzahlNachbar+1
              curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
              curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
              curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
              curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
              curr_nachb(AnzahlNachbar)%ib=ib1
              curr_nachb(AnzahlNachbar)%nTYPE='ie'
           END IF
         END IF
       END DO
     END IF         ! Floor(ib)%igx1==domain%ix1

!--------------------------------------------------------
!    South Boundary

     IF (Floor(ib)%igy0==domain%iy0) THEN
       IF (BCVel%South/='Period') THEN
         AnzahlNachbar=AnzahlNachbar+1
         curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
         curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
         curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
         curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
         curr_nachb(AnzahlNachbar)%ib=ib
         curr_nachb(AnzahlNachbar)%nTYPE='os'
       ELSE
         DO ib1=1,nb
           IF (domain%iy1==Floor(ib1)%igy1) THEN
             IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                 Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                 Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
                 Floor(ib)%igz1>Floor(ib1)%igz0) THEN
               AnzahlNachbar=AnzahlNachbar+1
               curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
               curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
               curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
               curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
               curr_nachb(AnzahlNachbar)%ib=ib1
               curr_nachb(AnzahlNachbar)%nTYPE='is'
             END IF
           END IF
         END DO
       END IF
     ELSE
       DO ib1=1,nb
         IF (Floor(ib)%igy0==Floor(ib1)%igy1) THEN
           IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
               Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
               Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
               Floor(ib)%igz1>Floor(ib1)%igz0) THEN
             AnzahlNachbar=AnzahlNachbar+1
             curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
             curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
             curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
             curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
             curr_nachb(AnzahlNachbar)%ib=ib1
             curr_nachb(AnzahlNachbar)%nTYPE='is'
           END IF
         END IF
       END DO
     END IF         ! Floor(ib)%igy0==domain%iy0

!--------------------------------------------------------
!    North Boundary

     IF (Floor(ib)%igy1==domain%iy1) THEN
       IF (BCVel%North/='Period') THEN
         AnzahlNachbar=AnzahlNachbar+1
         curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
         curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
         curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
         curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
         curr_nachb(AnzahlNachbar)%ib=ib
         curr_nachb(AnzahlNachbar)%nTYPE='on'
       ELSE
         DO ib1=1,nb
           IF (domain%iy0==Floor(ib1)%igy0) THEN
             IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                 Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                 Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
                 Floor(ib)%igz1>Floor(ib1)%igz0) THEN
               AnzahlNachbar=AnzahlNachbar+1
               curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
               curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
               curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
               curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
               curr_nachb(AnzahlNachbar)%ib=ib1
               curr_nachb(AnzahlNachbar)%nTYPE='in'
             END IF
           END IF
         END DO
       END IF
     ELSE
       DO ib1=1,nb
         IF (Floor(ib)%igy1==Floor(ib1)%igy0) THEN
           IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
               Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
               Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
               Floor(ib)%igz1>Floor(ib1)%igz0) THEN
             AnzahlNachbar=AnzahlNachbar+1
             curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
             curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
             curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
             curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
             curr_nachb(AnzahlNachbar)%ib=ib1
             curr_nachb(AnzahlNachbar)%nTYPE='in'
           END IF
         END IF
       END DO
     END IF         ! Floor(ib)%igy1==domain%iy1

!--------------------------------------------------------
!      Bottom Boundary

       IF (Floor(ib)%igz0==domain%iz0) THEN
         IF (BCVel%Bottom/='Period') THEN
           AnzahlNachbar=AnzahlNachbar+1
           curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
           curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
           curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
           curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
           curr_nachb(AnzahlNachbar)%ib=ib
           curr_nachb(AnzahlNachbar)%nTYPE='ob'
         ELSE
           DO ib1=1,nb
             IF (domain%iz1==Floor(ib1)%igz1) THEN
               IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                   Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                   Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                   Floor(ib)%igy1>Floor(ib1)%igy0) THEN
                 AnzahlNachbar=AnzahlNachbar+1
                 curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
                 curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
                 curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
                 curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
                 curr_nachb(AnzahlNachbar)%ib=ib1
                 curr_nachb(AnzahlNachbar)%nTYPE='ib'
               END IF
             END IF
           END DO
         END IF
       ELSE
          DO ib1=1,nb
             IF (Floor(ib)%igz0==Floor(ib1)%igz1) THEN
                IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                    Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                    Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                    Floor(ib)%igy1>Floor(ib1)%igy0) THEN
                   AnzahlNachbar=AnzahlNachbar+1
                   curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
                   curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
                   curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
                   curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
                   curr_nachb(AnzahlNachbar)%ib=ib1
                   curr_nachb(AnzahlNachbar)%nTYPE='ib'
                END IF
             END IF
           END DO
       END IF         ! Floor(ib)%igz0==domain%iz0

!--------------------------------------------------------
!      Top Boundary

       IF (Floor(ib)%igz1==domain%iz1) THEN
         IF (BCVel%Top/='Period') THEN
           AnzahlNachbar=AnzahlNachbar+1
           curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
           curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
           curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
           curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
           curr_nachb(AnzahlNachbar)%ib=ib
           curr_nachb(AnzahlNachbar)%nTYPE='ot'
         ELSE
           DO ib1=1,nb
             IF (domain%iz0==Floor(ib1)%igz0) THEN
               IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                   Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                   Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                   Floor(ib)%igy1>Floor(ib1)%igy0) THEN
                 AnzahlNachbar=AnzahlNachbar+1
                 curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
                 curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
                 curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
                 curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
                 curr_nachb(AnzahlNachbar)%ib=ib1
                 curr_nachb(AnzahlNachbar)%nTYPE='it'
               END IF
             END IF
           END DO
         END IF
       ELSE
          DO ib1=1,nb
             IF (Floor(ib)%igz1==Floor(ib1)%igz0) THEN
                IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                    Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                    Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                    Floor(ib)%igy1>Floor(ib1)%igy0) THEN
                   AnzahlNachbar=AnzahlNachbar+1
                   curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
                   curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
                   curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
                   curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
                   curr_nachb(AnzahlNachbar)%ib=ib1
                   curr_nachb(AnzahlNachbar)%nTYPE='it'
                END IF
             END IF
           END DO
       END IF         ! Floor(ib)%igz1==domain%iz1

       Floor(ib)%AnzahlNachbar=AnzahlNachbar
       ALLOCATE(Floor(ib)%Nachbars(AnzahlNachbar))
       Floor(ib)%Nachbars(1:AnzahlNachbar)= curr_nachb(1:AnzahlNachbar)

  END DO         ! ib

  DEALLOCATE(curr_nachb)

  DO ib=1,nb
    CALL Set(Floor(ib))
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      Nachbar%IncrX=ABS(-RefineX+RefineNachbarX)+1
      Nachbar%IncrY=ABS(-RefineY+RefineNachbarY)+1
      Nachbar%IncrZ=ABS(-RefineZ+RefineNachbarZ)+1
      IF (nType=='iw'.OR.nType=='ie') THEN
        Nachbar%CopyCase=2*(RefineY-RefineNachbarY+1)+ &
                           (RefineZ-RefineNachbarZ+1)
      ELSE IF (nType=='is'.OR.nType=='in') THEN
        Nachbar%CopyCase=2*(RefineX-RefineNachbarX+1)+ &
                           (RefineZ-RefineNachbarZ+1)
      ELSE IF (nType=='ib'.OR.nType=='it') THEN
        Nachbar%CopyCase=2*(RefineX-RefineNachbarX+1)+ &
                           (RefineY-RefineNachbarY+1)
      END IF
    END DO
  END DO

END SUBROUTINE GetNeighbor

SUBROUTINE NeighborBounds
!
!---------------------------------------------------------------
!---  Determine Neighbor boundary data
!---------------------------------------------------------------
    INTEGER :: ib, ib1, in

    DO ib=1,nb
      CALL Set(Floor(ib))
      DO in=1,AnzahlNachbar
        ib1=Nachbars(in)%ib
        Nachbars(in)%ix0=MAX(igx0,Floor(ib1)%igx0)*2.e0**RefineX
        Nachbars(in)%ix1=MIN(igx1,Floor(ib1)%igx1)*2.e0**RefineX
        Nachbars(in)%iy0=MAX(igy0,Floor(ib1)%igy0)*2.e0**RefineY
        Nachbars(in)%iy1=MIN(igy1,Floor(ib1)%igy1)*2.e0**RefineY
        Nachbars(in)%iz0=MAX(igz0,Floor(ib1)%igz0)*2.e0**RefineZ
        Nachbars(in)%iz1=MIN(igz1,Floor(ib1)%igz1)*2.e0**RefineZ
        IF (Nachbars(in)%nType=='iw') THEN
          Nachbars(in)%ix0=Floor(ib1)%igx1*2.e0**RefineX-1
          Nachbars(in)%ix1=Floor(ib1)%igx1*2.e0**RefineX
          Nachbars(in)%iNx0=Floor(ib1)%igx1*2.e0**Floor(ib1)%RefineX
          Nachbars(in)%iNx1=Floor(ib1)%igx1*2.e0**Floor(ib1)%RefineX
          Nachbars(in)%ixO=ix0
          Nachbars(in)%ixI=ix0+1
          Nachbars(in)%iNxO=Nachbars(in)%iNx1-1
          Nachbars(in)%iNxI=Nachbars(in)%iNx1
        ELSE IF (Nachbars(in)%nType=='ie') THEN
          Nachbars(in)%ix0=Floor(ib1)%igx0*2.e0**RefineX
          Nachbars(in)%ix1=Floor(ib1)%igx0*2.e0**RefineX+1
          Nachbars(in)%iNx0=Floor(ib1)%igx0*2.e0**Floor(ib1)%RefineX
          Nachbars(in)%iNx1=Floor(ib1)%igx0*2.e0**Floor(ib1)%RefineX
          Nachbars(in)%ixO=ix1+1
          Nachbars(in)%ixI=ix1
          Nachbars(in)%iNxO=Nachbars(in)%iNx0
          Nachbars(in)%iNxI=Nachbars(in)%iNx0+1
        ELSE IF (Nachbars(in)%nType=='is') THEN
          Nachbars(in)%iy0=Floor(ib1)%igy1*2.e0**RefineY-1
          Nachbars(in)%iy1=Floor(ib1)%igy1*2.e0**RefineY
          Nachbars(in)%iNy0=Floor(ib1)%igy1*2.e0**Floor(ib1)%RefineY-1
          Nachbars(in)%iNy1=Floor(ib1)%igy1*2.e0**Floor(ib1)%RefineY
          Nachbars(in)%iyO=iy0
          Nachbars(in)%iyI=iy0+1
          Nachbars(in)%iNyO=Nachbars(in)%iNy1-1
          Nachbars(in)%iNyI=Nachbars(in)%iNy1
        ELSE IF (Nachbars(in)%nType=='in') THEN
          Nachbars(in)%iy0=Floor(ib1)%igy0*2.e0**RefineY
          Nachbars(in)%iy1=Floor(ib1)%igy0*2.e0**RefineY+1
          Nachbars(in)%iNy0=Floor(ib1)%igy0*2.e0**Floor(ib1)%RefineY
          Nachbars(in)%iNy1=Floor(ib1)%igy0*2.e0**Floor(ib1)%RefineY
          Nachbars(in)%iyO=iy1+1
          Nachbars(in)%iyI=iy1
          Nachbars(in)%iNyO=Nachbars(in)%iNy0
          Nachbars(in)%iNyI=Nachbars(in)%iNy0+1
        ELSE IF (Nachbars(in)%nType=='ib') THEN
          Nachbars(in)%iz0=Floor(ib1)%igz1*2.e0**RefineZ-1
          Nachbars(in)%iz1=Floor(ib1)%igz1*2.e0**RefineZ
          Nachbars(in)%iNz0=Floor(ib1)%igz1*2.e0**Floor(ib1)%RefineZ
          Nachbars(in)%iNz1=Floor(ib1)%igz1*2.e0**Floor(ib1)%RefineZ
          Nachbars(in)%izO=iz0
          Nachbars(in)%izI=iz0+1
          Nachbars(in)%iNzO=Nachbars(in)%iNz1-1
          Nachbars(in)%iNzI=Nachbars(in)%iNz1
        ELSE IF (Nachbars(in)%nType=='it') THEN
          Nachbars(in)%iz0=Floor(ib1)%igz0*2.e0**RefineZ
          Nachbars(in)%iz1=Floor(ib1)%igz0*2.e0**RefineZ+1
          Nachbars(in)%iNz0=Floor(ib1)%igz0*2.e0**Floor(ib1)%RefineZ
          Nachbars(in)%iNz1=Floor(ib1)%igz0*2.e0**Floor(ib1)%RefineZ
          Nachbars(in)%izO=iz1+1
          Nachbars(in)%izI=iz1
          Nachbars(in)%iNzO=Nachbars(in)%iNz0
          Nachbars(in)%iNzI=Nachbars(in)%iNz0+1
        END IF
      END DO
    END DO

END SUBROUTINE NeighborBounds

SUBROUTINE ReadGrid(FileName)

!---------------------------------------------------------------
!---  Determine Multiblock Structure
!---------------------------------------------------------------

  CHARACTER(*) :: FileName
  INTEGER :: ib,i,j,erg,id_sb
  CHARACTER(300) :: Line
  REAL(8) :: breite,laenge,zwx
  INTEGER(IntKind) :: ierr_ingmv, ierr_inctrl
  INTEGER :: xRes,yRes,zRes,BX,BY,BZ,ix,iy,iz
  INTEGER :: xIncr,yIncr,zIncr
  LOGICAL :: SelfMultiBlock=.TRUE.
  LOGICAL :: MultiBlock=.TRUE.

!---------------------------------------------------------------
!   Read data from file
  out_area='n'
  nr_sb=0
  nr_soildef=0
  nr_landdef=0
  nr_lsoil=7
  struct_bound=std_bound

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'GridFileOut')>0) THEN
      BACKSPACE(InputUnit)
      CALL input_GridFileOut(InputUnit, ierr_ingmv)
    ELSE IF (INDEX(Line,'#Gitter')>0) THEN
    !-------------------------------------
       READ(InputUnit,*) indata_type
       READ(InputUnit,*) Domain%x0,Domain%x1,Domain%nx
       READ(InputUnit,*) Domain%y0,Domain%y1,Domain%ny
       READ(InputUnit,*) Domain%z0,Domain%z1,Domain%nz

    ELSE IF (INDEX(Line,'#NrSoilLayers')>0) THEN
    !-------------------------------------------
       READ(InputUnit,*) nr_lsoil

    ELSE IF (INDEX(Line,'&OutGMVControl')>0) THEN
    !--------------------------------------------
       BACKSPACE(InputUnit)
       CALL input_GMVControl(InputUnit, ierr_ingmv)
       WRITE(*,*)
       CALL Compute_Domain_Coord
       WRITE(*,*)
       WRITE(*,*) 'Preferences Output:' 
       IF((TRIM(indata_type)==geographical).AND.(invalue_to_out=='y')) THEN
       WRITE(*,*) '  Chosen Grid Simmulation :             ',out_wahlgrid," ", geographical
       ELSE
       WRITE(*,*) '  Chosen Grid Simmulation :             ',out_wahlgrid
       END IF
       WRITE(*,*) '  Chosen Output-Type b=binary/a=ascii : ',out_type
       WRITE(*,*)
       WRITE(*,*) '  Selected RadOutput : ',RadOutput, '  (numerator for parametrization)'
       WRITE(*,*) '  Selected ScaleRad  : ',ScaleRad,  '  (denominator for parametrization)'
       WRITE(*,*) '  Selected ScaleSoil : ',ScaleSoil, '  (scale dzi_soil in % dzi)'
 
    ELSE IF (INDEX(Line,'&GridDistsCtrl')>0) THEN
    !--------------------------------------------
       BACKSPACE(InputUnit)
       CALL input_dists_ctrl(InputUnit, ierr_inctrl)
       WRITE(*,*)
       WRITE(*,*) 'Preferences to Analyze Grid:'
       WRITE(*,*) '  Coefficient Distance Point(x) : ',distx_coeff
       WRITE(*,*) '  Coefficient Distance Point(y) : ',disty_coeff
       WRITE(*,*) '  Coefficient Distance Point(z) : ',distz_coeff
       WRITE(*,*)

    ELSE IF (INDEX(Line,'#OutputDomain')>0) THEN
    !-------------------------------------------
       READ(InputUnit,*) Domain%def_out_domain
       IF(TRIM(Domain%def_out_domain)==def_out_index) THEN
         READ(InputUnit,*) Domain%view_ixa,Domain%view_ixe
         READ(InputUnit,*) Domain%view_iya,Domain%view_iye
         READ(InputUnit,*) Domain%view_iza,Domain%view_ize
       ELSE IF (TRIM(Domain%def_out_domain)==def_out_coord) THEN
         READ(InputUnit,*) Domain%view_xa,Domain%view_xe
         READ(InputUnit,*) Domain%view_ya,Domain%view_ye
         READ(InputUnit,*) Domain%view_za,Domain%view_ze
       ELSE 
         Write(*,*) " Input-Struct \'#OutputDomain\' is incorrect, uncompleted! \n"
         Write(*,*) " A piece of data  XYZ_Number .or. XYZ_Coord ?"
         STOP ':  SUBROUTINE ReadGrid , #OutputDomain \n'
       END IF   
       out_area='y'

    ELSE IF (INDEX(Line,'#Multiblock')>0.AND.MultiBlock) THEN
       SelfMultiBlock=.FALSE.
    !-----------------------------------------
       !CALL Set_Domain_Fac_LonLat ! siehe Info in Subroutine
       Domain%ix0=0
       Domain%iy0=0
       Domain%iz0=0
       READ(InputUnit,*) Domain%ix0,Domain%ix1
       READ(InputUnit,*) Domain%iy0,Domain%iy1
       READ(InputUnit,*) Domain%iz0,Domain%iz1
       CALL Check_Multi_Gitter_Domain
       Domain%nx = Domain%ix1 - Domain%ix0   ! unter #Gitter gelesen
       Domain%ny = Domain%iy1 - Domain%iy0   ! ""
       Domain%nz = Domain%iz1 - Domain%iz0   ! ""
       Domain%nc = Domain%nx * Domain%ny * Domain%nz

 !     !Block structure
       !............... 
       READ(InputUnit,*)
       READ(InputUnit,*) nb
       ALLOCATE(Floor(nb))
 
       DO ib=1,nb
          READ(InputUnit,*)
          READ(InputUnit,*) Floor(ib)%igx0,Floor(ib)%igx1
          READ(InputUnit,*) Floor(ib)%igy0,Floor(ib)%igy1
          READ(InputUnit,*) Floor(ib)%igz0,Floor(ib)%igz1
          READ(InputUnit,*) Floor(ib)%Refine,Floor(ib)%RefineX,Floor(ib)%RefineY, &
                            Floor(ib)%RefineZ,Floor(ib)%RefLevel
          IF (Floor(ib)%Refine>=0) THEN
             Floor(ib)%ix0=Floor(ib)%igx0*2**Floor(ib)%RefineX
             Floor(ib)%ix1=Floor(ib)%igx1*2**Floor(ib)%RefineX
             Floor(ib)%iy0=Floor(ib)%igy0*2**Floor(ib)%RefineY
             Floor(ib)%iy1=Floor(ib)%igy1*2**Floor(ib)%RefineY
             Floor(ib)%iz0=Floor(ib)%igz0*2**Floor(ib)%RefineZ
             Floor(ib)%iz1=Floor(ib)%igz1*2**Floor(ib)%RefineZ
             Floor(ib)%xShift=0
             Floor(ib)%yShift=0
             Floor(ib)%zShift=0
          ELSE
             Floor(ib)%ix0=Floor(ib)%igx0/2**(-Floor(ib)%RefineX)
             Floor(ib)%ix1=Floor(ib)%igx1/2**(-Floor(ib)%RefineX)
             Floor(ib)%iy0=Floor(ib)%igy0/2**(-Floor(ib)%RefineY)
             Floor(ib)%iy1=Floor(ib)%igy1/2**(-Floor(ib)%RefineY)
             Floor(ib)%iz0=Floor(ib)%igz0/2**(-Floor(ib)%RefineZ)
             Floor(ib)%iz1=Floor(ib)%igz1/2**(-Floor(ib)%RefineZ)
             Floor(ib)%xShift=mod(Floor(ib)%igx0,2**(-Floor(ib)%RefineX))
             Floor(ib)%yShift=mod(Floor(ib)%igy0,2**(-Floor(ib)%RefineX))
             Floor(ib)%zShift=mod(Floor(ib)%igz0,2**(-Floor(ib)%RefineX))
          END IF
          Floor(ib)%ib=ib
 
          Floor(ib)%nx = Floor(ib)%ix1 - Floor(ib)%ix0
          Floor(ib)%ny = Floor(ib)%iy1 - Floor(ib)%iy0
          Floor(ib)%nz = Floor(ib)%iz1 - Floor(ib)%iz0
          Floor(ib)%nc = Floor(ib)%nx * Floor(ib)%ny * Floor(ib)%nz
 
          Floor(ib)%TypeW='iw'
          IF (Floor(ib)%igx0==domain%ix0.AND.BCVel%West/='Period') THEN
            Floor(ib)%TypeW='ow'
          ELSE IF (Floor(ib)%igx0==domain%ix0.AND.BCVel%West=='Period') THEN
             Floor(ib)%TypeW='pw'
          END IF
          Floor(ib)%TypeE='ie'
          IF (Floor(ib)%igx1==domain%ix1.AND.BCVel%East/='Period') THEN
            Floor(ib)%TypeE='oe'
          ELSE IF (Floor(ib)%igx1==domain%ix1.AND.BCVel%East=='Period') THEN
             Floor(ib)%TypeE='pe'
          END IF
          Floor(ib)%TypeS='is'
          IF (Floor(ib)%igy0==domain%iy0.AND.BCVel%South/='Period') THEN
            Floor(ib)%TypeS='os'
          ELSE IF (Floor(ib)%igy0==domain%iy0.AND.BCVel%South=='Period') THEN
             Floor(ib)%TypeS='ps'
          END IF
          Floor(ib)%TypeN='in'
          IF (Floor(ib)%igy1==domain%iy1.AND.BCVel%North/='Period') THEN
            Floor(ib)%TypeN='on'
          ELSE IF (Floor(ib)%igy1==domain%iy1.AND.BCVel%North=='Period') THEN
             Floor(ib)%TypeN='pn'
          END IF
          Floor(ib)%TypeB='ib'
          IF (Floor(ib)%igz0==domain%iz0.AND.BCVel%Bottom/='Period') THEN
            Floor(ib)%TypeB='ob'
           ELSE IF (Floor(ib)%igz0==domain%iz0.AND.BCVel%Bottom=='Period') THEN
             Floor(ib)%TypeB='pb'
          END IF
          Floor(ib)%TypeT='it'
          IF (Floor(ib)%igz1==domain%iz1.AND.BCVel%Top/='Period') THEN
            Floor(ib)%TypeT='ot'
          ELSE IF (Floor(ib)%igz1==domain%iz1.AND.BCVel%Top=='Period') THEN
             Floor(ib)%TypeT='pt'
          END IF
 
          Floor(ib)%NrW_FacesXY=0
          Floor(ib)%NrW_FacesZX=0
          Floor(ib)%NrW_FacesYZ=0
          Floor(ib)%NrRW_FacesXY=0
          Floor(ib)%NrRW_FacesZX=0
          Floor(ib)%NrRW_FacesYZ=0
          Floor(ib)%NrMP_FacesXY=0
          Floor(ib)%NrMP_FacesYZ=0
          Floor(ib)%NrMP_FacesZX=0
          Floor(ib)%NrW_Cells=0
          Floor(ib)%NrRW_Cells=0
          Floor(ib)%NrMP_Cells=0
          Floor(ib)%NrB_Cells=0
          Floor(ib)%nr_soildef=0
       END DO  ! ib



    ELSE IF (INDEX(Line,'#SelfMultiblock')>0.AND.SelfMultiBlock) THEN
      MultiBlock=.FALSE.
      Domain%ix0=0
      Domain%iy0=0
      Domain%iz0=0
      READ(InputUnit,*) Domain%ix0,Domain%ix1,BX
      READ(InputUnit,*) Domain%iy0,Domain%iy1,BY
      READ(InputUnit,*) Domain%iz0,Domain%iz1,BZ
      CALL Check_Multi_Gitter_Domain
      Domain%nx = Domain%ix1 - Domain%ix0
      Domain%ny = Domain%iy1 - Domain%iy0
      Domain%nz = Domain%iz1 - Domain%iz0
      Domain%nc = Domain%nx * Domain%ny * Domain%nz
      nb=BX*BY*BZ
      xRes=MOD(Domain%nx,BX)
      xIncr=Domain%nx/BX
      yRes=MOD(Domain%ny,BY)
      yIncr=Domain%ny/BY
      zRes=MOD(Domain%nz,BZ)
      zIncr=Domain%nz/BZ

      ALLOCATE(Floor(nb))
      iz0=0
      ib=1
      DO iz=1,BZ
        IF (iz<=zRes) THEN
          iz1=iz0+zIncr+1
        ELSE  
          iz1=iz0+zIncr
        END IF  
        iy0=0
        DO iy=1,BY
          IF (iy<=yRes) THEN
            iy1=iy0+yIncr+1
          ELSE  
            iy1=iy0+yIncr
          END IF  
          ix0=0 
          DO ix=1,BX
            IF (ix<=xRes) THEN
              ix1=ix0+xIncr+1
            ELSE  
              ix1=ix0+xIncr
            END IF  
            Floor(ib)%igx0=ix0
            Floor(ib)%igy0=iy0
            Floor(ib)%igz0=iz0
            Floor(ib)%igx1=ix1
            Floor(ib)%igy1=iy1
            Floor(ib)%igz1=iz1
            WRITE(*,*) 'block structure x0,x1,y0,y1,z0,z1',Floor(ib)%igx0,Floor(ib)%igx1, &
                        Floor(ib)%igy0,Floor(ib)%igy1,Floor(ib)%igz0,Floor(ib)%igz1

            Floor(ib)%Refine=0
            Floor(ib)%RefineX=0
            Floor(ib)%RefineY=0
            Floor(ib)%RefineZ=0
            Floor(ib)%RefLevel=4

            IF (Floor(ib)%Refine>=0) THEN
              Floor(ib)%ix0=Floor(ib)%igx0*2**Floor(ib)%RefineX
              Floor(ib)%ix1=Floor(ib)%igx1*2**Floor(ib)%RefineX
              Floor(ib)%iy0=Floor(ib)%igy0*2**Floor(ib)%RefineY
              Floor(ib)%iy1=Floor(ib)%igy1*2**Floor(ib)%RefineY
              Floor(ib)%iz0=Floor(ib)%igz0*2**Floor(ib)%RefineZ
              Floor(ib)%iz1=Floor(ib)%igz1*2**Floor(ib)%RefineZ
              Floor(ib)%xShift=0
              Floor(ib)%yShift=0
              Floor(ib)%zShift=0
            ELSE
              Floor(ib)%ix0=Floor(ib)%igx0/2**(-Floor(ib)%RefineX)
              Floor(ib)%ix1=Floor(ib)%igx1/2**(-Floor(ib)%RefineX)
              Floor(ib)%iy0=Floor(ib)%igy0/2**(-Floor(ib)%RefineY)
              Floor(ib)%iy1=Floor(ib)%igy1/2**(-Floor(ib)%RefineY)
              Floor(ib)%iz0=Floor(ib)%igz0/2**(-Floor(ib)%RefineZ)
              Floor(ib)%iz1=Floor(ib)%igz1/2**(-Floor(ib)%RefineZ)
              Floor(ib)%xShift=mod(Floor(ib)%igx0,2**(-Floor(ib)%RefineX))
              Floor(ib)%yShift=mod(Floor(ib)%igy0,2**(-Floor(ib)%RefineY)) ! Hinneburg (vorher RefineX)
              Floor(ib)%zShift=mod(Floor(ib)%igz0,2**(-Floor(ib)%RefineZ)) ! Hinneburg (vorher RefineX)
            END IF
            Floor(ib)%ib=ib

            Floor(ib)%nx = Floor(ib)%ix1 - Floor(ib)%ix0
            Floor(ib)%ny = Floor(ib)%iy1 - Floor(ib)%iy0
            Floor(ib)%nz = Floor(ib)%iz1 - Floor(ib)%iz0
            Floor(ib)%nc = Floor(ib)%nx * Floor(ib)%ny * Floor(ib)%nz

          Floor(ib)%TypeW='iw'
          IF (Floor(ib)%igx0==domain%ix0.AND.BCVel%West/='Period') THEN
            Floor(ib)%TypeW='ow'
          ELSE IF (Floor(ib)%igx0==domain%ix0.AND.BCVel%West=='Period') THEN
             Floor(ib)%TypeW='pw'
          END IF
          Floor(ib)%TypeE='ie'
          IF (Floor(ib)%igx1==domain%ix1.AND.BCVel%East/='Period') THEN
            Floor(ib)%TypeE='oe'
          ELSE IF (Floor(ib)%igx1==domain%ix1.AND.BCVel%East=='Period') THEN
             Floor(ib)%TypeE='pe'
          END IF
          Floor(ib)%TypeS='is'
          IF (Floor(ib)%igy0==domain%iy0.AND.BCVel%South/='Period') THEN
            Floor(ib)%TypeS='os'
          ELSE IF (Floor(ib)%igy0==domain%iy0.AND.BCVel%South=='Period') THEN
             Floor(ib)%TypeS='ps'
          END IF
          Floor(ib)%TypeN='in'
          IF (Floor(ib)%igy1==domain%iy1.AND.BCVel%North/='Period') THEN
            Floor(ib)%TypeN='on'
          ELSE IF (Floor(ib)%igy1==domain%iy1.AND.BCVel%North=='Period') THEN
             Floor(ib)%TypeN='pn'
          END IF
          Floor(ib)%TypeB='ib'
          IF (Floor(ib)%igz0==domain%iz0.AND.BCVel%Bottom/='Period') THEN
            Floor(ib)%TypeB='ob'
           ELSE IF (Floor(ib)%igz0==domain%iz0.AND.BCVel%Bottom=='Period') THEN
             Floor(ib)%TypeB='pb'
          END IF
          Floor(ib)%TypeT='it'
          IF (Floor(ib)%igz1==domain%iz1.AND.BCVel%Top/='Period') THEN
            Floor(ib)%TypeT='ot'
          ELSE IF (Floor(ib)%igz1==domain%iz1.AND.BCVel%Top=='Period') THEN
             Floor(ib)%TypeT='pt'
          END IF

            Floor(ib)%NrW_FacesXY=0
            Floor(ib)%NrW_FacesZX=0
            Floor(ib)%NrW_FacesYZ=0
            Floor(ib)%NrRW_FacesXY=0
            Floor(ib)%NrRW_FacesZX=0
            Floor(ib)%NrRW_FacesYZ=0
            Floor(ib)%NrMP_FacesXY=0
            Floor(ib)%NrMP_FacesYZ=0
            Floor(ib)%NrMP_FacesZX=0
            Floor(ib)%NrW_Cells=0
            Floor(ib)%NrRW_Cells=0
            Floor(ib)%NrMP_Cells=0
            Floor(ib)%NrB_Cells=0
            Floor(ib)%nr_soildef=0
            ix0=ix1
            ib=ib+1
          END DO  
          iy0=iy1
        END DO  
        iz0=iz1
      END DO  
    END IF

    IF (INDEX(Line,'#MultiLayerSoil')>0) THEN
      WRITE(*,*) '#MultiLayerSoil'
    !----------------------------------------
       READ(InputUnit,*) nr_sb    ! number blocks with soil_def
       DO ib=1,nr_sb
         READ(InputUnit,*)  ! comment or. empty line 
         READ(InputUnit,*) id_sb
         READ(InputUnit,*) nr_soildef
         Floor(id_sb)%nr_soildef=nr_soildef
         ALLOCATE(Floor(id_sb)%soil_type(1:nr_soildef,1:nr_lsoil))
         ALLOCATE(Floor(id_sb)%s_ixa(1:nr_soildef)) 
         ALLOCATE(Floor(id_sb)%s_ixe(1:nr_soildef)) 
         ALLOCATE(Floor(id_sb)%s_iya(1:nr_soildef)) 
         ALLOCATE(Floor(id_sb)%s_iye(1:nr_soildef)) 
         ALLOCATE(Floor(id_sb)%s_iza(1:nr_soildef)) 
         ALLOCATE(Floor(id_sb)%s_ize(1:nr_soildef)) 
         DO i=1,nr_soildef
            READ(InputUnit,*) Floor(id_sb)%s_ixa(i),Floor(id_sb)%s_ixe(i)
            READ(InputUnit,*) Floor(id_sb)%s_iya(i),Floor(id_sb)%s_iye(i)
            READ(InputUnit,*) Floor(id_sb)%s_iza(i),Floor(id_sb)%s_ize(i)
            READ(InputUnit,*) (Floor(id_sb)%soil_type(i,j), j=1,nr_lsoil)
         END DO
       END DO  ! nr_sb
       if(nr_sb>0) struct_bound=struct_bound+1
    ELSE IF (INDEX(Line,'#LandClassDef')>0) THEN
      WRITE(*,*) '#LandClassDef'
    !-------------------------------------------
       IF(nr_wahlfkt==nr_fkt_wgs_oro .OR. &
          nr_wahlfkt==nr_fkt_wgs_surf) THEN
           READ(InputUnit,*) str_corine
           READ(InputUnit,*) file_clc_list
           READ(InputUnit,*) file_ncorine  
           file_ncorine=ADJUSTL(file_ncorine)
           CALL Read_Corine_List()
           CALL Read_Corine_Data
           ifcorine=.TRUE.
           nr_landdef=Nr_Land_Cl  !symbolisch
           struct_bound=struct_bound+2 
       ELSE 
          READ(InputUnit,*) nr_landdef
          if( nr_landdef>0) then
             ALLOCATE(LandDef(1:2,1:3,1:nr_landdef))
             struct_bound=struct_bound+2
          end if
          DO i=1,nr_landdef
            READ(InputUnit,*) ! comment or. empty line
            READ(InputUnit,*) LandDef(1,1,i), LandDef(2,1,i)  !lc_ixa,lc_ixe
            READ(InputUnit,*) LandDef(1,2,i), LandDef(2,2,i)  !lc_iya,lc_iye
            READ(InputUnit,*) LandDef(1,3,i)                  !LandClass
          END DO
       END IF
       EXIT
    END IF

  END DO
  ! Domain_view --> out_gmvG
  IF(out_area=='n') THEN  
     Domain%view_ixa=Domain%ix0
     Domain%view_ixe=Domain%ix1
     Domain%view_iya=Domain%iy0
     Domain%view_iye=Domain%iy1
     Domain%view_iza=Domain%iz0
     Domain%view_ize=Domain%iz1
     !.........................
     Domain%view_xa=Domain%x0
     Domain%view_xe=Domain%x1
     Domain%view_ya=Domain%y0
     Domain%view_ye=Domain%y1
     Domain%view_za=Domain%z0
     Domain%view_ze=Domain%z1
  END IF

1 CONTINUE
  CLOSE(UNIT=InputUnit)

  CALL Prot_Result_ReadGrid()

END SUBROUTINE ReadGrid

SUBROUTINE Check_Multi_Gitter_Domain
    IF(TRIM(indata_type) == wgs84) THEN
       IF(((Domain%ix1-Domain%ix0)/=Domain%nx+1).OR. &
          ((Domain%iy1-Domain%iy0)/=Domain%ny+1).OR. &
          ((Domain%iz1-Domain%iz0)/=Domain%nz)) THEN
         Write(*,*) "Difference detected between,"
         Write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>> \n", &
                    "       Definitnon in Input-Struct \'#Gitter\'! \n \n", &
                    "       Definition in Input-Struct \'#Multiblock\' \n", &
                    "       Domain%ix1-Domain%ix0/=Domain%nx+1 \n", &
                    "       Domain%iy1-Domain%iy0/=Domain%ny+1 \n", &
                    "       Domain%iz1-Domain%iz0/=Domain%nz \n"
         Write(*,*) Domain%ix1-Domain%ix0,"/=",Domain%nx+1,", ",&
                    Domain%iy1-Domain%iy0,"/=",Domain%ny+1,", ",&
                    Domain%iz1-Domain%iz0,"/=",Domain%nz
         STOP ':  SUBROUTINE ReadGrid \n'
       END IF
      
    ELSE
       IF(((Domain%ix1-Domain%ix0)/=Domain%nx).OR. &
          ((Domain%iy1-Domain%iy0)/=Domain%ny).OR. &
          ((Domain%iz1-Domain%iz0)/=Domain%nz)) THEN
         Write(*,*) "Difference detected between,"
         Write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>> \n", &
                    "       Definitnon in Input-Struct \'#Gitter\'! \n \n", &
                    "       Definition in Input-Struct \'#Multiblock\' \n", &
                    "       Domain%ix1-Domain%ix0/=Domain%nx \n", &
                    "       Domain%iy1-Domain%iy0/=Domain%ny \n", &
                    "       Domain%iz1-Domain%iz0/=Domain%nz \n"
         Write(*,*) Domain%ix1-Domain%ix0,"/=",Domain%nx,", ",&
                    Domain%iy1-Domain%iy0,"/=",Domain%ny,", ",&
                    Domain%iz1-Domain%iz0,"/=",Domain%nz
         STOP ':  SUBROUTINE ReadGrid \n'
       END IF
    END IF
END SUBROUTINE Check_Multi_Gitter_Domain
SUBROUTINE Compute_Domain_Coord
  REAL(8) :: breite,laenge

       Write(*,*) "Coordinate Input:"
       IF(TRIM(indata_type) == gausskrueger ) THEN
         Write(*,*) "Gauss-Krueger:"
         Write(*,*) '  Domain%x0=', Domain%x0, '  Domain%y0=', Domain%y0
         Write(*,*) '  Domain%x1=', Domain%x1, '  Domain%y1=', Domain%y1
         IF(conv_ingrid=='spherical'.OR.out_wahlgrid=='G') THEN
           Write(*,*) 'Geographisch:'
           Call OGKGEO(Domain%y0/1000, Domain%x0/1000, breite, laenge)
           Domain%y0=(breite*4.0d0*ATAN(1.0d0))/180
           Domain%x0=(laenge*4.0d0*ATAN(1.0d0))/180
           Write(*,*) '  Domain%x0=',laenge, "  Domain%y0=",breite
           Call OGKGEO(Domain%y1/1000, Domain%x1/1000, breite, laenge)
           Domain%y1=(breite*4.0d0*ATAN(1.0d0))/180
           Domain%x1=(laenge*4.0d0*ATAN(1.0d0))/180
           Write(*,*) '  Domain%x1=',laenge, "  Domain%y1=",breite
           Write(*,*) 'Bogenmass:'
           Write(*,*) '  Domain%x0=', Domain%x0, '  Domain%y0=', Domain%y0
           Write(*,*) '  Domain%x1=', Domain%x1, '  Domain%y1=', Domain%y1
           IF (TRIM(Domain%def_out_domain)==def_out_coord) THEN
              Call OGKGEO(Domain%view_ya/1000, Domain%view_xa/1000, breite, laenge)
              Domain%view_ya=(breite*4.0d0*ATAN(1.0d0))/180
              Domain%view_xa=(laenge*4.0d0*ATAN(1.0d0))/180
              Call OGKGEO(Domain%view_ye/1000, Domain%view_xe/1000, breite, laenge)
              Domain%view_ye=(breite*4.0d0*ATAN(1.0d0))/180
              Domain%view_xe=(laenge*4.0d0*ATAN(1.0d0))/180
           END IF
         ELSE
           !Domain gausskrueger not to convert for computation
           !Nullpunkt verschieben, ermitteln
           gaussx_min=Domain%x0
           gaussy_min=Domain%y0
           !Domain%x1=Domain%x1-Domain%x0
           !Domain%x0=0.0
           !Domain%y1=Domain%y1-Domain%y0
           !Domain%y0=0.0
           Write(*,*) 'Domain gausskrueger not to convert for computation!'
           Write(*,*) 'zero point scale put off ---> lower edge domain'
         END IF
       ELSE IF (TRIM(indata_type)==geographical.OR.TRIM(indata_type)==wgs84) THEN
         Write(*,*) 'Geographisch:  '
         Write(*,*) '  Domain%x0=', Domain%x0, '  Domain%y0=', Domain%y0
         Write(*,*) '  Domain%x1=', Domain%x1, '  Domain%y1=', Domain%y1
         Domain%x0=(Domain%x0*4.0d0*ATAN(1.0d0))/180.0d0
         Domain%x1=(Domain%x1*4.0d0*ATAN(1.0d0))/180.0d0
         Domain%y0=(Domain%y0*4.0d0*ATAN(1.0d0))/180.0d0
         Domain%y1=(Domain%y1*4.0d0*ATAN(1.0d0))/180.0d0
         Write(*,*) 'Bogenmass:'
         Write(*,*) '  Domain%x0=', Domain%x0, '  Domain%y0=', Domain%y0
         Write(*,*) '  Domain%x1=', Domain%x1, '  Domain%y1=', Domain%y1
         Write(*,*) '  Domain%dx=', ABS((Domain%x1-Domain%x0)/Domain%nx)
         Write(*,*) '  Domain%dy=', ABS((Domain%y1-Domain%y0)/Domain%ny)
         Write(*,*) 'z-Achse  :'
         IF (TRIM(Domain%def_out_domain)==def_out_coord) THEN
            Domain%view_xa=x_GeoToRad(Domain%view_xa)
            Domain%view_xe=x_GeoToRad(Domain%view_xe)
            Domain%view_ya=y_GeoToRad(Domain%view_ya)
            Domain%view_ye=y_GeoToRad(Domain%view_ye)
         END IF
       ELSE IF (TRIM(indata_type)=='cyl') THEN
         Write(*,*) 'Zylindrisch:  '
         Domain%x0=(Domain%x0*4.0d0*ATAN(1.0d0))/180.0d0
         Domain%x1=(Domain%x1*4.0d0*ATAN(1.0d0))/180.0d0
         Write(*,*) '  Domain%x0=', Domain%x0, '  Domain%y0=', Domain%y0
         Write(*,*) '  Domain%x1=', Domain%x1, '  Domain%y1=', Domain%y1
         IF (TRIM(Domain%def_out_domain)==def_out_coord) THEN
           Domain%view_xa=(Domain%view_xa*4.0d0*ATAN(1.0d0))/180.0d0
           Domain%view_xe=(Domain%view_xe*4.0d0*ATAN(1.0d0))/180.0d0
         END IF
       ELSE
         Write(*,*) '  Domain%x0=', Domain%x0, '  Domain%y0=', Domain%y0
         Write(*,*) '  Domain%x1=', Domain%x1, '  Domain%y1=', Domain%y1
       END IF
       Write(*,*) '  Domain%z0=', Domain%z0
       Write(*,*) '  Domain%z1=', Domain%z1
       Write(*,*) '  Domain%dz=', (Domain%z1-Domain%z0)/Domain%nz

END SUBROUTINE Compute_Domain_Coord

SUBROUTINE  Prot_Result_ReadGrid

  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) trenn2 
  WRITE(OutUnitProt,*) leerzei3,'>>>>       Result  Read-Grid        <<<<'
  WRITE(OutUnitProt,*) trenn2 
  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) leerzei2,'Chosen Grid Simmulation :             ',out_wahlgrid
  WRITE(OutUnitProt,*) leerzei2,'Chosen Output-Type b=binary/a=ascii : ',out_type
  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) leerzei2,'Selected RadOutput : ',RadOutput, &
                       & '  (value to numerator for parametrization)'
  WRITE(OutUnitProt,*) leerzei2,'Selected ScaleRad  : ',ScaleRad, &
                       & '  (value to denominator for parametrization)'
  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) leerzei2,'Definitions of Analysis : '
  WRITE(OutUnitProt,*) leerzei2,'    Coefficient Distance Point(x) : ',distx_coeff
  WRITE(OutUnitProt,*) leerzei2,'    Coefficient Distance Point(y) : ',disty_coeff
  WRITE(OutUnitProt,*) leerzei2,'    Coefficient Distance Point(z) : ',distz_coeff
  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) leerzei2,'    Fine scaling distance to in_out-def: ', dist_fscv
  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) leerzei2,'    Distance Points coefficient of model border : '
  WRITE(OutUnitProt,*) leerzei2,'                       Point(x) :  ',dxViewLoc
  WRITE(OutUnitProt,*) leerzei2,'                       Point(y) :  ',dyViewLoc
  WRITE(OutUnitProt,*) leerzei2,'                       Point(z) :  ',dzViewLoc
  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) leerzei2,'    Number of cuts for volume-analysis : ',IncrVol
  WRITE(OutUnitProt,*)
  Write(OutUnitProt,*) leerzei2,"Coordinate Input:"
  Write(OutUnitProt,*) leerzei2,'  Domain%x0=', Domain%x0, '  Domain%y0=', Domain%y0
  Write(OutUnitProt,*) leerzei2,'  Domain%x1=', Domain%x1, '  Domain%y1=', Domain%y1
  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) leerzei2,"Function : ",TRIM(name_fkt)
  IF(nr_wahlfkt==9) THEN
     WRITE(OutUnitProt,*) leerzei13, TRIM(NameFun)
  END IF
  WRITE(OutUnitProt,*) leerzei13, TRIM(file_namefkt)
  WRITE(OutUnitProt,*)

  ! Ausgabe Grid-Input-Variablen: *.grid
  WRITE(OutUnitProt,*) &
    " ------------------------------------------------------------------\n", &
    " | Tabelle-Input-Variablen aus grid-File ->\n", &
    " ------------------------------------------------------------------\n", &
    " | &GridDistsCtrl:\n", & 
    " ------------------------------------------------------------------\n", &
    " |  dist_fscv      =",dist_fscv, &
        & " ! for fine scaling dist to in_out-def, CheckVertex\n", &
    " |  distx_coeff    =",distx_coeff, &
        & " ! Distance Point(x) coefficient for analyze\n", &
    " |  disty_coeff    =",disty_coeff, &
        & " ! Distance Point(y) coefficient for analyze\n", &
    " |  distz_coeff    =",distz_coeff, &
        & " ! Distance Point(z) coefficient for analyze\n", &
    " |  dxViewLoc      =",dxViewLoc, &
        & " ! Distance Point(x) coefficient of model border\n", &
    " |  dyViewLoc      =",dyViewLoc, &
        & " ! Distance Point(y) coefficient of model border\n", &
    " |  dzViewLoc      =",dzViewLoc, &
        & " ! Distance Point(z) coefficient of model border\n", &
    " |  IncrVol        =",IncrVol, &
        & " ! counter for cuts of volume-analysis (1 für Einheitswürfel)\n", &
    " |  dist_scMaxCell =",dist_scMaxCell, &
        & " ! adjustment value to the filter(screen) of cells with roughly max. Vol"
  Write(OutUnitProt,*) &
    " ------------------------------------------------------------------\n", &
    " | &OutGMVControl:\n", &
    " ------------------------------------------------------------------\n", &
    " |  out_wahlgrid   = ",out_wahlgrid, &
        & " ! Output-Grid: G-Global,C-Cartesian\n", &
    " |  out_type       = ",out_type, &
        & " ! Output-format: b->binary,a->ascii\n", &
    " |  invalue_to_out = ",invalue_to_out, &
        & " ! convert 'rad' to input value 'geo' for Output-GMV, only geo||wgs84 input\n", &
    " |  RadOutput      =",RadOutput, &
        & " ! RadOutput, value to numerator for parametrization\n", &
    " |  ScaleRad       =",ScaleRad, &
        & " ! ScaleRad, value to denominator for parametrization\n", &
    " ------------------------------------------------------------------"
  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) leerzei3,'>>>>     End  Result  Read-Grid        <<<<'
  WRITE(OutUnitProt,*) trenn2 
  WRITE(OutUnitProt,*)

END SUBROUTINE  Prot_Result_ReadGrid

SUBROUTINE ReadWeights(InputFile)

  CHARACTER*80 :: InputFile
  CHARACTER(300) :: Line
  INTEGER :: ib,i,j,k
  INTEGER :: ix,iy,iz
  INTEGER :: ixref, iyref, izref
  INTEGER :: levX,levY,levZ,shX,shY,shZ
  INTEGER :: levNX,levNY,levNZ,shNX,shNY,shNZ
  INTEGER :: InputUnit=10
  INTEGER :: nxx,nyy,nzz,nlines
  REAL(8) :: xw0,xw1,yw0,yw1,zw0,zw1
 
  OPEN(UNIT=InputUnit,FILE=TRIM(InputFile),STATUS='old')

  ALLOCATE(domain%dx(domain%ix0+1:domain%ix1))
  ALLOCATE(domain%dy(domain%iy0+1:domain%iy1))
  ALLOCATE(domain%dz(domain%iz0+1:domain%iz1))
  ALLOCATE(domain%xP(domain%ix0:domain%ix1))
  ALLOCATE(domain%yP(domain%iy0:domain%iy1))
  ALLOCATE(domain%zP(domain%iz0:domain%iz1))

  DO 
    READ(InputUnit,*,END=1) Line
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
1 CONTINUE
  REWIND(InputUnit)
  DO 
    READ(InputUnit,*,END=2) Line
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
2 CONTINUE
  REWIND(InputUnit)
  DO 
    READ(InputUnit,*,END=3) Line
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
!        WRITE(*,*) 'DOMAIN Z',zw0,zw1,domain%dz(iz),domain%zP(iz)
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
3 CONTINUE
  CLOSE(UNIT=InputUnit)
  !..........................................................................
  ! Verschoben-Kopie  InitAnalyzeArea
  domain%x0View=domain%xP(domain%ix0)
  domain%x1View=domain%xP(domain%ix1)
  domain%y0View=domain%yP(domain%iy0)
  domain%y1View=domain%yP(domain%iy1)
  domain%z0View=domain%zP(domain%iz0)
  domain%z1View=domain%zP(domain%iz1)
  !..........................................................................

  DO ib=1,nb
    CALL Set(Floor(ib))
    levX=-RefineX
    shX=MAX(1,2**levX)-1
    levY=-RefineY
    shY=MAX(1,2**levY)-1
    levZ=-RefineZ
    shZ=MAX(1,2**levZ)-1
    !......................................
    ix=ix0+1
    ixref = ix * 2.e0**levX
    xP(ix0)=domain%xP(ixref-shX-1)
    DO ix=ix0+1,ix1
      ixref = ix * 2.e0**levX
      dx(ix)=SUM(domain%dx(ixref-shX:ixref))
      dx(ix)=ABS(dx(ix))
      xP(ix)=xP(ix-1)+dx(ix)
    END DO
    !......................................
    iy=iy0+1
    iyref = iy * 2.e0**levY
    yP(iy0)=domain%yP(iyref-shY-1)
    DO iy=iy0+1,iy1
      iyref = iy * 2.e0**levY
      dy(iy)=SUM(domain%dy(iyref-shY:iyref))
      dy(iy)=ABS(dy(iy))
      yP(iy)=yP(iy-1)+dy(iy)
    END DO
    !......................................
    iz=iz0+1
    izref = iz * 2.e0**levZ
    zP(iz0)=domain%zP(izref-shZ-1)
    DO iz=iz0+1,iz1
      izref = iz * 2.e0**levZ
      dz(iz)=SUM(domain%dz(izref-shZ:izref))
      dz(iz)=ABS(dz(iz))
      zP(iz)=zP(iz-1)+dz(iz)
    END DO
  END DO

END SUBROUTINE ReadWeights

END MODULE GridInput_Mod
