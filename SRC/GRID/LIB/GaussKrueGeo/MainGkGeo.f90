PROGRAM MainProg
! Schmuecke:
! y         x        z
! 5619000.  4408000. 870.0
! 5609000.  4418000.

  IMPLICIT NONE
  REAL(8), PARAMETER :: RadOutput=1.0d4

  REAL(8) :: x,x1,x2,y,y1,y2,z, b, l,bb,lb
  INTEGER :: IREFS
  CHARACTER :: wahl
  y1=5619.0
  x1=4408.0
  y2=5609.0
  x2=4418.0

  x=5926.0
  y=3553.0
  z=870.0 

  WRITE (*,*) 
  WRITE (*,*) 'Program Start...'
  !::::::::::::::::::::::::::::::
  Write(*,*) "RadOutput:",RadOutput
   
  Write(*,*) "Gauss-Krueger: Bsp"
  Write(*,*) 'x=',x, "  y=", y
  Call OGKGEO(y, x, b, l)
  Write(*,*) 'Geographisch'
  Write(*,*) 'x=',l, "  y=",b 
  Write(*,*) 'Bogenmaß:'
  lb=(l*4.0d0*ATAN(1.0d0))/180
  bb=(b*4.0d0*ATAN(1.0d0))/180
  Write(*,*) 'x=',lb, "  y=",bb
!  Write(*,*) "Rueckrechnung Geographisch --> Gauss-Krueger: Bsp"
!  IREFS=?
!  Call OGEOGK(b, l, y, x, IREFS)
!  Write(*,*) 'x=',x, "  y=", y
  Write(*,*) ''
 

  Write(*,*) "Gauss-Krueger: Schmuecke-Ausschnitt"
  Write(*,*) 'x=',x1, "  y=", y1
  Call OGKGEO(y1, x1, b, l)
  Write(*,*) 'Geographisch'
  Write(*,*) 'x=',l, "  y=",b 
  Write(*,*) 'Bogenmaß:'
  lb=(l*4.0d0*ATAN(1.0d0))/180
  bb=(b*4.0d0*ATAN(1.0d0))/180
  Write(*,*) 'x=',lb, "  y=",bb
  Write(*,*) ''
  Write(*,*) 'x=',x2, "  y=", y2
  Call OGKGEO(y2, x2, b, l)
  Write(*,*) 'Geographisch'
  Write(*,*) 'x=',l, "  y=",b 
  Write(*,*) 'Bogenmaß:'
  lb=(l*4.0d0*ATAN(1.0d0))/180
  bb=(b*4.0d0*ATAN(1.0d0))/180
  Write(*,*) 'x=',lb, "  y=",bb
  Write(*,*) ''
   
  Write(*,*) "Bogenmaß 360 Grad =",(360*4.0d0*ATAN(1.0d0))/180  
  Write(*,*) "Bogenmaß 180 Grad =",(180*4.0d0*ATAN(1.0d0))/180
  Write(*,*) "Bogenmaß  90 Grad =",( 90*4.0d0*ATAN(1.0d0))/180
  Write(*,*) "Bogenmaß  45 Grad =",( 45*4.0d0*ATAN(1.0d0))/180
  Write(*,*) "Bogenmaß  30 Grad =",( 30*4.0d0*ATAN(1.0d0))/180
  WRITE (*,*) 'Program Ende'
 

  WRITE(*,*) "Would you like convert Gauss-Krueger --> Geo: y/n :"
  READ(*,*) wahl
  DO WHILE (wahl=='y')
    Write(*,*) "Gauss-Krueger Input : Bsp"
    Write(*,*) 'x=', "   y="
    READ(*,*) x,y
    Call OGKGEO(y, x, b, l)
    Write(*,*) 'Geographisch'
    Write(*,*) 'x=',l, "  y=",b
    Write(*,*) 'Bogenmaß:'
    lb=(l*4.0d0*ATAN(1.0d0))/180
    bb=(b*4.0d0*ATAN(1.0d0))/180
    Write(*,*) 'x=',lb, "  y=",bb
    Write(*,*) ''
 
  WRITE(*,*) "Would you like convert Gauss-Krueger --> Geo: y,n :"
  READ(*,*) wahl
  END DO

END PROGRAM MainProg
