MODULE IOControl_Mod

  USE Kind_Mod

  IMPLICIT NONE

   INTEGER :: OutUnitProt=24
   INTEGER :: OutUnitChL =25
   INTEGER :: OutUnitB   =91   !zu Bounds
   INTEGER :: OutUnitBpv =81   ! ""
   INTEGER :: TOutUnitB  =71
   INTEGER :: TOutUnitBpv=61

   CHARACTER(1)  :: leerzei1=' '
   CHARACTER(2)  :: leerzei2='  '
   CHARACTER(3)  :: leerzei3='   '
   CHARACTER(5)  :: leerzei5='     '
   CHARACTER(6)  :: leerzei6='      '
   CHARACTER(7)  :: leerzei7='       '
   CHARACTER(8)  :: leerzei8='        '
   CHARACTER(9)  :: leerzei9='         '
   CHARACTER(13) :: leerzei13="             "
   CHARACTER(14) :: leerzei14='              '
   CHARACTER(15) :: leerzei15='               '
   CHARACTER(16) :: leerzei16='                '
   CHARACTER(18) :: leerzei18='                  '
   CHARACTER(20) :: leerzei20='                    '
   CHARACTER(22) :: leerzei22='                      ' 
   CHARACTER(24) :: leerzei24='                        ' 
   CHARACTER(6)  :: pfeil=' ---> '
   CHARACTER(66) :: trenn1= &
        "------------------------------------------------------------------"
   CHARACTER(80) :: trenn2= &
        "==============================================================================="
   CHARACTER(75)  :: trenn_glob= &
    "......................................................................"
   CHARACTER(96)  :: version

CONTAINS

!-------------------------------------------------------------------------
! Progname_and_InputFile ; InfoProg 
!--------------------------------------------------------------------------
SUBROUTINE Prot_Progname_and_InputFile(ProgName,InputFile)
CHARACTER(100) :: ProgName
CHARACTER(80)  :: InputFile
  WRITE(OutUnitProt,*)
  Write(OutUnitProt,*) "###############################################################################"
  WRITE(OutUnitProt,*) TRIM(ProgName),"  ",TRIM(InputFile)
END SUBROUTINE Prot_Progname_and_InputFile

SUBROUTINE Prot_InfoProg(version)
  CHARACTER(96)  :: version
  Write(OutUnitProt,*) "###############################################################################"
  Write(OutUnitProt,*) "#       #                                                             #       #"
  Write(OutUnitProt,*) "#       #                         Program                             #       #"
  Write(OutUnitProt,*) "#       #         Creating: Weight-File for Model ASAM                #       #"
  Write(OutUnitProt,*) "#       #                           and                               #       #"
  Write(OutUnitProt,*) "#       #                   Output-File for GMV-Tool                  #       #"
  Write(OutUnitProt,*) "#       #                                                             #       #"
  Write(OutUnitProt,*) "###############################################################################"
  WRITE(OutUnitProt,*) version 
  WRITE(OutUnitProt,*)
END SUBROUTINE Prot_InfoProg

! Special beachten Linux,Mac, AIX
!SUBROUTINE Prot_ProgTime
!   INTEGER DATE_TIME(8)
!   CHARACTER (LEN =10) BIG_BEN(3)
!   CHARACTER(8) D, date
!   CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
!   ! date() unter Linux wenn gfortran nicht anzuwenden
!   !        wenn mac, kommentieren
!   !        wenn aix, okey         
!   D=date()
!   Write(OutUnitProt,'(a22,a8,a3,I2,a1,I2,a1,I2)') &
!                      leerzei22,D,leerzei3, &
!                      DATE_TIME(5),":",DATE_TIME(6),":",DATE_TIME(7)
!   Write(OutUnitProt,'(a22,a8,a3,a8)') leerzei22,"mm/dd/yy",leerzei3,"hh:mm:ss"
!   WRITE(OutUnitProt,*) trenn2
!END SUBROUTINE Prot_ProgTime

!program test_cpu_time
!             real :: start, finish
!             call cpu_time(start)
!                 ! put code to test here
!             call cpu_time(finish)
!             print '("Time = ",f6.3," seconds.")',finish-start
!end program test_cpu_time

!PROGRAM test_system_clock
!           INTEGER :: count, count_rate, count_max
!           CALL SYSTEM_CLOCK(count, count_rate, count_max)
!           WRITE(*,*) count, count_rate, count_max
!END PROGRAM

!program test_time_and_date  ! unter Linux okey, AIX sucht loader
!    character(8)  :: date
!    character(10) :: time
!    character(5)  :: zone
!    integer,dimension(8) :: values
!    ! using keyword arguments
!    call date_and_time(date,time,zone,values)
!    call date_and_time(DATE=date,ZONE=zone)
!    call date_and_time(TIME=time)
!    call date_and_time(VALUES=values)
!    print '(a,2x,a,2x,a)', date, time, zone
!    print '(8i5))', values
!end program test_time_and_date

!-------------------------------------------------------------------------
! Betreff   AnalyzeAllCells 
!--------------------------------------------------------------------------

SUBROUTINE WriteAuswAnalyzeAllCellsToProt
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) trenn2
   Write(OutUnitProt,*) leerzei3,">>>>     AnalyzeAllCells        <<<<"
   Write(OutUnitProt,*) trenn2
END SUBROUTINE WriteAuswAnalyzeAllCellsToProt

SUBROUTINE WriteEndAnalyzeAllCellsToProt
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) leerzei3,">>>>   Ende:  AnalyzeAllCells   <<<<"
   Write(OutUnitProt,*) trenn2
   WRITE(OutUnitProt,*)
END SUBROUTINE WriteEndAnalyzeAllCellsToProt
     
!-------------------------------------------------------------------------
! Betreff Weight2
!--------------------------------------------------------------------------
SUBROUTINE WriteAuswWeightToProt
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) trenn2
   Write(OutUnitProt,*) leerzei3, ">>>>     Auswertung  Weight2      <<<<"
   Write(OutUnitProt,*) trenn2
END SUBROUTINE WriteAuswWeightToProt

SUBROUTINE WriteEndWeightToProt
   Write(OutUnitProt,*) leerzei3
   Write(OutUnitProt,*) leerzei3, ">>>>   Ende Auswertung  Weight2   <<<<"
   Write(OutUnitProt,*) trenn2
   WRITE(OutUnitProt,*)
END SUBROUTINE WriteEndWeightToProt

SUBROUTINE WriteAuswWeightNullToProt
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) trenn2
   Write(OutUnitProt,*) leerzei6,"     >>>>   Auswertung  WNull   <<<<"
   Write(OutUnitProt,*) trenn2
END SUBROUTINE WriteAuswWeightNullToProt

SUBROUTINE WriteEndWeightNullToProt
   Write(OutUnitProt,*) leerzei6
   Write(OutUnitProt,*) leerzei6,"    >>>>   Ende Auswertung  WNull   <<<<"
   Write(OutUnitProt,*) trenn2
   WRITE(OutUnitProt,*)
END SUBROUTINE WriteEndWeightNullToProt

SUBROUTINE WriteBlkNrToProt(ib)
   INTEGER :: ib
   !Write(OutUnitProt,*) leerzei16,"Block :",ib
   Write(OutUnitProt,*) ""
   Write(OutUnitProt,*) leerzei7,"------------------------"
   Write(OutUnitProt,*) leerzei7," * Block :",ib
   Write(OutUnitProt,*) leerzei7,"-------------------------"
END SUBROUTINE WriteBlkNrToProt

SUBROUTINE WriteWarnNrFace(NrW,NrRW,NameFace)
   INTEGER :: NrW,NrRW
   CHARACTER(7) :: NameFace

   INTEGER :: erg

   Write(OutUnitProt,*) ""
   Write(OutUnitProt,*) leerzei9,pfeil,'Warning   NrW_',NameFace
   Write(OutUnitProt,*) leerzei6,'gezählte Faces  :',NrW,  '  NrW_',NameFace
   Write(OutUnitProt,*) leerzei6,'gelistete Faces :',NrRW, '  NrRW_',NameFace
   IF (NameFace=="FacesXY") THEN
       Write(OutUnitProt,*) leerzei6,"wenn 'ob' NrRN_FacesXY beachten !"
   END IF
   erg=NrW-NrRW
   IF(erg<0) THEN
     Write(OutUnitProt,*) leerzei9,erg," ",NameFace," zu wenig gezählt !"
   ELSE
     Write(OutUnitProt,*) leerzei9,erg," ",NameFace," zu viel gezählt !" 
   END IF
   Write(OutUnitProt,*) leerzei15, '! NrW_',NameFace,' angepasst !'
END SUBROUTINE WriteWarnNrFace

SUBROUTINE WriteWarnNrCells(NrW,NrRW)
   INTEGER :: NrW,NrRW

   INTEGER :: erg

   Write(OutUnitProt,*) ""
   Write(OutUnitProt,*) leerzei9,pfeil,"Warning   NrW_Cells"
   Write(OutUnitProt,*) leerzei6,"gezählte Cellen  :",NrW ,"   NrW_Cells"
   Write(OutUnitProt,*) leerzei6,"gelistete Cellen :",NrRW,"   NrRW_Cells"
   erg=NrW-NrRW
   IF(erg<0) THEN
     Write(OutUnitProt,*) leerzei6,erg," Cellen zu wenig gezählt"
   ELSE
     Write(OutUnitProt,*) leerzei6,erg," Cellen zu viel gezählt" 
   END IF
   Write(OutUnitProt,*) leerzei15,"! NrW_Cells angepasst !"
END SUBROUTINE WriteWarnNrCells


!-------------------------------------------------------------------------
! Betreff  SoilLayer 
!--------------------------------------------------------------------------

SUBROUTINE WriteAuswSoilLayerToProt
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) trenn2
   Write(OutUnitProt,*) leerzei3,">>>>     Auswertung  SoilLayer        <<<<"
   Write(OutUnitProt,*) trenn2
END SUBROUTINE WriteAuswSoilLayerToProt

SUBROUTINE WriteEndSoilLayerToProt
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) leerzei3,">>>>   Ende Auswertung  SoilLayer     <<<<"
   Write(OutUnitProt,*) trenn2
   WRITE(OutUnitProt,*)
END SUBROUTINE WriteEndSoilLayerToProt

!-------------------------------------------------------------------------
! Betreff   UpperBoundsLayer 
!--------------------------------------------------------------------------
SUBROUTINE WriteAuswUpperBoundsLayerProt
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) trenn2
   Write(OutUnitProt,*) leerzei3,">>>>     Auswertung  UpperBoundsLayer     <<<<"
   Write(OutUnitProt,*) trenn2
END SUBROUTINE WriteAuswUpperBoundsLayerProt

SUBROUTINE WriteEndUpperBoundsLayerProt
   Write(OutUnitProt,*) leerzei3
   Write(OutUnitProt,*) leerzei3,">>>>   Ende Auswertung  Output-ParaView   <<<<"
   Write(OutUnitProt,*) trenn2
   WRITE(OutUnitProt,*)
END SUBROUTINE WriteEndUpperBoundsLayerProt


!-------------------------------------------------------------------------
! Betreff   WriteAllCellsOroAsciiGMV (binary) 
!--------------------------------------------------------------------------

SUBROUTINE WriteAuswCellsOroGMVToProt
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) trenn2
   Write(OutUnitProt,*) leerzei3,">>>>     Auswertung  CellsOro...GMV      <<<<"
   Write(OutUnitProt,*) trenn2
END SUBROUTINE WriteAuswCellsOroGMVToProt

SUBROUTINE WriteEndCellsOroGMVToProt
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) leerzei3,">>>>  Ende Auswertung  CellsOro...GMV    <<<<"
   Write(OutUnitProt,*) trenn2
   WRITE(OutUnitProt,*)
END SUBROUTINE WriteEndCellsOroGMVToProt

!-------------------------------------------------------------------------
! Betreff   Output for ParaView 
!--------------------------------------------------------------------------
SUBROUTINE WriteAuswOutParaVToProt
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) trenn2
   Write(OutUnitProt,*) leerzei3,">>>>    Auswertung  Output-ParaView      <<<<"
   Write(OutUnitProt,*) trenn2
END SUBROUTINE WriteAuswOutParaVToProt

SUBROUTINE WriteEndOutParaVToProt
   Write(OutUnitProt,*) leerzei3
   Write(OutUnitProt,*) leerzei3,">>>>  Ende Auswertung  Output-ParaView   <<<<"
   Write(OutUnitProt,*) trenn2
   WRITE(OutUnitProt,*)
END SUBROUTINE WriteEndOutParaVToProt

!#########################################################################
!-------------------------------------------------------------------------
!  HP- Outputs Display 
!-------------------------------------------------------------------------

SUBROUTINE Display_Progname_and_InputFile(ProgName,InputFile)
CHARACTER(100) :: ProgName
CHARACTER(80)  :: InputFile
  WRITE (*,*)
  Write(*,*) "###############################################################################"
  WRITE (*,*) TRIM(ProgName),"  ",TRIM(InputFile) 
END SUBROUTINE Display_Progname_and_InputFile

SUBROUTINE Display_InfoProg(version)
CHARACTER(96)  :: version
  Write(*,*) "###############################################################################"
  Write(*,*) "#       #                                                             #       #"
  Write(*,*) "#       #                         Program                             #       #"
  Write(*,*) "#       #         Creating: Weight-File for Model ASAM                #       #"
  Write(*,*) "#       #                           and                               #       #"
  Write(*,*) "#       #                   Output-File for GMV-Tool                  #       #" 
  Write(*,*) "#       #                                                             #       #"
  Write(*,*) "###############################################################################"
  WRITE(*,*) version
  WRITE(*,*)

END SUBROUTINE Display_InfoProg


! Overview Display-Prints
!------------------------
! Write(*,*) ">....... Protocolls of Analysis   :  ", TRIM(OutputFile),'.prot'," ......."
! Write(*,*) ">....... Read Grid                :  ", TRIM(InputFileName)," ......."
! Write(*,*) ">....... Read Weights             :  ", TRIM(InputFileName)," ......."
! Write(*,*) ">....... Init Vertices            :  ", "........"
! Write(*,*) ">....... Analyze  Vertices        :  ", "........"
! Write(*,*) ">....... Init Edges               :  ", "........"
! Write(*,*) ">....... Analyze Edges            :  ", "........"
! Write(*,*) ">....... Init Faces               :  ", "........"
! Write(*,*) ">....... Analyze Faces            :  ", "........"
! Write(*,*) ">....... Init Cells               :  ", "........"
! Write(*,*) ">....... Analyze Cells            :  ", "........"
! Write(*,*) ">....... Analyze Cells of borders :  ", "........"
! Write(*,*) ">....... Output Weight            :  ", TRIM(OutputFile),'.Weights2'," ......."
! Write(*,*) ">....... Output Troposphere       :  ", "----", string_a, "----"
! Write(*,*) ">....... Output Troposphere       :  ", "----", string_b, "----"
! Write(*,*) ">....... Output GMV               :  ", TRIM(OutputFile),'.out.gmvG'," ......."
! Write(*,*) ">....... Output Cut GMV           :  ", TRIM(OutputFile),'.Cut.out.gmvG'," ......."
! Write(*,*) ">....... Output Cut2 GMV          :  ", TRIM(OutputFile),'.Cut2.out.gmvG'," ......."
! Write(*,*) ">....... Output Soil GMV          :  ", TRIM(OutputFile),'.Soil.out.gmvG'," ......."
! Write(*,*) ">....... SortVertAllFacesIn       :  ", "........"
! Write(*,*) ">....... SortVertCutAllCells      :  ", "........"
! Write(*,*) ">....... Output Orography GMV     :  ", TRIM(OutputFile),'.Oro.out.gmvG'," ......."
! Write(*,*) ">....... Protocolls of Analysis      "
! Write(*,*) ">....... Weight- and GMV-Output   :  ", TRIM(OutputFile),'.prot'," ......."
!
!Bsp: !WRITE(*,*) ("-",i=1,LEN_TRIM(name_fkt))
!--------------------------------------------------------------------------
SUBROUTINE Display_ProtFile(OutputFile)
  CHARACTER(80)  :: OutputFile
  Write(*,*) ">....... Protocolls of Analysis   :  ", TRIM(OutputFile),'.prot'
END SUBROUTINE Display_ProtFile

SUBROUTINE Display_GridFile(InputFileName)
  CHARACTER(80)  :: InputFileName
  Write(*,*) ">....... Read Grid                :  ", TRIM(InputFileName)
END SUBROUTINE Display_GridFile

SUBROUTINE Display_ProcWeights(InputFileName)
  CHARACTER(80)  :: InputFileName
  Write(*,*)
  Write(*,*) ">....... Read Weights             :  ", TRIM(InputFileName)
END SUBROUTINE Display_ProcWeights

SUBROUTINE Display_InitAllVertices
  Write(*,*) ">....... Init Vertices            :  ", "........"
END SUBROUTINE Display_InitAllVertices

SUBROUTINE Display_AnalyzeAllVertices
  Write(*,*) ">....... Analyze  Vertices        :  ", "........"
END SUBROUTINE Display_AnalyzeAllVertices

SUBROUTINE Display_InitAllEdges
  Write(*,*) ">....... Init Edges               :  ", "........"
END SUBROUTINE Display_InitAllEdges

SUBROUTINE Display_AnalyzeAllEdges
  Write(*,*) ">....... Analyze Edges            :  ", "........"
END SUBROUTINE Display_AnalyzeAllEdges

SUBROUTINE Display_InitAllFaces
  Write(*,*) ">....... Init Faces               :  ", "........"
END SUBROUTINE Display_InitAllFaces

SUBROUTINE Display_AnalyzeAllFaces
  Write(*,*) ">....... Analyze Faces            :  ", "........"
END SUBROUTINE Display_AnalyzeAllFaces

SUBROUTINE Display_InitAllCells
  Write(*,*) ">....... Init Cells               :  ", "........"
END SUBROUTINE Display_InitAllCells

SUBROUTINE Display_AnalyzeAllCells
  Write(*,*) ">....... Analyze Cells            :  ", "........"
END SUBROUTINE Display_AnalyzeAllCells

SUBROUTINE Display_AnalyzeAllCellsRand
  Write(*,*)
  Write(*,*) ">....... Analyze Cells of borders :  ", "........"
END SUBROUTINE Display_AnalyzeAllCellsRand


SUBROUTINE Display_OutWeightBlk(OutputFile)
  CHARACTER(80)  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Output Weight            :  ", TRIM(OutputFile),'.Weights2'
END SUBROUTINE Display_OutWeightBlk

SUBROUTINE Display_OutWeightNullBlk(OutputFile)
  CHARACTER(80)  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Output Weight            :  ", TRIM(OutputFile),'.WNull'
END SUBROUTINE Display_OutWeightNullBlk

SUBROUTINE Display_WriteOutputTopo(o_type)
  CHARACTER(1)  :: o_type 
  character(14) :: string_a=" ASCII-Format "
  character(15) :: string_b=" Binary-Format "
  if(o_type=='a') then
    Write(*,*)
    Write(*,*) "                                     ", "    --------------"
    Write(*,*) ">....... Output Troposphere       :  ", "----", string_a, "----"
    Write(*,*) "                                     ", "    --------------"
  else
    Write(*,*)
    Write(*,*) "                                     ", "    ---------------"
    Write(*,*) ">....... Output Troposphere       :  ", "----", string_b, "----"
    Write(*,*) "                                     ", "    ---------------"
  end if 
END SUBROUTINE Display_WriteOutputTopo

SUBROUTINE Display_OutGMVBlk(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Output GMV               :  ", TRIM(OutputFile),'.out.gmvG'
END SUBROUTINE Display_OutGMVBlk

SUBROUTINE Display_SetTroposOroParaV(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Set Troposphere and Orograhy for ParaView :  "
END SUBROUTINE Display_SetTroposOroParaV

SUBROUTINE Display_InitTropoOroOutStruct(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Init Tropo/Bounds/Oro    :  ", "Output-Struct"
END SUBROUTINE Display_InitTropoOroOutStruct

SUBROUTINE Display_OutTroposParaV(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Output ParaView          :  "
  Write(*,*) "          - Troposphere :            ", TRIM(OutputFile),'.pva.tropo'
END SUBROUTINE Display_OutTroposParaV

SUBROUTINE Display_OutTroposOroParaV(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*) "          - Troposphere-Orography :  ", TRIM(OutputFile),'.pva.gall'
END SUBROUTINE Display_OutTroposOroParaV

SUBROUTINE Display_OutUpperBoundsLayer(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Output BoundsLayer       :  ", TRIM(OutputFile),'.bound'
  Write(*,*) "                                     ", TRIM(OutputFile),'.pva.bound'
END SUBROUTINE Display_OutUpperBoundsLayer

SUBROUTINE Display_OutCutGMVBlk(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Output Cut GMV           :  ", TRIM(OutputFile),'.Cut.out.gmvG'
END SUBROUTINE Display_OutCutGMVBlk

SUBROUTINE Display_OutCut2GMVBlk(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Output Cut2 GMV          :  ", TRIM(OutputFile),'.Cut2.out.gmvG'
END SUBROUTINE Display_OutCut2GMVBlk

SUBROUTINE Display_OutSoilGMVBlk(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Output Soil GMV          :  ", TRIM(OutputFile),'.Soil.out.gmvG'
END SUBROUTINE Display_OutSoilGMVBlk

SUBROUTINE Display_SortVertAllFacesIn
  Write(*,*) ">....... SortVertAllFacesIn       :  ","........"
END SUBROUTINE Display_SortVertAllFacesIn

SUBROUTINE Display_SortVertCutAllCells
  Write(*,*) ">....... SortVertCutAllCells      :  ","........"
END SUBROUTINE Display_SortVertCutAllCells

SUBROUTINE Display_OutOroVolNullGMVBlk(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Output OrographyNull GMV :  ", TRIM(OutputFile),'.ONull.out.gmvG'
END SUBROUTINE Display_OutOroVolNullGMVBlk

SUBROUTINE Display_OutOroGMVBlk(OutputFile)
  CHARACTER*50  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Output Orography GMV     :  ", TRIM(OutputFile),'.Oro.out.gmvG'
END SUBROUTINE Display_OutOroGMVBlk

SUBROUTINE Display_InfoProtokollAnalysis(OutputFile)
CHARACTER(80)  :: OutputFile
  Write(*,*)
  Write(*,*) ">....... Protocol of all Analysis :  ", TRIM(OutputFile),'.prot'
END SUBROUTINE Display_InfoProtokollAnalysis

SUBROUTINE Display_Deallocation
  Write(*,*)
  Write(*,*) "                            ------------------------"
  Write(*,*) "                            | Deallocation storage |"
  Write(*,*) "                            ------------------------"
END SUBROUTINE Display_Deallocation

SUBROUTINE Display_InfoEndeProg(version)
CHARACTER(96)  :: version
  Write (*,*)
  WRITE (*,*) "         Program Ende!                           ", version
  WRITE (*,*) "###############################################################################" 
  Write (*,*)
END SUBROUTINE Display_InfoEndeProg

END MODULE IOControl_Mod

