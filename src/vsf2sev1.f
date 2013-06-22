C    Last change 2013 Iacopo
C    LAST CHANGE 07 MARCH 2006 RBT
C=======================================================================
C Create the new module for SEV and a new SEVTYPE variable: SEVDAT  
      Module  SEVMODULE
        INTEGER, SAVE, POINTER :: ISEVCB,ISEVOC,LCSEVOC,ISEVCL,ICSEVCL
        INTEGER, SAVE, POINTER :: LCPEV,LCHA,LCSRES,ICSEL
        REAL,    SAVE, DIMENSION(:,:), POINTER  :: PEV,SRES
        INTEGER, SAVE, DIMENSION(:,:), POINTER  :: SEL 
        REAL,    SAVE, DIMENSION(:,:), POINTER  :: SEVOC
        INTEGER, SAVE, DIMENSION(:,:), POINTER  :: SEVCL
        REAL,    SAVE, POINTER :: HA
        
        Type SEVTYPE
            INTEGER, POINTER :: ISEVCB,ISEVOC,LCSEVOC,ISEVCL,ICSEVCL
            INTEGER, POINTER :: LCPEV,LCHA,LCSRES,ICSEL
            REAL, DIMENSION(:,:),    POINTER  :: PEV,SRES
            INTEGER, DIMENSION(:,:), POINTER  :: SEL 
            REAL, DIMENSION(:,:),    POINTER  :: SEVOC
            INTEGER, DIMENSION(:,:), POINTER  :: SEVCL
            REAL, POINTER :: HA    
        End type SEVTYPE
        TYPE(SEVTYPE),SAVE :: SEVDAT(10)
      End module SEVMODULE
C================================
      Subroutine SVSF2SEV1PNT(IGRID)
C  Change SEV data to a different grid.
        use SEVMODULE
        
        ISEVCB=>SEVDAT(IGRID)%ISEVCB
        ISEVOC=>SEVDAT(IGRID)%ISEVOC
        LCSEVOC=>SEVDAT(IGRID)%LCSEVOC
        ISEVCL=>SEVDAT(IGRID)%ISEVCL
        ICSEVCL=>SEVDAT(IGRID)%ICSEVCL
        LCPEV=>SEVDAT(IGRID)%LCPEV
        LCHA=>SEVDAT(IGRID)%LCHA
        LCSRES=>SEVDAT(IGRID)%LCSRES
        ICSEL=>SEVDAT(IGRID)%ICSEL
        PEV=>SEVDAT(IGRID)%PEV
        SRES=>SEVDAT(IGRID)%SRES
        SEL=>SEVDAT(IGRID)%SEL 
        SEVOC=>SEVDAT(IGRID)%SEVOC
        SEVCL=>SEVDAT(IGRID)%SEVCL
        HA=>SEVDAT(IGRID)%HA
        return
      end
C================================
      Subroutine SVSF2SEV1PSV(IGRID)
C  Save global data for a grid.
        use SEVMODULE
        
        SEVDAT(IGRID)%ISEVCB=>ISEVCB
        SEVDAT(IGRID)%ISEVOC=>ISEVOC
        SEVDAT(IGRID)%LCSEVOC=>LCSEVOC
        SEVDAT(IGRID)%ISEVCL=>ISEVCL
        SEVDAT(IGRID)%ICSEVCL=>ICSEVCL
        SEVDAT(IGRID)%LCPEV=>LCPEV
        SEVDAT(IGRID)%LCHA=>LCHA
        SEVDAT(IGRID)%LCSRES=>LCSRES
        SEVDAT(IGRID)%ICSEL=>ICSEL
        SEVDAT(IGRID)%PEV=>PEV
        SEVDAT(IGRID)%SRES=>SRES
        SEVDAT(IGRID)%SEL=>SEL 
        SEVDAT(IGRID)%SEVOC=>SEVOC
        SEVDAT(IGRID)%SEVCL=>SEVCL
        SEVDAT(IGRID)%HA=>HA
        return
      end
C================================
      Subroutine SVSF2SEV1DA(IGRID)
C  Save global data for a grid.
        use SEVMODULE        
        deallocate(SEVDAT(IGRID)%ISEVCB)
        deallocate(SEVDAT(IGRID)%ISEVOC)
        deallocate(SEVDAT(IGRID)%LCSEVOC)
        deallocate(SEVDAT(IGRID)%ISEVCL)
        deallocate(SEVDAT(IGRID)%ICSEVCL)
        deallocate(SEVDAT(IGRID)%LCPEV)
        deallocate(SEVDAT(IGRID)%LCHA)
        deallocate(SEVDAT(IGRID)%LCSRES)
        deallocate(SEVDAT(IGRID)%ICSEL)
        deallocate(SEVDAT(IGRID)%PEV)
        deallocate(SEVDAT(IGRID)%SRES)
        deallocate(SEVDAT(IGRID)%SEL)
        deallocate(SEVDAT(IGRID)%SEVOC)
        deallocate(SEVDAT(IGRID)%SEVCL)
        deallocate(SEVDAT(IGRID)%HA)
        return
      end
C================================
C=======================================================================
      SUBROUTINE VSF2SEV1AR(IN,IGRID)
C-----Version 2011
C-----VERSION 07MARCH2006 RBT
C     ******************************************************************
C	ALLOCATE ARRAY SPACE FOR SURFACE EVAPORATION VARIABLES
C     ******************************************************************
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,IOUT
      use VSFMODULE,    ONLY:ISUM,ISUMIR,ZERO
      use SEVMODULE
C     ------------------------------------------------------------------
      INTEGER IN
      
      CHARACTER*200 LINE
C     ------------------------------------------------------------------
C Iac: 
C ------Allocate SEV scalar variables, which makes it possible 
C ------for multiple grids to be defined.
      allocate(ISEVCB,ISEVOC,LCSEVOC,ISEVCL,ICSEVCL)
      allocate(LCPEV,LCHA,LCSRES,ICSEL)
      allocate(HA)
C ------------      
  570 FORMAT(1X,'SURFACE EVAPORATION FLOW RATES WILL BE SAVED ON'
     &,' UNIT',I4)
  575 FORMAT(1X,'SURFACE EVAPORATION FLOW RATES WILL BE SAVED FOR'
     &,' ALL TIME STEPS IN EACH STRESS PERIOD')
  576 FORMAT(1X,'NO SURFACE EVAPORATION OUTPUT REQUESTED')
  577 FORMAT(1X,'SEV FLOW RATES WILL BE SAVED FOR ALL MODEL CELLS')
  578 FORMAT(1X,'SEV FLOW RATES WILL BE SAVED FOR ',1X,I5,1X'
     & USER-DEFINED CELL LOCATIONS')
  580 FORMAT(1X,'SURFACE EVAPORATION FLOW RATES WILL BE SAVED FOR',
     & 1X,I5,1X,'USER-DEFINED TIMES DURING THE SIMULATION')	
  585 FORMAT(1X,I10,' ELEMENTS IN RX ARRAY ARE USED BY SEV') 
  590 FORMAT(1X,I10,' ELEMENTS IN IR ARRAY ARE USED BY SEV')
C--------------------------------------------------------------------   	
C
C1-----READ ISEVCB, ISEVOC, AND ISEVCL FROM INPUT FILE
      LLOC=1
      CALL URDCOM(IN,0,LINE)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ISEVCB,R,IOUT,IN)
	IF (ISEVCB.NE.0) THEN
        LLOC=1
        CALL URDCOM(IN,0,LINE)
	  CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ISEVOC,R,IOUT,IN)
	  CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ISEVCL,R,IOUT,IN)
	ENDIF
C		
C2---- IF NO OUTPUT REQUESTED, WRITE NOTICE TO OUTPUT FILE AND SET FLAGS
      IF (ISEVCB.EQ.ZERO) THEN
	  WRITE(IOUT,576)
	  LCSEVOC=1
	  ISEVOC=0 
	ENDIF
C
C3---- PRINT UNIT NUMBER.TO MODFLOW OUTPUT FILE
	IF (ISEVCB.NE.0) WRITE(IOUT,570) ABS(ISEVCB)
C
C4--- PRINT OUTPUT CONTROL SETTING TO MODFLOW OUTPUT FILE
	IF (ISEVOC.LE.ZERO) THEN
	    WRITE(IOUT,575)
		LCSEVOC=1
	ELSE
	    WRITE(IOUT,580)ISEVOC
	ENDIF
	IF (ISEVCL.LE.ZERO) THEN
		WRITE (IOUT,577)
	    ICSEVCL=1
	ELSE
	   WRITE(IOUT,578)
	ENDIF
C
C Iac:  Allocate SEV arrays
      allocate(SEVOC(ISEVOC,2),SEVCL(ISEVCL,2))
      allocate(PEV(NCOL,NROW),SRES(NCOL,NROW),SEL(NCOL,NROW))
        
C5--- ALLOCATE SPACE IN RX ARRAY FOR EVAPORATION ARRAYS
      IRK=ISUM
	IIRK=ISUMIR
	IF (ISEVCB.NE.ZERO.AND.ISEVOC.GT.ZERO) THEN
		LCSEVOC=ISUM
		ISUM=ISUM+ISEVOC*2
	ENDIF
	IF (ISEVCB.NE.ZERO.AND.ISEVCL.GT.ZERO) THEN
		ICSEVCL=ISUMIR
		ISUMIR=ISUMIR+ISEVCL*2
	ENDIF
	LCPEV=ISUM
	ISUM=ISUM+NROW*NCOL
	LCSRES=ISUM
	ISUM=ISUM+NROW*NCOL
	ICSEL=ISUMIR
	ISUMIR=ISUMIR+NROW*NCOL
C 
C6---- CALCULATE & PRINT AMOUNT OF SPACE USED BY SPF PACKAGE.
      IRK=ISUM-IRK
	IIRK=ISUMIR-IIRK
      WRITE(IOUT,585)IRK
	WRITE(IOUT,590)IIRK
C
C Iac: -----SAVE POINTERS TO DATA AND RETURN.
      CALL SVSF2SEV1PSV(IGRID)
C7---- RETURN.
      RETURN
      END
C
C=======================================================================
      SUBROUTINE VSF2SEV1RPP(IN,FNAME,IGRID)
C-----Version 2011
C     VERSION 07MARCH2006 VSF1SEV1RPP
C     ******************************************************************
C     READ EVAPOTRANSPIRATION OUTPUT CONTROL
C     ******************************************************************
C
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,IOUT,DELR,DELC
      use VSFMODULE,    ONLY:ZERO
      use SEVMODULE
C     ------------------------------------------------------------------
C     ------------------------------------------------------------------
      INTEGER NUM
      CHARACTER(len=4) EXT
      CHARACTER(len=200) FNAME,OFILE,LINE
C
      EXT='.SEO'
C Iac: 
      CALL SVSF2SEV1PNT(IGRID)       
C     ------------------------------------------------------------------
C
C1---- READ SEVOC ARRAY IF REQUESTED (IF ISEVOC IS POSITIVE)
      IF (ISEVOC.GT.ZERO) THEN
        DO 1 K=1,ISEVOC
            LLOC=1
            CALL URDCOM(IN,0,LINE)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SEVOC(K,1),IOUT,IN)
	        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SEVOC(K,2),IOUT,IN)		  
    1   CONTINUE
      ENDIF 
C
C2---- READ SEVCL ARRAY IF REQUESTED (IF ISEVCL IS POSITIVE)
      IF (ISEVCL.GT.ZERO) THEN
        DO 2 K=1,ISEVCL
		  LLOC=1
          CALL URDCOM(IN,0,LINE)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SEVCL(K,1),R,IOUT,IN)
	      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SEVCL(K,2),R,IOUT,IN)
    2	CONTINUE
      ENDIF
C
C3---- OPEN OUTPUT FILE
      IF (ISEVCB.NE.0) THEN
	    OFILE=FNAME
	    NUM=INDEX(OFILE,' ')
        WRITE(OFILE(NUM-4:NUM),'(A4)')EXT
        IF(ISEVCB.LT.0)THEN
          OPEN (UNIT=ABS(ISEVCB),FILE=OFILE,FORM='UNFORMATTED',
     &	        STATUS='UNKNOWN')
        ELSE
	      OPEN (UNIT=ISEVCB,FILE=OFILE,FORM='FORMATTED',
     &	        STATUS='UNKNOWN')
        ENDIF
      ENDIF
C
C4--- RETURN
      RETURN
      END
C=======================================================================
      SUBROUTINE VSF2SEV1RP(IN,IGRID)
C-----Version 2011:
C     VERSION 07MARCH2006 VSF1SEV1RP
C     ******************************************************************
C     READ EVAPORATION DATA
C     ******************************************************************
C
C     SPECIFICATIONS:
C     IAC:
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,IOUT,DELR,DELC
      use VSFMODULE,    ONLY:ZERO
      use SEVMODULE
C     ------------------------------------------------------------------
      CHARACTER(len=24), dimension(4) :: ANAME
      CHARACTER(len=200) LINE
C     ------------------------------------------------------------------
C
      ANAME(1)= '    POTENTIAL EVAPORATION'
      ANAME(2)= '     ATMOSPHERIC PRESSURE'
      ANAME(3)= '       SURFACE RESISTANCE'
      ANAME(4)= '            SEV TOP LAYER'
C Iac: 
      CALL SVSF2SEV1PNT(IGRID)       
C     ------------------------------------------------------------------
C

C1------READ FLAGS SHOWING WHETHER DATA IS TO BE REUSED.
      CALL URDCOM(IN,0,LINE)
	LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,INPEV,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,INHA,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,INSRES,R,IOUT,IN)
	CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,INSEL,R,IOUT,IN)
C
C2------TEST INPEV TO SEE WHERE POTENTIAL EVAPORATION COMES FROM.
      IF(INPEV.GE.0)GO TO 32
C
C2A------IF INPEV<0 THEN REUSE POTENTIAL EVAPORATION ARRAY FROM LAST STRESS PERIOD
      WRITE(IOUT,3)
    3 FORMAT(1X,/1X,'REUSING PEV FROM LAST STRESS PERIOD')
      GO TO 35
C
C2B-------IF INPEV=>0 THEN CALL MODULE U2DREL TO READ POTENTIAL EVAPORATION.
   32 CALL U2DREL(PEV,ANAME(1),NROW,NCOL,0,IN,IOUT)
C
C3------MULTIPLY MAX ET RATE BY CELL AREA TO GET VOLUMETRIC RATE
      DO 40 IR=1,NROW
      DO 40 IC=1,NCOL
      IF (INPEV.GT.0) PEV(IC,IR)=PEV(IC,IR)*DELR(IC)*DELC(IR)
   40 CONTINUE
C
C4------TEST INHA TO SEE WHERE ATMOSPHERIC PRESSURE COMES FROM.
   35 IF(INHA.GE.0)GO TO 37
C
C4A-----IF INHA<0 THEN REUSE ATMOSPHERIC PRESSURE.
      WRITE(IOUT,4)
    4 FORMAT(1X,/1X,'REUSING HA FROM LAST STRESS PERIOD')
      GO TO 45
C
C4B------IF INHA=>0 READ ATMOSPHERIC PRESSURE.
   37 CALL URDCOM(IN,0,LINE)
	LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HA,IOUT,IN)
C
C5------TEST INSRES TO SEE WHERE SURFACE RESISTANCE COMES FROM
   45 IF(INSRES.GE.0)GO TO 47
C
C5A------IF INSRES<0 REUSE SURFACE RESISTANCE FROM LAST STRESS PERIOD
      WRITE(IOUT,5)
    5 FORMAT(1X,/1X,'REUSING SRES FROM LAST STRESS PERIOD')
      GO TO 100
C
C5B-------IF INSRES=>0 CALL MODULE U2DREL TO READ SURFACE RESISTANCE
   47 CALL U2DREL(SRES,ANAME(3),NROW,NCOL,0,IN,IOUT)
C
C6------TEST INSEL TO SEE WHERE TOP LAYER LOCATION COMES FROM
  100 IF(INSEL.GE.0)GO TO 57
C
C6A------IF INSEL<0 REUSE TOP LAYER LOCATION FROM LAST STRESS PERIOD
      WRITE(IOUT,6)
    6 FORMAT(1X,/1X,'REUSING SEL FROM LAST STRESS PERIOD')
      GO TO 110
C
C6B-------IF INSEL=>0 CALL MODULE U2DINT TO READ TOP LAYER LOCATION
   57 CALL U2DINT(SEL,ANAME(4),NROW,NCOL,0,IN,IOUT)
C
C7-----RETURN
  110 RETURN
      END
C
C====================================================================
C--------------------------------------------------------------------
      SUBROUTINE VSF2SEV1FM(IGRID)
C     Version 2011
C     VERSION 07MARCH2007 VSF1SEV1FM
C     ******************************************************************
C     COMPUTE SOIL EVAPORATION
C     ******************************************************************
C
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,HNEW,BOTM,NBOTM,IBOUND,
     &                       DELR,DELC,HCOF,RHS,LAYCBD,LBOTM,IOUT,BUFF
      use GWFLPFMODULE, ONLY:HK
      use VSFMODULE,    ONLY:ALPHA,VGN,PTAB,KTAB,STYPE,ISC,SLNUM,NTAB
      use SEVMODULE
C     ------------------------------------------------------------------
      REAL NODE,HD,BBOT,TTOP,PRESS,SR,THCK1,THCK2,THCK,VCOND,ALP1,
     &	   DLP,PTAB1,PTABN,ONE,KKREL
      INTEGER KK               
C	
C Iac: 
      CALL SVSF2SEV1PNT(IGRID) 
C     ------------------------------------------------------------------
	ONE=1.0
C
C1---- ESTABLISH BOUNDARY VALUES FOR INTERNAL SOIL CHARACTERISTICS TABLE
	ALP1=ALOG10(-PTAB(1))
	DLP=(ALOG10(-PTAB(NTAB))-ALP1)/(NTAB-1) 
	PTAB1=PTAB(1)
	PTABN=PTAB(NTAB)
C
C2------PROCESS EACH HORIZONTAL CELL LOCATION
      DO 10 IR=1,NROW
      DO 10 IC=1,NCOL
		IL=SEL(IC,IR)
		HD=HNEW(IC,IR,IL)
		BBOT=BOTM(IC,IR,LBOTM(IL))
	    TTOP=BOTM(IC,IR,LBOTM(IL)-1)
	    NODE=(TTOP+BBOT)/2		
	    PRESS=HD-NODE
		SR=SRES(IC,IR)
C
C3------IF THE CELL IS EXTERNAL IGNORE IT.
    4 IF(IBOUND(IC,IR,IL).LE.0)GO TO 10
C
C4---- CALCULATE CELL THICKNESS AND VERTICAL K
		M=STYPE(IC,IR,IL)
		IF (M.LE.0)CYCLE
C5--- CALCULATE RELATIVE PERMEABILITIES 
C
C5A--- IF CELL IS SATURATED, SET KKREL = 1.0
	    IF (HD.GE.NODE) THEN
		  KKREL=ONE
C
C5B--- OTHERWISE, LOOK UP K IN TABLE IF PRESSURE IS WITHIN SPECIFIED INTERVAL
	    ELSEIF (PRESS.GE.PTABN.AND.PRESS.LE.PTAB1)THEN
	      IT=INT((ALOG10(-PRESS)-ALP1)/DLP)+1
		  S1=(KTAB(IT+1,M)-KTAB(IT,M))/(PTAB(IT+1)-PTAB(IT))
		  CI=KTAB(IT,M)+S1*(PRESS-PTAB(IT))
		  KKREL=CI
C
C5C--- IF PRESSURE IS NOT WITHIN INTERVAL, COMPUTE DIRECTLY
	    ELSE
	      KKREL=SVSF1REF1KP(PRESS,ALPHA(M,1),VGN(M,1),ISC(M))	
	    ENDIF
		VCOND=HK(IC,IR,IL)*KKREL
C
C6--- COMPUTE ACTUAL EVAPORATION AND COMPARE IT TO POTENTIAL
C6--- NOTE THAT EVAPORATIVE RATE IS POSITIVE (MODFLOW CONVENTION)
		EV=VCOND*SR*(PRESS-HA)*DELC(IR)*DELR(IC)
		IF(EV.LT.0.0D0) THEN
		   EV=0.0D0
		   CYCLE
		ENDIF
C
		IF(EV.LT.PEV(IC,IR)) THEN
C6A--- SOIL EVAPORATION < POTENTIAL
			RHS(IC,IR,IL)=RHS(IC,IR,IL)+EV
		ELSE
C6B--- POTENTIAL EVAPORATION 
			RHS(IC,IR,IL)=RHS(IC,IR,IL)+PEV(IC,IR)
	    END IF
   10 CONTINUE
C
C7------RETURN
      RETURN
      END
C
C--------------------------------------------------------------------
      SUBROUTINE VSF2SEV1BD(KSTP,KPER,IGRID)
C-----Version 2011
C-----VERSION 07MARCH2006 VSF1SEV1BD
C     ******************************************************************
C     CALCULATE VOLUMETRIC BUDGET FOR EVAPOTRANSPIRATION
C     ******************************************************************
C
C        SPECIFICATIONS:
! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,HNEW,BOTM,NBOTM,BUFF,
     &                       DELR,DELC,LAYCBD,LBOTM,IBOUND,IOUT
      use GWFBASMODULE, ONLY:DELT,PERTIM,TOTIM,ICHFLG,HDRY,ICBCFL,
     &                       VBVL,VBNM,MSUM 
      use GWFLPFMODULE, ONLY:HK,LAYAVG,LAYTYP,WETFCT,
     &                       IWETIT,IHDWET,ILPFCB
      use GWFBCFMODULE, ONLY:IBCFCB
      use VSFMODULE
      use SEVMODULE
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      REAL    NODE,HD,BBOT,TTOP,PRESS,SR,THCK1,THCK2,THCK,VCOND,ALP1,
     &	      DLP,PTAB1,PTABN,ONE,KKREL

      DOUBLE PRECISION RATOUT,QQ,HH,SS,DD,XX
      DIMENSION UPER(ISEVOC),USTP(ISEVOC)
C
      DATA TEXT /'SOIL EVAPORATION'/
C     ------------------------------------------------------------------
C Iac: 
      CALL SVSF2SEV1PNT(IGRID)  
C
C1------CLEAR THE RATE ACCUMULATOR.
      RATOUT=ZERO
	  LIM=0.0
C
C2---- CLEAR THE BUFFER.
      DO 2 IL=1,NLAY
      DO 2 IR=1,NROW
      DO 2 IC=1,NCOL
      BUFF(IC,IR,IL)=ZERO
      
    2 CONTINUE
	ONE=1.0
C
C3---- ESTABLISH BOUNDARY VALUES FOR INTERNAL SOIL CHARACTERISTICS TABLE
      ALP1=ALOG10(-PTAB(1))
      DLP=(ALOG10(-PTAB(NTAB))-ALP1)/(NTAB-1) 
	  PTAB1=PTAB(1)
	  PTABN=PTAB(NTAB)
C
C4------PROCESS EACH HORIZONTAL CELL LOCATION
      DO 10 IR=1,NROW
      DO 10 IC=1,NCOL
C5------ ASSIGN TOP LAYER FROM SEL ARRAY
        IL=SEL(IC,IR)
		HD=HNEW(IC,IR,IL)
		BBOT=BOTM(IC,IR,LBOTM(IL))
	    TTOP=BOTM(IC,IR,LBOTM(IL)-1)
	    NODE=(TTOP+BBOT)/2		
	    PRESS=HD-NODE
		SR=SRES(IC,IR)
C
C6------IF THE CELL IS EXTERNAL IGNORE IT.
    4 IF(IBOUND(IC,IR,IL).LE.0) GO TO 10
        
C    
C7---- CALCULATE CELL THICKNESS AND VERTICAL K
        
        M=STYPE(IC,IR,IL)
        IF (M.LE.0) CYCLE
C
C8--- CALCULATE RELATIVE PERMEABILITIES 
C
C8A--- IF CELL IS SATURATED, SET KKREL = 1.0
        IF (HD.GE.NODE) THEN
		  KKREL=ONE
C
C8B--- OTHERWISE, LOOK UP K IN TABLE IF PRESSURE IS WITHIN SPECIFIED INTERVAL
        ELSE IF (PRESS.GE.PTABN.AND.PRESS.LE.PTAB1)THEN
	      IT=INT((ALOG10(-PRESS)-ALP1)/DLP)+1
		  S1=(KTAB(IT+1,M)-KTAB(IT,M))/(PTAB(IT+1)-PTAB(IT))
		  CI=KTAB(IT,M)+S1*(PRESS-PTAB(IT))
		  KKREL=CI
C
C8C--- IF PRESSURE IS NOT WITHIN INTERVAL, COMPUTE DIRECTLY
        ELSE
	      KKREL=SVSF1REF1KP(PRESS,ALPHA(M,1),VGN(M,1),ISC(M))	
        ENDIF
		VCOND=HK(IC,IR,IL)*KKREL
C
C9---- COMPUTE ACTUAL EVAPORATION AND COMPARE IT TO POTENTIAL
C9---- **NOTE THAT EVAPORATIVE RATE IS POSITIVE (MODFLOW CONVENTION)**
		EV=VCOND*SR*(PRESS-HA)*DELC(IR)*DELR(IC)
        IF(EV.LE.0.0D0) THEN
		   EV=0.0D0
		   CYCLE
        ENDIF
C
C10------IF ACTUAL EVAP >= PEV,SET Q=MAX ET RATE.
		 
        IF(EV.LT.PEV(IC,IR)) GOTO 7
		QQ=-PEV(IC,IR)
		GO TO 9
C
C10A------IF ACTUAL EVAP > PEV, Q=EV
    7		QQ=-EV
C
C10B-----ACCUMULATE TOTAL FLOW RATE.
    9		Q=QQ
		RATOUT=RATOUT-QQ
C
C11-----ADD Q TO BUFFER.
		BUFF(IC,IR,IL)=Q
C
   10 CONTINUE
       
C
C12------SET IBD TO INDICATE IF CELL-BY-CELL BUDGET VALUES WILL BE SAVED.
	IBD=0
C Iac: some problem here with the flag IBCFCB. However, this part is commented, 
C      beacuse of BCF is never used with the new version of VSF
      !IF(IBCFCB.LT.0) IBD=-1
      !IF(IBCFCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      !IF(ILPFCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      !IF(IBCFCB.GT.0) IBD=ICBCFL
      !IF(ILPFCB.GT.0) THEN
		!IBD=ICBCFL
		!IBCFCB=ILPFCB
      !ENDIF
C
C12A--- RECORD FLUXES IN CELL-BY-CELL FLOW FILE
	IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,
     1                       IBCFCB,SLM,NCOL,NROW,NLAY,IOUT)
      IF(IBD.EQ.2) CALL UBDSV1(KSTP,KPER,TEXT,IBCFCB,
     1            SLM,NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)

C
C13-----MOVE TOTAL ET RATE INTO VBVL FOR PRINTING BY BAS1OT.
      ROUT=RATOUT
      VBVL(3,MSUM)=ZERO
      VBVL(4,MSUM)=ROUT
C
C14-----ADD ET(ET_RATE TIMES STEP LENGTH) TO VBVL.
      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
C
C15-----MOVE BUDGET TERM LABELS TO VBNM FOR PRINT BY MODULE BAS1OT.
      VBNM(MSUM)=TEXT
C
C16-----INCREMENT BUDGET TERM COUNTER.
      MSUM=MSUM+1
C
C17-----PRINT THE FLOW FOR THE CELL IF REQUESTED.
C		
C17A--- DETERMINE REQUESTED STRESS PERIOD AND ELAPSED TIME
	IF (ISEVOC.GT.ZERO) THEN
	  DO 35 L=1,ISEVOC
		UPER(L)=INT(SEVOC(L,1))
		USTP(L)=SEVOC(L,2)
   35	  CONTINUE
   	  DO 42 L=1,ISEVOC
		IF (USTP(L).LT.ZERO) THEN
			STP=INT(ABS(USTP(L)))
			IF (KPER.EQ.UPER(L).AND.KSTP.EQ.STP) THEN
				GOTO 52
			ENDIF
		ELSE
		    DIFF=ABS(TOTIM-USTP(L))
		    IF (KPER.EQ.UPER(L).AND.DIFF.LE.LIM) THEN
				GOTO 52
			ENDIF
		ENDIF
   42	  CONTINUE
C18--- OUTPUT NOT REQUESTED FOR THIS TIMESTEP; RETURN
	  GOTO 901
	ENDIF
	IF (ISEVOC.EQ.ZERO) GOTO 901
C
C19---- WRITE OUTPUT 

   52 IF(ISEVCB.GT.0) THEN
C19A--- ASCII FILE
	  IBDLBL=0
C
	  DO 902 I=1,NROW
	  DO 902 J=1,NCOL
C19B------ ASSIGN TOP LAYER FROM SEL ARRAY
		IL=SEL(J,I)
C
C19C--- IF REQUESTED, CHECK CELL LOCATION FOR OUTPUT
		IF (ISEVCL.GT.ZERO) THEN
			ICL=0
			DO 903 L=1,ISEVCL
				IF (J.EQ.SEVCL(L,1).AND.I.EQ.SEVCL(L,2)) ICL=1
  903			CONTINUE
			IF (ICL.EQ.0) CYCLE
		ENDIF
C
	     if (IBDLBL.EQ.0)  then
            WRITE(ISEVCB,*) 'PERIOD ',KPER,'STEP ',KSTP,
     &                      'ELAPSED TIME ',TOTIM,
     &                      'TOTAL SEV FLOW FOR THIS TIME STEP =',ROUT
         end if

	    IBDLBL=1
C
C19D--- IF ISEVCL < 0, ONLY WRITE HEADER
	    IF (ISEVCL.LT.0) CYCLE
C
C19E--- WRITE FLOW RATE FOR CELL
          WRITE(ISEVCB,900) IL,I,J,BUFF(J,I,IL)
  900     FORMAT(1X,'LAYER',I4,'   ROW',I4,'   COL',I4,
     1       '   FLOW',G15.6)
  902   CONTINUE
	ELSEIF (ISEVCB.LT.0) THEN
C19F--- BINARY FILE
	  CALL UBDSV3(KSTP,KPER,TEXT,ISEVCB,BUFF,PFLG,1,
     1       NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
	ENDIF
C
C20-----RETURN.
  901 RETURN
      END


