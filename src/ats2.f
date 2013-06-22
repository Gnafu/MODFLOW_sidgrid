C     Last change  2013   Iacopo
C     Last change:  RBT   MARCH52006  12:30 PM
C================================
C Create the new module for ATS and a new ATSTYPE variable: ATSDAT  
      Module  ATSMODULE
        INTEGER, SAVE, POINTER ::  LCTIM,MXSTP,IMXSTP,ITIM
        REAL,    SAVE, POINTER ::  TMX,TMN,TSM,TSD
        REAL,    SAVE, DIMENSION(:),  POINTER  :: TIM 
        
        Type ATSTYPE
            INTEGER, POINTER ::  LCTIM,MXSTP,IMXSTP,ITIM
            REAL,    POINTER ::  TMX,TMN,TSM,TSD
            REAL,    DIMENSION(:),  POINTER  :: TIM 
        End type ATSTYPE
        
        TYPE(ATSTYPE),SAVE :: ATSDAT(10)
      
      End module ATSMODULE
C================================
      Subroutine SATS2PNT(IGRID)
C  Change ATS data to a different grid.
        use ATSMODULE
        
        LCTIM=>ATSDAT(IGRID)%LCTIM
        MXSTP=>ATSDAT(IGRID)%MXSTP
        IMXSTP=>ATSDAT(IGRID)%IMXSTP
        ITIM=>ATSDAT(IGRID)%ITIM
        TMX=>ATSDAT(IGRID)%TMX
        TMN=>ATSDAT(IGRID)%TMN
        TSM=>ATSDAT(IGRID)%TSM
        TSD=>ATSDAT(IGRID)%TSD
        TIM=>ATSDAT(IGRID)%TIM
                                        
        return
      end
C================================
      Subroutine SATS2PSV(IGRID)
C  Save global data for a grid.
        use ATSMODULE
                        
        ATSDAT(IGRID)%LCTIM=>LCTIM
        ATSDAT(IGRID)%MXSTP=>MXSTP
        ATSDAT(IGRID)%IMXSTP=>IMXSTP
        ATSDAT(IGRID)%ITIM=>ITIM
        ATSDAT(IGRID)%TMX=>TMX
        ATSDAT(IGRID)%TMN=>TMN
        ATSDAT(IGRID)%TSM=>TSM
        ATSDAT(IGRID)%TSD=>TSD
        ATSDAT(IGRID)%TIM=>TIM
        
        return
      end
C================================
      Subroutine SATS2DA(IGRID)
C  Save global data for a grid.
        use ATSMODULE       
C        ARGUMENTS

        INTEGER IGRID        
C ----------------------------------------------------------------------
        CALL SATS2PNT(IGRID)

        deallocate(LCTIM)
        deallocate(MXSTP)
        deallocate(IMXSTP)
        deallocate(ITIM)
        deallocate(TMX)
        deallocate(TMN)
        deallocate(TSM)
        deallocate(TSD)
        deallocate(TIM)
        
        return
      end
C================================
C=======================================================================
      SUBROUTINE ATS2AR(FNAME,IN,IGRID)
C ----Version 2011
C-----VERSION 28SEPT2004 RBT
C     ******************************************************************
C	READ OUTPUT FLAG AND OPEN TIME STEP INFO FILE
C     ALLOCATE ARRAY STORAGE FOR ADAPTIVE TIME STEPPING MODULE 
C     ******************************************************************
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:IOUT
      use VSFMODULE,    ONLY:ZERO,ISUM 
      use ATSMODULE
C     ------------------------------------------------------------------
      CHARACTER(len=200) LINE,EXT
      CHARACTER(len=200) OFILE,FNAME
      INTEGER            IFREFM,IRK,IN,NUM
      REAL               TSMULT,TSDIV
      
      EXT = '.tim'    
C     ------------------------------------------------------------------
C Iac: 
C Allocate ATS scalar variables:
      allocate(LCTIM,MXSTP,IMXSTP,ITIM,TMX,TMN,TSM,TSD)
C      
  570 FORMAT(1X,'TIME STEP INFORMATION WILL BE SAVED ON UNIT'
     &,I3)		
  585 FORMAT(1X,I10,' ELEMENTS IN X ARRAY ARE USED BY ATS')
  590 FORMAT(1X,'TIME STEP INFORMATION WILL NOT BE SAVED')
C
C
C1---- READ UNIT / FLAG FOR OUTPUT FILE AND TIME STEP ARRAY DIMENSION (MXSTP)
      LLOC=1
      CALL URDCOM(IN,IOUT,LINE)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ITIM,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXSTP,R,IOUT,IN)
C
C Iac:  Allocate ATS array TIM
      allocate(TIM(MXSTP))
C2----  IF TIME FILE REQUESTED (ITIM>0) OPEN UNFORMATTED OUTPUT FILE
      IF (ITIM.GT.0) THEN
C2A--- PRINT OUTPUT FILE NUMBER
       WRITE (IOUT,570) ITIM
C
	   OFILE=FNAME
	   NUM=INDEX(OFILE,' ')
	   WRITE(OFILE(NUM-4:NUM),'(A4)')EXT
	   OPEN (UNIT=ITIM,FILE=OFILE,FORM='UNFORMATTED',STATUS='UNKNOWN')
C 
C3---- ALLOCATE SPACE FOR THE ARRAY TIM(IF REQUESTED)
		IRK=ISUM
		LCTIM=ISUM
		ISUM=ISUM+MXSTP
C
C4---- CALCULATE & PRINT AMOUNT OF SPACE USED BY ATS PACKAGE.
	    IRK=ISUM-IRK
		WRITE(IOUT,585)IRK
      ELSE
		WRITE(IOUT,590)
		LCTIM=1 
      ENDIF
C
C5---- READ  MAX AND MIN TIMESTEPS AND TIMESTEP FACTORS
     
      CALL URDCOM(IN,IOUT,LINE)
      
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TMX,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TMN,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TSM,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,1,TSD,IOUT,IN)
      IMXSTP=MXSTP
C	
C Iac: -----SAVE POINTERS TO DATA AND RETURN.
      CALL SATS2PSV(IGRID)
C6---- RETURN.
      
      RETURN
      END
C
C=======================================================================
        SUBROUTINE ATS2ADV(PTIME,KSTP,DLT)
C-----Version 2011 
C-----VERSION 24SEP2004 RBT
C     ******************************************************************
C	ADVANCE TIME STEP BASED ON INITIAL USER-SPECIFIED DELT 
C     ******************************************************************
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
C        
        REAL DLT,RKSTP,DIFF
        INTEGER KSTP,IKSTP
        DOUBLE PRECISION PTIME
C     ------------------------------------------------------------------
C      
C1--- COMPUTE KSTP BASED ON ELAPSED TIME AND INTIAL DELTA T SIZE
	  IKSTP=INT(PTIME/DLT) 
	  RKSTP=PTIME/DLT
	  DIFF=RKSTP-IKSTP
      IF (DIFF.GT.0.9999)THEN 
	    KSTP=IKSTP+1
      ELSE
		KSTP=IKSTP
      ENDIF
      IF (KSTP.EQ.0) KSTP=1
C
C2---- RETURN
      RETURN
      END
C
C=======================================================================
      SUBROUTINE ATS2RST(IRST,ICNVG,AKSTP,PTIME,IGRID) 
C-----Version 2011
C-----VERSION 6JUNE2003 RBT
C     ******************************************************************
C	NON-CONVERGENCE: ADJUST TIME STEP ACCORDINGLY
C     ******************************************************************
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY
      use GWFBASMODULE, ONLY:DELT,PERTIM,TOTIM
      use ATSMODULE 
C     ------------------------------------------------------------------
      INTEGER IRST,ICNVG,AKSTP
      DOUBLE PRECISION PTIME
C     ------------------------------------------------------------------
C Iac: 
      CALL SATS2PNT(IGRID)      
C-----      
	IRST=0
C
C1---- DID NOT CONVERGE, CHECK DELT AND REDUCE IF POSSIBLE
	IF (ICNVG.NE.1) THEN
		IF (DELT.LE.TMN) THEN
			WRITE(*,5)
    5 FORMAT ('Simulation aborted: ATS - Did not converge with minimum
     & time step')
			STOP
		ELSEIF (DELT.GT.TMN) THEN
			PERTIM=PERTIM-DELT
			TOTIM=TOTIM-DELT
			PTIME=PTIME-DELT
			DELT=DELT/TSD
			IF(DELT.LT.TMN) DELT=TMN
		ENDIF
C
C2---- RESET PERIOD AND TOTAL SIMULATION TIME COUNTERS
		PERTIM=PERTIM+DELT
		TOTIM=TOTIM+DELT
		PTIME=PTIME+DELT
C
C3---- SET RESTART FLAG
		IRST=1
	ELSE
C
C4---- CONVERGED? STORE DELT VALUE IN TIME OUTPUT ARRAY (IF ACTIVATED)
		IF(ITIM.GT.0) TIM(AKSTP)=DELT
	ENDIF
C
C5---- RETURN
	RETURN
	END 
C
C
C=======================================================================
      SUBROUTINE ATS2ADJ(KKITER,NSTP,IDE4,ITMX,PERL,ENDPER,
     &                   KPER,KSTP,RST,IGRID)
C-----Version 2011
C-----VERSION 6MARCH2006 RBT
C     ******************************************************************
C	ADJUST TIME STEP BASED ON CURRENT TIME STEP'S ITERATION RATIO
C     ******************************************************************
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:MXITER
      use GWFBASMODULE, ONLY:DELT,TOTIM,IPEROC,ITSOC,PERTIM
      use ATSMODULE
C     ------------------------------------------------------------------
      REAL      MINCV,MAXCV,RATIO,KI,MI,PL,DT,OCTIME
      REAL      PERL  ! The PERLEN value (length of stress period) at this period
      INTEGER   KKITER,KPER,ENDPER,NUMSTP,AKSTP,FLG,RST
      integer   IDE4  ! The Unit number for input file *.DE4 (solver) 
C     ------------------------------------------------------------------     
C-----     
C
	ONE=1.0  
	ENDPER=0
	FLG=0
	RST=0
	KI=KKITER
	MI=MXITER
      IF (IDE4.GT.0) MI=ITMX
	MINCV=0.35
	MAXCV=0.65
	RATIO=KI/MI

C
C1---- IF CONVERGENCE IS WITHIN 35% OF MAXIMUM NUMBER OF ITERATIONS,
C1---- INCREASE BY TSMULT
      IF (RATIO.LT.MINCV) THEN
		DELT=DELT*TSM
C2---- IF CONVERGENCE IS ABOVE 65% OF MAXIMUM, REDUCE BY TSDIV
      ELSE IF (RATIO.GT.MAXCV)THEN
        DELT=DELT/TSD
      END IF

C3---- MAKE SURE TMN<DELT<TMX
10	   IF (DELT.GT.TMX) DELT=TMX
       IF (DELT.LT.TMN) DELT=TMN 
C4---- IF PERIOD LENGTH IS EXCEEDED, ADJUST DELT
      PL=DELT+PERTIM
      IF (PL.GT.PERLEN.AND.PERTIM.NE.PERL) THEN
		DELT=PERL-PERTIM
        !Iac: added the following check: 
        !IF (DELT.GT.TMX) DELT=TMX
        !IF (DELT.LT.TMN) DELT=TMN 
      END IF
C5---- CONSTRAIN TIME STEP TO REQUESTED OUTPUT CONTROL TIME
      IF (KPER.EQ.IPEROC)THEN
		OCTIME=ITSOC*(PERL/NSTP)
        IF(PL.GT.OCTIME) THEN
			DELT=OCTIME-PERTIM
        ENDIF
      ENDIF
C6---- IF DELT IS LAST STEP IN CURRENT STRESS PERIOD SET FLAG AND RETURN
	IF (PERTIM.GE.PERL)THEN
		ENDPER=1
		NUMSTP=KSTP
C7---- IF OUTPUT ACTIVATED, WRITE TIME INFORMATION FOR CURRENT STRESS PERIOD TO OUTPUT FILE
		IF (ITIM.GT.0) THEN
			WRITE (ITIM) KPER,PERL,NUMSTP
			WRITE (ITIM)(TIM(I),I=1,NUMSTP)
C8---- RE-INITIALIZE TIME STEP ARRAY
			TIM=0.
		ENDIF
	ENDIF
C
C9---- RETURN
100	RETURN
      END
