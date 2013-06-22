C Last change 2013 Iacopo
C Last change 07 March 2006 RBT 12:40 PM
C
C================================
C Create the new module for PND and a new PNDTYPE variable: PNDDAT  
      Module  PNDMODULE
        INTEGER, SAVE, POINTER ::  LCPOND,LCPNDOC,IPNDOC,IPNDCB
        INTEGER, SAVE, POINTER ::  IPOND
        REAL,    SAVE, DIMENSION(:,:),   POINTER  :: PFLG 
        REAL,    SAVE, DIMENSION(:,:),   POINTER  :: POND
        REAL,    SAVE, DIMENSION(:,:),   POINTER  :: PNDOC
        
        Type PNDTYPE
            INTEGER, POINTER ::  LCPOND,LCPNDOC,IPNDOC,IPNDCB
            INTEGER, POINTER ::  IPOND
            REAL,    DIMENSION(:,:),   POINTER  :: PFLG 
            REAL, DIMENSION(:,:),   POINTER  :: POND,PNDOC
        End type PNDTYPE
        
        TYPE(PNDTYPE),SAVE :: PNDDAT(10)
      
      End module PNDMODULE
C================================
      Subroutine SVSF2PND1PNT(IGRID)
C  Change PND data to a different grid.
        use PNDMODULE
        
        LCPOND=>PNDDAT(IGRID)%LCPOND
        LCPNDOC=>PNDDAT(IGRID)%LCPNDOC
        IPNDOC=>PNDDAT(IGRID)%IPNDOC
        IPNDCB=>PNDDAT(IGRID)%IPNDCB
        IPOND=>PNDDAT(IGRID)%IPOND
        PFLG=>PNDDAT(IGRID)%PFLG 
        POND=>PNDDAT(IGRID)%POND
        PNDOC=>PNDDAT(IGRID)%PNDOC
                       
        return
      end
C================================
      Subroutine SVSF2PND1PSV(IGRID)
C  Save global data for a grid.
        use PNDMODULE
        
        PNDDAT(IGRID)%LCPOND=>LCPOND
        PNDDAT(IGRID)%LCPNDOC=>LCPNDOC
        PNDDAT(IGRID)%IPNDOC=>IPNDOC
        PNDDAT(IGRID)%IPNDCB=>IPNDCB
        PNDDAT(IGRID)%IPOND=>IPOND
        PNDDAT(IGRID)%PFLG=>PFLG 
        PNDDAT(IGRID)%POND=>POND
        PNDDAT(IGRID)%PNDOC=>PNDOC
        
        return
      end
C================================
      Subroutine SVSF2PND1DA(IGRID)
C  Save global data for a grid.
        use PNDMODULE       
        deallocate(PNDDAT(IGRID)%LCPOND)
        deallocate(PNDDAT(IGRID)%LCPNDOC)
        deallocate(PNDDAT(IGRID)%IPNDOC)
        deallocate(PNDDAT(IGRID)%IPNDCB)
        deallocate(PNDDAT(IGRID)%IPOND)
        deallocate(PNDDAT(IGRID)%PFLG)
        deallocate(PNDDAT(IGRID)%POND)
        deallocate(PNDDAT(IGRID)%PNDOC)
        
        return
      end
C================================
C=======================================================================
      SUBROUTINE VSF2PND1AR(IN,IGRID)
C-----Version 2001
C-----VERSION 07MARCH2006 RBT
C     ******************************************************************
C	ALLOCATE ARRAY SPACE FOR SURFACE PONDING VARIABLES
C     ******************************************************************
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,BOTM,NBOTM,CC,CR,CV,
     &                       HOLD,HCOF,RHS,LAYCBD,LBOTM,IBOUND,
     &                       DELR,DELC,IOUT  
      use VSFMODULE,    ONLY:ZERO,ISUM,ISUMIR 
      use PNDMODULE  
C     ------------------------------------------------------------------
      INTEGER IN

      CHARACTER*200 LINE
C Iac: 
C Allocate PND scalar variable and PFLG, POND, arrays
      allocate(LCPOND,LCPNDOC,IPNDOC,IPNDCB)
      allocate(IPOND)
      allocate(PFLG(NCOL,NROW),POND(NCOL,NROW))
C     ------------------------------------------------------------------
  570 FORMAT(1X,'PONDED CELL LOCATIONS WILL BE SAVED ON'
     &,' UNIT',I4)
  575 FORMAT(1X,'PONDED CELL LOCATIONS WILL BE SAVED FOR'
     &,' ALL TIME STEPS IN EACH STRESS PERIOD')
  576 FORMAT(1X,'NO SURFACE PONDING OUTPUT REQUESTED')
  580 FORMAT(1X,'PONDED CELL LOCATIONS WILL BE SAVED FOR',
     &1X,I5,1X,'USER-DEFINED POINTS DURING THE SIMULATION')	
  585 FORMAT(1X,I10,' ELEMENTS IN RX ARRAY ARE USED BY PND') 
  590 FORMAT(1X,I10,' ELEMENTS IN IR ARRAY ARE USED BY PND')  	
C
C1-----READ IPNDCB AND IPNDOC FROM INPUT FILE
      LLOC=1
      CALL URDCOM(IN,0,LINE)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPNDCB,R,IOUT,IN)
	IF (IPNDCB.NE.0) THEN
        LLOC=1
        CALL URDCOM(IN,0,LINE)
	  CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPNDOC,R,IOUT,IN)
	ENDIF
C		
C2---- IF NO OUTPUT REQUESTED, WRITE NOTICE TO OUTPUT FILE AND SET FLAGS
      IF (IPNDCB.EQ.ZERO) THEN
	  WRITE(IOUT,576)
	  LCPNDOC=1
	  IPNDOC=0 
	ENDIF
C
C3---- PRINT UNIT NUMBER.TO MODFLOW OUTPUT FILE
	IF (IPNDCB.NE.0) WRITE(IOUT,570) ABS(IPNDCB)
C
C4--- PRINT OUTPUT CONTROL SETTING TO MODFLOW OUTPUT FILE
	IF (IPNDOC.LE.ZERO) THEN
	    WRITE(IOUT,575)
		LCPNDOC=1
	ELSE
	    WRITE(IOUT,580)IPNDOC
	ENDIF
C
C Iac:  Allocate PND array
      allocate(PNDOC(IPNDOC,2))
C5---- ALLOCATE SPACE IN RX AND IR ARRAY FOR POND, PFLG, PONDOC ARRAYS
      IRK=ISUM
      IIRK=ISUMIR
      IF (IPNDCB.NE.ZERO.AND.IPNDOC.GT.ZERO) THEN
		LCPNDOC=ISUM
		ISUM=ISUM+IPNDOC*2
      ENDIF
      LCPOND=ISUM
      ISUM=ISUM+NCOL*NROW
C
C Iac: canceled out variable ICPFLG:	ICPFLG=ISUMIR
      ISUMIR=ISUMIR+NCOL*NROW

C6---- CALCULATE & PRINT AMOUNT OF SPACE USED BY PND PACKAGE.
      IRK=ISUM-IRK
      IIRK=ISUMIR-IIRK
      WRITE(IOUT,585)IRK
      WRITE(IOUT,590)IIRK
C
C Iac: -----SAVE POINTERS TO DATA AND RETURN.
      CALL SVSF2PND1PSV(IGRID)
C7---- RETURN.
      RETURN
      END
C
C=======================================================================
      SUBROUTINE VSF2PND1RP(IN,FNAME,IGRID)
C-----Version 2001
C-----VERSION 07MARCH2006
C     ******************************************************************
C     READ PNDOC ARRAY (IF NEEDED) AND SURFACE PONDING DEPTHS
C
C     OPEN OUTPUT FILE
C     ******************************************************************
C     SPECIFICATIONS: 
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,IBOUND,IOUT
      use VSFMODULE,    ONLY:ZERO  
      use PNDMODULE   
C     ------------------------------------------------------------------
	
      INTEGER IN 
      CHARACTER(len=24)  ANAME
      CHARACTER(len=200) FNAME,OFILE,LINE
      CHARACTER(len=4)   EXT
C
      ANAME = 'SURFACE PONDING DEPTH'
      EXT   = '.PDO'
C Iac: 
      CALL SVSF2PND1PNT(IGRID) 
C     ------------------------------------------------------------------
C1---- READ PNDOC ARRAY IF REQUESTED (IF PNDOC IS POSITIVE)
	IF (IPNDOC.GT.ZERO) THEN
		DO 1 K=1,IPNDOC
	      LLOC=1
            CALL URDCOM(IN,0,LINE)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,PNDOC(K,1),IOUT,IN)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,PNDOC(K,2),IOUT,IN)		  
    1     CONTINUE
      ENDIF 
C
C2---- READ SURFACE PONDING DEPTHS
	CALL U2DREL(POND(1,1),ANAME,NROW,NCOL,1,IN,IOUT)
C
C Iac. 
C3---- SET ALL ACTIVE CELL FLAGS EQUAL TO NO-PONDING (1)
	PFLG=0
	DO 81 I=1,NROW
	DO 81 J=1,NCOL
		IF (IBOUND(J,I,1).EQ.1) PFLG(J,I)=1
   81	ENDDO		
C
C4---- OPEN OUTPUT FILE
	IF (IPNDCB.NE.0) THEN
	  OFILE=FNAME
	  NUM=INDEX(OFILE,' ')
	  WRITE(OFILE(NUM-4:NUM),'(A4)')EXT
	  IF(IPNDCB.LT.0)THEN
		OPEN (UNIT=ABS(IPNDCB),FILE=OFILE,FORM='UNFORMATTED',
     &	STATUS='UNKNOWN')
	  ELSE
	    OPEN (UNIT=IPNDCB,FILE=OFILE,FORM='FORMATTED',
     &	STATUS='UNKNOWN')
	  ENDIF
	ENDIF
C
C5---- RETURN
	RETURN
	END	
C
C=======================================================================
      SUBROUTINE VSF2PND1RST(IRST,IGRID)
C-----Version 2011: 
C-----VERSION 07MARCH2006 RBT 
C     ******************************************************************
C	CHECK INFILTRATION PONDING BOUNDARY CONDITION AND RESTART TIME
C	STEP IF NECESSARY
C     ******************************************************************
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,BOTM,NBOTM,CC,CR,CV,
     &                       HNEW,HOLD,LAYCBD,LBOTM,IBOUND
      use GWFRCHMODULE, ONLY:RECH
      use VSFMODULE,    ONLY:ZERO,CVS
      use PNDMODULE   
C     ------------------------------------------------------------------
C
C
      REAL VCOND,CHCH1,CHCH2,CHCH3,CHCH4,CHCH6,NODE
      integer IRST
C Iac: 
      CALL SVSF2PND1PNT(IGRID) 
C     ------------------------------------------------------------------
	IRST=0   
C1---- CHECK SURFACE PONDING AND SWITCH BOUNDARY CONDITION 
	DO 5 I=1,NROW
	DO 5 J=1,NCOL
	  IF (IBOUND(J,I,1).EQ.0) CYCLE
	  RRATE=RECH(J,I)
	  BBOT=BOTM(J,I,LBOTM(1))
	  TTOP=BOTM(J,I,LBOTM(1)-1)
	  NODE=(TTOP+BBOT)/2		
	  HDPOND=TTOP+POND(J,I)
	  HD1=HNEW(J,I,1)
	  THCK1=(BOTM(J,I,LBOTM(1)-1)-BOTM(J,I,LBOTM(2)-1))/2.
	  THCK2=(BOTM(J,I,LBOTM(2)-1)-BOTM(J,I,LBOTM(3)-1))/2.
	  THCK=THCK1+THCK2
 	  VCOND=CVS(J,I,1)*THCK
	  CHCH1=ZERO
	  CHCH2=ZERO
	  CHCH3=ZERO
	  CHCH4=ZERO
	  CHCH6=ZERO

C2---- CHECK PONDED CELLS  (PFLG = 2) IF RECHARGE IS LESS THAN K
      IF (PFLG(J,I).EQ.2.AND.RRATE.LE.VCOND) THEN
C
C2A--- CHECK FOR FLOW INTO DOMAIN FROM PONDED CELLS 
        IF (J.NE.1.AND.IBOUND(J-1,I,1).GT.0) THEN
C2B--- CALCULATE FLOW THROUGH LEFT FACE
            HDIFF=HNEW(J,I,1)-HNEW(J-1,I,1)
            CHCH1=HDIFF*CR(J-1,I,1)
        ENDIF
        IF (J.NE.NCOL.AND.IBOUND(J+1,I,1).GT.0) THEN
C2C--- CALCULATE FLOW THROUGH RIGHT FACE.
            HDIFF=HNEW(J,I,1)-HNEW(J+1,I,1)
            CHCH2=HDIFF*CR(J,I,1)
        ENDIF
        IF (I.NE.1.AND.IBOUND(J,I-1,1).GT.0) THEN
C2D--- CHECK FLOW THROUGH BACK FACE
            HDIFF=HNEW(J,I,1)-HNEW(J,I-1,1)
            CHCH3=HDIFF*CC(J,I-1,1)
        ENDIF
        IF (I.NE.NROW.AND.IBOUND(J,I+1,1).GT.0) THEN
C2E--- CHECK FLOW THROUGH FRONT FACE
            HDIFF=HNEW(J,I,1)-HNEW(J,I+1,1)
            CHCH4=HDIFF*CC(J,I,1)
        ENDIF			
        IF (IBOUND(J,I,2).GT.0) THEN
C2G--- CHECK FLOW THROUGH LOWER FACE
            HDIFF=HNEW(J,I,1)-HNEW(J,I,2)
            CHCH6=HDIFF*CV(J,I,1)
        ENDIF
C2H--- SUM FLOWS THROUGH ALL FACES:RATE= total flow out 
        RATE = CHCH1+CHCH2+CHCH3+CHCH4+CHCH6

C2I--- SWITCH BOUNDARY CONDITION IF TOTAL FLOW OUT IS 1% GREATER THAN APPLIED FLUX
        IF (RATE.GT.1.01*RRATE) THEN
			IBOUND(J,I,1) = 1	   					
			PFLG(J,I)=1
			IRST = 1
        ENDIF
C
C3---- IF HEAD IN SURFACE CELL EXCEEDS PONDED ELEVATION OR SURFACE FLUX > K
C3---- SWITCH TO CONSTANT HEAD (PONDED) CELL	   
      ELSEIF (PFLG(J,I).EQ.1)THEN
        IF (RRATE.GT.VCOND.OR.HD1.GT.HDPOND) THEN
			IBOUND(J,I,1)=-1
			HNEW(J,I,1)  = HDPOND
C
C3A---- SWITCH FLAG TO PONDED (PFLG = 2)
			PFLG(J,I)=2
			IRST = 1
        ENDIF
      ENDIF
    5	ENDDO
C
C4---- INCREASE BOUNDARY CONDITION SWITCHING FLAG (IPOND)
      IF (IRST.EQ.1) IPOND = IPOND+1
C4A--- STOP SIMULATION IF > 30 SWITCHES OCCUR
      IF (IPOND.GT.30) THEN
        WRITE (*,101)
        STOP
      ENDIF
  101 FORMAT ('ERROR: EXCESSIVE OSCILLATION OF PONDING BOUNDARY',
     &' CONDITION.',/1X,'SIMULATION ABORTED.')
C
C5---- RETURN 
      RETURN
      END
C
C=======================================================================
      SUBROUTINE VSF2PND1BD(KSTP,KPER,IGRID)
C-----Version 2011
C-----VERSION 07MARCH2006 VSF1PND1BD RBT
C     ******************************************************************
C     COMPUTE FLOW OUT FROM PONDED CELLS
C     ******************************************************************
C
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,BOTM,NBOTM,CC,CR,CV,NPER,
     &                       HNEW,HOLD,LAYCBD,LBOTM,IBOUND,IOUT,BUFF
      use GWFBASMODULE, ONLY:DELT,PERTIM,TOTIM,ICHFLG,ICBCFL
      use VSFMODULE,    ONLY:ZERO
      use PNDMODULE
            
C     ------------------------------------------------------------------
      CHARACTER(len=16) TEXT, TEXT1
      DOUBLE PRECISION HD,CHIN,CHOUT,XX1,XX2,XX3,XX4,XX5,XX6
      INTEGER  KSTP,COUNT
      REAL     LIM

C
      DIMENSION UPER(IPNDOC), USTP(IPNDOC)
      PARAMETER (LIM=0.005)
C
      TEXT = 'PONDED FLOW     '
      
        
C Iac: 
      CALL SVSF2PND1PNT(IGRID) 

C     ------------------------------------------------------------------

C
C1-----CLEAR BUDGET ACCUMULATORS.      
      CHIN=ZERO
      CHOUT=ZERO
      IBDLBL=0
C
C2------CLEAR BUFFER.
      DO 5 K=1,NLAY
      DO 5 I=1,NROW
      DO 5 J=1,NCOL
      BUFF(J,I,K)=ZERO
5     CONTINUE
C
C3-----COUNT PONDED CELLS
	NCH=0
      DO 7 I=1,NROW
      DO 7 J=1,NCOL
       IF(IBOUND(J,I,1).LT.0.AND.PFLG(J,I).EQ.2) NCH=NCH+1
7     CONTINUE
	IF (ICBCFL.NE.0) THEN
	  WRITE(IOUT,8)NCH,KSTP,KPER
8	FORMAT(2X,I4,1X,'CELLS ARE PONDED DURING TIMESTEP',I4,
     &' PERIOD ',I4)
	ENDIF
C
C3A--- RETURN IF NO OUTPUT REQUESTED
	IF (IPNDCB.EQ.ZERO) RETURN
C
C4------LOOP THROUGH EACH CELL AND CALCULATE FLOW IN EACH
C4------PONDED CELL.
      DO 200 I=1,NROW
      DO 200 J=1,NCOL
	K=1
C
C5------IF CELL IS NOT PONDED SKIP IT & GO ON TO NEXT CELL.
      IF (IBOUND(J,I,K).NE.-1.OR.PFLG(J,I).NE.2)GO TO 200
C
C6------CLEAR VALUES FOR FLOW RATE THROUGH EACH FACE OF CELL.
      X1=ZERO
      X2=ZERO
      X3=ZERO
      X4=ZERO
      X5=ZERO
      X6=ZERO
      CHCH1=ZERO
      CHCH2=ZERO
      CHCH3=ZERO
      CHCH4=ZERO
      CHCH5=ZERO
      CHCH6=ZERO
C
C7------CALCULATE FLOW THROUGH THE LEFT FACE.
C7A-----IF THERE IS NO FLOW TO CALCULATE THROUGH THIS FACE, THEN GO ON
C7A-----TO NEXT FACE.  NO FLOW OCCURS AT THE EDGE OF THE GRID, TO AN
C7A-----ADJACENT NO-FLOW CELL, OR TO AN ADJACENT CONSTANT-HEAD CELL.
      IF(J.EQ.1) GO TO 30
      IF(IBOUND(J-1,I,K).EQ.0) GO TO 30
      IF(IBOUND(J-1,I,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 30
C
C7B-----CALCULATE FLOW THROUGH THIS FACE INTO THE ADJACENT CELL.
      HDIFF=HNEW(J,I,K)-HNEW(J-1,I,K)
      CHCH1=HDIFF*CR(J-1,I,K)
      IF(IBOUND(J-1,I,K).LT.0) GO TO 30
      X1=CHCH1
      XX1=X1
C
C7C-----ACCUMULATE POSITIVE AND NEGATIVE FLOW.
      IF (X1) 10,30,20
   10 CHOUT=CHOUT-XX1
      GO TO 30
   20 CHIN=CHIN+XX1
C
C8------CALCULATE FLOW THROUGH THE RIGHT FACE.
   30 IF(J.EQ.NCOL) GO TO 60
      IF(IBOUND(J+1,I,K).EQ.0) GO TO 60
      IF(IBOUND(J+1,I,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 60
      HDIFF=HNEW(J,I,K)-HNEW(J+1,I,K)
      CHCH2=HDIFF*CR(J,I,K)
      IF(IBOUND(J+1,I,K).LT.0) GO TO 60
      X2=CHCH2
      XX2=X2
      IF(X2)40,60,50
   40 CHOUT=CHOUT-XX2
      GO TO 60
   50 CHIN=CHIN+XX2
C
C9------CALCULATE FLOW THROUGH THE BACK FACE.
   60 IF(I.EQ.1) GO TO 90
      IF(IBOUND(J,I-1,K).EQ.0) GO TO 90
      IF(IBOUND(J,I-1,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 90
      HDIFF=HNEW(J,I,K)-HNEW(J,I-1,K)
      CHCH3=HDIFF*CC(J,I-1,K)
      IF(IBOUND(J,I-1,K).LT.0) GO TO 90
      X3=CHCH3
      XX3=X3
      IF(X3) 70,90,80
   70 CHOUT=CHOUT-XX3
      GO TO 90
   80 CHIN=CHIN+XX3
C
C10-----CALCULATE FLOW THROUGH THE FRONT FACE.
   90 IF(I.EQ.NROW) GO TO 150
      IF(IBOUND(J,I+1,K).EQ.0) GO TO 150
      IF(IBOUND(J,I+1,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 150
      HDIFF=HNEW(J,I,K)-HNEW(J,I+1,K)
      CHCH4=HDIFF*CC(J,I,K)
      IF(IBOUND(J,I+1,K).LT.0) GO TO 150
      X4=CHCH4
      XX4=X4
      IF (X4) 100,150,110
  100 CHOUT=CHOUT-XX4
      GO TO 150
  110 CHIN=CHIN+XX4
C
C11-----CALCULATE FLOW THROUGH THE LOWER FACE.
  150 IF(K.EQ.NLAY) GO TO 180
      IF(IBOUND(J,I,K+1).EQ.0) GO TO 180
      IF(IBOUND(J,I,K+1).LT.0 .AND. ICHFLG.EQ.0) GO TO 180
      HD=HNEW(J,I,K+1)
C
  152 HDIFF=HNEW(J,I,K)-HD
      CHCH6=HDIFF*CV(J,I,K)
      IF(IBOUND(J,I,K+1).LT.0) GO TO 180
      X6=CHCH6
      XX6=X6
      IF(X6) 160,180,170
  160 CHOUT=CHOUT-XX6
      GO TO 180
  170 CHIN=CHIN+XX6
C
C12-----SUM THE FLOWS THROUGH FIVE FACES OF CONSTANT HEAD CELL, AND
C12-----STORE SUM IN BUFFER.
 180  RATE=CHCH1+CHCH2+CHCH3+CHCH4+CHCH6
      BUFF(J,I,K)=RATE
C		
  200 CONTINUE
C
C13-----PRINT THE FLOW FOR THE CELL IF REQUESTED.
C		
C13A--- DETERMINE REQUESTED STRESS PERIOD AND ELAPSED TIME 
	IF (IPNDOC.GT.ZERO) THEN
	  DO 35 L=1,IPNDOC
		UPER(L)=INT(PNDOC(L,1))
		USTP(L)=PNDOC(L,2)
   35	  CONTINUE
   	  DO 42 L=1,IPNDOC
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
C14--- OUTPUT NOT REQUESTED FOR THIS TIMESTEP; RETURN
	  GOTO 901
	ENDIF
C
C15--- WRITE OUTPUT 

   52 IF(IPNDCB.GT.0) THEN
C15A-- ASCII FILE
	  IBDLBL=0
	  IF (NCH.EQ.0)THEN 
	    WRITE(IPNDCB,899) TEXT,KPER,KSTP,TOTIM,NCH
		GOTO 901
	  ENDIF
	  DO 902 I=1,NROW
	  DO 902 J=1,NCOL
	    IF (IBOUND(J,I,K).NE.-1.OR.PFLG(J,I).NE.2) CYCLE
	    IF(IBDLBL.EQ.0)  WRITE(IPNDCB,899) TEXT,KPER,KSTP,TOTIM,NCH
  899     FORMAT(A,' PERIOD',I4,'   STEP',I4,
     &     '  SIMULATION TIME',1X,G15.6,'  # OF PONDED CELLS:',I4)
          WRITE(IPNDCB,900) 1,I,J,BUFF(J,I,1)
  900     FORMAT(1X,'LAYER',I4,'   ROW',I4,'   COL',I4,
     1       '   FLOW',G15.6)
          IBDLBL=1
  902	  CONTINUE
	ELSEIF (IPNDCB.LT.0) THEN
C15B-- BINARY FILE
	  CALL UBDSV3(KSTP,KPER,TEXT,IPNDCB,BUFF,PFLG,1,
     1       NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
	ENDIF
C

C16--- RETURN.                            
  901 RETURN
      END
C
C=======================================================================
      SUBROUTINE VSF2PND1FL(IRST,KSTP,KITER,IGRID,Kkper,casc_flg)
C-----Version 2011 The subroutine depends also on the stress period: Kkper
C                  cascadingflow flag: casc_flg  
C-----VERSION 07MARCH2006 VSF1PND1FL RBT
C     ******************************************************************
C     RESET FLAGS FOR PONDED CELLS
C     ******************************************************************
C
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,BOTM,NBOTM,
     &                       HNEW,LAYCBD,LBOTM,IBOUND
      use GWFRCHMODULE, ONLY:RECH     
      use VSFMODULE,    ONLY:CVS
      use PNDMODULE
C -- Iac:      
      use cascading_module
      use GWFUZFMODULE, ONLY:EXCESPP
      implicit none
C -- !!         
C     ------------------------------------------------------------------
      INTEGER KSTP,KITER,IRST
      REAL    HDPOND,VCOND,THCK,THCK1,THCK2,RRATE,TTOP
C -- Iac:      
      real    IncomingFlux, OutgoingFlux
      integer casc_flg !! This is the unit number for the CASCADING input file
      integer Kkper, IGRID, I, J
C -- !!      

C Iac: 
      CALL SVSF2PND1PNT(IGRID) 
C     ------------------------------------------------------------------
C
C1--- ONLY RESET FLAGS AT BEGINING OF STRESS PERIOD
      IF (KSTP.NE.1.OR.KITER.NE.1.OR.IRST.NE.0) RETURN
C2--- RESET IPOND FLAG
	IPOND = 1
      DO 300 I=1,NROW	
      DO 300 J=1,NCOL
		RRATE=RECH(J,I)
		THCK1=(BOTM(J,I,LBOTM(1)-1)-BOTM(J,I,LBOTM(2)-1))/2.
		THCK2=(BOTM(J,I,LBOTM(2)-1)-BOTM(J,I,LBOTM(3)-1))/2.
		THCK=(THCK1+THCK2)/2
		VCOND=CVS(J,I,1)*THCK
		TTOP=BOTM(J,I,LBOTM(1)-1)
		HDPOND=TTOP+POND(J,I)
C IAC: 
C      ---> IF the new algorithm is activate, then call the subroutines
C      ---> cascade_in, computing the IncomingFlux
C      ---> cascade_out, computing the OutgoingFlux
C      ---> ELSE, IncomingFlux = 0  and  OutgoingFlux =0
C      ---> Add (IncomingFlux - OutgoingFlux) to rainfall                       
       if (casc_flg > 0 .AND. Kkper >= 2 ) then
          call CASCADE_IN(IncomingFlux,i, j, Kkper,0,Igrid)
          call CASCADE_OUT(OutgoingFlux,i, j, Kkper,0,Igrid)                             
       else 
          IncomingFlux = 0
          OutgoingFlux = 0
       end if  
       ! Update the Input Infiltration, adding/subtracting
       !        the cascading flow
       RRATE = RRATE + (IncomingFlux - OutgoingFlux)    
   
        IF (PFLG(J,I).GT.0.AND.RRATE.GT.VCOND) THEN
C2A--- IF RECHARGE RATE EXCEEDS SURFACE K, SWITCH TO PONDED	   
				IBOUND(J,I,1)=-1
				HNEW(J,I,1)=HDPOND
				PFLG(J,I)=2
C2Abis - Iac -- Compute Excesspp() if applied flux is greater than infiltration:
                if (casc_flg > 0 .AND. KKPER >=2) then
                  EXCESPP(J,I) = RRATE - VCOND
                end if 
C2B--- OTHERWISE, RESET ALL FLAGS
        ELSEIF (PFLG(J,I).GT.0) THEN 
                PFLG(J,I)=1
				IBOUND(J,I,1)=1
		ENDIF
300   END DO
C
C3-----RETURN
      RETURN
      END
