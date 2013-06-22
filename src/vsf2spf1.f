C Last change 2013  Iacopo
C    LAST CHANGE  07 MARCH 2006 RBT
C
C ================================
C Create the new module for SPF and a new SPFTYPE variable: SPFDAT  
      Module  SPFMODULE
        INTEGER, SAVE, POINTER ::     LCSEEP,LCSPFOC,ISPFOC,ISPFCB
        REAL,    SAVE, DIMENSION(:,:),   POINTER  :: SPFOC
        INTEGER, SAVE, DIMENSION(:,:,:), POINTER  :: SEEP
        
        Type SPFTYPE
            INTEGER, POINTER ::  LCSEEP,LCSPFOC,ISPFOC,ISPFCB
            REAL,    DIMENSION(:,:),   POINTER  :: SPFOC
            INTEGER, DIMENSION(:,:,:), POINTER  :: SEEP    
        End type SPFTYPE
        
        TYPE(SPFTYPE),SAVE :: SPFDAT(10)
      
      End module SPFMODULE
C================================
      Subroutine SVSF2SPF1PNT(IGRID)
C  Change SEV data to a different grid.
        use SPFMODULE
        
        ISPFCB=>SPFDAT(IGRID)%ISPFCB
        ISPFOC=>SPFDAT(IGRID)%ISPFOC
        LCSEEP=>SPFDAT(IGRID)%LCSEEP 
        LCSPFOC=>SPFDAT(IGRID)%LCSPFOC
        SPFOC=>SPFDAT(IGRID)%SPFOC
        SEEP=>SPFDAT(IGRID)%SEEP
        
        return
      end
C================================
      Subroutine SVSF2SPF1PSV(IGRID)
C  Save global data for a grid.
        use SPFMODULE
        
        SPFDAT(IGRID)%ISPFCB=>ISPFCB
        SPFDAT(IGRID)%ISPFOC=>ISPFOC
        SPFDAT(IGRID)%LCSEEP=>LCSEEP 
        SPFDAT(IGRID)%LCSPFOC=>LCSPFOC
        SPFDAT(IGRID)%SPFOC=>SPFOC
        SPFDAT(IGRID)%SEEP=>SEEP
        
        return
      end
C================================
      Subroutine SVSF2SPF1DA(IGRID)
C  Save global data for a grid.
        use SPFMODULE       
        deallocate(SPFDAT(IGRID)%ISPFCB)
        
        deallocate(SPFDAT(IGRID)%ISPFOC)
        deallocate(SPFDAT(IGRID)%LCSEEP)
        deallocate(SPFDAT(IGRID)%LCSPFOC)
        deallocate(SPFDAT(IGRID)%SPFOC)
        deallocate(SPFDAT(IGRID)%SEEP)
        
        return
      end
C================================

      
C=======================================================================
      SUBROUTINE VSF2SPF1AR(IN,IGRID)     
C-----Version 2011
C-----VERSION 07MARCH2006 RBT
C     ******************************************************************
C	ALLOCATE ARRAY SPACE FOR SEEPAGE FACE VARIABLES
C     ******************************************************************
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,IOUT
      use VSFMODULE,    ONLY:ISUM,ISUMIR,ZERO
      use GWFUZFMODULE, ONLY:EXCESPP
      use SPFMODULE 
C     ------------------------------------------------------------------
      INTEGER  IN
      CHARACTER(len=200) LINE
C Iac: 
      CALL SVSF2SPF1PNT(IGRID) 
C     ------------------------------------------------------------------
C ------Allocate SPF scalar variables, which makes it possible 
C ------for multiple grids to be defined.
      allocate(LCSEEP,LCSPFOC,ISPFOC,ISPFCB)
    
C-----      
  570 FORMAT(1X,'SEEPAGE FACE LOCATIONS WILL BE SAVED ON'
     &,' UNIT',I4)
  575 FORMAT(1X,'SEEPAGE FACE LOCATIONS WILL BE SAVED FOR'
     &,' ALL TIME STEPS IN EACH STRESS PERIOD')
  576 FORMAT(1X,'NO SEEPAGE FACE OUTPUT REQUESTED')
  580 FORMAT(1X,'SEEPAGE FACE LOCATIONS WILL BE SAVED FOR',
     &1X,I5,1X,'USER-DEFINED POINTS DURING THE SIMULATION')	
  585 FORMAT(1X,I10,' ELEMENTS IN RX ARRAY ARE USED BY SPF') 
  590 FORMAT(1X,I10,' ELEMENTS IN IR ARRAY ARE USED BY SPF') 
	  	
C
C1-----READ ISPFCB AND ISPFOC FROM INPUT FILE
      LLOC=1
      CALL URDCOM(IN,0,LINE)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ISPFCB,R,IOUT,IN)
	IF (ISPFCB.NE.0) THEN
        LLOC=1
        CALL URDCOM(IN,0,LINE)
	  CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ISPFOC,R,IOUT,IN)
	ENDIF
C		
C2---- IF NO OUTPUT REQUESTED, WRITE NOTICE TO OUTPUT FILE AND SET FLAGS
      IF (ISPFCB.EQ.ZERO) THEN
	  WRITE(IOUT,576)
	  LCSPFOC=1
	  ISPFOC=0 
	ENDIF
C
C3---- PRINT UNIT NUMBER.TO MODFLOW OUTPUT FILE
	IF (ISPFCB.NE.0) WRITE(IOUT,570) ABS(ISPFCB)
C
C4--- PRINT OUTPUT CONTROL SETTING TO MODFLOW OUTPUT FILE
	IF (ISPFOC.LE.ZERO) THEN
	    WRITE(IOUT,575)
		LCSPFOC=1
	ELSE
	    WRITE(IOUT,580)ISPFOC
	ENDIF
C
C Iac:  Allocate SPF arrays
      allocate(SPFOC(ISPFOC,2),SEEP(NCOL,NROW,NLAY))
      
C5---- ALLOCATE SPACE IN RX AND IR ARRAY FOR SEEP AND SPFOC ARRAYS
      IRK=ISUM
	IIRK=ISUMIR
	IF (ISPFCB.NE.ZERO.AND.ISPFOC.GT.ZERO) THEN
		LCSPFOC=ISUM
		ISUM=ISUM+ISPFOC*2
	ENDIF
	LCSEEP=ISUMIR
	ISUMIR=ISUMIR+NCOL*NROW*NLAY

C6---- CALCULATE & PRINT AMOUNT OF SPACE USED BY SPF PACKAGE.
      IRK=ISUM-IRK
	IIRK=ISUMIR-IIRK
      WRITE(IOUT,585)IRK
	WRITE(IOUT,590)IIRK
C
C Iac: -----SAVE POINTERS TO DATA AND RETURN.
      CALL SVSF2SPF1PSV(IGRID)
C7---- RETURN.
      RETURN
      END
C
C=======================================================================
      SUBROUTINE VSF2SPF1RP(IN,FNAME,IGRID)
C-----Version 2011
C-----VERSION 07MARCH2006 VSF1SPF1RP
C     ******************************************************************
C     READ SPFOC ARRAY (IF NEEDED) AND POTENTIAL SEEPAGE FACE LOCATIONS
C
C     OPEN OUTPUT FILE
C     ******************************************************************
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,IOUT,LAYCBD,LBOTM,IBOUND
      use VSFMODULE,    ONLY:ZERO
      use SPFMODULE 
C     ------------------------------------------------------------------
C
      REAL HD,PD
      INTEGER IN 
C
      CHARACTER(len=200) FNAME,OFILE,LINE
      CHARACTER(len=4)   EXT
      CHARACTER(len=24)  ANAME
C
      EXT   ='.SFO'
      ANAME = 'SEEPAGE LOCATIONS'
C Iac: 
      CALL SVSF2SPF1PNT(IGRID)       
C     ------------------------------------------------------------------
C1---- READ SPFOC ARRAY IF REQUESTED (IF ISPFOC IS POSITIVE)
	IF (ISPFOC.GT.ZERO) THEN
		DO 1 K=1,ISPFOC
	      LLOC=1
            CALL URDCOM(IN,0,LINE)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SPFOC(K,1),IOUT,IN)
	      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SPFOC(K,2),IOUT,IN)		  
    1     CONTINUE
      ENDIF 
C
C2---- READ POTENTIAL SEEPAGE FACE LOCATIONS
	DO 2 K=1,NLAY
	  CALL U2DINT(SEEP(1,1,K),ANAME,NROW,NCOL,1,IN,IOUT)
    2	CONTINUE
c
C3---- OPEN OUTPUT FILE
	IF (ISPFCB.NE.0) THEN
	  OFILE=FNAME
	  NUM=INDEX(OFILE,' ')
	  WRITE(OFILE(NUM-4:NUM),'(A4)')EXT
	  IF(ISPFCB.LT.0)THEN
		OPEN (UNIT=ABS(ISPFCB),FILE=OFILE,FORM='UNFORMATTED',
     &	STATUS='UNKNOWN')
	  ELSEIF (ISPFCB.GT.ZERO) THEN
	    OPEN (UNIT=ISPFCB,FILE=OFILE,FORM='FORMATTED',
     &	STATUS='UNKNOWN')
	  ENDIF
	ENDIF
C
C4---- RETURN
      RETURN
      END	
C
C=======================================================================
      SUBROUTINE VSF2SPF1FM(IGRID,Kkper,casc_flg)
C-----Version 2011 The subroutine depends also on the stress period: Kkper
C                  cascadingflow flag: casc_flg 
C-----VERSION 07MARCH2006 VSF1SPF1FM RBT 
C     ******************************************************************
C	  CHECK SEEPAGE FACE BOUNDARY CONDITION AND ADJUST POSITION
C     IF NECESSARY
C     ******************************************************************
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,IOUT,LAYCBD,LBOTM,IBOUND,
     &                       HNEW,BOTM,NBOTM,CR,CC,CV,DELC,DELR
      use VSFMODULE,    ONLY:ZERO
      use GWFUZFMODULE, ONLY:SEEPOUT
      use SPFMODULE
      use cascading_module
       

C     ------------------------------------------------------------------
      REAL NODE,CHCH1,CHCH2,CHCH3,CHCH4,CHCH6,RATE
C Iac: 
      integer casc_flg !! This is the unit number for the CASCADING input file
      integer Kkper, IGRID
      real    IncomingFlux, OutgoingFlux, NewRate
      
C Iac: 
      CALL SVSF2SPF1PNT(IGRID)       
C1---- CHECK SEEPAGE FACE POSITION AND ADJUST AS NECESSARY
	DO 80 K = 1,NLAY
	DO 80 I = 1,NROW
	DO 80 J = 1,NCOL
		IF (SEEP(J,I,K).LE.ZERO) CYCLE
		HD=HNEW(J,I,K)
	    BBOT=BOTM(J,I,LBOTM(K))
	    TTOP=BOTM(J,I,LBOTM(K)-1)
C --        
		CHCH1=ZERO
		CHCH2=ZERO
		CHCH3=ZERO
		CHCH4=ZERO
		CHCH5=ZERO
		CHCH6=ZERO
		IF (IBOUND(J,I,K).EQ.-1) THEN
C
C2A--- CHECK FOR FLOW INTO DOMAIN FROM SEEPAGE FACE CELLS
			IF (J.NE.1.AND.IBOUND(J-1,I,K).GT.0) THEN
C2B--- CALCULATE FLOW THROUGH LEFT FACE
				HDIFF=HNEW(J,I,K)-HNEW(J-1,I,K)
				CHCH1=HDIFF*CR(J-1,I,K)
			ENDIF
			IF (J.NE.NCOL.AND.IBOUND(J+1,I,K).GT.0) THEN
C2C--- CALCULATE FLOW THROUGH RIGHT FACE.
				HDIFF=HNEW(J,I,K)-HNEW(J+1,I,K)
				CHCH2=HDIFF*CR(J,I,K)
			ENDIF
			IF (I.NE.1.AND.IBOUND(J,I-1,K).GT.0) THEN
C2D--- CHECK FLOW THROUGH BACK FACE
				HDIFF=HNEW(J,I,K)-HNEW(J,I-1,K)
				CHCH3=HDIFF*CC(J,I-1,K)
			ENDIF
			IF (I.NE.NROW.AND.IBOUND(J,I+1,K).GT.0) THEN
C2E--- CHECK FLOW THROUGH FRONT FACE
				HDIFF=HNEW(J,I,K)-HNEW(J,I+1,K)
				CHCH4=HDIFF*CC(J,I,K)
			ENDIF			
			IF (K.NE.1.AND.IBOUND(J,I,K-1).GT.0) THEN
C2F--- CHECK FLOW THROUGH UPPER FACE
				HDIFF=HNEW(J,I,K)-HNEW(J,I,K-1)
				CHCH5=HDIFF*CV(J,I,K-1)
			ENDIF
			IF (K.NE.NLAY.AND.IBOUND(J,I,K+1).GT.0) THEN
C2G--- CHECK FLOW THROUGH LOWER FACE
				HDIFF=HNEW(J,I,K)-HNEW(J,I,K+1)
				CHCH6=HDIFF*CV(J,I,K)
			ENDIF
C2H--- SUM FLOWS THROUGH ALL SIX FACES
			RATE = CHCH1+CHCH2+CHCH3+CHCH4+CHCH5+CHCH6  

C
C--Iac 2012-
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
       !        the cascading flow to top cells
       if (K == 1) then 
            NewRate = RATE + (IncomingFlux - OutgoingFlux)
            RATE = NewRate
       end if
C--
C2I--- SWITCH BOUNDARY CONDITION IF TOTAL FLOW IS POSITIVE (FLOW OUT OF CELL IN THE DOMAIN)
       IF (RATE.GT.ZERO)  IBOUND(J,I,K) = 1
C2Ibis  Iac:
C           IF cascading flow is activated,  
C           if  there is a seepage flux (OUT from the domain), 
C           then store it to the variable SEEPOUT:          
            if (casc_flg > 0 .AND. K == 1 ) then
               if (RATE < 0) then            
                    SEEPOUT(J,I) = -RATE
               end if
            end if
		ELSE           
             
C2J--- CHECK FOR POSITIVE PRESSURES ABOVE SEEPAGE FACE		
			PD=NODE+0.01*(TTOP-BBOT)
            IF (HD.GT.PD) THEN
				IBOUND(J,I,K) = -1
				HNEW(J,I,K) = NODE
	        ENDIF
		ENDIF
  80	END DO
C
C3-----RETURN
	RETURN
	END
C
C-----------------------------------------------------------------------
      SUBROUTINE VSF2SPF1BD(KSTP,KPER,IGRID,casc_flg)
C-----Version 2013
C-----VERSION 07MARCH2006 VSF1SPF1BD RBT
C     ******************************************************************
C     COMPUTE FLOW OUT OF SEEPAGE FACE CELLS
C     ******************************************************************
C
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,HNEW,BOTM,NBOTM,CC,CR,CV,
     &                       LAYCBD,LBOTM,IOUT,BUFF,IBOUND,DELR,DELC
      use GWFBASMODULE, ONLY:DELT,PERTIM,TOTIM,ICHFLG,HDRY,
     &                       ICBCFL,VBVL,VBNM,MSUM
      use GWFLPFMODULE, ONLY:HK,VKA,VKCB,HANI,SC1,SC2,WETDRY,LAYFLG,
     &                       LAYWET,LAYVKA,CHANI,LAYAVG,LAYTYP,WETFCT,
     &                       IWETIT,IHDWET,ILPFCB 
      use GWFBCFMODULE, ONLY:LAYCON
      use VSFMODULE,    ONLY:ZERO
      use SPFMODULE
      use GWFUZFMODULE, ONLY:SEEPOUT
      use cascading_module
      
            
C-----------------------------------------------------------------------
      CHARACTER*16   TEXT
      DOUBLE PRECISION HD,CHIN,CHOUT,XX1,XX2,XX3,XX4,XX5,XX6
CIac
      integer casc_flg !! This is the unit number for the CASCADING input file
C -- Iac:      
      real    IncomingFlux, OutgoingFlux
      integer Kkper, IGRID, I, J
C
      INTEGER UPER
      REAL  USTP
C
      DIMENSION USTP(ISPFOC),UPER(ISPFOC)
C
      TEXT = 'SEEPAGE FACE   '
C Iac: 
      CALL SVSF2SPF1PNT(IGRID)       
C     ------------------------------------------------------------------
	  LIM=0.005
C
C1----- CLEAR BUDGET ACCUMULATORS.
      CHIN=ZERO
      CHOUT=ZERO
      IBDLBL=0
C
C2----- CLEAR BUFFER.
      DO 5 K=1,NLAY
      DO 5 I=1,NROW
      DO 5 J=1,NCOL
      BUFF(J,I,K)=ZERO
5     CONTINUE
C
C3---- COUNT SEEPAGE FACE CELLS
      NCH=0
      DO 7 K=1,NLAY
      DO 7 I=1,NROW
      DO 7 J=1,NCOL
        IF(IBOUND(J,I,K).LT.0.AND.SEEP(J,I,K).GT.ZERO) NCH=NCH+1
7     CONTINUE
      IF (ICBCFL.NE.0) THEN
         WRITE(IOUT,8)NCH,KSTP,KPER
8	    FORMAT(2X,I4,1X,'SEEPAGE FACE CELLS DURING TIMESTEP',I4,
     & ' PERIOD ',I4)
      ENDIF
C      
C4--- RETURN IF NO OUTPUT REQUESTED
	IF (ISPFCB.EQ.ZERO) RETURN
C
C5------LOOP THROUGH EACH CELL AND CALCULATE FLOW IN EACH
C5------SEEPAGE FACE CELL.
      DO 200 K=1,NLAY
C Iac 2012: cancel LAYCON(N), it seems not useful.
        !LC=LAYCON(K)
      DO 200 I=1,NROW
      DO 200 J=1,NCOL
C
C5A-----IF CELL IS NOT SEEPAGE FACE SKIP IT & GO ON TO NEXT CELL.
      IF (IBOUND(J,I,K).NE.-1.OR.SEEP(J,I,K).LE.ZERO)GO TO 200
C
C5B-----CLEAR VALUES FOR FLOW RATE THROUGH EACH FACE OF CELL.
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
C6------CALCULATE FLOW THROUGH THE LEFT FACE.
C6-----IF THERE IS NO FLOW TO CALCULATE THROUGH THIS FACE, THEN GO ON
C6-----TO NEXT FACE.  NO FLOW OCCURS AT THE EDGE OF THE GRID, TO AN
C6-----ADJACENT NO-FLOW CELL, OR TO AN ADJACENT CONSTANT-HEAD CELL.
      IF(J.EQ.1) GO TO 30
      IF(IBOUND(J-1,I,K).EQ.0) GO TO 30
      IF(IBOUND(J-1,I,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 30
C
C6A-----CALCULATE FLOW THROUGH THIS FACE INTO THE ADJACENT CELL.
      HDIFF=HNEW(J,I,K)-HNEW(J-1,I,K)
      CHCH1=HDIFF*CR(J-1,I,K)
      IF(IBOUND(J-1,I,K).LT.0) GO TO 30
      X1=CHCH1
      XX1=X1
C
C6B-----ACCUMULATE POSITIVE AND NEGATIVE FLOW.
      IF (X1) 10,30,20
   10 CHOUT=CHOUT-XX1
      GO TO 30
   20 CHIN=CHIN+XX1
C
C6C------CALCULATE FLOW THROUGH THE RIGHT FACE.
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
C6D------CALCULATE FLOW THROUGH THE BACK FACE.
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
C6E-----CALCULATE FLOW THROUGH THE FRONT FACE.
   90 IF(I.EQ.NROW) GO TO 120
      IF(IBOUND(J,I+1,K).EQ.0) GO TO 120
      IF(IBOUND(J,I+1,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 120
      HDIFF=HNEW(J,I,K)-HNEW(J,I+1,K)
      CHCH4=HDIFF*CC(J,I,K)
      IF(IBOUND(J,I+1,K).LT.0) GO TO 120
      X4=CHCH4
      XX4=X4
      IF (X4) 100,120,110
  100 CHOUT=CHOUT-XX4
      GO TO 120
  110 CHIN=CHIN+XX4
C
C6F-----CALCULATE FLOW THROUGH THE UPPER FACE.
  120 IF(K.EQ.1) GO TO 150
      IF(IBOUND(J,I,K-1).EQ.0) GO TO 150
      IF(IBOUND(J,I,K-1).LT.0 .AND. ICHFLG.EQ.0) GO TO 150
      HD=HNEW(J,I,K)
C
  122 HDIFF=HD-HNEW(J,I,K-1)
      CHCH5=HDIFF*CV(J,I,K-1)
      IF(IBOUND(J,I,K-1).LT.0) GO TO 150
      X5=CHCH5
      XX5=X5
      IF(X5) 130,150,140
  130 CHOUT=CHOUT-XX5
      GO TO 150
  140 CHIN=CHIN+XX5
C
C6G-----CALCULATE FLOW THROUGH THE LOWER FACE.
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
C6H-----SUM THE FLOWS THROUGH SIX FACES OF CONSTANT HEAD CELL, AND
C6H-----STORE SUM IN BUFFER.
 180  RATE=CHCH1+CHCH2+CHCH3+CHCH4+CHCH5+CHCH6
      BUFF(J,I,K)=RATE
    
C2Itris  Iac: 
C     IF cascading flow is activated,  
C     if  there is a seepage flux (OUT from the domain), 
C     then store it to the variable SEEPOUT:          
C      if (casc_flg > 0 .AND. RATE < 0) then
C         SEEPOUT(J,I) = -RATE
C      end if
      
C
  200 CONTINUE
C		
C7--- DETERMINE REQUESTED STRESS PERIOD AND ELAPSED TIME 
	IF (ISPFOC.GT.ZERO) THEN
	  DO 35 L=1,ISPFOC
		UPER(L)=INT(SPFOC(L,1))
		USTP(L)=SPFOC(L,2)
   35	  CONTINUE
   	  DO 42 L=1,ISPFOC
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
C7A--- OUTPUT NOT REQUESTED FOR THIS TIMESTEP; RETURN
	  GOTO 901
	ENDIF
	IF (ISPFOC.EQ.ZERO) GOTO 901
C
        
C8---- WRITE OUTPUT 

   52 IF(ISPFCB.GT.0) THEN
C8A--- ASCII FILE
	  IF (NCH.EQ.0)THEN 
	    WRITE(ISPFCB,899) TEXT,KPER,KSTP,TOTIM,NCH
		GOTO 901
	  ENDIF
	  IBDLBL=0
	  DO 902 K=1,NLAY
	  DO 902 J=1,NCOL
	  DO 902 I=1,NROW
	    IF (IBOUND(J,I,K).NE.-1.OR.SEEP(J,I,K).LE.ZERO)CYCLE
	    IF(IBDLBL.EQ.0) WRITE(ISPFCB,899) TEXT,KPER,KSTP,TOTIM,NCH
  899     FORMAT(A,' PERIOD',I4,'   STEP',I4,
     &     '  SIMULATION TIME',1x,G15.6,'  # OF SEEPAGE CELLS:',I4)
          WRITE(ISPFCB,900) K,I,J,BUFF(J,I,K)
  900     FORMAT(1X,'LAYER',I4,'   ROW',I4,'   COL',I4,
     1       '   FLOW',G15.6)
          IBDLBL=1
  902	  CONTINUE
	ELSEIF (ISPFCB.LT.0) THEN
C8B--- BINARY FILE
        CALL UBUDSV(KSTP,KPER,TEXT,
     1              ISPFCB,BUFF,NCOL,NROW,NLAY,IOUT)
	ENDIF
C
C9-----RETURN.
  901 RETURN
      END
C
C     ------------------------------------------------------------------
      SUBROUTINE VSF2SPF1FL(KSTP,IGRID)
C-----Version 2011
C-----VERSION 07MARCH2006 VSF1SPF1FL RBT
C     ******************************************************************
C     RESET FLAGS FOR PONDED CELLS
C     ******************************************************************
C
C     SPECIFICATIONS:
!! IAC
      use GLOBAL,       ONLY:NCOL,NROW,NLAY,IBOUND
      use VSFMODULE,    ONLY:ZERO
      use SPFMODULE
C     ------------------------------------------------------------------
      INTEGER KSTP
C Iac: 
      CALL SVSF2SPF1PNT(IGRID)       
C     ------------------------------------------------------------------
C 
C1---- SEEPAGE FACE CELLS
      DO 200 K=1,NLAY   
      DO 200 I=1,NROW
      DO 200 J=1,NCOL
C2----- IF BEGINING OF STRESS PERIOD, RESET FLAGS
        IF (KSTP.EQ.1.AND.SEEP(J,I,K).GT.ZERO) THEN
		 IBOUND(J,I,K)=1
        END IF
200   END DO
      
C3---- RETURN
      RETURN
      END
