C Last version: Iacopo Borsi, 2013
C     ******************************************************************
C     MAIN CODE FOR U.S. GEOLOGICAL SURVEY MODULAR MODEL -- MODFLOW-2005
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C1------USE package modules.
      USE GLOBAL
      USE GWFBASMODULE
      USE GWFHUFMODULE, ONLY:IOHUFHDS,IOHUFFLWS
      USE GWFEVTMODULE, ONLY:NEVTOP
      USE GWFRCHMODULE, ONLY:NRCHOP
      USE GWFLAKMODULE, ONLY:NLAKESAR,THETA,STGOLD,STGNEW,VOL
      USE GWFUZFMODULE, ONLY: IUZFBND, FINF, VKS
      USE PCGMODULE
c      USE LMGMODULE
      USE SIPMODULE
      USE DE4MODULE
      USE GMGMODULE
      USE LGRMODULE
C Iac
      use cascading_module
      use VSFMODULE
      use ATSMODULE,  ONLY: MXSTP, IMXSTP
      
      INCLUDE 'openspec.inc'
C
C-------ASSIGN VERSION NUMBER AND DATE
      CHARACTER*40 VERSION
      CHARACTER*10 MFVNAM
      PARAMETER (VERSION='SID&GRID Project  --  2013')
      PARAMETER (MFVNAM='-SIDGRID')
C
      CHARACTER*80 HEADNG(2)
      CHARACTER*200 FNAME
      INTEGER IBDT(8)

CATS-- ADAPTIVE TIME STEPPING VARIABLES ADDED 
C
       DOUBLE PRECISION PTIME      
C Iac:
CVSF -- adding the CUNIT relative to ATS & VSF packages (from 64 to 69)
C Iac: adding CUNIT = CFL for cascading process
      CHARACTER(len=4) CUNIT(NIUNIT)
      CUNIT(1:71)=(/ 'BCF6','WEL ','DRN ','RIV ','EVT ', '    ', 'GHB ',  !  7
     &           'RCH ', 'SIP ', 'DE4 ', '    ', 'OC  ', 'PCG ', 'lmg ',  ! 14
     &           'gwt ', 'FHB ', 'RES ', 'STR ', 'IBS ', 'CHD ', 'HFB6',  ! 21
     &           'LAK ', 'LPF ', 'DIS ', '    ', 'PVAL', '    ', 'HOB ',  ! 28
     &           '    ', '    ', 'ZONE', 'MULT', 'DROB', 'RVOB', 'GBOB',  ! 35
     &           'STOB', 'HUF2', 'CHOB', 'ETS ', 'DRT ', '    ', 'GMG ',  ! 42
     &           'HYD ', 'SFR ', '    ', 'GAGE', 'LVDA', 'BFH ', 'LMT6',  ! 49
     &           'MNW2', 'MNWI', 'MNW1', 'KDEP', 'SUB ', 'UZF ', 'gwm ',  ! 56
     &           'SWT ', 'cfp ', '    ', '    ', '    ', '    ', 'nrs ',  ! 63
     &           'REF ', 'ATS ', 'PND ', 'SPF ', 'SEV ', 'RZE ', 'CFL ',  ! 70
     &           '    '/)
C     ------------------------------------------------------------------
c ---
C2------WRITE BANNER TO SCREEN AND DEFINE CONSTANTS.
      WRITE (*,1) MFVNAM,VERSION
    1 FORMAT (/,34X,'MODFLOW',A,/,
     &4X,'U.S. GEOLOGICAL SURVEY MODULAR FINITE-DIFFERENCE',
     &' GROUND-WATER FLOW MODEL',/,29X,'Version ',A/)
      INUNIT = 99
      INUNIT = 99
      ILUNIT = 98
      NGRIDS = 1 
      ILGR = 0
C
C3------GET THE NAME OF THE NAME FILE
      CALL GETNAMFIL(ILGR,NGRIDS,FNAME,ILUNIT)
      DO IGRID = 1, NGRIDS
        MAXUNIT= INUNIT
C
C4A------IF USING LGR, READ NAMES FROM LGR NAME FILE
        IF(ILGR .NE. 0) CALL GETNAMFILLGR(ILUNIT,FNAME,IGRID)
C4B-----OPEN NAME FILE.
        OPEN (UNIT=INUNIT,FILE=FNAME,STATUS='OLD',ACTION=ACTION(1))
        NC=INDEX(FNAME,' ')
        WRITE(*,490)' Using NAME file: ',FNAME(1:NC)
  490   FORMAT(A,A)
C
C5------Get current date and time, assign to IBDT, and write to screen
        IF(IGRID .EQ. 1)THEN
          CALL DATE_AND_TIME(VALUES=IBDT)
          WRITE(*,2) (IBDT(I),I=1,3),(IBDT(I),I=5,7)
    2     FORMAT(1X,'Run start date and time (yyyy/mm/dd hh:mm:ss): ',
     &    I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2,/)
        ENDIF
C
C6------ALLOCATE AND READ (AR) PROCEDURE
        NSOL=1
        CALL GWF2BAS7AR(INUNIT,CUNIT,VERSION,24,31,32,MAXUNIT,IGRID,12,
     1                  HEADNG,26,MFVNAM)
        IF(IUNIT(50).GT.0 .AND. IUNIT(52).GT.0) THEN
          WRITE(IOUT,'(1X,/,1X,A)')
     1    'MNW1 and MNW2 cannot both be active in the same grid'
          CALL USTOP(' ')
        ENDIF
        IF(ILGR .NE.0) CALL GWF2LGR1AR(ILUNIT,FNAME,NGRIDS,IGRID)
        
CVSF-- CALL REF ALLOCATE SUBROUTINE, if active:
        IF (IUNIT(64).GT.0) CALL VSF2REF1AR(IUNIT(64),IGRID)        
C        
        IF(IUNIT(1).GT.0)  CALL GWF2BCF7AR(IUNIT(1),IGRID)
C -- VSF         
        if(IUNIT(23).GT.0.AND.IUNIT(64).LE.0) then
            CALL GWF2LPF7AR(IUNIT(23),IGRID)
        else if (IUNIT(23).GT.0.AND.IUNIT(64).GT.0) then
            CALL VSF2LPF7AR(IUNIT(23),IGRID)
        end if
C --    
CATS
C Iac:
C-- IF ATS IS ACTIVE, CALL ATS ALLOCATE SUBROUTINE
C-- ELSE, allocate ATS variable, and set them to zero, since they are not 
C         really necessary, but they are named anyway in the code
C         even if ATS is not activated
        if (IUNIT(65).GT.0) then
            CALL ATS2AR(FNAME,IUNIT(65),IGRID)
        else 
             allocate(MXSTP, IMXSTP)
             MXSTP=0
             IMXSTP=0
        end if
C --
        IF(IUNIT(37).GT.0) CALL GWF2HUF7AR(IUNIT(37),IUNIT(47),
     1                                     IUNIT(53),IGRID)
        IF(IUNIT(2).GT.0)  CALL GWF2WEL7AR(IUNIT(2),IGRID)
        IF(IUNIT(3).GT.0)  CALL GWF2DRN7AR(IUNIT(3),IGRID)
        IF(IUNIT(4).GT.0)  CALL GWF2RIV7AR(IUNIT(4),IGRID)
        IF(IUNIT(5).GT.0)  CALL GWF2EVT7AR(IUNIT(5),IGRID)
        IF(IUNIT(7).GT.0)  CALL GWF2GHB7AR(IUNIT(7),IGRID)
        IF(IUNIT(8).GT.0)  CALL GWF2RCH7AR(IUNIT(8),IGRID)
        IF(IUNIT(16).GT.0) CALL GWF2FHB7AR(IUNIT(16),IGRID)
        IF(IUNIT(17).GT.0) CALL GWF2RES7AR(IUNIT(17),IGRID)
        IF(IUNIT(18).GT.0) CALL GWF2STR7AR(IUNIT(18),IGRID)
        IF(IUNIT(19).GT.0) CALL GWF2IBS7AR(IUNIT(19),IUNIT(54),IGRID)
        IF(IUNIT(20).GT.0) CALL GWF2CHD7AR(IUNIT(20),IGRID)
        IF(IUNIT(21).GT.0) CALL GWF2HFB7AR(IUNIT(21),IGRID)
        IF(IUNIT(44).GT.0) CALL GWF2SFR7AR(IUNIT(44),IUNIT(1),IUNIT(23),
     1                     IUNIT(37),IUNIT(15),NSOL,IOUTS,IGRID,Iroute)
      
        IF(IUNIT(55).GT.0) CALL GWF2UZF1AR(IUNIT(55),IUNIT(1),
     1                                     IUNIT(23),IUNIT(37),IGRID)
        IF(IUNIT(22).GT.0 .OR. IUNIT(44).GT.0) CALL GWF2LAK7AR(
     1               IUNIT(22),IUNIT(44),IUNIT(15),IUNIT(55),NSOL,IGRID)
        IF(IUNIT(46).GT.0) CALL GWF2GAG7AR(IUNIT(46),IUNIT(44),
     1                                     IUNIT(22),IGRID)
        IF(IUNIT(39).GT.0) CALL GWF2ETS7AR(IUNIT(39),IGRID)
        IF(IUNIT(40).GT.0) CALL GWF2DRT7AR(IUNIT(40),IGRID)
        IF(IUNIT(54).GT.0) CALL GWF2SUB7AR(IUNIT(54),IGRID)
        IF(IUNIT(48).GT.0) CALL GWF2BFH1AR(IUNIT(48),ILGR,IGRID)
        IF(IUNIT(9).GT.0) CALL SIP7AR(IUNIT(9),MXITER,IGRID)
        IF(IUNIT(10).GT.0) CALL DE47AR(IUNIT(10),MXITER,IGRID)
        IF(IUNIT(13).GT.0) CALL PCG7AR(IUNIT(13),MXITER,IGRID)
c        IF(IUNIT(14).GT.0) CALL LMG7AR(IUNIT(14),MXITER,IGRID)
        IF(IUNIT(42).GT.0) CALL GMG7AR(IUNIT(42),MXITER,IGRID)
        IF(IUNIT(50).GT.0) CALL GWF2MNW27AR(IUNIT(50),IGRID)
        IF(IUNIT(51).GT.0) CALL GWF2MNW2I7AR(IUNIT(51),IUNIT(50),IGRID)
        IF(IUNIT(52).GT.0) CALL GWF2MNW17AR(IUNIT(52),IUNIT(9),
     1                       IUNIT(10),0,IUNIT(13),
     2                       0,IUNIT(42),FNAME,IGRID)
        IF(IUNIT(57).GT.0) CALL GWF2SWT7AR(IUNIT(57),IGRID)
        IF(IUNIT(43).GT.0) CALL GWF2HYD7BAS7AR(IUNIT(43),IGRID)
        IF(IUNIT(43).GT.0 .AND. IUNIT(19).GT.0)
     1                     CALL GWF2HYD7IBS7AR(IUNIT(43),IGRID)
        IF(IUNIT(43).GT.0 .AND. IUNIT(54).GT.0)
     1                     CALL GWF2HYD7SUB7AR(IUNIT(43),IGRID)
        IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0)
     1                     CALL GWF2HYD7STR7AR(IUNIT(43),IGRID)
        IF(IUNIT(43).GT.0 .AND. IUNIT(44).GT.0)
     1                     CALL GWF2HYD7SFR7AR(IUNIT(43),IGRID)
C Iac: allocate & read procedure for cascading flow
        if (IUNIT(70).GT.0) then 
        call CASCADEAR(IUNIT(70),IUNIT(55),IGRID,IUNIT(44),IUNIT(22))
        end if  
C
C VSF: allocate & read procedure for stress packages:    
CVSF-- CALL PND ALLOCATE SUBROUTINE 
        IF (IUNIT(66).GT.0) CALL VSF2PND1AR(IUNIT(66),IGRID) 
CVSF-- CALL SPF ALLOCATE SUBROUTINE  
        IF (IUNIT(67).GT.0) CALL VSF2SPF1AR(IUNIT(67),IGRID) 
CVSF-- CALL SEV ALLOCATE SUBROUTINE    
        IF (IUNIT(68).GT.0) CALL VSF2SEV1AR(IUNIT(68),IGRID)
CVSF-- CALL RZE ALLOCATE SUBROUTINE            
        IF (IUNIT(69).GT.0) CALL VSF2RZE1AR(IUNIT(69),IGRID)
C-------end of VSF allocate procedure
C
C  Observation allocate and read
        CALL OBS2BAS7AR(IUNIT(28),IGRID)
        IF(IUNIT(33).GT.0) CALL OBS2DRN7AR(IUNIT(33),IUNIT(3),IGRID)
        IF(IUNIT(34).GT.0) CALL OBS2RIV7AR(IUNIT(34),IUNIT(4),IGRID)
        IF(IUNIT(35).GT.0) CALL OBS2GHB7AR(IUNIT(35),IUNIT(7),IGRID)
        IF(IUNIT(36).GT.0) CALL OBS2STR7AR(IUNIT(36),IUNIT(18),IGRID)
        IF(IUNIT(38).GT.0) CALL OBS2CHD7AR(IUNIT(38),IGRID)

CVSF -- Read and prepare for VSF packages not depending on stress period
C    -- In particular:
C       READ VSF OUTPUT CONTROL ARRAY,SOIL CHARATERISTIC PARAMETERS (ALPHA AND VGN)
C       AND RESIDUAL SOIL WATER SATURATIONS (RSWC), AND OPEN OUTPUT FILE (IVSFCB)   
C-- CALL REF read & prepare
        IF (IUNIT(64).GT.0) CALL VSF2REF1RP(IUNIT(64),FNAME,IGRID)
C-- CALL PND read & prepare
        IF (IUNIT(66).GT.0) CALL VSF2PND1RP(IUNIT(66),FNAME,IGRID)
C-- CALL SPF read & prepare
        IF (IUNIT(67).GT.0) CALL VSF2SPF1RP(IUNIT(67),FNAME,IGRID)
C-- CALL SEV read & prepare not depending on stress period 
        IF (IUNIT(68).GT.0) CALL VSF2SEV1RPP(IUNIT(68),FNAME,IGRID)
C-- CALL REV read & prepare not depending on stress period
        IF (IUNIT(69).GT.0) CALL VSF2RZE1RPP(IUNIT(69),FNAME,IGRID)
  
C END LOOP FOR ALLOCATING AND READING DATA FOR EACH GRID
      ENDDO
      CLOSE(ILUNIT)
C

C7------SIMULATE EACH STRESS PERIOD.    
      DO 100 KPER = 1, NPER
        KKPER = KPER
        DO IGRID = 1, NGRIDS
          CALL GWF2BAS7ST(KKPER,IGRID)
          IF(IUNIT(19).GT.0) CALL GWF2IBS7ST(KKPER,IGRID)
          IF(IUNIT(54).GT.0) CALL GWF2SUB7ST(KKPER,IGRID)
          IF(IUNIT(57).GT.0) CALL GWF2SWT7ST(KKPER,IGRID)
C
C7B-----READ AND PREPARE INFORMATION FOR STRESS PERIOD.
C----------READ USING PACKAGE READ AND PREPARE MODULES.
          IF(ILGR .NE. 0) CALL GWF2LGR1RP(KKPER,IGRID)
          IF(IUNIT(2).GT.0) CALL GWF2WEL7RP(IUNIT(2),IGRID)
          IF(IUNIT(3).GT.0) CALL GWF2DRN7RP(IUNIT(3),IGRID)
          IF(IUNIT(4).GT.0) CALL GWF2RIV7RP(IUNIT(4),IGRID)
          IF(IUNIT(5).GT.0) CALL GWF2EVT7RP(IUNIT(5),IGRID)
          IF(IUNIT(7).GT.0) CALL GWF2GHB7RP(IUNIT(7),IGRID)
          IF(IUNIT(8).GT.0) CALL GWF2RCH7RP(IUNIT(8),IGRID)
          IF(IUNIT(17).GT.0) CALL GWF2RES7RP(IUNIT(17),IGRID)
          IF(IUNIT(18).GT.0) CALL GWF2STR7RP(IUNIT(18),IGRID)
          IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0)
     1                       CALL GWF2HYD7STR7RP(IUNIT(43),KKPER,IGRID)
          IF(IUNIT(20).GT.0) CALL GWF2CHD7RP(IUNIT(20),IGRID)
          IF(IUNIT(44).GT.0) CALL GWF2SFR7RP(IUNIT(44),IUNIT(15),
     1                                       IUNIT(22),KKPER,NSOL,
     2                                       IOUTS,IGRID)
          IF(IUNIT(43).GT.0 .AND. IUNIT(44).GT.0)
     1                       CALL GWF2HYD7SFR7RP(IUNIT(43),KKPER,IGRID)
C iac: adding cascading flow capability:     
          IF(IUNIT(55).GT.0) 
     1                  CALL GWF2UZF1RP(IUNIT(55),KKPER,IGRID,IUNIT(70))
C ---
C iac: call Recharge Update Subroutine if VSF, RCH and CFL are active:      
          IF(KKPER>0.AND.IUNIT(8)>0 .AND. IUNIT(64)>0 .AND. IUNIT(70)>0)
     1           call CASCADERFU(KKPER,IGRID)
C --                  
          IF(IUNIT(22).GT.0) CALL GWF2LAK7RP(IUNIT(22),IUNIT(1),
     1                 IUNIT(15),IUNIT(23),IUNIT(37),IUNIT(44),
     2                 IUNIT(55),KKPER,NSOL,IOUTS,IGRID)
          IF(IUNIT(46).GT.0.AND.KKPER.EQ.1) CALL GWF2GAG7RP(IUNIT(15),
     1               IUNIT(22),IUNIT(55),NSOL,IGRID)
          IF(IUNIT(39).GT.0) CALL GWF2ETS7RP(IUNIT(39),IGRID)
          IF(IUNIT(40).GT.0) CALL GWF2DRT7RP(IUNIT(40),IGRID)
          IF(IUNIT(50).GT.0) CALL GWF2MNW27RP(IUNIT(50),KKPER,IUNIT(9),
     1                       IUNIT(10),0,IUNIT(13),0,IUNIT(42),0,IGRID)
          IF(IUNIT(51).GT.0.AND.KKPER.EQ.1) CALL GWF2MNW2I7RP(IUNIT(51),
     1                       0,IGRID)
          IF(IUNIT(52).GT.0) CALL GWF2MNW17RP(IUNIT(52),IUNIT(1),
     1                              IUNIT(23),IUNIT(37),KKPER,IGRID)
          IF(IUNIT(48).GT.0) CALL GWF2BFH1RP(IUNIT(48),KKPER,IUNIT(1),
     1                                       IUNIT(23),IUNIT(37),IGRID) 
CVSF -- Read and prepare for VSF packages depending on stress period     
C-- CALL SEV read & prepare depending on stress period 
        IF (IUNIT(68).GT.0) CALL VSF2SEV1RP(IUNIT(68),IGRID)
C-- CALL REV read & prepare depending on stress period
        IF (IUNIT(69).GT.0) CALL VSF2RZE1RP(IUNIT(69),IGRID)
        
        ENDDO 
C
C Iac: adding ATS capability (if active)
C
C-ATS- RESET PERIOD COUNTERS
        PTIME=0.0
        KSTP=0
        KKSTP=0
        MXSTP=IMXSTP
C-ATS- IF ATS PACKAGE INACTIVE OR STEADY STATE, SET MXSTP = NSTP
        IF(IUNIT(65).LE.0.OR.ISSFLG(KPER).EQ.1) MXSTP=NSTP(KKPER)
C-ATS- SET INITIAL TIME STEP SIZE AND TIME STEP MULTIPLIER IF TRANSIENT
        IF (IUNIT(65).GT.0.AND.ISSFLG(KPER).EQ.0) THEN
			DLT=PERLEN(KKPER)/NSTP(KKPER)
            IF (TSMULT(KPER).NE.1.0) DELT=DLT
			TSMULT(KPER)=1.0
        ENDIF
        
C
C-ATS- TIME STEP DO LOOP CHANGED TO A DO WHILE LOOP FOR ADAPTIVE TIME STEPPING
C-ATS- SIMULATE EACH TIME STEP.
 
        DO 90 WHILE (PTIME.LT.PERLEN(KKPER) .AND. KSTP.LT.MXSTP)
C
C7C No ATS Version: DO 90 KSTP = 1, NSTP(KPER)
C Iac:

C-ATS- ADVANCE TIME COUNTERS (KSTP,PTIME)
         IF (KSTP.GT.0) THEN 
		    PTIME=PTIME+DELT*TSMULT(KPER)
         ELSE
			PTIME=PTIME+DELT
         ENDIF
         IF (IUNIT(65).GT.0.AND.ISSFLG(KPER).EQ.0) THEN 
            CALL ATS2ADV(PTIME,KSTP,DLT)
         ELSE
            KSTP = KSTP+1	  
         ENDIF
C--------	end of ATS advance     


          KKSTP = KSTP

   
          DO IGRID = 1, NGRIDS
C
C7C1----CALCULATE TIME STEP LENGTH. SET HOLD=HNEW.
            CALL GWF2BAS7AD(KKPER,KKSTP,IGRID)
            IF(IUNIT(20).GT.0) CALL GWF2CHD7AD(KKPER,IGRID)
            IF(IUNIT(1).GT.0) CALL GWF2BCF7AD(KKPER,IGRID)
            IF(IUNIT(17).GT.0) CALL GWF2RES7AD(KKSTP,KKPER,IGRID)
            IF(IUNIT(23).GT.0) CALL GWF2LPF7AD(KKPER,IGRID)
            IF(IUNIT(37).GT.0) CALL GWF2HUF7AD(KKPER,IGRID)
            IF(IUNIT(16).GT.0) CALL GWF2FHB7AD(IGRID)
            IF(IUNIT(22).GT.0) CALL GWF2LAK7AD(KKPER,KKSTP,IUNIT(15),
     1                                             IGRID)
            IF(IUNIT(50).GT.0) THEN
              IF (IUNIT(1).GT.0) THEN
                CALL GWF2MNW27BCF(KPER,IGRID)
              ELSE IF (IUNIT(23).GT.0) THEN
                CALL GWF2MNW27LPF(KPER,IGRID)
              ELSE IF(IUNIT(37).GT.0) THEN
                CALL GWF2MNW27HUF(KPER,IGRID)
              ELSE
                WRITE(IOUT,1000)
 1000           FORMAT(/1X,
     &        '***ERROR: MNW2 PACKAGE DOES NOT SUPPORT',/,
     &        ' SELECTED FLOW PACKAGE',/,
     &        ' (MNW2 DOES FULLY SUPPORT BCF, LPF, AND HUF PACKAGES)',/,
     &      ' -- STOP EXECUTION')
                CALL USTOP('MNW2 error-flow package')
              END IF
              CALL GWF2MNW27AD(KKSTP,KKPER,IGRID)
            END IF
            IF(IUNIT(52).GT.0) CALL GWF2MNW17AD(IUNIT(1),IUNIT(23),
     1                                         IUNIT(37),IGRID)
            IF(IUNIT(48).GT.0) CALL GWF2BFH1AD(IUNIT(48),IGRID)
C--Iac            
C--VSF-- CONSTANT HEAD FLAGS RESET FOR PONDED AND SEEPAGE FACE CELLS
            IF(IUNIT(64).GT.0.AND.IUNIT(67).GT.0) then
                CALL VSF2SPF1FL(KKSTP,IGRID)
            end if           

C
C---------INDICATE IN PRINTOUT THAT SOLUTION IS FOR HEADS
            CALL UMESPR('SOLVING FOR HEAD',' ',IOUT)
          ENDDO
C           
C7C1A-VSF-RESTARTING POINT FOR ADAPTIVE TIME STEPPING
   99 CONTINUE   
                      
          WRITE(*,25) KPER,KSTP
   25     FORMAT(' Solving:  Stress period: ',i5,4x,
     &       'Time step: ',i5,4x,'Ground-Water Flow Eqn.')  
          
C          
C---------BEGIN LOOP FOR ITERATING BETWEEN GRIDS (LGR ITERATIONS)
          LGRCNVG = 0            
          LGRITER = 0            
          ISHELL = 1
          DO WHILE (LGRCNVG .EQ. 0)
            LGRITER = LGRITER + 1

C7C2----ITERATIVELY FORMULATE AND SOLVE THE FLOW EQUATIONS FOR EACH GRID
            DO IGRID = 1, NGRIDS
              CALL SGWF2BAS7PNT(IGRID)
C-------------CHECK IF LGR IS ACTIVE
              IF(ILGR .NE. 0)THEN
                CALL SGWF2LGR1PNT(IGRID)
C---------------CHECK IF PARENT OR CHILD GRID 
                IF(ISCHILD .EQ. -1)THEN
                  ISHELL = 1
                ELSE
C-----------------PUT PARENT HEAD TO CHILD GRID SHARED NODES AND RELAX
                  CALL GWF2LGR1BH(KKPER,KKSTP,LGRITER,
     1                            GLOBALDAT(1)%HNEW,  
     2                            GLOBALDAT(1)%NCOL,
     3                            GLOBALDAT(1)%NROW,
     4                            GLOBALDAT(1)%NLAY) 
                  ISHELL = 3
                ENDIF          
C---------------ADJUST STORAGE IF LGR IS ACTIVE AND TR SIMULATION
                IF(ITRSS .NE. 0)THEN 
                  IF(ISCHILD .EQ. -1)THEN
                    DO LG =2,NGRIDS
                      IF(LGRDAT(IGRID)%IBPFLG(LG) .NE. 0)
     1                  CALL GWF2LGR1FMPBS(KKPER,KKSTP,LGRITER,IUNIT(1),
     2                                   IUNIT(23),IUNIT(37),LG,IGRID)
                    ENDDO
                  ELSE
                    CALL GWF2LGR1FMCBS(KKPER,KKSTP,LGRITER,IUNIT(1),
     1                                 IUNIT(23),IUNIT(37),IGRID)
                  ENDIF
                ENDIF
              ENDIF
C-------------BEGIN CAGE-SHELL INTERPOLATION LOOP
              IBSKIP2 = 0
              DO IBSOLV =1, ISHELL             
                IF(ILGR.NE.0 .AND. IGRID.GT.1) CALL GWF2LGR1FMIB(IBSOLV,
     1                                                          IBSKIP2)
                IF(IGRID.GT.1 .AND. IBSKIP2.EQ.1) CYCLE
C
C7C2----ITERATIVELY FORMULATE AND SOLVE THE FLOW EQUATIONS.
C- 
              
                DO 30 KITER = 1, MXITER
                  KKITER = KITER
C                 
C7C2A---FORMULATE THE FINITE DIFFERENCE EQUATIONS.
                  CALL GWF2BAS7FM(IGRID)
C VSF: change the CALL statement, choosing between Standard GWF or VSF:                  
                  !BFC - VSF case IS NOT AVAILABLE IN THIS RELEASE !!!
                  if(IUNIT(1).GT.0.0.AND.IUNIT(64).LE.0) then
                    call GWF2BCF7FM(KKITER,KKSTP,KKPER,IGRID)
                  end if
                  !LPF case:
                  if(IUNIT(23).GT.0.AND.IUNIT(64).LE.0) then
                    call GWF2LPF7FM(KKITER,KKSTP,KKPER,IGRID)
                  else if (IUNIT(23).GT.0.AND.IUNIT(64).GT.0) then
                    call VSF2LPF7FM(KKITER,KKSTP,KKPER,IGRID)
                  end if
C ----                           
                  IF(IUNIT(37).GT.0) CALL GWF2HUF7FM(KKITER,
     1                                   KKSTP,KKPER,IUNIT(47),IGRID)
CVSF-- ADJUST CONDUCTANCES FOR UNSATURATED CELLS 
CVSF-- AND CALCULATE UNSATURATED STORAGE
                  if (IUNIT(64).GT.0.AND.IUNIT(1).GT.0) then
                     call VSF2REF1FM(KKSTP,KKPER,KKITER,IGRID,IUNIT(70),
     &                               IUNIT(44),IUNIT(22))
	              else if (IUNIT(64).GT.0.AND.IUNIT(23).GT.0) then
                    call VSF2REF1FM(KKSTP,KKPER,KKITER,IGRID,IUNIT(70),
     &                               IUNIT(44),IUNIT(22))
                  end if
                  if (IUNIT(64).GT.0.AND.IUNIT(66).GT.0) 
     &		   call  VSF2PND1FL(IRST,KKSTP,KKITER,IGRID,KKPER,IUNIT(70))
C      
                  if (IUNIT(64).GT.0.AND.IUNIT(67).GT.0)                     
     &		         call  VSF2SPF1FM(IGRID,Kkper,IUNIT(70))

C-----------------ADJUST CONDUCTANCES IF LGR IS ACTIVE
                  IF(ILGR .NE. 0)THEN
                    IF(IGRID .EQ. 1)THEN  
                      DO LG =2,NGRIDS
                        IF(LGRDAT(IGRID)%IBPFLG(LG) .NE. 0)
     1                    CALL GWF2LGR1FMPBC(KKPER,KKSTP,KKITER,LGRITER,
     2                       IUNIT(1),LGRDAT(LG)%NPCBEG,
     3                       LGRDAT(LG)%NPRBEG,LGRDAT(LG)%NPLBEG,
     4                       LGRDAT(LG)%NPCEND,LGRDAT(LG)%NPREND,
     5                       LGRDAT(LG)%NPLEND,LGRDAT(LG)%IBOTFLG,
     6                       LGRDAT(LG)%PFLUX,LG)
                      ENDDO
                    ELSEIF(ISCHILD .GE. 0)THEN    
                      CALL GWF2LGR1FMCBC(KKPER,KKSTP,KKITER,LGRITER,
     1                                   IUNIT(1))
                    ENDIF
                  ENDIF
C VSF: add CALL statements for VSF stress packages                  
                  ! SEV package
                  if (IUNIT(68).GT.0) call VSF2SEV1FM(IGRID) 
                  ! RZE package
                  if (IUNIT(69).GT.0) call VSF2RZE1FM(IGRID)
C --                 
                  IF(IUNIT(21).GT.0) CALL GWF2HFB7FM(IGRID)
                  IF(IUNIT(2).GT.0)  CALL GWF2WEL7FM(IGRID)
                  IF(IUNIT(3).GT.0)  CALL GWF2DRN7FM(IGRID)
                  IF(IUNIT(4).GT.0)  CALL GWF2RIV7FM(IGRID)
                  IF(IUNIT(5).GT.0)  THEN
                    IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) CALL GWF2LAK7ST(
     1                                                          0,IGRID)
                    CALL GWF2EVT7FM(IGRID)
                    IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) CALL GWF2LAK7ST(
     1                                                          1,IGRID)
                  END IF
                  IF(IUNIT(7).GT.0) CALL GWF2GHB7FM(IGRID)
                  IF(IUNIT(8).GT.0) THEN
                     IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(
     1                                                          0,IGRID)
                     CALL GWF2RCH7FM(IGRID)
                     IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(
     1                                                          1,IGRID)
                  END IF
                  IF(IUNIT(16).GT.0) CALL GWF2FHB7FM(IGRID)
                  IF(IUNIT(17).GT.0) CALL GWF2RES7FM(IGRID)
                  IF(IUNIT(18).GT.0) CALL GWF2STR7FM(IGRID)
                  IF(IUNIT(19).GT.0) CALL GWF2IBS7FM(KKPER,IGRID)
                  IF(IUNIT(39).GT.0) CALL GWF2ETS7FM(IGRID)
                  IF(IUNIT(40).GT.0) CALL GWF2DRT7FM(IGRID)
                  IF(IUNIT(55).GT.0) CALL GWF2UZF1FM(KKPER,KKSTP,KKITER,
     1                                 IUNIT(44),IUNIT(22),IUNIT(58),
     2                                 IGRID)
                  IF(IUNIT(44).GT.0) CALL GWF2SFR7FM(KKITER,KKPER,KKSTP,
     1                                    IUNIT(22),ILGR,LGRITER,
     2                                    NGRIDS,IGRID)
                  IF(IUNIT(22).GT.0) CALL GWF2LAK7FM(KKITER,KKPER,KKSTP,
     1                                           IUNIT(44),IUNIT(55),
     2                                           IGRID)
                  IF(IUNIT(50).GT.0) THEN
                     IF (IUNIT(1).GT.0) THEN
                       CALL GWF2MNW27BCF(KPER,IGRID)
                    ELSE IF (IUNIT(23).GT.0) THEN
                       CALL GWF2MNW27LPF(KPER,IGRID)
                    ELSE IF(IUNIT(37).GT.0) THEN
                       CALL GWF2MNW27HUF(KPER,IGRID)
                    END IF
                    CALL GWF2MNW27FM(KKITER,kkstp,kkper,IGRID)
                  END IF
                  IF(IUNIT(52).GT.0) CALL GWF2MNW17FM(KKITER,IUNIT(1),
     1                                     IUNIT(23),IUNIT(37),IGRID)
                  IF(IUNIT(54).GT.0) CALL GWF2SUB7FM(KKPER,KKITER,
     1                                               IUNIT(9),IGRID)
                  IF(IUNIT(57).GT.0) CALL GWF2SWT7FM(KKPER,IGRID)
                  IF(IUNIT(48).GT.0) CALL GWF2BFH1FM(KKPER,KKSTP,KKITER,
     1                                              IUNIT(1),IGRID) 
C
C7C2B---MAKE ONE CUT AT AN APPROXIMATE SOLUTION.
                  IERR=0
                  IF (IUNIT(9).GT.0) THEN
                         CALL SIP7PNT(IGRID)
                         CALL SIP7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,EL,
     1                     FL,GL,V,W,HDCG,LRCH,NPARM,KKITER,HCLOSE,ACCL,
     2                     ICNVG,KKSTP,KKPER,IPCALC,IPRSIP,MXITER,
     3                     NSTP(KKPER),NCOL,NROW,NLAY,NODES,IOUT,0,IERR)
                  END IF
                  IF (IUNIT(10).GT.0) THEN
                         CALL DE47PNT(IGRID)
                         CALL DE47AP(HNEW,IBOUND,AU,AL,IUPPNT,IEQPNT,
     1                     D4B,MXUP,MXLOW,MXEQ,MXBW,CR,CC,CV,HCOF,RHS,
     2                     ACCLDE4,KITER,ITMX,MXITER,NITERDE4,HCLOSEDE4,
     3                     IPRD4,ICNVG,NCOL,NROW,NLAY,IOUT,LRCHDE4,
     4                     HDCGDE4,IFREQ,KKSTP,KKPER,DELT,NSTP(KKPER),
     5                     ID4DIR,ID4DIM,MUTD4,DELTL,NBWL,NUPL,NLOWL,
     6                     NLOW,NEQ,NUP,NBW,IERR)  
                  END IF
                  IF (IUNIT(13).GT.0) THEN
                         CALL PCG7PNT(IGRID)                         
                         CALL PCG7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,VPCG,
     1                     SS,P,CD,HCHG,LHCH,RCHG,LRCHPCG,KKITER,NITER,
     2                     HCLOSEPCG,RCLOSEPCG,ICNVG,KKSTP,KKPER,IPRPCG,
     3                     MXITER,ITER1,NPCOND,NBPOL,NSTP(KKPER),NCOL,
     4                     NROW,NLAY,NODES,RELAXPCG,IOUT,MUTPCG,IT1,
     5                     DAMPPCG,BUFF,HCSV,IERR,HPCG,DAMPPCGT,
     6                     ISSFLG(KKPER),HDRY)
                        
                  END IF
c                  IF (IUNIT(14).GT.0) THEN
c                    CALL LMG7PNT(IGRID)
c                    CALL LMG7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,A,IA,JA,
c     1                 U1,FRHS,IG,ISIZ1,ISIZ2,ISIZ3,ISIZ4,KKITER,BCLOSE,
c     2                 DAMPLMG,ICNVG,KKSTP,KKPER,MXITER,MXCYC,NCOL,NROW,
c     3                 NLAY,NODES,HNOFLO,IOUT,IOUTAMG,ICG,IADAMPLMG,
c     4                 DUPLMG,DLOWLMG) 
c                  END IF
                  IF (IUNIT(42).GT.0) THEN
                         CALL GMG7PNT(IGRID)
                         CALL GMG7AP(HNEW,RHS,CR,CC,CV,HCOF,HNOFLO,
     1                               IBOUND,IITER,MXITER,RCLOSEGMG,
     2                               HCLOSEGMG,KKITER,KKSTP,KKPER,NCOL,
     3                               NROW,NLAY,ICNVG,SITER,TSITER,
     4                               DAMPGMG,IADAMPGMG,IOUTGMG,IOUT,
     5                               GMGID,IUNITMHC,DUP,DLOW,CHGLIMIT,
     6                               BIGHEADCHG,HNEWLAST)
                  ENDIF                  
                  IF(IERR.EQ.1) CALL USTOP(' ')
C
C7C2C---IF CONVERGENCE CRITERION HAS BEEN MET STOP ITERATING.
                  IF (ICNVG.EQ.1) GOTO 33  
  30            CONTINUE
                
                KITER = MXITER
C
  33            CONTINUE
             
C7C2D-VSF-- NONCONVERGENCE? ADJUST TIMESTEP AND REITERATE IF ADAPTIVE TIME STEPPING ACTIVE
           if (IUNIT(65).GT.0.AND.ISSFLG(KPER).EQ.0) then               
              call ATS2RST(IRST,ICNVG,KKSTP,PTIME,IGRID)             
              if (IRST.EQ.1)  GOTO 99
           end if 
C
C7C2D-VSF-- IF INFILTRATION PONDING ACTIVE, CHECK BOUNDARY CONDITION AND ADJUST IF NECESSARY
           if (IUNIT(64).GT.0.AND.IUNIT(66).GT.0.AND.ICNVG.EQ.1) then 
              call VSF2PND1RST(IRST,IGRID) 
              if (IRST.EQ.1)  GOTO 99
           end if 
C
C-------------END CAGE-SHELL INTERPOLATION LOOOP
              ENDDO
C-------------PREPARE THE NEXT GRID FOR LGR ITERATION
              IF(ILGR.NE.0)THEN 
                IF(ISCHILD .EQ. -1)THEN
                  DO LG =2,NGRIDS
                    CALL GWF2LGR1INITP(KKPER,KKSTP,LGRITER,
     1                      LGRDAT(LG)%NPCBEG,LGRDAT(LG)%NPRBEG,
     2                      LGRDAT(LG)%NPLBEG,LGRDAT(LG)%NPCEND,
     3                      LGRDAT(LG)%NPREND,LGRDAT(LG)%NPLEND,
     4                      LGRDAT(LG)%IBOTFLG,LGRDAT(LG)%ISHFLG,
     5                      LGRDAT(LG)%MXLGRITER,LG,IGRID) 
                  ENDDO
                ELSEIF(IGRID.NE.1)THEN
C-----------------CALCULATE FLUX ENTERING THE CHILD INTERFACE AND RELAX
                  CALL GWF2LGR1FMBF(LGRITER,KKPER,KKSTP,IUNIT(1),
     1                              IUNIT(23),IUNIT(37)) 
                ENDIF
              ENDIF
C-----------END GRID LOOP
            ENDDO 
C-----------CHECK CONVEGENCE OF LGR IF LGR IS ACTIVE
            IF(ILGR .EQ. 0)THEN
              LGRCNVG = 1
            ELSE
              CALL GWF2LGR1CNVG(IGRID,NGRIDS,LGRCNVG,LGRITER,KKPER,
     1                          KKSTP)
            ENDIF
C---------END LGR ITERATION LOOP
          ENDDO       
C
C7C3----DETERMINE WHICH OUTPUT IS NEEDED FOR EACH GRID
          DO IGRID = 1, NGRIDS
            CALL GWF2BAS7OC(KKSTP,KKPER,ICNVG,
     1                      GLOBALDAT(IGRID)%IUNIT(12),IGRID)

C---------SWAP POINTERS FOR LGR DATA
            IF(ILGR.NE.0) CALL SGWF2LGR1PNT(IGRID)
C          
C7C4----CALCULATE BUDGET TERMS. SAVE CELL-BY-CELL FLOW TERMS.
            MSUM = 1
C             
            IF (IUNIT(1).GT.0.AND.IUNIT(64).LE.0) THEN
              CALL GWF2BCF7BDS(KKSTP,KKPER,IGRID)
              CALL GWF2BCF7BDCH(KKSTP,KKPER,IGRID)
              IBDRET=0
              IC1=1
              IC2=NCOL
              IR1=1
              IR2=NROW
              IL1=1
              IL2=NLAY
              DO 37 IDIR = 1, 3
                CALL GWF2BCF7BDADJ(KKSTP,KKPER,IDIR,IBDRET,
     1                            IC1,IC2,IR1,IR2,IL1,IL2,IGRID)
C   
   37         CONTINUE
            ENDIF
C VSF -- Changed the CALL statement                        
            IF(IUNIT(23).GT.0.0.AND.IUNIT(64).LE.0) THEN
              CALL GWF2LPF7BDS(KKSTP,KKPER,IGRID)
              CALL GWF2LPF7BDCH(KKSTP,KKPER,IGRID)
              IBDRET=0
              IC1=1
              IC2=NCOL
              IR1=1
              IR2=NROW
              IL1=1
              IL2=NLAY
              DO 157 IDIR=1,3
                CALL GWF2LPF7BDADJ(KKSTP,KKPER,IDIR,IBDRET,
     &                          IC1,IC2,IR1,IR2,IL1,IL2,IGRID)
  157         CONTINUE     
C VSF:
            else if (IUNIT(23).GT.0.AND.IUNIT(64).GT.0) then 
C7C4B-VSF- CALL STATEMENT ALTERED FOR VSF CAPABILITY              
              call VSF2LPF7BDS(KKSTP,KKPER,IGRID)                
              call VSF2LPF7BDCH(KKSTP,KKPER,IGRID)              

              IBDRET=0
              IC1=1
              IC2=NCOL
              IR1=1
              IR2=NROW
              IL1=1
              IL2=NLAY
              do 158 IDIR=1,3
              
C7C4B-VSF- CALL STATEMENT ALTERED FOR VSF CAPABILITY
              call  VSF2LPF7BDADJ(KKSTP,KKPER,IDIR,IBDRET,
     &                          IC1,IC2,IR1,IR2,IL1,IL2,IGRID)
     
  158           CONTINUE
            ENDIF
            IF(IUNIT(37).GT.0) THEN
              CALL GWF2HUF7BDS(KKSTP,KKPER,IGRID)
              CALL GWF2HUF7BDCH(KKSTP,KKPER,IUNIT(47),IGRID)
              IBDRET=0
              IC1=1
              IC2=NCOL
              IR1=1
              IR2=NROW
              IL1=1
              IL2=NLAY
              DO 159 IDIR=1,3
                CALL GWF2HUF7BDADJ(KKSTP,KKPER,IDIR,IBDRET,
     &                          IC1,IC2,IR1,IR2,IL1,IL2,IUNIT(47),IGRID)
159           CONTINUE
            ENDIF
            IF(IUNIT(2).GT.0) CALL GWF2WEL7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(3).GT.0) CALL GWF2DRN7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(4).GT.0) CALL GWF2RIV7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(5).GT.0) THEN
               IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) CALL GWF2LAK7ST(
     1                                                       0,IGRID)
               CALL GWF2EVT7BD(KKSTP,KKPER,IGRID)
               IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) CALL GWF2LAK7ST(
     1                                                       1,IGRID)
            END IF
C VSF Stress packages ---
C VSF---SOIL EVAPORATION
		    if(IUNIT(64).GT.0.AND.IUNIT(68).GT.0) then
               call VSF2SEV1BD(KKSTP,KKPER,IGRID)
            end if
C VSF EVAPOTRANSPIRATION
            if (IUNIT(64).GT.0.AND.IUNIT(69).GT.0) then
               call VSF2RZE1BD(KKSTP,KKPER,IGRID)               
            end if
C VSF--SEEPAGE FACE BUDGET 
		    if (IUNIT(64).GT.0.AND.IUNIT(67).GT.0) then 
               call VSF2SPF1BD(KKSTP,KKPER,IGRID,IUNIT(70))
            end if
C VSF--PONDING BUDGET 
		    if (IUNIT(64).GT.0.AND.IUNIT(66).GT.0) then 
     		   call VSF2PND1BD(KKSTP,KKPER,IGRID)
            end if 
C--- end of VSF Stress packages             
            IF(IUNIT(7).GT.0) CALL GWF2GHB7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(8).GT.0) THEN
               IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(
     1                                                       0,IGRID)
               CALL GWF2RCH7BD(KKSTP,KKPER,IGRID)
               IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(
     1                                                       1,IGRID)
            END IF
            IF(IUNIT(16).GT.0) CALL GWF2FHB7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(17).GT.0) CALL GWF2RES7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(18).GT.0) CALL GWF2STR7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(19).GT.0) CALL GWF2IBS7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(39).GT.0) CALL GWF2ETS7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(40).GT.0) CALL GWF2DRT7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(55).GT.0) CALL GWF2UZF1BD(KKSTP,KKPER,IUNIT(22),
     1                               IGRID)
            IF(IUNIT(44).GT.0) CALL GWF2SFR7BD(KKSTP,KKPER,IUNIT(15),
     1                         IUNIT(22),IUNIT(46),IUNIT(55),NSOL,ILGR,
     2                         LGRITER,NGRIDS,IGRID)

            IF(IUNIT(22).GT.0) CALL GWF2LAK7BD(KKSTP,KKPER,IUNIT(15),
     1                         IUNIT(46),IUNIT(44),IUNIT(55),NSOL,IGRID)
            IF(IUNIT(50).GT.0) CALL GWF2MNW27BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(52).GT.0) CALL GWF2MNW17BD(NSTP(KPER),KKSTP,KKPER,
     1                        IGRID)
            IF(IUNIT(54).GT.0) CALL GWF2SUB7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(57).GT.0) CALL GWF2SWT7BD(KKSTP,KKPER,IGRID)
            IF(IUNIT(48).GT.0) CALL GWF2BFH1BD(KKSTP,KKPER,IUNIT(1),
     1                                        IUNIT(23),IUNIT(37),IGRID)
CVSF-- VARIABLY SATURATED FLOW BUDGET
		    if (IUNIT(64).GT.0) then           
                call VSF2REF1BD(KKSTP,KKPER,IBDRET,IGRID)
            end if
C --
C Iac: IF  cascading flow is active AND Ponding and/or SeepageFlow of VSF
C          or UZF are active, then compute the total runoff budget and print out the results.      
      !if (IUNIT(70).GT.0.AND.KKSTEP==NSTP(KKPER)) then 
            if (IUNIT(70)>0) then
            if (IUNIT(55)>0 .OR. IUNIT(64)>0) then
                call CASCADEBD(KSTP,KPER,IGRID)
            end if 
            end if
C - end Iac
            IF(ILGR .NE.0 .AND. ISCHILD .LE. 0) CALL GWF2LGR1PBD(KKSTP,
     1                                        KKPER,NGRIDS)
     
CLMT
CLMT----CALL LINK-MT3DMS SUBROUTINES TO SAVE FLOW-TRANSPORT LINK FILE
CLMT----FOR USE BY MT3DMS FOR TRANSPORT SIMULATION
CLMT
            INCLUDE 'lmt7.inc'
CLMT                              
C
C  Observation and hydrograph simulated equivalents
            CALL OBS2BAS7SE(IUNIT(28),IGRID)
            IF(IUNIT(33).GT.0) CALL OBS2DRN7SE(IGRID)
            IF(IUNIT(34).GT.0) CALL OBS2RIV7SE(IGRID)
            IF(IUNIT(35).GT.0) CALL OBS2GHB7SE(IGRID)
            IF(IUNIT(36).GT.0) CALL OBS2STR7SE(IGRID)
            IF(IUNIT(38).GT.0) CALL OBS2CHD7SE(KKPER,IGRID)
            IF(IUNIT(43).GT.0) CALL GWF2HYD7BAS7SE(1,IGRID)
            IF(IUNIT(43).GT.0 .AND. IUNIT(19).GT.0)
     1                                CALL GWF2HYD7IBS7SE(1,IGRID)
            IF(IUNIT(43).GT.0 .AND. IUNIT(54).GT.0)
     1                                CALL GWF2HYD7SUB7SE(1,IGRID)
            IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0)
     1                                CALL GWF2HYD7STR7SE(1,IGRID)
            IF(IUNIT(43).GT.0 .AND. IUNIT(44).GT.0)
     1                                CALL GWF2HYD7SFR7SE(1,IGRID)
C
C7C5---PRINT AND/OR SAVE DATA.
            CALL GWF2BAS7OT(KKSTP,KKPER,ICNVG,1,IGRID)
            IF(IUNIT(19).GT.0) CALL GWF2IBS7OT(KKSTP,KKPER,IUNIT(19),
     1                                         IGRID)
            IF(IUNIT(37).GT.0)THEN
              IF(IOHUFHDS .NE.0 .OR.IOHUFFLWS .NE.0)
     1           CALL GWF2HUF7OT(KKSTP,KKPER,ICNVG,1,IGRID)
            ENDIF
            IF(IUNIT(51).NE.0) CALL GWF2MNW2I7OT(NSTP(KKPER),KKSTP,
     1                         KKPER,IGRID)
            IF(IUNIT(54).GT.0) CALL GWF2SUB7OT(KKSTP,KKPER,IUNIT(54),
     1                                         IGRID)
            IF(IUNIT(57).GT.0) CALL GWF2SWT7OT(KKSTP,KKPER,IGRID)
            IF(IUNIT(43).GT.0) CALL GWF2HYD7BAS7OT(KKSTP,KKPER,IGRID)
C
C------PRINT FLUXES FROM CHILD SPECIFIED HEAD BOUNDARY CONDITIONS 
            IF(ILGR.NE.0 .AND. ISCHILD .GE. 0) CALL GWF2LGR1CBD(KKSTP,
     1                                            KKPER,IGRID,   
     2                                            GLOBALDAT(1)%NCOL,
     3                                            GLOBALDAT(1)%NROW,
     4                                            GLOBALDAT(1)%NLAY, 
     5                                            GLOBALDAT(1)%IBOUND) 
C
C------CHECK FOR CHANGES IN HEAD AND FLUX BOUNDARY CONDITIONS 
            IF(IUNIT(48).GT.0) CALL GWF2BFH1OT(KKSTP,KKPER,IGRID)
C
C7C6---JUMP TO END OF PROGRAM IF CONVERGENCE WAS NOT ACHIEVED.
            IF(ICNVG.EQ.0 .AND. ILGR .EQ. 0) GO TO 110
C---------END GRID OT GRID LOOP
          ENDDO
C
C
CVSF-- ADJUST TIME STEP BASED ON CONVERGENCE IF ADAPTIVE TIME STEPPING IS ACTIVE
        IF (IUNIT(65).GT.0.AND.ISSFLG(KPER).EQ.0) THEN
           
           CALL ATS2ADJ(KKITER,NSTP(KKPER),IUNIT(10),ITMX,
     &                  PERLEN(KKPER),ENDPER,KPER,KKSTP,RST,IGRID) 
CVSF-- IF STRESS PERIOD ENDED AND TIMESTEP WAS ADJUSTED,  
CVSF-- RESTART ITERATION LOOP
           IF (RST.GT.0) GOTO 99

CVSF-- IF STRESS PERIOD ENDED, EXIT TIMESTEP LOOP
           IF (ENDPER.GT.0) EXIT
        END IF 
C
C-----END OF TIME STEP (KSTP) LOOP:   
   90      CONTINUE 
C
C-ATS- IF MAXIMUM TIME STEPS REACHED, CHECK PERIOD TIME AND STOP IF ERROR
      IF (IUNIT(65).GT.0.AND.ISSFLG(KPER).EQ.0) THEN
      
      IF (KSTP.GE.MXSTP.AND.ABS(PTIME/PERLEN(KKPER)-1.0).GT.1.0E-4)THEN
         WRITE (*,101)
         STOP
      END IF
  101 FORMAT(2X,'MAXIMUM TIME STEPS SPECIFIED BY ATS MODULE EXCEEDED.'
     &,/,2X,'SIMULATION ABORTED.')
      END IF

C-ATS
C-----END OF STRESS PERIOD (KPER) LOOP:    
  100 CONTINUE
C
C8------END OF SIMULATION
C
  110 DO IGRID = 1, NGRIDS
        CALL SGWF2BAS7PNT(IGRID)
        IF(IUNIT(52).NE.0 .AND. ICNVG.NE.0)
     1     CALL GWF2MNW17OT(IGRID)
C
C-------SAVE RESTART RECORDS FOR SUB PACKAGE
         IF(IUNIT(54).GT.0) CALL GWF2SUB7SV(IGRID)
C
C  Observation output
        IF(IUNIT(28).GT.0) CALL OBS2BAS7OT(IUNIT(28),IGRID)
        IF(IUNIT(33).GT.0) CALL OBS2DRN7OT(IGRID)
        IF(IUNIT(34).GT.0) CALL OBS2RIV7OT(IGRID)
        IF(IUNIT(35).GT.0) CALL OBS2GHB7OT(IGRID)
        IF(IUNIT(36).GT.0) CALL OBS2STR7OT(IGRID)
        IF(IUNIT(38).GT.0) CALL OBS2CHD7OT(IGRID)
      ENDDO
C
        CALL GLO1BAS6ET(IOUT,IBDT,1)
C
      DO IGRID = 1, NGRIDS
C9------CLOSE FILES AND DEALLOCATE MEMORY.  GWF2BAS7DA MUST BE CALLED
C9------LAST BECAUSE IT DEALLOCATES IUNIT.
        CALL SGWF2BAS7PNT(IGRID)
C ATS -- Deallocate ATS:
C Necessary?
        if(IUNIT(65).GT.0) call SATS2DA(IGRID)
C --
        IF(IUNIT(1).GT.0) CALL GWF2BCF7DA(IGRID)
        IF(IUNIT(2).GT.0) CALL GWF2WEL7DA(IGRID)
        IF(IUNIT(3).GT.0) CALL GWF2DRN7DA(IGRID)
        IF(IUNIT(4).GT.0) CALL GWF2RIV7DA(IGRID)
        IF(IUNIT(5).GT.0) CALL GWF2EVT7DA(IGRID)
        IF(IUNIT(7).GT.0) CALL GWF2GHB7DA(IGRID)
        IF(IUNIT(8).GT.0) CALL GWF2RCH7DA(IGRID)
        IF(IUNIT(9).GT.0) CALL SIP7DA(IGRID)
        IF(IUNIT(10).GT.0) CALL DE47DA(IGRID)
        IF(IUNIT(13).GT.0) CALL PCG7DA(IGRID)
c        IF(IUNIT(14).GT.0) CALL LMG7DA(IGRID)
        IF(IUNIT(16).GT.0) CALL GWF2FHB7DA(IGRID)
        IF(IUNIT(17).GT.0) CALL GWF2RES7DA(IGRID)
        IF(IUNIT(18).GT.0) CALL GWF2STR7DA(IGRID)
        IF(IUNIT(19).GT.0) CALL GWF2IBS7DA(IGRID)
        IF(IUNIT(20).GT.0) CALL GWF2CHD7DA(IGRID)
        IF(IUNIT(21).GT.0) CALL GWF2HFB7DA(IGRID)
        IF(IUNIT(22).GT.0 .OR. IUNIT(44).GT.0)CALL GWF2LAK7DA(IUNIT(22),
     1                                              IGRID)
        IF(IUNIT(23).GT.0) CALL GWF2LPF7DA(IGRID)
        IF(IUNIT(37).GT.0) CALL GWF2HUF7DA(IGRID)
        IF(IUNIT(39).GT.0) CALL GWF2ETS7DA(IGRID)
        IF(IUNIT(40).GT.0) CALL GWF2DRT7DA(IGRID)
        IF(IUNIT(42).GT.0) CALL GMG7DA(IGRID)
        IF(IUNIT(44).GT.0) CALL GWF2SFR7DA(IGRID)
        IF(IUNIT(46).GT.0) CALL GWF2GAG7DA(IGRID)
        IF(IUNIT(48).GT.0) CALL GWF2BFH1DA(IGRID)
        IF(IUNIT(50).GT.0) CALL GWF2MNW27DA(IGRID)
        IF(IUNIT(51).GT.0) CALL GWF2MNW2I7DA(IGRID)
        IF(IUNIT(52).GT.0) CALL GWF2MNW17DA(IGRID)
        IF(IUNIT(54).GT.0) CALL GWF2SUB7DA(IGRID)
        IF(IUNIT(55).GT.0) CALL GWF2UZF1DA(IGRID)
        IF(IUNIT(57).GT.0) CALL GWF2SWT7DA(IGRID)
C --Iac: Deallocate VSF Data:              
        if(IUNIT(64).GT.0) call SVSF2VSF1DA(IGRID)
        
        if(IUNIT(66).GT.0) call SVSF2PND1DA(IGRID)
        if(IUNIT(67).GT.0) call SVSF2SPF1DA(IGRID)
        if(IUNIT(68).GT.0) call SVSF2SEV1DA(IGRID)
        if(IUNIT(69).GT.0) call SVSF2RZE1DA(IGRID)
C --Iac: Deallocate CASCADE Data:
        if(IUNIT(70).GT.0) call CASCADEDA(IUNIT(55),IGRID)
C ---         
        CALL OBS2BAS7DA(IUNIT(28),IGRID)
        IF(IUNIT(33).GT.0) CALL OBS2DRN7DA(IGRID)
        IF(IUNIT(34).GT.0) CALL OBS2RIV7DA(IGRID)
        IF(IUNIT(35).GT.0) CALL OBS2GHB7DA(IGRID)
        IF(IUNIT(36).GT.0) CALL OBS2STR7DA(IGRID)
        IF(IUNIT(38).GT.0) CALL OBS2CHD7DA(IGRID)
        IF(IUNIT(43).GT.0) CALL GWF2HYD7DA(IGRID)
        CALL GWF2BAS7DA(IGRID)
      ENDDO
C
C10-----END OF PROGRAM.
      IF(ICNVG.EQ.0) THEN
        WRITE(*,*) 'FAILED TO MEET SOLVER CONVERGENCE CRITERIA'
      ELSE
        WRITE(*,*) ' Normal termination of simulation'
      END IF
C
      CALL USTOP(' ')
C
      END
      SUBROUTINE GETNAMFIL(ILGR,NGRIDS,FNAME,ILUNIT)
C     ******************************************************************
C     GET THE NAME OF THE NAME FILE
C     ******************************************************************
C        SPECIFICATIONS:
C
C     ------------------------------------------------------------------
      CHARACTER*(*) FNAME
      CHARACTER*200 COMLIN,LINE
      LOGICAL EXISTS
      INCLUDE 'openspec.inc'
C     ------------------------------------------------------------------
C
C Get name file from command line or user interaction.
        FNAME=' '
        COMLIN=' '
C *** Subroutines GETARG and GETCL are extensions to Fortran 90/95 that
C *** allow a program to retrieve command-line arguments.  To enable
C *** Modflow-2000 to read the name of a Name file from the command
C *** line, either GETARG or GETCL must be called, but not both.  As
C *** distributed, the call to GETARG is uncommented.  For compilers
C *** that support GETCL but not GETARG, comment out the call to GETARG
C *** and uncomment the call to GETCL.  The calls to both GETARG and
C *** GETCL may be commented out for compilers that do not support
C *** either extension.
        CALL GETARG(1,COMLIN)
C        CALL GETCL(COMLIN)
        ICOL = 1
        IF(COMLIN.NE.' ') THEN
          FNAME=COMLIN
        ELSE
   15   WRITE (*,*) ' Enter the name of the NAME FILE or LGR CONTROL ',
     &              'FILE:'
          READ (*,'(A)') FNAME
          CALL URWORD(FNAME,ICOL,ISTART,ISTOP,0,N,R,0,0)
          FNAME=FNAME(ISTART:ISTOP)
          IF (FNAME.EQ.' ') GOTO 15
        ENDIF
        INQUIRE (FILE=FNAME,EXIST=EXISTS)
        IF(.NOT.EXISTS) THEN
          NC=INDEX(FNAME,' ')
          FNAME(NC:NC+3)='.nam'
          INQUIRE (FILE=FNAME,EXIST=EXISTS)
          IF(.NOT.EXISTS) THEN
            WRITE (*,480) FNAME(1:NC-1),FNAME(1:NC+3)
  480       FORMAT(1X,'Can''t find name file ',A,' or ',A)
            CALL USTOP(' ')
          ENDIF
        ENDIF
C1A-----CHECK FOR LGR KEYWORD.  IF AN LGR SIMULATION, THEN READ THE 
C1A-----NUMBER OF GRIDS AND LEAVE FILE OPEN FOR PARSING.
C1A-----IF NOT LGR, THEN CLOSE FILE AND CONTINUE AS NORMAL.
      OPEN(UNIT=ILUNIT,FILE=FNAME,STATUS='OLD',ACTION=ACTION(1))
      CALL URDCOM(ILUNIT,0,LINE)
      ICOL=1
      CALL URWORD(LINE,ICOL,ISTART,ISTOP,1,N,R,0,ILUNIT)
      IF(LINE(ISTART:ISTOP) .EQ. 'LGR') THEN
        ILGR = 1
        WRITE(*,*) ' RUNNING MODFLOW WITH LGR '
        CALL URDCOM(ILUNIT,0,LINE)
        ICOL=1
        CALL URWORD(LINE,ICOL,ISTART,ISTOP,2,NGRIDS,R,0,ILUNIT)
        WRITE(*,*) 'NGRIDS = ', NGRIDS 
      ELSE
        CLOSE(ILUNIT)
      ENDIF
C
      RETURN
      END
      SUBROUTINE GLO1BAS6ET(IOUT,IBDT,IPRTIM)
C     ******************************************************************
C     Get end time and calculate elapsed time
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      INTEGER IBDT(8), IEDT(8), IDPM(12)
      DATA IDPM/31,28,31,30,31,30,31,31,30,31,30,31/ ! Days per month
      DATA NSPD/86400/  ! Seconds per day
C     ------------------------------------------------------------------
C
C     Get current date and time, assign to IEDT, and write.
      CALL DATE_AND_TIME(VALUES=IEDT)
      WRITE(*,1000) (IEDT(I),I=1,3),(IEDT(I),I=5,7)
 1000 FORMAT(1X,'Run end date and time (yyyy/mm/dd hh:mm:ss): ',
     &I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2)
      IF(IPRTIM.GT.0) THEN
        WRITE(IOUT,'(1X)')
        WRITE(IOUT,1000) (IEDT(I),I=1,3),(IEDT(I),I=5,7)
      END IF
C
C     Calculate elapsed time in days and seconds
      NDAYS=0
      LEAP=0
      IF (MOD(IEDT(1),4).EQ.0) LEAP = 1
      IBD = IBDT(3)            ! BEGIN DAY
      IED = IEDT(3)            ! END DAY
C     FIND DAYS
      IF (IBDT(2).NE.IEDT(2)) THEN
C       MONTHS DIFFER
        MB = IBDT(2)             ! BEGIN MONTH
        ME = IEDT(2)             ! END MONTH
        NM = ME-MB+1             ! NUMBER OF MONTHS TO LOOK AT
        IF (MB.GT.ME) NM = NM+12
        MC=MB-1
        DO 10 M=1,NM
          MC=MC+1                ! MC IS CURRENT MONTH
          IF (MC.EQ.13) MC = 1
          IF (MC.EQ.MB) THEN
            NDAYS = NDAYS+IDPM(MC)-IBD
            IF (MC.EQ.2) NDAYS = NDAYS + LEAP
          ELSEIF (MC.EQ.ME) THEN
            NDAYS = NDAYS+IED
          ELSE
            NDAYS = NDAYS+IDPM(MC)
            IF (MC.EQ.2) NDAYS = NDAYS + LEAP
          ENDIF
   10   CONTINUE
      ELSEIF (IBD.LT.IED) THEN
C       START AND END IN SAME MONTH, ONLY ACCOUNT FOR DAYS
        NDAYS = IED-IBD
      ENDIF
      ELSEC=NDAYS*NSPD
C
C     ADD OR SUBTRACT SECONDS
      ELSEC = ELSEC+(IEDT(5)-IBDT(5))*3600.0
      ELSEC = ELSEC+(IEDT(6)-IBDT(6))*60.0
      ELSEC = ELSEC+(IEDT(7)-IBDT(7))
      ELSEC = ELSEC+(IEDT(8)-IBDT(8))*0.001
C
C     CONVERT SECONDS TO DAYS, HOURS, MINUTES, AND SECONDS
      NDAYS = ELSEC/NSPD
      RSECS = MOD(ELSEC,86400.0)
      NHOURS = RSECS/3600.0
      RSECS = MOD(RSECS,3600.0)
      NMINS = RSECS/60.0
      RSECS = MOD(RSECS,60.0)
      NSECS = RSECS
      RSECS = MOD(RSECS,1.0)
      MSECS = NINT(RSECS*1000.0)
      NRSECS = NSECS
      IF (RSECS.GE.0.5) NRSECS=NRSECS+1
C
C     Write elapsed time to screen
        IF (NDAYS.GT.0) THEN
          WRITE(*,1010) NDAYS,NHOURS,NMINS,NRSECS
 1010     FORMAT(1X,'Elapsed run time: ',I3,' Days, ',I2,' Hours, ',I2,
     &      ' Minutes, ',I2,' Seconds',/)
        ELSEIF (NHOURS.GT.0) THEN
          WRITE(*,1020) NHOURS,NMINS,NRSECS
 1020     FORMAT(1X,'Elapsed run time: ',I2,' Hours, ',I2,
     &      ' Minutes, ',I2,' Seconds',/)
        ELSEIF (NMINS.GT.0) THEN
          WRITE(*,1030) NMINS,NSECS,MSECS
 1030     FORMAT(1X,'Elapsed run time: ',I2,' Minutes, ',
     &      I2,'.',I3.3,' Seconds',/)
        ELSE
          WRITE(*,1040) NSECS,MSECS
 1040     FORMAT(1X,'Elapsed run time: ',I2,'.',I3.3,' Seconds',/)
        ENDIF
C
C     Write times to file if requested
      IF(IPRTIM.GT.0) THEN
        IF (NDAYS.GT.0) THEN
          WRITE(IOUT,1010) NDAYS,NHOURS,NMINS,NRSECS
        ELSEIF (NHOURS.GT.0) THEN
          WRITE(IOUT,1020) NHOURS,NMINS,NRSECS
        ELSEIF (NMINS.GT.0) THEN
          WRITE(IOUT,1030) NMINS,NSECS,MSECS
        ELSE
          WRITE(IOUT,1040) NSECS,MSECS
        ENDIF
      ENDIF
C
      RETURN
      END
