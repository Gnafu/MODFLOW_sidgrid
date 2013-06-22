C    LAST CHANGE:  07MARCH2006 RBT
C
C----NOTE:	ISC FLAG ARRAY DETERMINES THE FORM OF THE SOIL CHARACTERISTIC
C----			FUNCTION.
C----			ISC > 0:  VAN GENUCHTEN / MUALEM FUNCTION (VAN GENUCHTEN, 1980)
C----			ISC =-2:  VAUCLIN / HAVERKAMP FUNCTION (VAUCLIN ET AL., 1975)
C----		    ISC =-3:  HAVERKAMP LOG FUNCTION (HAVERKAMP ET AL., 1977)
C=======================================================================
      FUNCTION SVSF1REF1PS(EFSAT,ALPHA,VGN,ISC)
C-----VERSION 07MARCH2006 RBT
C     ******************************************************************
C	CALCULATE PORE-WATER PRESSURE AS A FUNCTION OF EFFECTIVE 
C	SATURATION USING VAN GENUCHTEN-MUALEM FUNCTION (ISC>1), 
C	VAUCLIN FUNCTION (ISC=-2) OR HAVERKAMP LOG FUNCTION (ISC=-3)
C     ******************************************************************
C     SPECIFICATIONS:
      REAL EFSAT,ALPHA,VGN,M,PRESS
      INTEGER ISC 
C     ******************************************************************
C
	IF (ISC.GT.0) THEN
	  M=1-(1/VGN)
	  PRESS=(1./ALPHA)*(1./(EFSAT**(1/M))-1)**(1./VGN)
	  SVSF1REF1PS=-PRESS    
	ELSEIF(ISC.EQ.-2) THEN
	  PRESS=((ALPHA/EFSAT)-ALPHA)**(1/VGN)   
	  SVSF1REF1PS=-PRESS
	ELSE
	  P=((ALPHA/EFSAT)-ALPHA)**(1/VGN)
	  SVSF1REF1PS=-EXP(P)
	ENDIF
C
      END 
C
C=======================================================================
	FUNCTION SVSF1REF1KP(PRESS,ALPHA,VGN,ISC) 
C-----VERSION 07MARCH2006 RBT
C     ******************************************************************
C	CALCULATE RELATIVE PERMEABILITY AS A FUNCTION OF 
C	PORE-WATER PRESSURE USING VAN GENUCHTEN FUNCTION (ISC>1) 
C	OR VAUCLIN FUNCTION (ISC<0),
C     ******************************************************************
C     SPECIFICATIONS:
	REAL PRESS,ALPHA,VGN,M
	INTEGER ISC
C     ******************************************************************
C	
	IF (ISC.GT.0) THEN
	  M=1.-(1./VGN)
	  PRESS=-PRESS
    	  TOP=(1-(ALPHA*PRESS)**(VGN-1)*(1+(ALPHA*PRESS)**VGN)**(-M))**2
	  SVSF1REF1KP=TOP/((1+(ALPHA*PRESS)**VGN)**(M/2))
	  PRESS=-PRESS
	ELSE
  	  PRESS=-PRESS 
        SVSF1REF1KP=ALPHA/(ALPHA+(PRESS)**VGN)  
	  PRESS=-PRESS 
	ENDIF
C
	END
C
C=======================================================================
	FUNCTION SVSF1REF1SP(PRESS,RSAT,ALPHA,VGN,ISC)
C-----VERSION 07MARCH2006 RBT
C     ******************************************************************
C	CALCULATE SOIL WATER SATURATION AS A FUNCTION OF PORE-WATER 
C	PRESSURE USING VAN GENUCHTEN-MUALEM FUNCTION (ISC>1), 
C	VAUCLIN FUNCTION (ISC=-2), OR HAVERKAMP LOG FUNCTION (ISC=-3)
C     ******************************************************************
C     SPECIFICATIONS:
	REAL PRESS,RSAT,ALPHA,VGN,M
	INTEGER ISC
C     ******************************************************************
C	
	IF (ISC.GT.0) THEN
	  M=1.-(1./VGN)
	  PRESS=-PRESS
	  SVSF1REF1SP=RSAT+(1.-RSAT)/((1.+(ALPHA*PRESS)**VGN)**M)
	  PRESS=-PRESS
	ELSEIF (ISC.EQ.-2) THEN
	  PRESS=-PRESS
	  SVSF1REF1SP=RSAT+(ALPHA*(1-RSAT))/(ALPHA+PRESS**VGN)
	  PRESS=-PRESS
	ELSEIF (ISC.EQ.-3) THEN
	  PRESS=-PRESS
	  SVSF1REF1SP=RSAT+(ALPHA*(1-RSAT))/(ALPHA+LOG(PRESS)**VGN)
	  PRESS=-PRESS
	ENDIF
	END
C
C======================================================================
	FUNCTION SVSF1REF1CP(PRESS,RSAT,POR,ALPHA,VGN,ISC)
C-----VERSION 07MARCH2006 RBT
C	*****************************************************************
C	CALCULATE SPECIFIC MOISTURE CAPACITY AS A FUNCTION OF PORE-WATER
C	PRESSURE USING VAN GENUCHTEN FUNCTION (ISC>1),
C	VAUCLIN FUNCTION (ISC=-2), OR HAVERKAMP LOG FUNCTION (ISC=-3)
C	*****************************************************************
C	SPECIFICATIONS
	REAL SAT,RSAT,POR,ALPHA,VGN,M,EFSAT,BUFF,PRESS,RSWC
	INTEGER ISC
C     ******************************************************************
C
	IF (ISC.GT.0) THEN
	  M=1.-(1./VGN) 
	  PRESS= -PRESS
	  EFSAT=(1./(1.+(PRESS*ALPHA)**VGN))**M
	  BUFF=(1./(1.-M))*ALPHA*M*(POR-RSAT)*(EFSAT**(1./M))
	  SVSF1REF1CP=BUFF*((1.-EFSAT**(1./M))**M)
	  PRESS=-PRESS
	ELSEIF (ISC.EQ.-2) THEN
	  PRESS=-PRESS
	  RSWC=RSAT*POR
	  BUFF=ALPHA*VGN*(POR-RSWC)*PRESS**(VGN-1)
        SVSF1REF1CP=BUFF/((ALPHA+PRESS**VGN)**2.)
	  PRESS=-PRESS
	ELSEIF (ISC.EQ.-3) THEN
	  PRESS=-PRESS
	  RSWC=RSAT*POR
	  BUFF=ALPHA*VGN*(POR-RSWC)*(LOG(PRESS))**(VGN-1)
	  SVSF1REF1CP=BUFF/(PRESS*((ALPHA+(LOG(PRESS))**VGN)**2.))
	  PRESS=-PRESS
	ENDIF
C
	END  
