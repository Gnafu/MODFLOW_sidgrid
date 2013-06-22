MODULE cascading_module
! Define the variables:
integer , SAVE, DIMENSION(:,:),         POINTER :: direction  
real, SAVE, DIMENSION(:,:), POINTER :: slope
real, SAVE, DIMENSION(:,:), POINTER :: n_Mann
real, save, DIMENSION(:,:), POINTER :: CFLFLOW
integer, SAVE, DIMENSION(:,:), POINTER :: IRUNVSF

! Define the new type (used in case of LGR):
TYPE cascadingtype
	integer , DIMENSION(:,:),         POINTER :: direction  
	real, DIMENSION(:,:), POINTER :: slope
	real, DIMENSION(:,:), POINTER :: n_Mann
	real, DIMENSION(:,:), POINTER :: CFLFLOW
	integer, DIMENSION(:,:), POINTER :: IRUNVSF
	
END TYPE

TYPE(cascadingtype), SAVE:: CASCADINGDAT(10)
END MODULE
! ==============================================================
SUBROUTINE CASCADEDA(IUZF,Igrid)
!   Deallocate CASCADING DATA. 
    use cascading_module
    USE GWFUZFMODULE
    ! ---
    integer Igrid, IUZF
    deallocate (CASCADINGDAT(Igrid)%direction)
    deallocate (CASCADINGDAT(Igrid)%slope)
    deallocate (CASCADINGDAT(Igrid)%n_Mann)
    deallocate (CASCADINGDAT(Igrid)%CFLFLOW)
    ! if UZF not active, deallocate variables:
    if (IUZF<=0) then 
		deallocate (GWFUZFDAT(Igrid)%EXCESPP)
		deallocate (GWFUZFDAT(Igrid)%SEEPOUT)
		deallocate (CASCADINGDAT(Igrid)%IRUNVSF)
	end if  
    
end subroutine CASCADEDA
! ==============================================================
SUBROUTINE CASCADEPNT(IUZF,Igrid)
    ! SET POINTERS FOR THE CASCADING PART:
    use cascading_module
    USE GWFUZFMODULE
    ! ---
    integer Igrid, IUZF
    !
    direction=>CASCADINGDAT(Igrid)%direction
    slope=>CASCADINGDAT(Igrid)%slope
    n_Mann=>CASCADINGDAT(Igrid)%n_Mann
    CFLFLOW=>CASCADINGDAT(Igrid)%CFLFLOW
    ! if UZF not active, set UZF pointers:
    if (IUZF<=0) then 
		SEEPOUT=>GWFUZFDAT(Igrid)%SEEPOUT
        EXCESPP=>GWFUZFDAT(Igrid)%EXCESPP
        IRUNVSF=>CASCADINGDAT(Igrid)%IRUNVSF
	end if
end subroutine CASCADEPNT

! ==============================================================
SUBROUTINE CASCADESV(IUZF,Igrid)
    ! SAVE THE CASCADING POINTER for grid:
    use cascading_module
    USE GWFUZFMODULE
    ! ---
    integer Igrid,IUZF
	CASCADINGDAT(Igrid)%direction=>direction
	CASCADINGDAT(Igrid)%slope=>slope
	CASCADINGDAT(Igrid)%n_Mann=>n_Mann
	CASCADINGDAT(Igrid)%CFLFLOW=>CFLFLOW
	! if UZF not active, save UZF pointers:
    if (IUZF<=0) then 
		GWFUZFDAT(Igrid)%SEEPOUT=>SEEPOUT 
        GWFUZFDAT(Igrid)%EXCESPP=>EXCESPP 
        CASCADINGDAT(Igrid)%IRUNVSF=>IRUNVSF
	end if
end subroutine CASCADESV


! ==============================================================
!! Read and allocate procedure for Cascading Flow process
subroutine CASCADEAR(IN,IUZF,IGRID,Iunitsfr,Iunitlak)
	! Input:
	! IN  is the unit number of the cascading input file
	! IUZF is a flag. IUZF>0 or <=0 if UZF is active or not, respectively.
	! IGRID, the integer of the grid (LGR method)
	! Output:
	! allocation of cascading variables and reading in input file
	use cascading_module
	use GWFUZFMODULE, only: EXCESPP, SEEPOUT
	use GLOBAL,       only: NCOL, NROW, IOUT
	implicit none
	
	! Declare global variables:
	integer :: IN,IUZF,IGRID
	integer :: Iunitsfr,Iunitlak  ! Unit file for streams and lakes
	
	! Declare inner variable
	integer :: ir, ic
	
	! IF UZF is not active (IUFZ<=0) then 
	! provide allocation and save for UZF variables needed for cascading
	! provide allocation and save variable IRUNVSF needed for cascading
	if (IUZF<=0) then 
		allocate (SEEPOUT(NCOL,NROW), EXCESPP(NCOL,NROW))
		allocate (IRUNVSF(NCOL,NROW))
	end if 
    ! write in the main output file
	WRITE (IOUT, *) '                                    '
    WRITE (IOUT, *) ' --- ************************** --- '
    WRITE (IOUT, *) ' --- USING NEW CASCADING MODULE --- '
    WRITE (IOUT, *) ' --- ************************** --- '
    WRITE (IOUT, *) '                                    '

    !Allocate space for Cascading variables
    allocate (direction(NROW,NCOL))
    allocate (slope(NROW,NCOL))
    allocate (n_Mann(NROW,NCOL))
    allocate (CFLFLOW(NROW,NCOL))
        
    ! Read data from the input file
    do ir=1, NROW
		do ic=1, NCOL
            READ (IN,*) slope(ir,ic), direction(ir,ic), n_Mann(ir,ic)
        end do
    end do
    ! IF UZF NOT active, then proceede and read also the IRUNVSF variable
    !                    IF Streams or lakes are specified in the model:
    if (IUZF<=0) then
      if (Iunitsfr>0.OR.Iunitlak>0) then  
		do ir=1, NROW
		   do ic=1, NCOL
               READ (IN,*) IRUNVSF(ic,ir)
           end do
        end do
      else 
        IRUNVSF=0
      end if 
	end if 
    
! Save POINTERS FOR GRID AND RETURN.
    CALL CASCADESV(IUZF,Igrid)
    RETURN
end subroutine CASCADEAR	


! ==============================================================
! Compute the outgoing  flow from the cell (irow,jcol)
!
subroutine CASCADE_OUT(OutgoingFlux,irow,jcol,KPER,IUZF,IGRID)
	
	use GLOBAL,       only: DELC, DELR, PERLEN, ITMUNI,LENUNI         	
	use GWFUZFMODULE, only: EXCESPP
	use cascading_module
	implicit none
		
	! Declare global variables:
	integer :: irow,jcol, KPER, IUZF,IGRID
	real :: OutgoingFlux
	
	! Declare local variables
	real  :: L, w, Const 
	real  :: discharge_out, Qout
	
	! SET POINTERS FOR THE CURRENT GRID.
	CALL CASCADEPNT(IUZF,IGRID)
	
	! Select the Manning Constant, according to System of Units
	Const = 1.0 ! in m^(1/3)/s
	if (ITMUNI == 1) then ! Time in sec
		if (LENUNI == 1) Const = 1.486         ! Length in feet
		if (LENUNI == 3) Const = 4.64          ! Length in cm
	end if
	if (ITMUNI == 2) then ! Time in minutes
		if (LENUNI == 1) Const = 1.486*60.0    ! Length in feet
		if (LENUNI == 2) Const = 60.0          ! Length in m  
		if (LENUNI == 3) Const = 4.64*60.0     ! Length in cm 
	end if
	if (ITMUNI == 3) then ! Time in hours
		if (LENUNI == 1) Const = 1.486*3600    ! Length in feet
		if (LENUNI == 2) Const = 3600.0        ! Length in m
		if (LENUNI == 3) Const = 4.64*3600 ! Length in cm 
	end if
	if (ITMUNI == 4) then ! Time in days
		if (LENUNI == 1) Const = 1.486*86400
		if (LENUNI == 2) Const = 86400.0
		if (LENUNI == 3) Const = 86400/0.000001
	end if
	! -------
			
	!! IF the cell slope is not zero, then route the excess flow outside
	IF (slope(irow,jcol)/=0.AND.EXCESPP(jcol,irow)/=0) THEN
		! compute length (L) and width (w) of the cell with respect to its slope:
	 call slope_geometry(L,w,direction(irow,jcol),DELC(irow),DELR(jcol))
		Qout = discharge_out(EXCESPP(jcol,irow), PERLEN(KPER-1), & 
  &		       PERLEN(KPER),slope(irow,jcol), w,L,n_Mann(irow,jcol),Const)
		
		OutgoingFlux = Qout/(DELC(irow)*DELR(jcol))
		
	ELSE
		OutgoingFlux = 0
	END IF
	
end subroutine CASCADE_OUT	

! =====================================================================
! Compute the total incoming flux for the cell (irow,jcol) for the KPER-th stress period
!
subroutine CASCADE_IN(IncomingFlux, irow,jcol, KPER, IUZF, IGRID)
	    	
	use GLOBAL,       only: DELC, DELR, PERLEN, NROW, NCOL,ITMUNI,LENUNI 
	use GWFUZFMODULE, only: EXCESPP
	use cascading_module
	implicit none
	
	! Declare local variables
	real :: discharge_out, Flux_temp
	integer :: i,j,MR,MC
	real :: DeltaTold, DeltaT, Area, L, w
	real, dimension(:,:), allocatable :: outflow
	
	! Declare global variables:
	integer :: irow,jcol, KPER, IUZF, IGRID
	real :: IncomingFlux, Const
		
	!!! ---
	! SET POINTERS FOR THE CURRENT GRID.
	CALL CASCADEPNT(IUZF,IGRID)
	
	
	! Select the Manning COnstant, according to System of Units
	Const = 1.0 ! in m^3/s
	if (ITMUNI == 1) then ! Time in sec
		if (LENUNI == 1) Const = 1.486
		if (LENUNI == 3) Const = 1.0/0.000001
	end if
	if (ITMUNI == 2) then ! Time in minutes
		if (LENUNI == 1) Const = 1.486*60.0
		if (LENUNI == 2) Const = 60.0
		if (LENUNI == 3) Const = 60.0/0.000001
	end if
	if (ITMUNI == 3) then ! Time in hours
		if (LENUNI == 1) Const = 1.486*3600
		if (LENUNI == 2) Const = 3600.0
		if (LENUNI == 3) Const = 3600/0.000001
	end if
	if (ITMUNI == 4) then ! Time in days
		if (LENUNI == 1) Const = 1.486*86400
		if (LENUNI == 2) Const = 86400.0
		if (LENUNI == 3) Const = 86400/0.000001
	end if
	! -------
	
	! Introduce an internal (faster) notation:
	i=irow
	j=jcol
	MR=NROW !size(EXCESPP(1,:))! MAX_row
	MC=NCOL !size(EXCESPP(:,1))! MAX_col
	
	
	DeltaTold = PERLEN(KPER-1) ! previous stress period length
	DeltaT    = PERLEN (KPER)  ! stress period length
	
	ALLOCATE (outflow(MR,MC))
	
	! Set the index of direction to 0 for cells with slope==0:
	do i=1,MR
		do j=1,MC
			if (slope(i,j).LE.0) then 
				direction(i,j) = 0
			end if
		end do
	end do
	! Update cell position	 
	i=irow
	j=jcol	
	

	
	! Compute, for any cell, the outgoing flow (outflow)
	do i = 1, MR
		do j =1, MC
		    ! compute length (L) and width (w) of the cell with respect to its slope:
		    call slope_geometry(L,w,direction,DELC(i),DELR(j))
			! compute the flow calling the function "discharge_out"
			if (EXCESPP(j,i)>=0.AND.slope(i,j) > 0) then
				outflow(i,j) = discharge_out(EXCESPP(j,i), DeltaTold, &
				               & DeltaT, slope(i,j), w, L, n_Mann(i,j),Const) 
                             
			else
				outflow(i,j)= 0
			end if
		end do
	end do
	
	
	! Update cell position	 
	i=irow
	j=jcol	

	Area = DELC(i)*DELR(j)

	! Inizialize the incoming flux:
	Flux_temp = 0 
	
	if (i /= 1.AND.i/=MR.AND.j/=1.AND.j/= MC) then
		if (direction(i-1,j)==5) then 
		 Flux_temp = Flux_temp + outflow(i-1,j)/Area
		end if 
		if (direction(i-1,j+1)==6)  Flux_temp = Flux_temp + outflow(i-1,j+1)/Area
		if (direction(i,j+1)==7)  Flux_temp = Flux_temp + outflow(i,j+1)/Area
		if (direction(i+1,j+1)==8)  Flux_temp = Flux_temp + outflow(i+1,j+1)/Area
		if (direction(i+1,j)==1)  Flux_temp = Flux_temp + outflow(i+1,j)/Area
		if (direction(i+1,j-1)==2)  Flux_temp = Flux_temp + outflow(i+1,j-1)/Area
		if (direction(i,j-1)==3)  Flux_temp = Flux_temp + outflow(i,j-1)/Area
		if (direction(i-1,j-1)==4)  Flux_temp = Flux_temp + outflow(i-1,j-1)/Area
	end if 
	
	if (i==1) then
		if (j==1) then
			if (direction(1,2)==7) Flux_temp = Flux_temp + outflow(1,2)/Area
			if (direction(2,2)==8) Flux_temp = Flux_temp + outflow(2,2)/Area
			if (direction(2,1)==1) Flux_temp = Flux_temp + outflow(2,1)/Area
		else if (j==MC) then
			if (direction(1,MC-1)==3) Flux_temp = Flux_temp + outflow(1,MC-1)/Area
			if (direction(2,MC-1)==2) Flux_temp = Flux_temp + outflow(2,MC-1)/Area
			if (direction(1,MC)==1) Flux_temp = Flux_temp + outflow(1,MC)/Area
		else 
			if (direction(i,j+1)==7)  Flux_temp = Flux_temp + outflow(i,j+1)/Area
			if (direction(i+1,j+1)==8)  Flux_temp = Flux_temp + outflow(i+1,j+1)/Area
			if (direction(i+1,j)==1)  Flux_temp = Flux_temp + outflow(i+1,j)/Area
			if (direction(i+1,j-1)==2)  Flux_temp = Flux_temp + outflow(i+1,j-1)/Area
			if (direction(i,j-1)==3)  Flux_temp = Flux_temp + outflow(i,j-1)/Area
		end if
	end if
	
	if (i==MR) then
		if (j==1) then
			if (direction(MR-1,1)==5) Flux_temp = Flux_temp + outflow(MR-1,1)/Area
			if (direction(MR-1,2)==6) Flux_temp = Flux_temp + outflow(MR-1,2)/Area
			if (direction(MR,2)==7) Flux_temp = Flux_temp + outflow(MR,2)/Area
		else if (j==MC) then
			if (direction(MR-1,MC)==5) Flux_temp = Flux_temp + outflow(MR-1,MC)/Area
			if (direction(MR,MC-1)==3) Flux_temp = Flux_temp + outflow(MR,MC-1)/Area
			if (direction(MR-1,MC-1)==4) Flux_temp = Flux_temp + outflow(MR-1,MC-1)/Area
		else 
			if (direction(i-1,j)==5)  Flux_temp = Flux_temp + outflow(i-1,j)/Area
			if (direction(i-1,j+1)==6)  Flux_temp = Flux_temp + outflow(i-1,j+1)/Area
			if (direction(i,j+1)==7)  Flux_temp = Flux_temp + outflow(i,j+1)/Area
			if (direction(i,j-1)==3)  Flux_temp = Flux_temp + outflow(i,j-1)/Area
			if (direction(i-1,j-1)==4)  Flux_temp = Flux_temp + outflow(i-1,j-1)/Area
		end if
	end if	
	
	if (j==1.AND.i/=1.AND.i/=MR) then
		if (direction(i-1,1)==5) Flux_temp = Flux_temp + outflow(i-1,1)/Area
		if (direction(i-1,2)==6) Flux_temp = Flux_temp + outflow(i-1,2)/Area
		if (direction(i,2)==7)   Flux_temp = Flux_temp + outflow(i,2)/Area
		if (direction(i+1,2)==8) Flux_temp = Flux_temp + outflow(i+1,2)/Area
		if (direction(i+1,1)==1) Flux_temp = Flux_temp + outflow(i+1,1)/Area
	end if
	
	if (j==MC.AND.i/=1.AND.i/=MR) then
		if (direction(i-1,MC)==5) Flux_temp = Flux_temp + outflow(i-1,MC)/Area
		if (direction(i+1,MC)==1) Flux_temp = Flux_temp + outflow(i+1,MC)/Area
		if (direction(i+1,MC-1)==2) Flux_temp = Flux_temp + outflow(i+1,MC-1)/Area
		if (direction(i,MC-1)==3) Flux_temp = Flux_temp + outflow(i,MC-1)/Area
		if (direction(i-1,MC-1)==4) Flux_temp = Flux_temp + outflow(i-1,MC-1)/Area
	end if
	
	!The output is Flux_temp 
	IncomingFlux = Flux_temp 
	IF (IncomingFlux.LT.0) IncomingFlux=0
	
end subroutine CASCADE_IN

real function discharge_out(excessflow, DTold, DT, So, w, L, n,Const)
	
	implicit none
	
	! Global variables
	real :: Qout
	real :: excessflow
	real :: So, n, Const
	real :: DTold, DT, L, w
	
	! Inner variables
	real ::  Hin, Qin, appoggio
	double precision :: A
	
	!! The output is the volumetric flow (dims = L^3/T) which discharges from the cell:
	! calculate the parameters
	A = DT*( (n/(w*sqrt(So)))**(5/3)  )/(w*L)
	Hin = 2*DTold*excessflow/(w*L)
	Qin = (Hin**(5/3))*Const*w*sqrt(So)/n
	if (A>0) then 
	    discharge_out = (( sqrt(1+4*A*(Qin**(5/3))) -1)/(2*A) )**2 !(DTold/DT)*(L*w)*excessflow !
	else 
	    discharge_out = 0
	end if 
	 
	! check and correct:
	if (discharge_out > excessflow) discharge_out = excessflow
	
end function discharge_out
!! =====================================================
SUBROUTINE VSF2OLF(Iunitsfr, Iunitlak)
!     ******************************************************************
!     ASSIGN OVERLAND RUNOFF AS INFLOW TO STREAMS AND LAKES
!     VERSION SID&GRID 2001
!     === This subroutine corresponds to the subroutine SGWF2UZF1OLF
!         of package UZF. Minor revisions are needed to adapt it 
!         to the new version of VSF
!     ******************************************************************
      USE GWFUZFMODULE, ONLY: SEEPOUT, EXCESPP
      USE GLOBAL,       ONLY: NCOL, NROW
      USE GWFSFRMODULE, ONLY: NSS, NSTRM, ISTRM, SEG, STRM
      USE GWFLAKMODULE, ONLY: NLAKES, OVRLNDRNF
      use cascading_module, ONLY: IRUNVSF
      IMPLICIT NONE
!     -----------------------------------------------------------------
!     SPECIFICATIONS:
!     -----------------------------------------------------------------
!     ARGUMENTS
!     -----------------------------------------------------------------
      INTEGER Iunitsfr, Iunitlak
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
      INTEGER i, ic, ir, nseg, irun
      REAL rchlen, seglen, gwrunof, seepout1, TOTRUNOFF
!     -----------------------------------------------------------------
!C1-----INITIALIZE OVERLAND RUNOFF TO STREAMS AND LAKES TO ZERO.
      gwrunof = 0.0
      TOTRUNOFF = 0.0
!C
!C2-----OVERLAND RUNOFF TO STREAMS FROM UZF SEPARATE FROM 
!C        OVERLAND RUNOFF SPECIFIED FOR STREAMS.
      IF ( Iunitsfr.GT.0 ) THEN
        DO i = 1, NSS
          SEG(26, i) = 0.0
        END DO
        DO i = 1, NSTRM
          STRM(24, i) = 0.0
        END DO
      END IF
      IF ( Iunitlak.GT.0 ) THEN
!
!C3------OVERLAND RUNOFF TO LAKES SET TO ZERO EVEN IF
!C        VALUE IS SPECIFIED FOR A LAKE.
        DO i = 1, NLAKES
          OVRLNDRNF(i) = 0.0
        END DO
      END IF
!C     

!C4------LOOP THROUGH IRUNVSF ARRAY AND ADD SEEPOUT PLUS EXCESSPP TO
!C         CORRECT STREAM SEGMENT OR LAKE.
      DO ir = 1, NROW
        DO ic = 1, NCOL
          seepout1 = SEEPOUT(ic, ir) + EXCESPP(ic, ir) 
          TOTRUNOFF = TOTRUNOFF + seepout1
          IF ( seepout1.GT.0.0 ) THEN
            irun = IRUNVSF(ic, ir)
            IF ( irun.GT.0 .AND. irun.LE.NSS .AND. Iunitsfr.GT.0 ) THEN
              SEG(26, irun) = SEG(26, irun) + seepout1
            ELSE IF (irun<=0.AND.ABS(irun)<=NLAKES .AND.Iunitlak>0) THEN
              OVRLNDRNF(ABS(irun)) = OVRLNDRNF(ABS(irun)) + seepout1
            END IF
          END IF
          SEEPOUT(ic, ir) = 0.0
        END DO
      END DO
!
!C5------PROPORTION RUNOFF TO REACHES ON BASIS OF STREAM LENGTH.
      IF ( Iunitsfr.GT.0 ) THEN
        DO i = 1, NSTRM
          nseg = ISTRM(4, i)
          seglen = SEG(1, nseg)
          gwrunof = SEG(26, nseg)
          rchlen = STRM(1, i)
          STRM(24, i) = gwrunof*(rchlen/seglen)
        END DO
      END IF
!
!C6-----RETURN.
      RETURN
END SUBROUTINE VSF2OLF





!! =======================================================
subroutine  slope_geometry(L,w,direction,DC,DR)
	implicit none
	
	! global variables
	real :: L,w,DC,DR
	integer :: direction
	
	! inner variables
	real :: diag
	
	! -- define the diagonal of the cell
	diag = sqrt(DC**2 + DR**2)
	
	! define (L,w) according to the direction of maximum slope:
	
	if (direction == 1 .OR. direction == 5) then
		L = DC
		w = DR
	else if (direction == 3 .OR. direction == 7) then
		L = DR
		w = DC
	else
		L = diag
		w = 0.66*diag
	end if
	
end subroutine slope_geometry

! ==============================================================
! If VSF is active (=> UZF not active):
! after any stress period, compute the total runoff flow (Hortonian and Dunnian)
! and write the specific output files.

subroutine CASCADEBD(KSTEP,KPER,IGRID)
    use GWFUZFMODULE, only: SEEPOUT, EXCESPP
    use GLOBAL,       only: NCOL, NROW, NPER
    use cascading_module
    
    implicit none

    real :: tot_Horton, tot_Dunnian, tot_Cascading
    integer :: I,J, KSTEP,KPER, IGRID
    ! ---
    if (KSTEP >= 2) return
	tot_Horton=0
	tot_Dunnian=0
	tot_Cascading = 0
	do I=1,NROW
	    do  J=1,NCOL
            tot_Horton = tot_Horton + EXCESPP(J,I)
            tot_Dunnian = tot_Dunnian + SEEPOUT(J,I)
            tot_Cascading = tot_Cascading + ABS(CFLFLOW(i,j))
        end do
    end do
    
    ! At the First Stress Period open the files
    if (KPER==1) then 
		open(unit=900,file="tot_Horton_runoff.dat")
		open(unit=910,file="tot_Dunnian_runoff.dat")
		open(unit=920,file="Cascading_flow.dat")
	end if
	
    ! Write the output in files
    if (tot_Horton >0) then
      write (900,*) '## Stress Period ', 'Time Step ', 'Total Hortonian Runoff (M^3/L)'
      write (900,*) KPER, KSTEP, tot_Horton
      write (900,*) '-- Array of top layer cells and relative Hortonian flow'
	  do i =1,NROW 
	  	  write (900,'(11E13.5)') (EXCESPP(j,i),j=1,NCOL)
	  end do
    end if
      
    if (tot_Dunnian  >0) then
	  write (910,*) '## Stress Period ', 'Time Step ', '  Total Dunnian Runoff (M^3/L)'	
	  write (910,*) KPER, KSTEP, tot_Dunnian
	  write (910,*) '-- Array of top layer cells and relative Dunnian flow'
	  do i =1,NROW 
	  	  write (910,'(11E13.5)') (SEEPOUT(j,i),j=1,NCOL)
	  end do 
    end if
    
    if (tot_Cascading .NE. 0) then
	  write (920,*) '## ---- Stress Period ', 'Time Step ', '  Total Cascading flow (M^3/L)'	
	  write (920,*) KPER, KSTEP, tot_Cascading
	  write (920,*) '-- Array of top layer cells and relative cascading flow'
	  do i =1,NROW 
	  	  write (920,'(11E13.5)') (CFLFLOW(i,j),j=1,NCOL)
	  end do
    end if
    ! At the Last Stress Period close the files
    if (KPER==NPER) then
		close(900)
		close(910)
		close(920)
	end if

end subroutine CASCADEBD


!C=======================================================================
SUBROUTINE CASCADERFU(KKPER,IGRID)
!C-----Version 2012  -  Iac
!C     ******************************************************************
!C     RECHARGE FLOW UPDATE 
!C     --------------------
!C     CALLED by MAIN FILE ONLY IF CFL Package IS ACTIVE
!C     ADD to EXCESSPP variable the (possible) EXCESSFLOW
!C 
      use GLOBAL, ONLY:NCOL,NROW,NLAY,BOTM,NBOTM,LAYCBD,LBOTM,DELR,DELC,PERLEN
      use GWFLPFMODULE, ONLY:VKA,LAYVKA,HK

      use GWFRCHMODULE, ONLY:RECH     
      use VSFMODULE,    ONLY:CVS
      use cascading_module
      use GWFUZFMODULE, ONLY:EXCESPP
      implicit none
!C -- !!         
!C     -----------------------------------------------------------------
!C     LOCAL VARIABLES
!C     ------------------------------------------------------------------
      INTEGER i, j, KKPER, KSTP,KITER,IGRID
      REAL    HDPOND,VCOND,THCK,THCK1,THCK2,RRATE,TTOP,VertKS
      real    surfinf, IncomingFlux, OutgoingFlux

!------SET POINTERS FOR THE CURRENT GRID.
	  !CALL SGWF2UZF1PNT(Igrid)
      !CALL SGWF2RCH7PNT(Igrid)
      !CALL CASCADEPNT(0,Igrid)
      
!C     ------------------------------------------------------------------	        
      if (kkper < 2) return
      
      do i = 1, NROW
         do j  = 1, NCOL
            THCK1=(BOTM(j,i,LBOTM(1)-1)-BOTM(j,i,LBOTM(2)-1))/2.
		    THCK2=(BOTM(j,i,LBOTM(2)-1)-BOTM(j,i,LBOTM(3)-1))/2.
		    THCK=(THCK1+THCK2)/2
		    
		    ! Get Vertical Hydraulyc Conductivity from LPF variables
		    if (LAYVKA(1) == 0) THEN
				VertKS = VKA(j,i,1)
			else
				VertKS = HK(J,I,1)/VKA(J,I,1)
			end if		    
		     
!C Get Recharge rate (made volumetric by GWF2RCH7RP ) and return the 
!C     flow rate for unit area:
			surfinf = RECH(j,i)/(DELR(j)*DELC(i))
!C      ---> call the subroutines
!C      ---> cascade_in, computing the IncomingFlux
!C      ---> cascade_out, computing the OutgoingFlux
!C      ---> and add (IncomingFlux - OutgoingFlux) to rainfall 
            call CASCADE_IN(IncomingFlux,i,j,Kkper,0,Igrid)
            call CASCADE_OUT(OutgoingFlux,i,j,Kkper,0,Igrid) 
            ! Save in buffer CFLFLOW the Volumetric flow 
            CFLFLOW(i,j) = (IncomingFlux - OutgoingFlux)*DELR(j)*DELC(i) 
            
            surfinf = surfinf  +  (IncomingFlux - OutgoingFlux)
!C --- check for negative indiltration rate and, just in case, set to 0         
            if (surfinf < 0.0) then
				surfinf = 0.0
!C5------SET INFILTRATION RATE TO SATURATED VERTICAL K WHEN RATE IS
!        GREATER THAN K AND ROUTE EXCESS WATER TO STREAM IF 
!        IRUNFLG IS NOT EQUAL TO ZERO.
! -----------------
            else if (surfinf.GT.VertKS) then
				EXCESPP(j,i) = (surfinf-VertKS)*DELC(i)*DELR(j)
                surfinf  = VertKS 
            end if       
!
! ------Update rainfall data (and MULTIPLY BY CELL AREA TO GET VOLUMETRIC RATE). 
            RECH(j,i) = surfinf*DELR(j)*DELC(i)  
          end do
        end do 
RETURN
END SUBROUTINE CASCADERFU
!=======================================      
