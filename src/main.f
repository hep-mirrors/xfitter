      program HERAFitter 
C--------------------------------------------------------
C
C> HERA PDF Fit Program
C
C-------------------------------------------------------

      implicit none
      external fcn

      include 'steering.inc'
      include 'thresholds.inc'
      include 'couplings.inc'
      include 'for_debug.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'systematics.inc'
      include 'g_offset.inc'	
      include 'covar.inc'

      integer icond
      integer nOffset
      ! logical doOffset  ! defined in 'systematics.inc'
      double precision www(NTot)
      double precision test(4,4),anui(4,4),test2(4,4)
      integer ii,jj,kk,narg
      integer command_argument_count
      character*80 command
C-----------------------------------------------------
*     ------------------------------------------------
*     Print info message
*     ------------------------------------------------
      call herafitterInfo

      narg = command_argument_count()
      if (narg.gt.0) then
         call get_command_argument(1,command)
         if (index(command,'--convert-covar').ne.0) then
            if (index(command,'=').ne.0) then
               command = command(index(command,'=')+1:len(command))
            else
               command = 'covar.in'
            endif
            call CovMatrixConverter(command)
            goto 36
         endif
      endif

*     ------------------------------------------------
*     Read the steering file steering.txt
*     ------------------------------------------------ 
      call read_steer

* Init random numbers 
      call init_rnd_seeds()

      call hf_errlog(12020501,
     +     'I: steering.txt has been read successfully') 

*     ------------------------------------------------
*     Read the measured data points
*     ------------------------------------------------
      call read_data
      call hf_errlog(12020502,
     +     'I: data tables have been read successfully') 

*     ------------------------------------------------
*     Initialise theory modules
*     ------------------------------------------------
      call init_theory_modules
      call hf_errlog(12020503,
     +     'I: theory modules initialised successfully') 


*     ------------------------------------------------
*     Do NNPDF if initialized
*     ------------------------------------------------
      if (FLAGRW) then 
         call pdfreweighting
         if (DORWONLY) goto 36
      endif

      if (LHAPDFErrors) then  ! PDF errors
         call get_lhapdferrors
         goto 36
      endif

*     ------------------------------------------------
*     Do the fit
*     ------------------------------------------------
      
c Modifications for the Offset method by WS & JT
c WS: 2012-10-28 subroutine Do_Fit defined in 'minuit_ini.f'
c ..........................................................

      call flush(6)
      nOffset = ProbeOffset(SysForm, nSys)
      doOffset = nOffset.gt.0
      ! print *,' --- Offset: ',nOffset,doOffset
      ! --- ignore CorSysIndex if no Offset type errors
      if (doOffset) then
        print *,' '
        print *,' --- Offset mode: '
        print *,'  ',nOffset,' offset type errors'
        if (CorSysIndex .gt. NSYSMAX) then
          print *,'  ALL offset fits selected.'
        else
          print *,'  SINGLE offset fit selected:  CorSysIndex = ',CorSysIndex
        endif
        print *,'  UsePrevFit = ',UsePrevFit
      else
        CorSysIndex = 0
      endif
      
      if (doOffset .and. CorSysIndex .gt. NSYSMAX) then
        do CorSysIndex = 0,nOffset
          call Do_Fit
        enddo
        do CorSysIndex = -nOffset,-1
          call Do_Fit
        enddo
      else
        call Do_Fit
      endif
      
      
      if (doOffset) then
        call Offset_Finalize(icond)
        call flush(6)
        if(icond .ne. 0) goto 36
        Call RecovCentrPars
        if (DOBANDS) then
           write(6,*) ' --- Calculating error bands from Offset errors...'
           ! --- set the scale by defining delta chi2 value
           ! --- for MNCOMD(fcn,'ITERATE 10',...) it is set by SET ERRDEF dchi2
           ! --- with default value of 5
           Call DecorDiag(5.d0)
           call Error_Bands_Pumplin
           ! Error_Bands_Pumplin calls GetUmat(i,j) which needs only umat from common /umatco/
           ! See minuit/src/iterate.F
        endif
      else
        if (ControlFitSplit) then
           Call FindBestFCN3  !> Overfitting protection.
        else
*
* Write out central parameters
*
           call write_pars(0)
        endif

        if (DOBANDS) then
          write(6,*) ' --- Calculate error bands ...'
          lprint = .false.    
          call hf_errlog(12020506, 'I: Calculation of error bands required')
          call MNCOMD(fcn,'ITERATE 10',icond,0)
          call MNCOMD(fcn,'MYSTUFF 1000',icond,0)
          call MNCOMD(fcn,'MYSTUFF 2000',icond,0)
          call Error_Bands_Pumplin
        endif
      
        close (24)
        close (25)
      endif

 36   continue

      call close_theor_eval

*     ------------------------------------------------
*     Print error log summary
*     ------------------------------------------------
      call HF_errsum(6)

*     ------------------------------------------------
*     Done
*     ------------------------------------------------
      stop
      end

C-----------------------------------------------------



      Subroutine HERAfitterinfo
*     ------------------------------------------------
      print *,' '
      print *,' '
      print *,'----------------------------------------------------------------------------------------'
      print *,'                                                                                        '
      print *,'  HH   HH  EEEEEEE  RRRRR      AAA    FFFFFFF  II  TTTTTTTT  TTTTTTTT  EEEEEEE  RRRRR   '
      print *,'  HH   HH  EE       RR   RR  AA   AA  FF       II     TT        TT     EE       RR   RR '
      print *,'  HHHHHHH  EEEEE    RR   RR  AA   AA  FFFFF    II     TT        TT     EEEEEE   RR   RR '
      print *,'  HH   HH  EE       RRRRR    AAAAAAA  FF       II     TT        TT     EE       RRRRR   '
      print *,'  HH   HH  EE       RR  RR   AA   AA  FF       II     TT        TT     EE       RR  RR  '
      print *,'  HH   HH  EEEEEEE  RR   RR  AA   AA  FF       II     TT        TT     EEEEEEE  RR   RR '
      print *,'                                                                                        '
      print *,'  Version trunk                                                                         '
      print *,'  http://herafitter.org                             herafitter-help@desy.de    '
      print *,'----------------------------------------------------------------------------------------'
      print *,' '
      print *,' '
      end
