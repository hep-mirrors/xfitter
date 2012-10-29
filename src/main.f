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

      integer icond
C-----------------------------------------------------
*     ------------------------------------------------
*     Print info message
*     ------------------------------------------------
      call h1fitterInfo

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
      
c WS: 2012-10-28 subroutine Do_Fit defined in 'minuit_ini.f'
c ..........................................................

      if (CorrSystByOffset .and. CorSysIndex .gt. NSYSMAX) then
        do CorSysIndex = 0,nSys
          call Do_Fit
        enddo
        do CorSysIndex = -nSys,-1
          call Do_Fit
        enddo
      else
        call Do_Fit
      endif
      
      
      if (CorrSystByOffset) then
        call Offset_Finalize
        goto 36
        ! --- no Bands in the current version
      else
        if (ControlFitSplit) then
           Call FindBestFCN3  !> Overfitting protection.
        else
*
* Write out central parameters
*
           call write_pars(0)
        endif
      endif

      
      if (DOBANDS) then
         write(6,*) ' --- Calculate error bands ...'
         lprint = .false.    
         call hf_errlog(12020506,
     +     'I: Calculation of error bands required') 
         call MNCOMD(fcn,'ITERATE 10',icond,0)
         call MNCOMD(fcn,'MYSTUFF 1000',icond,0)
         call MNCOMD(fcn,'MYSTUFF 2000',icond,0)
         call Error_Bands_Pumplin
         ! Error_Bands_Pumplin calls GetUmat(i,j),
         ! which needs only umat from common /umatco/
         ! See minuit/src/iterate.F
      endif
      close (24)
      close (25)

 36   continue

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



      Subroutine H1fitterinfo
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
      print *,'  Version 0.2.0                                                                         '
      print *,'  http://herafitter.hepforge.org                             herafitter-help@desy.de    '
      print *,'----------------------------------------------------------------------------------------'
      print *,' '
      print *,' '
      end
