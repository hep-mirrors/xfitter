      subroutine xFitter 
C--------------------------------------------------------
C
C> xFitter PDF Fit Program
C
C-------------------------------------------------------

      implicit none
      external fcn

#include "steering.inc"
#include "thresholds.inc"
#include "couplings.inc"
#include "for_debug.inc"
#include "ntot.inc"
#include "indata.inc"
#include "systematics.inc"
#include "g_offset.inc"
#include "covar.inc"
#include "theorexpr.inc"
#include "chi2scan.inc"
#include "extrapars.inc"

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
*     Print HFitter banner
*     ------------------------------------------------
C      call hfbanner

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

*
*  Read parameters:
*
      nExtraParam = 0
      call parse_params() !read parameters.yaml

*     ------------------------------------------------
*     Read the steering file steering.txt
*     ------------------------------------------------ 
      call read_steer

      call hf_errlog(12020501,
     +  'I: steering.txt has been read successfully')

* Init random numbers 
      call init_rnd_seeds()

*     ------------------------------------------------
*     Read the measured data points
*     ------------------------------------------------
      call read_data
      call hf_errlog(12020502,
     +     'I: data tables have been read successfully') 

      if (SCAN) then            ! chi2 scan
         call chi2_scan
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

      call init_minimizer()

      call run_minimizer()

      call run_error_analysis()
      
      
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
      
*Offset is broken in 2.2.1
*      if (doOffset) then
*        call Offset_Finalize(icond)
*        call flush(6)
*        if(icond .ne. 0) goto 36
*        Call RecovCentrPars
*      else
      if (ControlFitSplit) then
         Call FindBestFCN3  !> Overfitting protection.
      else
*
* Write out central parameters
*
        call write_pars(0)
      endif
      close(24)
      close(25)

      call report_convergence_status()

 36   continue

      call close_theor_eval

*     ------------------------------------------------
*     Print error log summary
*     ------------------------------------------------
      call HF_errsum(6)
      end
