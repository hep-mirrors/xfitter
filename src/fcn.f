C---------------------------------------------------------
C> @brief  Main minimization subroutine for MINUIT
C>
C> @param  npar      number of currently variable parameters
C> @param  g_dummy   the (optional) vector of first derivatives
C> @param  chi2out   the calculated function value
C> @param  parminuit vector of (constant and variable) parameters from minuit.in.txt
C> @param  iflag     flag set by minuit (1-init, 2-iteration, 3-finalisation)
C> @param  futil     name of utilitary routine
C---------------------------------------------------------
      subroutine  fcn(npar,g_dummy,chi2out,parminuit,iflag,futil)

      implicit none

*     ---------------------------------------------------------
*     declaration related to minuit
*     ---------------------------------------------------------

      integer npar,iflag
      double precision g_dummy(*),parminuit(*),chi2out,futil
      external futil

#include "ntot.inc"
#include "fcn.inc"
#include "endmini.inc"
#include "for_debug.inc"
#include "systematics.inc"
      integer i
      double precision chi2data_theory !function

C Store FCN flag in a common block:
      IFlagFCN = IFlag

      NparFCN  = npar

C Store params in a common block:
      do i=1,MNE
         parminuitsave(i) = parminuit(i)
      enddo

C Count number of FCN calls:
      if (iflag.eq.3) then
         nfcn3 = nfcn3 + 1
         if (nfcn3.gt.MaxFCN3) then
            print *,'Fatal error: too many FCN 3 call'
            print *,'Increase number of MaxFCN3 calls in endmini.inc'
            print *,'Or reduce number of FNC 3 calls'
            print *,'Stop'
            call hf_stop
         endif
      endif

C Store only if IFlag eq 3:
      if (iflag.eq.3) then
         do i=1,MNE
            pkeep(i) = parminuit(i)
C !> Also store for each fcn=3 call:
            pkeep3(i,nfcn3) = parminuit(i)
         enddo
      endif

      call HF_errlog(12020515,'I: FCN is called')

C     Print MINUIT extra parameters
c which are actually all parameters
      call printminuitextrapars(iflag)
C Copy new parameter values from MINUIT to whereever parameterisations
c will take them from
      call copy_minuit_extrapars(parminuit)

#ifdef TRACE_CHISQ
      call MntInpGetparams ! calls MInput.GetMinuitParams();
#endif

C Count the free POIs (always: that count is the c_CI denominator and the GoF
C dof, whether or not Bartlett is enabled), then take the theory derivatives
C w.r.t. them. Both at the MLE, BEFORE the chi2 call, so the Bartlett factors
C are available at iflag=3 for Results.txt, XRANGE and the error bands.
C Bartlett_ComputeD restores the parameters and the theory before returning.
      if (iflag .eq. 3) then
         call Bartlett_CountFreePOI
         call Bartlett_ComputeD(parminuit)
      endif

*
* Evaluate the chi2:
*
      chi2out = chi2data_theory(iflag)

#ifdef TRACE_CHISQ
      if (iflag.eq.1) then
        ! print *,'INIT'
        call MntInpWritepar('minuit.all_in.txt')
        call MntShowVNames(ndfMINI)
      endif
      call MntShowVValues(chi2out)
#endif
      call flush
      return
      end
C------------------------------------------------------
C> @brief Helper for C++
C------------------------------------------------------
      subroutine update_theory_iteration
      implicit none
#include "ntot.inc"
#include "datasets.inc"
      integer idataset
      character*128 Msg

      call init_at_iteration
      do idataset=1,Ndatasets
         if(NDATAPOINTS(idataset).gt.0) then
            call GetTheoryForDataset(idataset)
         else
             write (Msg,
     $           '(''W: Data set '',i2
     $,'' contains no data points, will be ignored'')')
     $           idataset
           call hf_errlog(29052013,Msg)
         endif
      enddo

      end

C------------------------------------------------------------------------------
C> @brief     Calculate predictions for the data samples and return total chi2.
C> @details   Created by splitting original fcn() function
C> @param[in] iflag minuit flag indicating minimisation stage
C> @authors   Sasha Glazov, Voica Radescu
C> @date      22.03.2012
C------------------------------------------------------------------------------
      double precision function chi2data_theory(iflag)

      implicit none
C--------------------------------------------------------------
      integer iflag

#include "steering.inc"
#include "for_debug.inc"
#include "ntot.inc"
#include "datasets.inc"
#include "systematics.inc"
#include "theo.inc"
#include "indata.inc"
#include "thresholds.inc"
#include "fcn.inc"
#include "polarity.inc"
#include "endmini.inc"
#include "fractal.inc"
#include "bartlett_ci.inc"
#include "bartlett_fd.inc"

*     ---------------------------------------------------------
*     declaration related to chisquare
*     ---------------------------------------------------------
      double precision chi2out
      double precision fchi2, fcorchi2
!     double precision DeltaLength
      double precision BSYS(NSYSMax), RSYS(NSYSMax)
      double precision EBSYS(NSYSMax),ERSYS(NSYSMax)
      double precision pchi2(nset),chi2_log
      double precision pchi2offs(nset)
      
      double precision chi2out_print, chi2off_print, fcorchi2_print
      double precision chi2_log_print, pchi2_print, chi2_poi_print
      double precision c_bart_chi2, c_bart_ci
      double precision ersys_raw, bart_weight, bart_weight_lr
      double precision eps_i, b_theta_i, btil_i, r_i

*     ---------------------------------------------------------
*     declaration related to code flow/debug
*     ---------------------------------------------------------

      double precision time1,time2,time3
      logical od !to check if a unit is open with INQUIRE

*     ---------------------------------------------------------
*     declaration related to others
*     ---------------------------------------------------------
      double precision x
      double precision quv,qdv, qus, qds, qst, qch, qbt, qgl
      integer iq, ix, nndi, ndi,ndi2
      character*300 base_pdfname
      integer npts(nset)
      integer nPOI, nExtSyst, ndf_bart
      double precision f2SM,f1SM,flSM
      integer i,j,jsys,ndf,n0,h1iset,jflag,k,pr,nwds
      logical refresh
      integer isys,ipoint,jpoint
      integer idataset
      double precision TempChi2
      double precision GetTempChi2   ! Temperature penalty for D, E... params.
      double precision OffsDchi2   ! correction for final Offset calculation

C  x-dependent fs:
      double precision fs0
      double precision fshermes

c updf stuff
      logical firsth
      double precision auh
      common/f2fit/auh(50),firsth
      Logical Firstd,Fccfm1,Fccfm2
      Common/ myfirst/Firstd,Fccfm1,Fccfm2
      Integer IGLU
      Common/CAGLUON/Iglu
      character filename*132
      Common/updfout/filename
      character CCFMfile*132
      Common/CCFMout/CCFMfile
      CHARACTER   evolfNAME*132
      Common/gludatf/evolfname
      Integer idx

      character*2 TypeC, FormC, TypeD
      character*64 Msg

      double precision rmass,rmassp,rcharge
      COMMON /MASSES/ rmass(150),rmassp(50),rcharge(150)

C Penalty from MINUIT extra parameters constraints
      double precision extraparsconstrchi2
C---------------------------------------------------

      ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx to be removed
      if (kmuc .eq. 0) then
         kmuc = 1.
         kmub = 1.
         kmut = 1.
      endif

C--OZ 21.04.2016 Increment IfcnCount here instead of fcn routine
      IfcnCount=IfcnCount+1
      if (lprint) then
        write(6,*) ' ===========  Calls to fcn= IfcnCount ',IfcnCount
      endif
C--------------------------------------------------------------
*     ---------------------------------------------------------
*     initialise variables
*     ---------------------------------------------------------
      chi2out = 0.d0
      fchi2 = 0.d0

      iflagfcn = iflag

      do jsys=1,nsys
         bsys(jsys) = 0.d0
         rsys(jsys) = 0.d0
         ebsys(jsys) = 0.d0
         ersys(jsys) = 0.d0
      enddo


      do i=1,ntot !why for both used and unused points? --Ivan
         THEO(i) = 0.d0
         THEO_MOD(i) = 0.d0
      enddo ! on second thought, why clear these anyway?

      if (iflag.eq.1) then
         open(87,file=TRIM(OutDirName)//'/pulls.first.txt')
      endif

      if ((iflag.eq.3).or.(iflag.eq.4)) then
         open(88,file=TRIM(OutDirName)//'/pulls.last.txt')
         do i=1,nset
            npts(i) = 0
         enddo
      endif

*     ---------------------------------------------------------
*     Initialise various c++ code per iteration
*     ---------------------------------------------------------
      call init_at_iteration
*     ---------------------------------------------------------
*     Calculate theory for datasets:
*     ---------------------------------------------------------
      do idataset=1,NDATASETS
         if(NDATAPOINTS(idataset).gt.0) then
            call GetTheoryForDataset(idataset)
         else
             write (Msg,
     $  '(''W: Data set '',i2,
     $           '' contains no data points, will be ignored'')') idataset
           call hf_errlog(29052013,Msg)
         endif
      enddo

      if (Debug) then
         print*,'after GetTheoryfordataset'
      endif

      call cpu_time(time1)

      !Count datapoints in each dataset?
      if(iflag.eq.3)then
        do i=1,npoints
          h1iset = JSET(i)
          npts(h1iset)=npts(h1iset)+1
        enddo
      endif
      ndf=npoints-nparFCN !degrees of freedom
      n0 =npoints
      if(iflag.eq.1)then !at first iteration
        if(lrand.and.DataToTheo)then
          call hf_errlog(2019033100,
     $'F: DataToTheo and MCError cannot be used at the same time')
        endif
        if(lrand)then !MCErrors
          call MC_Method()
        endif

        NSysData = 0
        do i=1,NSys
          if(ISystType(i).eq.iDataSyst)then
            NSysData=NSysData+1
          endif
        enddo

        if(DataToTheo)then !Copy theory to data
          do i=1,npoints
            daten(i)=theo(i)
            !Update total uncorrelated uncertainty
            alpha(i)=daten(i)*sqrt(
     &      e_stat_poisson(i)**2+
     &      e_stat_const(i)**2+   !Should I use e_stat_const or e_sta_const or e_sta? I am not sure... --Ivan
     &      e_uncor_poisson(i)**2+
     &      e_uncor_const(i)**2+  !or e_unc_const?
     &      e_uncor_mult(i)**2+
     &      e_uncor_logNorm(i)**2)
          enddo
        endif
      endif

*     ---------------------------------------------------------
*     calculate chisquare
*     ---------------------------------------------------------
      OffsDchi2 = 0.d0
      if (doOffset .and. iflag.eq.3) then
        Chi2OffsRecalc = .true.
        Chi2OffsFinal = .true.
        call GetNewChisquare(iflag,n0,OffsDchi2,rsys,ersys,
     $       pchi2offs,fcorchi2)
      else
        Chi2OffsRecalc = .false.
      endif
      Chi2OffsFinal = .false.
      call GetNewChisquare(iflag,n0,fchi2,rsys,ersys,pchi2,fcorchi2)
      if (doOffset .and. iflag.eq.3) then
        Chi2OffsRecalc = .false.
        OffsDchi2 = OffsDchi2 - fchi2
      endif

      if (ControlFitSplit) then
         print '(''Fit     chi2/Npoint = '',F10.4,I4,F10.4)',chi2_fit
     $        , NFitPoints,chi2_fit/NFitPoints
         print '(''Control chi2/Npoint = '',F10.4,I4,F10.4)',chi2_cont
     $        , NControlPoints,chi2_cont/NControlPoints
      endif

*
* Save NDF
*
      ndfMINI = ndf

      if (iflag.eq.1) close(87)

      if (iflag.eq.3) then
C Broken since 2.2.0
!        if (dobands) then
!           print *,'SAVE PDF values'
!        endif

         TheoFCN3 = Theo  ! save
         TheoModFCN3 = Theo_Mod
         ALphaModFCN3 = ALPHA_Mod

*     ---------------------------------------------------------
*     write out data points with fitted theory
*     ---------------------------------------------------------
         call writefittedpoints

*        --------------------------------
*        Temporary output for HERAverager
*        --------------------------------
         if (Debug) then
           call WriteCSforAverager
         endif

      endif

C Broken since 2.2.0
! Temperature regularisation:
c     if (Temperature.ne.0) then
c        TempChi2 = GetTempChi2()
c        print *,'Temperature chi2=',TempChi2
c        fchi2 = fchi2 + TempChi2
c     endif


c Penalty from MINUIT extra parameters constraints (only for fits)
C However when/if LHAPDFErrors mode will be combined with minuit, this will need modification.
      if (.not. LHAPDFErrors) then
         call getextraparsconstrchi2(extraparsconstrchi2)
         fchi2 = fchi2 + extraparsconstrchi2
      endif

      chi2out = fchi2+
     $     shift_polRHp**2+shift_polRHm**2+
     $     shift_polLHp**2+shift_polLHm**2+
     $     shift_polL**2+shift_polT**2
c If for any reason we got chi2==NaN, set it to +inf so that that
c a minimizer would treat it as very bad
      if(chi2out/=chi2out)then !if chi2out is NaN
        chi2out=1e10 !set it to a very large (but finite) number, so
        !that the minimizer would move away from this point
        !We used to use +Infinity, but that breaks MINUIT
      endif

c Print time, number of calls, chi2
         call cpu_time(time3)
         print '(''cpu_time'',3F10.2)', time1, time3, time3-time1
         write(6,'(A20,i6,F12.2,i6,F12.2)') '
     $        xfitter chi2out,ndf,chi2out/ndf ',ifcncount, chi2out,
     $        ndf, chi2out/ndf
      call flush

! ----------------  RESULTS OUTPUT ---------------------------------
! Reopen "Results.txt" file if it is not open
! It does not get opened by this point when using CERES
      INQUIRE(85,OPENED=od)
      if(.not. od)then
        call IOFileNamesMini()
      endif
      if (iflag.eq.1) then
         write(85,*) 'First iteration ',chi2out,ndf,chi2out/ndf
         call flush
      endif

      if (iflag.eq.3) then

         ! ----- Bartlett-only-for-printing scaling -----
         nExtSyst = 0
         do i=1,nsys
            if ( SysForm(i) .eq. isExternal) then
               nExtSyst = nExtSyst + 1
            endif
         enddo

         ! nparFCN - nExtSyst = freePOI - fixedExt, so it under-counts the POIs
         ! whenever an external systematic is a FIXED MINUIT parameter. Use the
         ! free-POI count that Bartlett_CountFreePOI actually established.
         nPOI = nparFCN - nExtSyst
         if (BartlettHaveNPOI) nPOI = BartlettNPOI
         print *, 'External systs = ', nExtSyst
         c_bart_chi2 = 1.0D0
         c_bart_ci = 1.0D0
         if (BartlettEnabled .and. EoEEnabled) then
            ! External NPs are penalty-constrained: each adds one unit to
            ! E[chi2], so the dof of the GoF statistic is npoints - nPOI
            ! (Eq. 38 of arXiv:2407.05322). The ndf+nExtSyst form is only
            ! equivalent when every external NP is free.
            ndf_bart = ndf + nExtSyst
            if (BartlettHaveNPOI) ndf_bart = npoints - nPOI
            if (ndf_bart .gt. 0) then
               c_bart_chi2 = 1.0D0 / ( 1.0D0 + BartlettGoFFactor / dble(ndf_bart) )
               ! b_q >= 0 for r in [0,1], so c_bart_chi2 > 1 signals a
               ! pathological per-source factor (J_ss > 2*sigma_u_hat^2)
               if (c_bart_chi2 .gt. 1.0D0) then
                  call hf_errlog(25100101,
     $'W: c_bart_chi2 > 1; mathematically should be <=1, check EoE inputs or set Enable_Bartlett = .false.')
               endif
            endif
            if (nPOI .gt. 0) then
               c_bart_ci = sqrt(1.0D0 + BartlettLRFactor / dble(nPOI))
            endif
         endif
C        Store CI Bartlett factor for use by XRANGE (minuit.out.txt / xfitter-draw)
         xfitter_bart_ci = c_bart_ci

         ! Print RAW chi2 (unscaled minimization output).  Preserve the
         ! historical heading for ordinary fits; EoE output keeps the
         ! explicit Chi2 prefix used by the EoE result parsers.
         if (EoEEnabled) then
            write(85,'(''Chi2 after minimisation      '',F10.2,I6,F10.3)')
     $           chi2out,ndf,chi2out/ndf
         else
            write(85,'(''After minimisation '',F10.2,I6,F10.3)')
     $           chi2out,ndf,chi2out/ndf
         endif
         if (doOffset) then
            write(85,'(''  Offset corrected (raw)     '',F10.2,I6,F10.3)')
     $           chi2out+OffsDchi2,ndf,(chi2out+OffsDchi2)/ndf
         endif

         ! Print Bartlett-corrected chi2 if enabled
         if (BartlettEnabled .and. EoEEnabled) then
            chi2out_print = chi2out * c_bart_chi2
            write(85,'(''Corrected Chi2               '',F10.2,I6,F10.3)')
     $           chi2out_print,ndf,chi2out_print/ndf
            if (doOffset) then
               chi2off_print = (chi2out + OffsDchi2) * c_bart_chi2
               write(85,'(''  Offset corrected (Bartlett)'',F10.2,I6,F10.3)')
     $              chi2off_print,ndf,chi2off_print/ndf
            endif
         else
            chi2out_print = chi2out
            if (doOffset) chi2off_print = chi2out + OffsDchi2
         endif
         write(85,*)

         if (BartlettEnabled .and. EoEEnabled) then
            write(85,'(A,F12.6)') 'Chi2 Bartlett factor (rescales the chi2) = ', 
     $        c_bart_chi2
            write(85,'(A,F12.6)') 'Confidence Intervals Bartlett factor (rescales PDF uncertainties) = ', 
     $        c_bart_ci
         endif

         write(6,*)
         if (EoEEnabled) then
            write(6,'(''Chi2 after minimisation      '',F10.2,I6,F10.3)')
     $           chi2out,ndf,chi2out/ndf
         else
            write(6,'(''After minimisation '',F10.2,I6,F10.3)')
     $           chi2out,ndf,chi2out/ndf
         endif
         if (doOffset) then
            write(6,'(''  Offset corrected (raw)     '',F10.2,I6,F10.3)')
     $           chi2out+OffsDchi2,ndf,(chi2out+OffsDchi2)/ndf
         endif
         if (BartlettEnabled .and. EoEEnabled) then
            write(6,'(''Corrected Chi2               '',F10.2,I6,F10.3)')
     $           chi2out_print,ndf,chi2out_print/ndf
            if (doOffset) then
               write(6,'(''  Offset corrected (Bartlett)'',F10.2,I6,F10.3)')
     $              chi2off_print,ndf,chi2off_print/ndf
            endif
            write(6,*)
            write(6,'(A,F12.6)')
     $        'Chi2 Bartlett factor (rescales the chi2) = ', c_bart_chi2
            write(6,'(A,F12.6)')
     $        'CI Bartlett factor (rescales PDF uncertainties) = ',c_bart_ci
         endif
         write(6,*)

! ----------------  END OF RESULTS OUTPUT ---------------------------------

         ! Store minuit parameters
         call write_pars(nfcn3)

         if (ControlFitSplit) then
            print
     $     '(''Fit     chi2/Npoint, after fit = '',F10.4,I4,F10.4)'
     $           ,chi2_fit
     $           , NFitPoints,chi2_fit/NFitPoints
            print
     $     '(''Control chi2/Npoint, after fit = '',F10.4,I4,F10.4)'
     $           ,chi2_cont
     $           , NControlPoints,chi2_cont/NControlPoints

c            write (71,'(4F10.4)')
c     $           paruval(4),paruval(5),chi2_fit/NFitPoints
c     $           ,chi2_cont/NControlPoints

            ! Store chi2 per fcn3 call values:
            chi2cont3(nfcn3) = chi2_cont
            chi2fit3(nfcn3)  = chi2_fit
         endif

      endif



      if (iflag.eq.3) then

         if (doOffset) then
            fcorchi2 = 0d0
            do h1iset=1,nset
               pchi2(h1iset) = pchi2offs(h1iset)
               fcorchi2 = fcorchi2 + pchi2offs(h1iset)
            enddo
!     fcorchi2 = chi2out+OffsDchi2
         endif

! ----------------  RESULTS OUTPUT ---------------------------------
         write(85,*) ' Partial chi2s '
         chi2_log = 0
         do h1iset=1,nset
            if ( Chi2PoissonCorr ) then
               if (npts(h1iset).gt.0) then
                  chi2_log = chi2_log + chi2_poi(h1iset)
                  
                  pchi2_print    = pchi2(h1iset)    * c_bart_chi2
                  chi2_poi_print = chi2_poi(h1iset) * c_bart_chi2

                  write(6,'(''Dataset '',i4,F10.2,''('',SP,F6.2,SS,'')'',
     $                 i6,''  '',A48)')
     $               h1iset, pchi2_print, chi2_poi_print, npts(h1iset),
     $               datasetlabel(h1iset)

                  write(85,'(''Dataset '',i4,F10.2,''('',SP,F6.2,SS,'')'',
     $                 i6,''  '',A48)')
     $               h1iset, pchi2_print, chi2_poi_print, npts(h1iset),
     $               datasetlabel(h1iset)
               endif
            else
               if (npts(h1iset).gt.0) then
                  pchi2_print = pchi2(h1iset) * c_bart_chi2
                  write(6,'(''Dataset '',i4,F10.2,i6,''  '',A48)')
     $               h1iset, pchi2_print, npts(h1iset), datasetlabel(h1iset)
                  write(85,'(''Dataset '',i4,F10.2,i6,''  '',A48)')
     $               h1iset, pchi2_print, npts(h1iset), datasetlabel(h1iset)
               endif
            endif
         enddo
         write(85,*)

         ! Print raw Correlated Chi2
         write(85,*) 'Correlated Chi2 ', fcorchi2
         write(6,*)  'Correlated Chi2 ', fcorchi2

         ! Print Bartlett-corrected Correlated Chi2 if enabled
         if (BartlettEnabled .and. EoEEnabled) then
            fcorchi2_print = fcorchi2 * c_bart_chi2
            write(85,*) 'Corrected Correlated Chi2 ', fcorchi2_print
            write(6,*)  'Corrected Correlated Chi2 ', fcorchi2_print
         endif
! ----------------  END OF RESULTS OUTPUT ---------------------------------

         if (Chi2PoissonCorr) then
            ! Print raw Log penalty Chi2
            write(85,*) 'Log penalty Chi2 ', chi2_log
            write(6,*)  'Log penalty Chi2 ', chi2_log

            ! Print Bartlett-corrected Log penalty Chi2 if enabled
            if (BartlettEnabled .and. EoEEnabled) then
               chi2_log_print = chi2_log * c_bart_chi2
               write(85,*) 'Corrected Log penalty Chi2 ',chi2_log_print
               write(6,*)  'Corrected Log penalty Chi2 ',chi2_log_print
            endif
         endif

         base_pdfname = TRIM(OutDirName)//'/pdfs_q2val_'
         if (CorSysIndex.eq.0) then
            open (76,file=TRIM(OutDirName)//'/lhapdf.block.txt',status='unknown')

            call store_pdfs(base_pdfname)
            call print_lhapdf6
         endif

c WS: print NSYS --- needed for batch Offset runs
         write(85,*) 'Systematic shifts ',NSYS
         write(85,*) ' '
         if (BartlettEnabled .and. EoEEnabled) then
            ! Extended header with all EoE/Bartlett columns
            write(85,'(A5,1X,A55,1X,A9,1X,A9,1X,A9,1X,A7,1X,A9,1X,A10,
     $           1X,A10,1X,A12,1X,A12,1X,A6)')
     $        'Index','Name','Shift','Error','Corr Err','Epsilon',
     $        'r','b_tilde','b_mutheta','GoF Contrib','LR Contrib','Type'
         else
            ! Preserve the historical non-EoE header.
            write(85,'(A5,'' '',A55,'' '',A9,''   +/-'',A9,A10,A4)')
     $        ' ', 'Name     ', 'Shift','Error',' ','Type'
         endif
         do jsys=1,nsys
C     !> Store also type of systematic source info
            if ( SysForm(jsys) .eq. isNuisance ) then
               FormC = ':N'
            elseif ( SysForm(jsys) .eq. isMatrix ) then
               FormC = ':C'
            elseif ( SysForm(jsys) .eq. isOffset ) then
               FormC = ':O'
            elseif ( SysForm(jsys) .eq. isExternal ) then
               FormC = ':E'
            endif

            if ( SysScalingType(jsys) .eq. isPoisson ) then
               TypeC = ':P'
            elseif ( SysScalingType(jsys) .eq. isNoRescale ) then
               TypeC = ':A'
            elseif ( SysScalingType(jsys) .eq. isLinear ) then
               TypeC = ':M'
            endif

            if (ISystType(jsys).eq. iDataSyst) then
               TypeD = ':D'
            elseif (ISystType(jsys).eq. iTheorySyst ) then
               TypeD = ':T'
            endif

            if (BartlettEnabled .and. EoEEnabled) then
               ! Extended format with EoE columns
               eps_i     = EoEEpsilon(jsys)
               btil_i    = BartlettBTilde(jsys)     ! b~_theta   (Eq. 35)
               b_theta_i = BartlettSysFactor(jsys)  ! b_mutheta  (Eq. 34)
               r_i       = BartlettRatio(jsys)      ! j_ss / sigma_u_hat^2

               ! Raw error: divide out the sqrt(1+b~_theta) inflation
               if (1.0D0 + btil_i .gt. 0.0D0) then
                  ersys_raw = ersys(jsys) / sqrt(1.0D0 + btil_i)
               else
                  ersys_raw = ersys(jsys)
               endif

               ! Per-source contributions, one formula for :N and :E alike.
               ! Both are >= 0 by construction (0 <= b~ <= b_mutheta <= 3eps^2).
               ! A FIXED external parameter is excluded from the factors, so
               ! its printed contributions are zero as well.
               bart_weight    = 0.0D0
               bart_weight_lr = 0.0D0
               if (EoEActive(jsys) .and.
     $             (SysForm(jsys).eq.isNuisance .or.
     $              (SysForm(jsys).eq.isExternal .and.
     $               .not. SysExtFixed(jsys)))) then
                  bart_weight    = 3.0D0*eps_i*eps_i - b_theta_i
                  bart_weight_lr = b_theta_i - btil_i
               endif

               write(85,'(I5,1X,A55,1X,F9.4,1X,F9.4,1X,F9.4,1X,F7.4,
     $              1X,F9.5,1X,F10.5,1X,F10.5,1X,F12.5,1X,F12.5,1X,3A2)')
     $           jsys, SYSTEM(jsys), rsys(jsys), ersys_raw, ersys(jsys),
     $           eps_i, r_i, btil_i, b_theta_i,
     $           bart_weight, bart_weight_lr,
     $           FormC, TypeC, TypeD

            else
               ! Preserve the historical non-EoE row format.
               write(85,'(I5,''  '',A55,'' '',F9.4,''   +/-'',F9.4,
     $              A8,3A2)')
     $           jsys, SYSTEM(jsys), rsys(jsys), ersys(jsys),
     $           ' ', FormC, TypeC, TypeD
            endif
         enddo

C Print legend for Bartlett corrections
         if (BartlettEnabled .and. EoEEnabled) then
            write(85,*)
            write(85,'(A)') '--- Bartlett Correction Legend ---'
            write(85,'(A)')
     $        'Error     : Raw uncertainty sqrt(j~_ss) from the extended NP Hessian'
            write(85,'(A)')
     $        '            (quadratic convention, same statistic for :N and :E)'
            write(85,'(A)')
     $        'Corr Err  : Corrected = Error * sqrt(1 + b_tilde)'
            write(85,'(A)')
     $        'Epsilon   : Error-on-error parameter for this source'
            write(85,'(A)')
     $        'r         : j_ss / sigma_u_hat^2, marginal over the POIs. r <= 1 is'
            write(85,'(A)')
     $        '            a theorem; r > 1 flags non-converged EoE profiling.'
            write(85,'(A)')
     $        'b_tilde   : (4r~-r~^2)*eps^2 with r~ = [A^-1]_ss/sigma_u_hat^2  (Eq. 35)'
            write(85,'(A)')
     $        'b_mutheta : (4r -r ^2)*eps^2 with r  = j_ss   /sigma_u_hat^2    (Eq. 34)'
            write(85,'(A)')
     $        '            j = (A - G K^-1 G^T)^-1, K = D^T W D, G = Gamma^T W D,'
            write(85,'(A)')
     $        '            D = d theo/d mu. Same formula for :N and :E sources.'
            write(85,'(A)')
     $        'GoF Contrib: 3*eps^2 - b_mutheta   (Eq. 39), >= 0 by construction'
            write(85,'(A)')
     $        'LR Contrib : b_mutheta - b_tilde   (Eq. 32), >= 0 by construction'
            write(85,'(A)') ''
            write(85,'(A)') 'Global factors:'
            write(85,'(A)')
     $        '  Chi2 Bartlett = 1 / (1 + sum(GoF Contrib)/ndf_bart)'
            write(85,'(A)')
     $        '  CI Bartlett   = sqrt(1 + sum(LR Contrib)/nPOI)'
            write(85,'(A,I6,A,I4)')
     $        '  with ndf_bart = npoints - nPOI = ', ndf_bart,
     $        ' ,  nPOI = ', nPOI
            if (BartlettFailed) then
               write(85,'(A)')
     $ '  WARNING: Bartlett factors UNAVAILABLE (singular POI Hessian'
               write(85,'(A)')
     $ '  or missing theory derivatives); corrections are set to zero.'
            endif
            write(85,'(A)')
     $        '  (Eq. refs: arXiv:2407.05322)'
            write(85,'(A)')
     $        'All factors are exact: no midpoint approximation, no dependence on'
            write(85,'(A)')
     $        'the error-band method, and identical for :N and :E treatments of'
            write(85,'(A)')
     $        'the same physical source.'
            write(85,'(A)') ''
            write(85,'(A)')
     $        'Type: :N=Nuisance :C=Covar :O=Offset :E=External'
            write(85,'(A)')
     $        '      :P=Poisson :A=Additive :M=Mult :D=Data :T=Theory'
            write(85,*)
         endif

C Trigger reactions:
         call fcn3action

         call cpu_time(time2)
         print '(''cpu_time'',3F10.2)', time1, time2, time2-time1


      endif

C Return the chi2 value:
      chi2data_theory = chi2out
      end

C-----------------------------------------------------------------------
C> @brief Central finite differences of the THEORY w.r.t. the free POIs.
C>
C> Fills BartlettD(i,m) = d theo_i / d mu_m at the MLE, used by
C> chi2_calc_syst_shifts to build the POI blocks of the quadratic Hessian
C>    K = D^T W D ,   G = Gamma^T W D
C> and hence the marginal NP covariance j = (A - G K^-1 G^T)^-1 (Eq. 34 of
C> arXiv:2407.05322) identically for external and nuisance sources.
C>
C> Two things are essential here:
C>  * we difference the THEORY, not the chi2. Differencing the chi2 (which
C>    re-profiles the nuisances) would return the Schur complement
C>    K - G A^-1 G^T instead of K, and the construction would collapse.
C>  * external systematic NPs are NOT POIs: theo does not depend on them,
C>    they enter the residual through Gamma*theta. They are skipped.
C>
C> Cost: 2*nPOI+1 theory-only evaluations (no chi2, no profiling).
C>
C> @param pcen central MINUIT parameter vector (the MLE)
C-----------------------------------------------------------------------
      subroutine Bartlett_CountFreePOI

      implicit none
#include "ntot.inc"
#include "endmini.inc"
#include "systematics.inc"
#include "extrapars.inc"
#include "bartlett_fd.inc"

      double precision val, err, xlo, xhi
      integer i, ipar, idx, isys
      integer GetParameterIndex   !function
      character*80 parname
      logical isExtNP(MNE)
C-----------------------------------------------------------------------
      BartlettHaveNPOI = .false.
      BartlettNPOI     = 0

C Mark the MINUIT parameters that are external systematic NPs, and flag the
C FIXED ones: a fixed parameter is a constant, not a fitted NP, so it is
C excluded from the extended NP Hessian and from the Bartlett factors (its
C Gamma*theta contribution to the residual stays, through ShiftExternal), and
C its reported error remains MINUIT's (zero for a fixed parameter).
      do i=1,MNE
         isExtNP(i) = .false.
      enddo
      do isys=1,nsys
         SysExtFixed(isys) = .false.
      enddo
      do isys=1,nsys
         if ( SysForm(isys) .eq. isExternal ) then
            idx = GetParameterIndex(system(isys))
            if (idx .gt. 0) then
               i = iExtraParamMinuit(idx)
               if (i.ge.1 .and. i.le.MNE) then
                  isExtNP(i) = .true.
                  call mnpout(i, parname, val, err, xlo, xhi, ipar)
                  SysExtFixed(isys) = ( ipar .le. 0 )
                  if ( ipar .le. 0 .and. EoEActive(isys) ) then
                     call hf_errlog(25100109,
     $'I: Bartlett: EoE-active external systematic is FIXED in MINUIT; excluded from the Bartlett factors')
                  endif
               endif
            endif
         endif
      enddo

C Free POIs = variable MINUIT parameters that are not external NPs. DEPENDENT
C parameters are not MINUIT parameters; the sum rules resolve them inside
C init_at_iteration.
      do i=1,MNE
         call mnpout(i, parname, val, err, xlo, xhi, ipar)
         if ( ipar .gt. 0 .and. .not. isExtNP(i) ) then
            if (BartlettNPOI .ge. BartMaxPOI) then
               call hf_errlog(25100105,
     $'W: Bartlett: nPOI exceeds BartMaxPOI, exact Bartlett factors disabled')
               BartlettNPOI     = 0
               BartlettHaveNPOI = .false.
               return
            endif
            BartlettNPOI = BartlettNPOI + 1
            BartlettPOIMinuit(BartlettNPOI) = i
         endif
      enddo

      BartlettHaveNPOI = .true.

      end

C-----------------------------------------------------------------------
C> @brief Central finite differences of the THEORY w.r.t. the free POIs.
C> Requires Bartlett_CountFreePOI to have run.
C-----------------------------------------------------------------------
      subroutine Bartlett_ComputeD(pcen)

      implicit none
#include "ntot.inc"
#include "endmini.inc"
#include "systematics.inc"
#include "indata.inc"
#include "theo.inc"
#include "bartlett_fd.inc"

      double precision pcen(*)
      double precision a(MNE)
      double precision, allocatable :: tp(:)
      double precision h, val, err, xlo, xhi
      integer i, m, ipoi, ipar
      character*80 parname
C-----------------------------------------------------------------------
      BartlettHaveD = .false.

      if (.not. (BartlettEnabled .and. EoEEnabled)) return
      if (.not. BartlettHaveNPOI) return

C nPOI = 0 is a legitimate limit: G = 0, so j = A^-1 = j-tilde and every LR
C contribution vanishes. Nothing to difference.
      if (BartlettNPOI .eq. 0) then
         BartlettHaveD = .true.
         return
      endif

      allocate(tp(NTOT))

      do ipoi = 1, BartlettNPOI
         m = BartlettPOIMinuit(ipoi)
         call mnpout(m, parname, val, err, xlo, xhi, ipar)
         h = 0.1D0 * abs(err)
         if (h .le. 0.0D0) h = 1.0D-4 * max( abs(pcen(m)), 1.0D0 )

         do i=1,MNE
            a(i) = pcen(i)
         enddo

C set_scan_parameters keeps ExtraParamValue and parminuitsave consistent. Only
C the POI entries move here, so the external-NP entries that
C Chi2_calc_readExternal reads afterwards are untouched; but a POI read through
C XParValueByName (which goes via parminuitsave) would otherwise silently miss
C the displacement and corrupt this column of D.
         a(m) = pcen(m) + h
         call set_scan_parameters(a)
         call update_theory_iteration
         do i=1,npoints
            tp(i) = theo(i)
         enddo

         a(m) = pcen(m) - h
         call set_scan_parameters(a)
         call update_theory_iteration
         do i=1,npoints
            BartlettD(i,ipoi) = ( tp(i) - theo(i) ) / ( 2.0D0*h )
         enddo
      enddo

      deallocate(tp)

C Restore the central parameters (both commons) and the theory/evolution grids.
      do i=1,MNE
         a(i) = pcen(i)
      enddo
      call set_scan_parameters(a)
      call update_theory_iteration

      BartlettHaveD = .true.

      end

C copy parameters from minuit to wherever c++ components will read them from
C @param[in] p - vector of parameter values (I am not sure in what order)
C this should be called whenever minuit parameters change
C this replaces old subroutine PDF_param_iteration
      subroutine copy_minuit_extrapars(p)
      implicit none
      double precision p(*)
#include "extrapars.inc"
      integer i
      double precision getfittedparamfromarrayd
      do i=1,nExtraParam
        ExtraParamValue(i)=getfittedparamfromarrayd
     $       (ExtraParamNames(i),p)
      enddo
      end

C
C Set scan parameters consistently in BOTH commons read by the chi2:
C ExtraParamValue (via copy_minuit_extrapars) and parminuitsave, which
C Chi2_calc_readExternal and XParValueByName read for external systematics.
C parminuitsave is normally refreshed inside fcn(), which a direct call to
C chi2data_theory bypasses — without this, scans over displaced parameters
C are blind to the external-NP components of the displacement.
C
      subroutine set_scan_parameters(p)
      implicit none
      double precision p(*)
#include "endmini.inc"
      integer i
      do i=1,MNE
         parminuitsave(i) = p(i)
      enddo
      call copy_minuit_extrapars(p)
      end

C
C Reset extra parameters
C
      subroutine reset_extra_parameters
      implicit none
#include "extrapars.inc"
      nextraparam = 0
      end
