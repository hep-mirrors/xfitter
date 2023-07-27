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

#include "fcn.inc"
#include "endmini.inc"
#include "for_debug.inc"
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
*     ---------------------------------------------------------
      call ACOTindexRESET()
*     ---------------------------------------------------------
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
! ----------------  RESULTS OUTPUT ---------------------------------
! Reopen "Results.txt" file if it is not open
! It does not get opened by this point when using CERES
      INQUIRE(85,OPENED=od)
      if(.not. od)then
        call IOFileNamesMini()
      endif
      if (iflag.eq.1) then
         write(85,*) 'First iteration ',chi2out,ndf,chi2out/ndf
      endif

      if (iflag.eq.3) then
!          write(85,*),'NFCN3 ',nfcn3
         write(85,'(''After minimisation '',F10.2,I6,F10.3)') chi2out,ndf,chi2out/ndf
!          if (doOffset .and. iflag.eq.3)
         if (doOffset)
     $    write(85,'(''  Offset corrected '',F10.2,I6,F10.3)') chi2out+OffsDchi2,ndf,(chi2out+OffsDchi2)/ndf
         write(85,*)

         write(6,*)
         write(6,'(''After minimisation '',F10.2,I6,F10.3)') chi2out,ndf,chi2out/ndf
!          if (doOffset .and. iflag.eq.3)
         if (doOffset)
     $    write(6,'(''  Offset corrected '',F10.2,I6,F10.3)') chi2out+OffsDchi2,ndf,(chi2out+OffsDchi2)/ndf
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
                  write(6,'(''Dataset '',i4,F10.2,''('',SP,F6.2,SS,'')'',
     $                 i6,''  '',A48)')
     $                  h1iset,pchi2(h1iset),chi2_poi(h1iset),npts(h1iset)
     $                 ,datasetlabel(h1iset)
                  write(85,'(''Dataset '',i4,F10.2,
     $                 ''('',SP,F6.2,SS,'')'',
     $                 i6,''  '',A48)')
     $                  h1iset,pchi2(h1iset),chi2_poi(h1iset),npts(h1iset)
     $                 ,datasetlabel(h1iset)
               endif
            else
               if (npts(h1iset).gt.0) then
                  write(6,'(''Dataset '',i4,F10.2,
     $                 i6,''  '',A48)')
     $                  h1iset,pchi2(h1iset)
     $                 ,npts(h1iset)
     $                 ,datasetlabel(h1iset)
                  write(85,'(''Dataset '',i4,F10.2,i6,''  '',A48)')
     $                  h1iset,pchi2(h1iset),npts(h1iset)
     $                 ,datasetlabel(h1iset)
               endif
            endif
         enddo
         write(85,*)
         write(85,*) 'Correlated Chi2 ', fcorchi2
! ----------------  END OF RESULTS OUTPUT ---------------------------------

         write(6,*) 'Correlated Chi2 ', fcorchi2

         if (Chi2PoissonCorr) then
            write(6,*) 'Log penalty Chi2 ', chi2_log
            write(85,*) 'Log penalty Chi2 ', chi2_log
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
         write(85,'(A5,'' '',A35,'' '',A9,''   +/-'',A9,A10,A4)')
     $        ' ', 'Name     ', 'Shift','Error',' ','Type'
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

            write(85,'(I5,''  '',A35,'' '',F9.4,''   +/-'',F9.4,A8,3A2)')
     $           jsys,SYSTEM(jsys),rsys(jsys),ersys(jsys),' ',FormC,
     $           TypeC,TypeD
         enddo

C Trigger reactions:
         call fcn3action

         call cpu_time(time2)
         print '(''cpu_time'',3F10.2)', time1, time2, time2-time1


      endif

C Return the chi2 value:
      chi2data_theory = chi2out
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
      integer GetParameterIndex  !function
      do i=1,nExtraParam
        ExtraParamValue(i)=p(iExtraParamMinuit(GetParameterIndex(trim(ExtraParamNames(i)))))
      enddo
      end
c ==========================================
c ==========================================
c ==========================================
c     THIS IS NOT THE RIGHT PLACE FOR THIS FUNCTION
c ==========================================
c ==========================================
c ==========================================
      subroutine ACOTindexRESET()
      implicit none
      integer index
      common /acotIndex/  index

      index=0
      return
      end
