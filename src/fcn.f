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

C function:
      double precision chi2data_theory
! [--- WS 2015-10-10
      double precision XParValueByName
! ---]

C-----------------------------------------------------------------
      
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

C Print MINUIT extra parameters
      call printminuitextrapars

*     ---------------------------------------------------------
*     PDF parameterisation at the starting scale
*     ---------------------------------------------------------

      call PDF_Param_Iteration(parminuit,iflag)

! [--- WS 2015-10-10
#ifdef TRACE_CHISQ
      call MntInpGetparams ! calls MInput.GetMinuitParams();
#endif
!       if(doHiTwist) then
!         print '(''HiTwist:'',6F10.4)',XParValueByName('HT_x0'), XParValueByName('HT_sig0'), XParValueByName('HT_lambda') 
!       endif
! ---]
      
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
#include "pdfparam.inc"
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
      double precision DeltaLength
      double precision BSYS(NSYSMax), RSYS(NSYSMax)
      double precision EBSYS(NSYSMax),ERSYS(NSYSMax)
      double precision pchi2(nset),chi2_log
      double precision pchi2offs(nset)

*     ---------------------------------------------------------
*     declaration related to code flow/debug
*     ---------------------------------------------------------

      double precision time1,time2,time3


*     ---------------------------------------------------------
*     declaration related to others
*     ---------------------------------------------------------
      double precision x
      double precision quv,qdv, qus, qds, qst, qch, qbt, qgl
      integer iq, ix, nndi, ndi,ndi2
      character*300 base_pdfname
      integer npts(nset)
      double precision f2SM,f1SM,flSM
      integer i,j,kflag,jsys,ndf,n0,h1iset,jflag,k,pr,nwds
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
      Integer Itheory_ca
      Common/theory/Itheory_ca
      Integer idx

      character*2 TypeC, FormC, TypeD
      character*64 Msg
      
      double precision rmass,rmassp,rcharge
      COMMON /MASSES/ rmass(150),rmassp(50),rcharge(150)

C Penalty from MINUIT extra parameters constraints
      double precision extraparsconstrchi2


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
      ndf = -nparFCN
      n0 = 0

      iflagfcn = iflag

      itheory_ca = itheory      

      do jsys=1,nsys
         bsys(jsys) = 0.d0
         rsys(jsys) = 0.d0
         ebsys(jsys) = 0.d0
         ersys(jsys) = 0.d0
      enddo


      do i=1,ntot
         THEO(i) = 0.d0
         THEO_MOD(i) = 0.d0
      enddo



      if (Itheory.eq.3) then
        Itheory = 0
      endif

      if(Itheory.ge.100) then
        auh(1) = parminuitsave(1)
        auh(2) = parminuitsave(2)
        auh(3) = parminuitsave(3)
        auh(4) = parminuitsave(4)
        auh(5) = parminuitsave(5)
        auh(6) = parminuitsave(6)
        auh(7) = parminuitsave(7)
        auh(8) = parminuitsave(8)
        auh(9) = parminuitsave(9)
c        write(6,*) ' fcn npoint ',npoints
         firsth=.true.
         Fccfm1=.true.
        
      endif

*     ---------------------------------------------------------
*     Extra constraints on input PDF due to momentum and quark 
*     counting sum rules:
*     ---------------------------------------------------------

      kflag=0
      if (Itheory.eq.0.or.Itheory.eq.10.or.itheory.eq.11
     $.or.itheory.eq.35)  then 
         call SumRules(kflag)
      endif
      if (kflag.eq.1) then
         write(6,*) ' --- problem in SumRules, kflag = 1'
         call HF_errlog(12020516,
     +        'F: FCN - problem in SumRules, kflag = 1')
      endif

      if (iflag.eq.1) then
         open(87,file=TRIM(OutDirName)//'/pulls.first.txt')
      endif

      if ((iflag.eq.3).or.(iflag.eq.4)) then
         open(88,file=TRIM(OutDirName)//'/pulls.last.txt')
         do i=1,nset
            npts(i) = 0
         enddo
      endif

C
C     Call a subrotine which vanishes nonvalence DGLAP contribution
C     for dipole model fits.
C
      if (DipoleModel.eq.3.or.DipoleModel.eq.4) then
         call LeaveOnlyValenceQuarks
      endif

*     ---------------------------------------------------------  	 
*     Call evolution
*     ---------------------------------------------------------  	 
      if (Debug) then
         print*,'before evolution'
      endif
      if (itheory.eq.0.or.itheory.eq.10.or.itheory.eq.11.or.
     1    itheory.eq.35) then         
         call Evolution
      elseif(Itheory.ge.100) then
          if(itheory.eq.101) then 
             Iglu=1112 
             filename='ccfm-test.dat'
             call calcglu
c iglu is set in sigcalc to iglu=1111
          elseif(itheory.eq.102) then 
c iglu in sigcalc is set by steer-ep 
          elseif(itheory.eq.103) then 
c iglu is set in sigcalc to iglu=1113
             Iglu = 1113
          elseif(itheory.eq.104) then 
             Iglu=1112            
c call evolution to generate grid file
             evolfname='Cascade/updfGrids/ccfm-grid.dat'
             if(iflag.eq.1) call evolve_tmd
c iglu is set in sigcalc to iglu=1111
          elseif(itheory.eq.105) then 
c call evolution to generate grid file
             evolfname='ccfm-test.dat'
             if(iflag.eq.1) call evolve_tmd
c iglu is set in sigcalc to iglu=1111
          endif 
          firsth=.false.
      endif

c fill PDFs for ABKM scheme
      if (mod(HFSCHEME,10).eq.4) then 
c update heavy-quark masses (for FF ABM and FF ABM RUNM schemes)
        if (HFSCHEME.eq.4.or.HFSCHEME.eq.444) then
          call UpdateHQMasses()
          rmass(8)  = HF_MASS(1)
          rmass(10) = HF_MASS(2)
        endif
        call PDFFILLGRID
c this will need to BMSN             
c             call fillvfngrid
      endif

      if (Debug) then
         print*,'after evolution'
      endif


	
*     ---------------------------------------------------------  	 
*     Initialise theory calculation per iteration
*     ---------------------------------------------------------  	 
      call GetTheoryIteration
      if(itheory.ge.103) call GetGridkt


      if (Debug) then
         print*,'after GetTheoryIteration'
      endif
		 
*     ---------------------------------------------------------  	       
*     Calculate theory for datasets:
*     ---------------------------------------------------------  	 
      do idataset=1,NDATASETS
         if(NDATAPOINTS(idataset).gt.0) then
            call GetTheoryForDataset(idataset)
         else 
             write (Msg,
     $  '(''W: Data set '',i2,'' contains no data points, will be ignored'')') idataset
           call hf_errlog(29052013,Msg)
         endif
      enddo

      if (Debug) then
         print*,'after GetTheoryfordataset'
      endif


      call cpu_time(time1)

*     ---------------------------------------------------------  	 
*     Start of loop over data points: 	 
*     ---------------------------------------------------------

      do 100 i=1,npoints
         
         h1iset = JSET(i)
         
         if (iflag.eq.3) npts(h1iset) = npts(h1iset) + 1
         
         n0 = n0 + 1
         ndf = ndf + 1

         
         
 100  continue
*     ---------------------------------------------------------
*     end of data loop
*     ---------------------------------------------------------
      if (Debug) then
         print*,'after data loop'
      endif


* -----------------------------------------------------------
*     Toy MC samples:
* -----------------------------------------------------------

      if (IFlag.eq.1 .and. lrand) then
         call MC_Method()         
      endif 

      if (IFlag.eq.1) then
         NSysData = 0
         do i=1,NSys
            if ( ISystType(i) .eq. iDataSyst) then
               NSysData = NSysData + 1
            endif
         enddo
      endif

      if ( (IFlag.eq.1).and.(DataToTheo)) then
         do i=1,npoints
            daten(i) = theo(i)
         enddo
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
         if (dobands) then
            print *,'SAVE PDF values'
         endif
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

*     ---------------------------------------------------------
*     pdf lenght term -- used for Chebyshev Polynomial
*     ---------------------------------------------------------
      if (ILenPdf.gt.0) then
         call  PDFLength(DeltaLength)      
         print *,'Chi2 from PDF length:',DeltaLength,fchi2
      else
         DeltaLength = 0.
      endif

      fchi2 = fchi2 + DeltaLength
      
! Temperature regularisation:
      if (Temperature.ne.0) then
         TempChi2 = GetTempChi2()
         print *,'Temperature chi2=',TempChi2
         fchi2 = fchi2 + TempChi2
      endif
      
      
c Penalty from MINUIT extra parameters constraints
      call getextraparsconstrchi2(extraparsconstrchi2)
      fchi2 = fchi2 + extraparsconstrchi2


      chi2out = fchi2+
     $     shift_polRHp**2+shift_polRHm**2+
     $     shift_polLHp**2+shift_polLHm**2+
     $     shift_polL**2+shift_polT**2

      


      if (lprint) then
         call cpu_time(time3)
         print '(''cpu_time'',3F10.2)', time1, time3, time3-time1 
         write(6,'(A20,i6,F12.2,i6,F12.2)') '
     $        xfitter chi2out,ndf,chi2out/ndf ',ifcncount, chi2out, 
     $        ndf, chi2out/ndf



      endif ! end  lprint

! ----------------  RESULTS OUTPUT ---------------------------------
      if (iflag.eq.1) then
         write(85,*) 'First iteration ',chi2out,ndf,chi2out/ndf
      endif

      if (iflag.eq.3) then
!          write(85,*),'NFCN3 ',nfcn3
         write(85,'(''After minimisation '',F10.2,I6,F10.3)'),chi2out,ndf,chi2out/ndf
!          if (doOffset .and. iflag.eq.3)
         if (doOffset)
     $    write(85,'(''  Offset corrected '',F10.2,I6,F10.3)'),chi2out+OffsDchi2,ndf,(chi2out+OffsDchi2)/ndf
         write(85,*)
         
         write(6,*)
         write(6,'(''After minimisation '',F10.2,I6,F10.3)'),chi2out,ndf,chi2out/ndf
!          if (doOffset .and. iflag.eq.3)
         if (doOffset)
     $    write(6,'(''  Offset corrected '',F10.2,I6,F10.3)'),chi2out+OffsDchi2,ndf,(chi2out+OffsDchi2)/ndf
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
!          if(itheory.eq.0) then      
!          endif
         
      if (doOffset) then
        fcorchi2 = 0d0
        do h1iset=1,nset
          pchi2(h1iset) = pchi2offs(h1iset)
          fcorchi2 = fcorchi2 + pchi2offs(h1iset)
        enddo
!         fcorchi2 = chi2out+OffsDchi2
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
     $                 ,h1iset,pchi2(h1iset),chi2_poi(h1iset),npts(h1iset)
     $                 ,datasetlabel(h1iset)
                  write(85,'(''Dataset '',i4,F10.2,
     $                 ''('',SP,F6.2,SS,'')'',
     $                 i6,''  '',A48)')
     $                 ,h1iset,pchi2(h1iset),chi2_poi(h1iset),npts(h1iset)
     $                 ,datasetlabel(h1iset)
               endif
            else
               if (npts(h1iset).gt.0) then
                  write(6,'(''Dataset '',i4,F10.2,
     $                 i6,''  '',A48)')
     $                 ,h1iset,pchi2(h1iset)
     $                 ,npts(h1iset)
     $                 ,datasetlabel(h1iset)
                  write(85,'(''Dataset '',i4,F10.2,i6,''  '',A48)')
     $                 ,h1iset,pchi2(h1iset),npts(h1iset)
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

       if (ITheory.ne.2) then

          IF((Itheory.ge.100).or.(Itheory.eq.50)) then
              auh(1) = parminuitsave(1)
              auh(2) = parminuitsave(2)
              auh(3) = parminuitsave(3)
              auh(4) = parminuitsave(4)
              auh(5) = parminuitsave(5)
              auh(6) = parminuitsave(6)
              auh(7) = parminuitsave(7)
              auh(8) = parminuitsave(8)
              auh(9) = parminuitsave(9)
              if(Itheory.ge.103) then
              idx = index(ccfmfile,' ')-1
                 if(idx.gt.2) then
                   filename=ccfmfile(1:idx)//'.dat'
                   firsth=.true.
                   call calcglu
                 endif
              endif         
              open(91,file=TRIM(OutDirName)//'/params.txt')
              write(91,*) auh(1),auh(2),auh(3),auh(4),auh(5),auh(6),auh(7),auh(8),auh(9)

             

          else
C  Hardwire:
             if ( ReadParsFromFile .and. DoBandsSym) then
                call ReadPars(ParsFileName, pkeep)
                call PDF_param_iteration(pkeep,2)
                kflag = 0
                call SumRules(kflag)
             endif
             
             call Evolution

C LHAPDF output:
c WS: for the Offset method save central fit only
            if (CorSysIndex.eq.0) then
              open (76,file=TRIM(OutDirName)//'/lhapdf.block.txt',status='unknown')

              call store_pdfs(base_pdfname)
              call fill_c_common
              call print_lhapdf6
            endif
          endif
       endif

c WS: print NSYS --- needed for batch Offset runs
       write(85,*) 'Systematic shifts ',NSYS
       write(85,*) ' '
       write(85,'(A5,'' '',A35,'' '',A9,''   +/-'',A9,A10,A4)')
     $         ' ', 'Name     ', 'Shift','Error',' ','Type'
       do jsys=1,nsys
C !> Store also type of systematic source info
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
     $         jsys,SYSTEM(jsys),rsys(jsys),ersys(jsys),' ',FormC, 
     $         TypeC,TypeD
       enddo

c AS release applgrids
c AS applgrid example
       if ( useapplg ) then
          call ag_releasegrids
       endif
       call cpu_time(time2)
       print '(''cpu_time'',3F10.2)', time1, time2, time2-time1

       
      endif

      
C Return the chi2 value:
      chi2data_theory = chi2out

      end


C---------------------------------------------------------------------
!> @brief   Calculate penalty term for higher oder parameters using "temperature" 
!> @details Currently works only for standard param-types (10p-13p-like)
C---------------------------------------------------------------------
      double precision function GetTempChi2()

      implicit none
#include "pdfparam.inc"
      integer i
      double precision chi2
      double precision xscale(3)
      data xscale/0.01,0.01,0.01/
      chi2 = 0.

C Over d,e and F
      do i=1,3
         chi2 = chi2 + (paruval(i+3)*xscale(i))**2 
         chi2 = chi2 + (pardval(i+3)*xscale(i))**2 
         chi2 = chi2 + (parubar(i+3)*xscale(i))**2 
         chi2 = chi2 + (pardbar(i+3)*xscale(i))**2 
         if (i.le.2) then
            chi2 = chi2 + (parglue(i+3)*xscale(i))**2 
         endif
      enddo


      GetTempChi2 = chi2*Temperature
C---------------------------------------------------------------------
      end


C---------------------------------------------------------
C Update heavy-quark mass values for FF ABM RUNM or FF ABM
C---------------------------------------------------------
      subroutine UpdateHQMasses
      implicit none
#include "couplings.inc"
#include "steering.inc"
#include "extrapars.inc"
      integer imc,imb,imt
      data imc,imb,imt /0,0,0/
      save imc,imb,imt
      integer GetParameterIndex
      character*256 parname
      double precision par,unc,ll,ul
      integer st
      ! check scheme
      if (HFSCHEME.ne.4.and.HFSCHEME.ne.444) then
        continue
      endif
      ! mc
      if(imc.eq.0) then
        imc=GetParameterIndex('mc')
        if(imc.eq.0) then
          imc=-1
        else
          imc=iExtraParamMinuit(imc)
        endif
      endif
      if(imc.gt.0) then
        call MNPOUT(imc,parname,par,unc,ll,ul,st)
        HF_MASS(1)=par
        mch=par
      endif
      ! mb
      if(imb.eq.0) then
        imb=GetParameterIndex('mb')
        if(imb.eq.0) then
          imb=-1
        else
          imb=iExtraParamMinuit(imb)
        endif
      endif
      if(imb.gt.0) then
        call MNPOUT(imb,parname,par,unc,ll,ul,st)
        HF_MASS(2)=par
        mbt=par
      endif
      ! mt
      if(imt.eq.0) then
        imt=GetParameterIndex('mt')
        if(imt.eq.0) then
          imt=-1
        else
          imt=iExtraParamMinuit(imt)
        endif
      endif
      if(imt.gt.0) then
        call MNPOUT(imt,parname,par,unc,ll,ul,st)
        HF_MASS(3)=par
        mtp=par
      endif
      end
      
