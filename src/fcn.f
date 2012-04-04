c----------------------------------------------------------------------
      subroutine  fcn(npar,g_dummy,chi2out,parminuit,iflag,futil)
*     ---------------------------------------------------------
*     Main minimization subroutine for MINUIT
*     ---------------------------------------------------------
      implicit none

*     ---------------------------------------------------------
*     declaration related to minuit
*     ---------------------------------------------------------

      integer npar,iflag
      double precision g_dummy(*),parminuit(*),chi2out,futil
      external futil

      include 'fcn.inc'
      include 'endmini.inc'
      
      integer i

C function:
      double precision chi2data_theory

C-----------------------------------------------------------------
      
C Store FCN flag in a common block:
      IFlagFCN = IFlag

      NparFCN  = npar

C Store params in a common block:
      do i=1,MNE
         parminuitsave(i) = parminuit(i)
      enddo

C Store only if IFlag eq 3:
      if (iflag.eq.3) then
         do i=1,MNE
            pkeep(i) = parminuit(i)
         enddo
      endif


      IfcnCount=IfcnCount+1
      write(6,*) ' ===========  Calls to fcn= IfcnCount ',IfcnCount

      call HF_errlog(12020515,'I: FCN is called')

*     ---------------------------------------------------------
*     PDF parameterisation at the starting scale
*     ---------------------------------------------------------
      call PDF_Param_Iteration(parminuit,iflag)

*
* Evaluate the chi2:
*     
      chi2out = chi2data_theory(iflag)
* 
      
      return
      end


      double precision function chi2data_theory(iflag)
C-------------------------------------------------------------
C Created 22/03/2012 by SG and VR by splitting original fcn() function
C
C Calculate predictions for the data samples and return total chi2.
C  Input: minuit flag iflag, indicating minimisation stage
C--------------------------------------------------------------
      implicit none
C--------------------------------------------------------------
      integer iflag

      include 'steering.inc'
      include 'pdfparam.inc'
      include 'alphas.inc'
      include 'for_debug.inc'
      include 'couplings.inc'
      include 'ntot.inc'
      include 'datasets.inc'
      include 'systematics.inc'
      INCLUDE 'theo.inc'
      INCLUDE 'indata.inc'
      INCLUDE 'thresholds.inc'
      include 'fcn.inc'
      include 'polarity.inc'

*     ---------------------------------------------------------
*     declaration related to alphas
*     for RT code, transfer alpha S
*     ---------------------------------------------------------
      double precision alphaszero
      double precision hf_get_alphas

*     ---------------------------------------------------------
*     declaration related to chisquare
*     ---------------------------------------------------------
      double precision chi2out
      double precision fchi2, fcorchi2
      double precision DeltaLength
      double precision BSYS(NSYSMax), RSYS(NSYSMax)
      double precision EBSYS(NSYSMax),ERSYS(NSYSMax)
      double precision pchi2(nset)

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
      character*25 base_pdfname
      integer npts(nset)
      double precision f2SM,f1SM,flSM
      integer i,j,kflag,jsys,ndf,n0,h1iset,jflag,k,pr,nwds
      logical refresh
      integer isys,ipoint,jpoint
      integer idataset
C  x-dependent fs:
      double precision fs0,epsi
      double precision fshermes
      external LHAPDFsubr

C--------------------------------------------------------------
*     ---------------------------------------------------------
*     initilise variables 
*     ---------------------------------------------------------
      chi2out = 0.d0
      fchi2 = 0.d0
      ndf = -nparFCN
      n0 = 0

      iflagfcn = iflag

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

*     ---------------------------------------------------------
*     Extra constraints on input PDF due to momentum and quark 
*     counting sum rules:
*     ---------------------------------------------------------

      kflag=0
      if (Itheory.eq.0)  call SumRules(kflag)
      if (kflag.eq.1) then
         write(6,*) ' --- problem in SumRules, kflag = 1'
         call HF_errlog(12020516,
     +        'F: FCN - problem in SumRules, kflag = 1')
      endif

*     -----------------------------------------------------
*      set alphas
*     -----------------------------------------------------

      if(itheory.eq.0) then 
         call setalf(dble(alphas),Mz*Mz)
         alphaSzero= hf_get_alphas(1.0D0)
         call RT_SetAlphaS(alphaSzero)
         if(IPDFSET.eq.5) then
            call PDFINP(LHAPDFsubr, IPDFSET, dble(0.001), epsi, nwds)
         endif
      endif 

      
      if (iflag.eq.1) then
         open(87,file='output/pulls.first.txt')
      endif

      if ((iflag.eq.3).or.(iflag.eq.4)) then
         open(88,file='output/pulls.last.txt')
         do i=1,nset
            npts(i) = 0
         enddo
      endif


C
C     Call a subrotine which vanishes nonvalence DGLAP contribution
C      for dipole model fits.

      if (DipoleModel.gt.2) then
         call LeaveOnlyValenceQuarks
      endif


*     ---------------------------------------------------------  	 
*     Call evolution
*     ---------------------------------------------------------  	 
      if (Debug) then
         print*,'before evolution'
      endif

      if (Itheory.eq.0) then         
         call Evolution
      else
         print*,'ITHEORY ne 0'
      endif

c fill PDFs for ABKM scheme
      if (mod(HFSCHEME,10).eq.4) then 
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


      if (Debug) then
         print*,'after GetTheoryIteration'
      endif
		 
*     ---------------------------------------------------------  	       
*     Calculate theory for datasets:
*     ---------------------------------------------------------  	 
      do idataset=1,NDATASETS
         call GetTheoryForDataset(idataset)
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
*     Toy MC samples based on predictions:
* -----------------------------------------------------------

      if (IFlag.eq.1 .and. lrand .and. .not. LRandData) then
         call MC_Method()
      endif 

*     ---------------------------------------------------------
*     calculate chisquare
*     ---------------------------------------------------------
      call GetChisquare(iflag,n0,fchi2,rsys,ersys,pchi2,fcorchi2)

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

      chi2out = fchi2 + DeltaLength
      

      chi2out = fchi2+
     $     shift_polRHp**2+shift_polRHm**2+
     $     shift_polLHp**2+shift_polLHm**2+
     $     shift_polL**2+shift_polT**2

      


      if (lprint) then
         call cpu_time(time3)
         print '(''cpu_time'',3F10.2)', time1, time3, time3-time1 
         write(6,'(A20,i6,F12.2,i4,F12.2)') '
     $        FitPDF chi2out,ndf,chi2out/ndf ',ifcncount, chi2out, 
     $        ndf, chi2out/ndf



      endif ! end  lprint


      if (iflag.eq.1) then
         write(85,*) 'First iteration ',chi2out,ndf,chi2out/ndf
      endif

      if (iflag.eq.3) then
         write(85,'(''After minimisation '',F10.2,I6,F10.3)'),chi2out,ndf,chi2out/ndf
         write(85,*)
         write(6,*)
         write(6,'(''After minimisation '',F10.2,I6,F10.3)'),chi2out,ndf,chi2out/ndf
         write(6,*)
      endif



      if (iflag.eq.3) then
         if(itheory.eq.0) then      
         endif
         write(85,*) ' Partial chi2s '
         do h1iset=1,nset
            if (npts(h1iset).gt.0) then
               write(6,'(''Dataset '',i4,F10.2,i6)')
     $              ,h1iset,pchi2(h1iset),npts(h1iset)
               write(85,'(''Dataset '',i4,F10.2,i6)')
     $              ,h1iset,pchi2(h1iset),npts(h1iset)
            endif
         enddo
         write(85,*)

         write(6,*) 'Correlated Chi2 ', fcorchi2
         write(85,*) 'Correlated Chi2 ', fcorchi2

       base_pdfname = 'output/pdfs_q2val_'

       if (ITheory.ne.2) then
          call Evolution
C LHAPDF output:
          open (76,file='output/lhapdf.block.txt',status='unknown')
          call store_pdfs(base_pdfname)
       endif

       write(85,*) 'Systematic shifts '
c       open(unit=77,file='output/systematics_polar.txt')
       do jsys=1,nsys
c          write(77,*)jsys,' ', SYSTEM(jsys),rsys(jsys),' +/- ',ersys(jsys)
          write(85,'(I5,'' '',A20,'' '',F9.4,''   +/-'',F9.4)')
     $         jsys,SYSTEM(jsys),rsys(jsys),ersys(jsys)
       enddo
c       close(77)

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

