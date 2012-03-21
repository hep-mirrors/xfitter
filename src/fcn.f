c----------------------------------------------------------------------
      subroutine  fcn(npar,g,f,parminuit,iflag,futil)
*     ---------------------------------------------------------
*     Main minimization subroutine for MINUIT
*     ---------------------------------------------------------
      implicit none

      include 'steering.inc'
      include 'pdfparam.inc'
      include 'alphas.inc'
      include 'endmini.inc'
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
*     declaration related to minuit
*     ---------------------------------------------------------
      integer npar,iflag
      double precision g(*),parminuit(*),f,futil
      external futil

*     ---------------------------------------------------------
*     declaration related to alphas
*     for RT code, transfer alpha S
*     ---------------------------------------------------------
      double precision alphaszero
      double precision asfunc

*     ---------------------------------------------------------
*     declaration related to chisquare
*     ---------------------------------------------------------
      double precision fchi2, fcorchi2
      double precision DeltaLength
      double precision BSYS(NSYSMax), RSYS(NSYSMax)
      double precision EBSYS(NSYSMax),ERSYS(NSYSMax)
      double precision pchi2(nset)



*     ---------------------------------------------------------
*     declaration related to code flow/debug
*     ---------------------------------------------------------

      double precision time1,time2,time3

      Integer Icount
      data icount/0/

!*** count FCN calls,  Iprint=0 suppress printing
      integer IfcnCount, Iprint  
      data IfcnCount,Iprint  /0,0/ 
      save IfcnCount,iPrint

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
      integer ipoints_jet(NSET)
      logical refresh, refresh_DIS
      integer isys,ipoint,jpoint
      integer nflav , ierr
      integer idataset
C  x-dependent fs:
      double precision fs0,epsi
      double precision fshermes
      external LHAPDFsubr
C-----------------------------------------------------------------
      
C Store FCN flag in a common block:
      IFlagFCN = IFlag


      IfcnCount=IfcnCount+1
      write(6,*) ' ===========  Calls to fcn= IfcnCount ',IfcnCount

      call HF_errlog(12020515,'I: FCN is called')

*     ---------------------------------------------------------
*     initilise variables 
*     ---------------------------------------------------------
      f = 0.d0
      fchi2 = 0.d0
      ndf = -npar
      n0 = 0


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
*     PDF parameterisation at the starting scale
*     ---------------------------------------------------------
      call PDF_Param_Iteration(parminuit,iflag)

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
         alphaSzero=asfunc(1.0D0,nflav,ierr)
         call RT_SetAlphaS(alphaSzero)
         if(IPDFSET.eq.5) then
            call PDFINP(LHAPDFsubr, IPDFSET, dble(0.001), epsi, nwds)
         endif
      endif 

      
*     -----------------------------------------------------
*      prntouts for debug is set in steering.txt
*     -----------------------------------------------------      
      if (debug) then
         if (Itheory.ne.2) then
            if (itheory.eq.1) then
            elseif (iparam.eq.1) then
c               write(6,*) 'couplings ',cvu,cau,cvd,cad
               write(6,*) 'ag bg cg dg ',ag,bg,cg,dg
               write(6,*) 'au ad aubar adbar fu ',
     +              au,ad,aubar,adbar,fu
               write(6,*) 'dd du bu cu cd cub cdb ',
     +              dd,du,bu,cu,cd,Cubar,Cdbar
               write(6,*) 'alphas_s(M_Z) ',alphas
               
            elseif (iparam.eq.2.or.iparam.eq.22.or.iparam.eq.225
     $              .or.iparam.eq.221.or.iparam.eq.222
     $              .or.iparam.eq.229.or.iparam.eq.227) then
               
               write(6,*) 'iparam alphas couplings ',iparam,alphas!,cvu,cau,cvd,cad
               write(6,*) 'ag bg cg dg apg bpg cpg ',ag,bg,cg,dg,apg,bpg,cpg
               write(6,*) 'aUv,bUv,cUv,dUv,eUv,fUv ',
     +              aUv,bUv,cUv,dUv,eUv,fUv
               write(6,*) 'aDv,bDv,cDv,dDv,fDv ',
     +              aDv,bDv,cDv,dDv,fDv
               write(6,*) 'aUb,bUb,cUb,dUb ',
     +              aUbar,bUbar,cUbar,dUbar
               write(6,*) 'aDb,bDb,cDb,dDb ',
     +              aDbar,bDbar,cDbar,dDbar
            elseif (iparam.eq.4) then
               write(6,*) 'iparam alphas couplings ',iparam,alphas!,cvu,cau,cvd,cad
               write(6,*) 'ag bg cg dg apg bpg cpg ',ag,bg,cg,dg,apg,bpg,cpg
               write(6,*) 'auv,buv,cuv,duv,euv,fuv ',
     +              aUv,bUv,cUv,dUv,eUv,fUv
               write(6,*) 'adv,bdv,cdv,ddv,fdv ',
     +              aDv,bDv,cDv,dDv,fDv
               write(6,*) 'asea,bsea,csea',
     +              asea,bsea,csea
               
            elseif (iparam.eq.301) then
               write(6,*) 'iparam alphas couplings ',iparam,alphas
               write(6,*) 'ag,bg,cg, aUv,bUv,cUv ',ag,bg,cg,aUv,bUv,cUv
 
 
            endif               ! end itheory =1 
            
            
         endif                  ! Itheory<>2
         

      endif                     ! (debug)


*     --------------------------------------------------

      if (iflag.eq.1) then
         open(87,file='output/pulls.first.txt')
      endif

      if ((iflag.eq.3).or.(iflag.eq.4)) then
         open(88,file='output/pulls.last.txt')
         do i=1,nset
            npts(i) = 0
         enddo
      endif



      do i=31,42
         ipoints_jet(i) = 0
      enddo
      refresh_DIS = .true.


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


      if (iflag.eq.1) close(87)


C Store params in a common block:
      do i=1,MNE
         parminuitsave(i) = parminuit(i)
      enddo


      if (iflag.eq.3) then
         if (dobands) then
            print *,'SAVE PDF values'
         endif
         do i=1,MNE
            pkeep(i) = parminuit (i)
         enddo

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

      f = fchi2 + DeltaLength
      

      f=fchi2+
     $     shift_polRHp**2+shift_polRHm**2+
     $     shift_polLHp**2+shift_polLHm**2+
     $     shift_polL**2+shift_polT**2

      icount = icount + 1

      if (lprint) then
         call cpu_time(time3)
         print '(''cpu_time'',3F10.2)', time1, time3, time3-time1 
         write(6,'(A20,i6,F12.2,i4,F12.2)') ' FitPDF f,ndf,f/ndf ',icount, f, ndf, f/ndf



      endif ! end  lprint


      if (iflag.eq.1) then
         write(85,*) 'First iteration ',f,ndf,f/ndf
      endif

      if (iflag.eq.3) then
         write(85,'(''After minimisation '',F10.2,I6,F10.3)'),f,ndf,f/ndf
         write(85,*)
         write(6,*)
         write(6,'(''After minimisation '',F10.2,I6,F10.3)'),f,ndf,f/ndf
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

 999  continue
      return
      end
