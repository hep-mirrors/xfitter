c----------------------------------------------------------------------
      subroutine  fcn(npar,g,f,p,iflag,futil)
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


*     ---------------------------------------------------------
*     declaration related to minuit
*     ---------------------------------------------------------
      integer npar,iflag
      double precision g(*),p(*),f,futil
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
      integer i,j,kflag,jsys,ndf,n0,h1iset,jflag,k,pr
      integer ipoints_jet(NSET)
      logical refresh, refresh_DIS
      integer isys,ipoint,jpoint
      integer nflav , ierr
      integer idataset
C  x-dependent fs:
      double precision fs0
      double precision fshermes
C-----------------------------------------------------------------
      
C Store FCN flag in a common block:
      IFlagFCN = IFlag


      IfcnCount=IfcnCount+1
      write(6,*) ' ===========  Calls to fcn= IfcnCount ',IfcnCount


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
      call PDF_Param_Iteration(p,iflag)


*     ---------------------------------------------------------
*     Extra constraints on input PDF due to momentum and quark 
*     counting sum rules:
*     ---------------------------------------------------------
      kflag=0
      if (Itheory.eq.0)  call SumRules(kflag)
      if (kflag.eq.1) then
         write(6,*) ' --- problem in SumRules, kflag = 1'
         stop
      endif

*     -----------------------------------------------------
*      set alphas
*     -----------------------------------------------------

      if(itheory.eq.0) then 
         call setalf(dble(alphas),Mz*Mz)
         alphaSzero=asfunc(1.0D0,nflav,ierr)
         call RT_SetAlphaS(alphaSzero)
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
            npts(h1iset) = 0
         enddo
      endif



      do i=31,42
         ipoints_jet(i) = 0
      enddo
      refresh_DIS = .true.


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

*     ---------------------------------------------------------  	       
*     Calculate theory for datasets:
*     ---------------------------------------------------------  	 
      do idataset=1,NDATASETS
         call GetTheoryForDataset(idataset)
      enddo



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


*     ---------------------------------------------------------
*     calculate chisquare
*     ---------------------------------------------------------
      call GetChisquare(iflag,n0,fchi2,rsys,ersys,pchi2,fcorchi2)



      if (iflag.eq.1) close(87)


      if (iflag.eq.3) then
         mpar = npar
         do i=1,mpar0
            pkeep(i) = p(i)
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



*     ---------------------------------------------------------  	 
*     Calculate chisquare:
*      - first get error matrix
*      - invert error matric to get errors
*      - calculate chisquare
*     ---------------------------------------------------------


      subroutine GetChisquare(flag_in,n0_in,fchi2_in,rsys_in,ersys_in,pchi2_in,fcorchi2_in)

      implicit none

      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'for_debug.inc'
      include 'systematics.inc'
      INCLUDE 'theo.inc'


      integer ir(nsysmax),n0_in 
      integer isys,ipoint,jpoint,ifail,flag_in
      double precision chisq,fchi2_in
      double precision d,t,error,errorunc
      double precision errorsta, fac, fcorchi2_in
      double precision d_i, t_i, error_i, d_j, error_j, t_j         
      integer i,j,jsys,h1iset
      double precision pchi2_in(nset)
      double precision EBSYS_in(NSYSMax),ERSYS_in(NSYSMax)
      double precision BSYS_in(NSYSMax),RSYS_in(NSYSMax)
      double precision sub


*     ----------------------------------------------------------
*     Initialise
*     ----------------------------------------------------------
      chisq=0.d0
      fchi2_in = 0.d0

      sub = 0.d0

      do jsys=1,nsys
         ir(nsys)=0.d0
         bsys_in(jsys) = 0.d0
         rsys_in(jsys) = 0.d0
         ebsys_in(jsys) = 0.d0
         ersys_in(jsys) = 0.d0
      enddo

      do i=1,nset
         pchi2_in(i)=0.d0
      enddo
*     -- for ICHI2=2, the matrix sysa is filled already in Systematics

      if (ICHI2.ne.2) then
         do isys=1,nsys
            do jsys=1,nsys
               sysa(isys,jsys) = 0.d0
               if (jsys.eq.isys) sysa(isys,jsys) = 1.d0
            enddo
         enddo
      endif


*     ----------------------------------------------------------
*     Pascaud-like chi2 plus error scaled modifications
*     ----------------------------------------------------------
      if (mod(ICHI2, 10).eq.1) then



*     ---------------------------------------------------------
*     now calculate the bsys(nsys) and the sysa matrix
*     ---------------------------------------------------------

         do isys = 1,nsys
            
            do 200 ipoint=1,n0_in
               
               d = DATEN(ipoint)
               t = THEO(ipoint)
               error = ALPHA(ipoint)
                  

***   scale errors for chi2 calculation
***   in principle:  unc*(t/d),  sta*dsqrt(t/d)
               if (ICHI2.eq.11) then
***   mixed scaling - decompose - scale - recombine
                  errorunc = E_UNC(ipoint)*d/100.
                  if (errorunc.gt.error) then
                     errorsta = 0.
                  else
                     errorsta = dsqrt(error**2-errorunc**2)
                  endif
                  if (t.gt.0) then
                     errorsta = errorsta*dsqrt(abs(t/d))
                     errorunc = errorunc*(abs(t/d))
                  endif
                  error = dsqrt(errorsta**2+errorunc**2)
                  
                  
c                    if ((h1iset.eq.101).or.(h1iset.eq.102)
c     $                    .or.(h1iset.eq.103).or.(h1iset.eq.104)) then
c                        error = alpha(ipoint)
c                     endif



               else if (ICHI2.eq.21) then
***   linear scaling
                  error = error*(abs(t/d))
               else if (ICHI2.eq.31) then
***   sqrt scaling
                  error = error*dsqrt(abs(t/d))
               endif
                  
               bsys_in(isys) = bsys_in(isys) 
     +              + t*(d-t)*BETA(isys,ipoint)/error**2 
               
               ebsys_in(isys) = ebsys_in(isys)
     +              + t * BETA(isys,ipoint)/error
               
               do 300 jsys=1,nsys
                  sysa(isys,jsys) = sysa(isys,jsys)
     +                 + beta(isys,ipoint)*beta(jsys,ipoint)*t**2/error**2
 300           continue
               
 200        continue

         enddo



*     ---------------------------------------------------------
*     inverse sysa and find the shifts
*     ---------------------------------------------------------

         
         if (nsys.gt.0) then
            CALL DINV (NSys,sysa,NSYSMAX,IR,IFAIL)
         endif

         do isys=1,nsys
            do jsys=1,nsys
               rsys_in(isys) = rsys_in(isys)-sysa(isys,jsys)*bsys_in(jsys)
            enddo
         enddo

         if (DEBUG.and.flag_in.eq.1) then
            write(78,*)
            do isys=1,nsys
               write(78,*) 'isys rsys ',isys,rsys_in(isys)
            enddo
            write(78,*)
         endif




*     ---------------------------------------------------------
*     now calculate the chi2
*     ---------------------------------------------------------

         do ipoint=1,n0_in

            h1iset = JSET(ipoint)
            d = daten(ipoint)
            t = theo(ipoint)

            error = alpha(ipoint)

***   scale errors for chi2 calculation - as above!
            if (ICHI2.eq.11) then
***   mixed scaling - decompose - scale - recombine
               errorunc = E_UNC(ipoint)*d/100.
               if (errorunc.gt.error) then
                  errorsta = 0.
               else
                  errorsta = dsqrt(error**2-errorunc**2)
               endif
               if (t.gt.0) then
                  errorsta = errorsta*dsqrt(abs(t/d))
                  errorunc = errorunc*(abs(t/d))
               endif
               error = dsqrt(errorsta**2+errorunc**2)

               if ((h1iset.eq.101).or.(h1iset.eq.102)
     $              .or.(h1iset.eq.103).or.(h1iset.eq.104)) then
                  error = alpha(ipoint)
               endif



            else if (ICHI2.eq.21) then
***   linear scaling
               error = error*(abs(t/d))
            else if (ICHI2.eq.31) then
***   sqrt scaling
               error = error*dsqrt(abs(t/d))
            endif

            fac = 1.d0
            do isys=1,nsys
               fac = fac - rsys_in(isys)*beta(isys,ipoint)                
            enddo 

            if (DEBUG.and.flag_in.eq.1) then
               write(78,*) 'ipoint fac ',ipoint,fac
            endif


            t = t*fac
            THEO_MOD(ipoint)=t
            ALPHA_MOD(ipoint)=error
            chisq = (d-t)**2/error**2
            fchi2_in = fchi2_in + chisq

            if (flag_in.eq.3) then
               pchi2_in(h1iset) = pchi2_in(h1iset)+chisq
            endif

            if (flag_in.eq.1) then
               write(87,880) h1iset, 
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
            endif
            if (flag_in.eq.3) then
               write(88,880) h1iset,
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
 880           format(1x, i2, 2x, G12.6, 2x, G12.4, 2x, G12.6, 3(2x, G12.4))
            endif

*     -- errors on the shifts 
            do isys=1,nsys
               ersys_in(isys) = sqrt(sysa(isys,isys))
            enddo

         enddo

         do isys=1,nsys
            fchi2_in = fchi2_in + rsys_in(isys)**2
         enddo

         

c....print out the correlated chi2
         if (flag_in.eq.3) then
            fcorchi2_in=0.0d0
            do isys=1,nsys
               fcorchi2_in= fcorchi2_in +rsys_in(isys)**2
            enddo
         endif


c...........................


      elseif (ICHI2.eq.2) then  !  CTEQ-like chi2

         do ipoint=1,n0_in

            d_i = daten(ipoint)
            t_i = theo(ipoint)
            error_i = alpha(ipoint)

            chisq = (d_i-t_i)**2 / error_i**2
            fchi2_in = fchi2_in + chisq

            do jsys=1,nsys
               bsys_in(jsys) = bsys_in(jsys)
     +              + beta(jsys,ipoint)*(d_i-t_i)/error_i**2
            enddo

            if (flag_in.eq.3) then
               h1iset = JSET(ipoint)
               pchi2_in(h1iset) = pchi2_in(h1iset)+chisq
            endif


         enddo

*     
*     - Now calculate the term to subtract from chi2 (cf CTEQ)

         do i=1,nsys
            do j=1,nsys
               sub = sub + bsys_in(i) * sysa(i,j) * bsys_in(j)
            enddo
         enddo


         print *,'haha',fchi2_in,sub
         fchi2_In = fchi2_in - sub


*     -- and get the systematic shifts :

         do i=1,nsys
            RSYS_in(i) = 0.
            do j=1,nsys
               rsys_in(i) = rsys_in(i) + sysa(i,j) * bsys_in(j)
            enddo
         enddo

c....print out the correlated chi2
         if (flag_in.eq.3) then
            fcorchi2_in=0.0d0
            do isys=1,nsys
               fcorchi2_in= fcorchi2_in +rsys_in(isys)**2
            enddo
         endif
c...........................

      endif

      return
      end

