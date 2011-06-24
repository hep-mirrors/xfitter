c     ---------------------------------------------------------------------
      subroutine  fcn(npar,g,f,p,iflag,futil)
C----------------------------------------------------------------------
C
C Main minimization subroutine for MINUIT
C
c---------------------------------------------------------------------
      implicit none

      include 'steering.inc'
      include 'pdfparam.inc'
      include 'alphas.inc'
      include 'endmini.inc'
      include 'for_debug.inc'
      include 'couplings.inc'
      include 'lprint.inc'
      include 'ntot.inc'
      include 'covar_chi2.inc'
      include 'datasets.inc'
      include 'systematics.inc'
      INCLUDE 'theo.inc'
      INCLUDE 'indata.inc'
      INCLUDE 'thresholds.inc'
      include 'hadcor.inc'
      include 'fcn.inc'

      integer npar,iflag
      double precision g(*),p(*),f,futil

      external futil


      double precision chisq,fchi2

      double precision d,t,error,errorunc
      double precision errorsta, fac, fcorchi2

      double precision d_i, t_i, error_i, d_j, error_j, t_j
      double precision DeltaLength

      double precision time1,time2,time3

      double precision x

      double precision quv,qdv, qus, qds, qst, qch, qbt, qgl
      double precision sub

      double precision BSYS(NSYSMax), RSYS(NSYSMax)
      double precision EBSYS(NSYSMax),ERSYS(NSYSMax)

      integer kpr,cdebug
      double precision xpr(9)
      data xpr/0.00008d0,0.0008d0,0.008d0,0.08d0,
     $     0.13d0,0.18d0,0.25d0,0.4d0,0.65d0/
      data kpr,cdebug/9,1/
      Integer Icount
      data icount/0/
      integer iq, ix, nndi, ndi,ndi2
      character*25 base_pdfname
      integer ir(nsysmax)
      integer npts
      dimension npts(nset)
cv for saving sm predictions
      double precision f2SM,f1SM,flSM
      double precision pchi2(nset)
      integer i,j,kflag,jsys,ndf,n0,h1iset,jflag,k,pr
      integer ipoints_jet(NSET)
      logical refresh, refresh_DIS
      integer isys,ipoint,jpoint,ifail

      integer nflav , ierr

      integer idataset

cv voica
      double precision alphaS0
      COMMON/INPUT/alphaS0


C  x-dependent fs:
      double precision fs0
      double precision fshermes


CV ======================================================
CV ======================================================
CV passing info from ACOT codes to h1fitter
CV
CV 
      integer Isch, Iset, Iflg, Ihad
      Common /Ischeme/ Isch, Iset, Iflg, Ihad  !*** pass info out to Fnc123 and Fcc123
C----------------------------------------------------
c     First time through ACOT, we'll compute and store K-factors 
c     After that, we'll use a massive LO calc
      integer IfirstACOT, Idata  !*** 
c     IfirstACOT = 1  Calculate K-Factor Table (First pass)
c     IfirstACOT = 0  Use K-Factor Table  (2nd Pass) 
c     IfirstACOT =-1  Compute Full NLO each time (Over-ride: VERY SLOW)
      Common /acotBLK/ IfirstACOT, Idata  !*** pass switch & data point # to ACOT 


      data   IfirstACOT /1/  !*** We'll set to zero after first pass
      integer ifirstGRID
      data   IfirstGRID /1/  !*** We'll set to zero after first pass
C----------------------------------------------------
      integer IfcnCount, Iprint  !*** count FCN calls,  Iprint=0 suppress printing
      data IfcnCount,Iprint  /0,0/ 
      save IfcnCount,iPrint

C
C Functions:
C
      double precision gluon, singlet, Uminus,Dminus,Sea
      double precision asfunc

C------------------------------------------------------------------------

C
C Store FCN parameters:
C
      IFlagFCN = IFlag

      Iset =11  !** pdf set-> force by hand to use external pdfs
      Iflg=0    !*** dummy: not yet used 
      Ihad=1    !*** proton
C     Count function calls and print:
      IfcnCount=IfcnCount+1
      write(6,*) ' ===========  Calls to fcn= IfcnCount ',IfcnCount
      open(95,file='output/kfactors.dat', status='unknown')

CV ======================================================
CV===========================================================
      f = 0.d0

C-1- 29/07/2010: add Itheory=3 condition ----------
      if (Itheory.eq.3) then
        Itheory = 0
      endif
C-2- 29/07/2010: end of the addition --------------


      if (iflag.eq.3) then 
         ifirstACOT=3
         ifirstGRID=1
      endif


      if ((iflag.eq.1).or.(iflag.ge.10)) then

         if (iflag.ge.10) ifirstACOT=iflag
         
         mpar0 = npar
         
         if (q0.ge.qc) then
            fcharm=charm_frac
         else
            fcharm = 0.
         endif

*     set fstrange via steering
         fstrange = strange_frac
         
C 31 Oct 2009:  add x-dependent strange 
C
         if (ifsttype.eq.0) then
            fs0 = fstrange 
         else
            fs0 = fshermes(0.D0)  ! strange fraction at x=0
         endif


         if (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
*     set fcharm via steering
            fcharm =charm_frac
         endif
      endif


C
C PDF parameterisation at the starting scale
C

      call PDF_Param_Iteration(p,iflag)


C 
C Extra constraints on input PDF due to momentum and quark counting sum rules:
C
      kflag=0
      if (Itheory.eq.0)  call SumRules(kflag)
      if (kflag.eq.1) then
         write(6,*) ' --- problem in SumRules, kflag = 1'
         stop
      endif

      do jsys=1,nsys
         bsys(jsys) = 0.d0
         rsys(jsys) = 0.d0
         ebsys(jsys) = 0.d0
         ersys(jsys) = 0.d0
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

      do i=1,ntot
         THEO(i) = 0.d0
      enddo


      if(itheory.eq.0) then 
         call setalf(dble(alphas),Mz*Mz)
         alphas0=asfunc(1.0D0,nflav,ierr)
      endif 

*     -- set alphas

      if (lprint.and.ONLINE) then

C-1- 27/07/2010: added condition 'Itheory<>2' -----
        if (Itheory.ne.2) then
C-2- 27/10/2010: end of the addition --------------

         if (itheory.eq.1) then

         elseif (iparam.eq.1) then
	    write(6,*) 'couplings ',cvu,cau,cvd,cad
            write(6,*) 'ag bg cg dg ',ag,bg,cg,dg
            write(6,*) 'au ad aubar adbar fu ',
     +           au,ad,aubar,adbar,fu
            write(6,*) 'dd du bu cu cd cub cdb ',
     +           dd,du,bu,cu,cd,Cubar,Cdbar
            write(6,*) 'alphas_s(M_Z) ',alphas

         elseif (iparam.eq.2.or.iparam.eq.22.or.iparam.eq.225
     $           .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then

c            write(6,*) 'iparam alphas couplings ',iparam,alphas,cvu,cau,cvd,cad
c            write(6,*) 'ag bg cg dg apg bpg cpg ',ag,bg,cg,dg,apg,bpg,cpg
c            write(6,*) 'aUv,bUv,cUv,dUv,eUv,fUv ',
c     +           aUv,bUv,cUv,dUv,eUv,fUv
c            write(6,*) 'aDv,bDv,cDv,dDv,fDv ',
c     +           aDv,bDv,cDv,dDv,fDv
c            write(6,*) 'aUb,bUb,cUb,dUb ',
c     +           aUbar,bUbar,cUbar,dUbar
c            write(6,*) 'aDb,bDb,cDb,dDb ',
c     +           aDbar,bDbar,cDbar,dDbar
c            write(6,*) 'Rfudge, afudge, F2ht1, f2ht2 ',
c     +           Rfudge, afudge, f2ht1, f2ht2
         elseif (iparam.eq.4) then
            write(6,*) 'iparam alphas couplings ',iparam,alphas,cvu,cau,cvd,cad
            write(6,*) 'ag bg cg dg apg bpg cpg ',ag,bg,cg,dg,apg,bpg,cpg
            write(6,*) 'auv,buv,cuv,duv,euv,fuv ',
     +           aUv,bUv,cUv,dUv,eUv,fUv
            write(6,*) 'adv,bdv,cdv,ddv,fdv ',
     +           aDv,bDv,cDv,dDv,fDv
            write(6,*) 'asea,bsea,csea',
     +           asea,bsea,csea

         endif

C-1- 27/07/2010: closing of added condition -------
        endif   ! Itheory<>2
C-2- 27/07/2010: end of the addition --------------

      endif                     ! (lprint.and.ONLINE)


      chisq=0.d0
      fchi2 = 0.d0
      ndf = -npar
      n0 = 0
*     ---------------------------------------------------------
*     -- a first loop to fill in the theoretical x-sections
*     ---------------------------------------------------------

      if (iflag.eq.1) then
         open(87,file='output/pulls.first.txt')
      endif

      if ((iflag.eq.3).or.(iflag.eq.4)) then
         open(88,file='output/pulls.last.txt')
         do i=1,nset
            npts(h1iset) = 0
         enddo
      endif


*     --------------------------------------------------


      do i=31,42
         ipoints_jet(i) = 0
      enddo
      refresh_DIS = .true.


*     ---------------------------------------------------------  	 
*     -- start of loop over data points: 	 
*     ---------------------------------------------------------
      if (Debug) then
         print*,'before evolution'
      endif

      if (Itheory.eq.0) then
         
         call Evolution
      else
         print*,'ITHEORY ne 0'
      endif

      if ((HFSCHEME.eq.4).or.(vfnsINDX.eq.4)) then
c         call omegridini(xbmin,nxpgrid,nsmgrid)
c         call apeqsgrid(q2min,q2max,xbmin)
c         call fillvfngrid
      endif

      if (Debug) then
         print*,'after evolution'
      endif

C Initialise theory calculation per iteration
      call GetTheoryIteration
      
C Calculate theory for datasets:
      do idataset=1,NDATASETS
         call GetTheoryForDataset(idataset)
      enddo


      call cpu_time(time1)
      do 100 i=1,npoints
         idata=i                !*** Pass this to ACOT for K-Factor use
         
         h1iset = JSET(i)
         
         if (iflag.eq.3) npts(h1iset) = npts(h1iset) + 1
         
         n0 = n0 + 1
         ndf = ndf + 1

         
         
 100  continue

      if (Debug) then
         print*,'after data loop'
      endif

*     ---------------------------------------------------------------------------
*     Pascaud-like chi2 plus error scaled modifications
      if (mod(ICHI2, 10).eq.1) then


         if (.not.lTHEO) then

*     ---------------------------------------------------------
*     -- now calculate the bsys(nsys) and the sysa matrix
*     ---------------------------------------------------------

            do isys = 1,nsys

               do 200 ipoint=1,n0

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


                    if ((h1iset.eq.101).or.(h1iset.eq.102)
     $                    .or.(h1iset.eq.103).or.(h1iset.eq.104)) then
                        error = alpha(ipoint)
                     endif



                  else if (ICHI2.eq.21) then
***   linear scaling
                     error = error*(abs(t/d))
                  else if (ICHI2.eq.31) then
***   sqrt scaling
                     error = error*dsqrt(abs(t/d))
                  endif
                  
                  bsys(isys) = bsys(isys) 
     +                 + t*(d-t)*BETA(isys,ipoint)/error**2 

                  ebsys(isys) = ebsys(isys)
     +                 + t * BETA(isys,ipoint)/error

                  do 300 jsys=1,nsys
                     sysa(isys,jsys) = sysa(isys,jsys)
     +               + beta(isys,ipoint)*beta(jsys,ipoint)*t**2/error**2
 300              continue

 200           continue

            enddo



*     ---------------------------------------------------------
*     --  inverse sysa and find the shifts
*     ---------------------------------------------------------


*            do isys=1,nsys
*              write(6,*) 'isys rsys ',isys,rsys(isys)
*            enddo

            CALL DINV  (NSys,sysa,NSYSMAX,IR,IFAIL)
            do isys=1,nsys
               do jsys=1,nsys
                  rsys(isys) = rsys(isys)-sysa(isys,jsys)*bsys(jsys)
               enddo
            enddo

            if (DEBUG.and.iflag.eq.1) then
               write(78,*)
               do isys=1,nsys
                  write(78,*) 'isys rsys ',isys,rsys(isys)
               enddo
               write(78,*)
            endif


         endif


*     ---------------------------------------------------------
*     -- now calculate the chi2
*     ---------------------------------------------------------


         if (Debug) then
            write(6,*) 'iflag ICHI2 lTHEO',iflag, ICHI2, lTHEO
         endif

         do ipoint=1,n0

            h1iset = JSET(ipoint)
            d = daten(ipoint)
            t = theo(ipoint)

            error = alpha(ipoint)

c            if (h1iset.eq.102) print*,'voica atlas', ipoint,d,t,error

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
               fac = fac - rsys(isys)*beta(isys,ipoint)                
            enddo 

            if (DEBUG.and.iflag.eq.1) then
               write(78,*) 'ipoint fac ',ipoint,fac
            endif


            t = t*fac
            chisq = (d-t)**2/error**2
            fchi2 = fchi2 + chisq

c              write(6,*) 'ipoint d t error/d*100', ipoint, d, t, error/d*100 
c              write(6,*) 'ipoint fac chisq fchi2',ipoint,fac, chisq, fchi2

            if (iflag.eq.3) then
               pchi2(h1iset) = pchi2(h1iset)+chisq
            endif

            if (iflag.eq.1) then
               write(87,880) h1iset, 
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
            endif
            if (iflag.eq.3) then
               write(88,880) h1iset,
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
 880           format(1x, i2, 2x, F9.6, 2x, F10.4, 2x, F9.6, 3(2x, F9.4))
            endif

*     -- errors on the shifts (correct ??)
            do isys=1,nsys
               ersys(isys) = ersys(isys) + 2.*theo(ipoint)**2
     +              *beta(isys,ipoint)**2 / error**2
            enddo

cv            print*,'testvvv', ipoint, d,t,t/fac,chisq, fchi2
         enddo

         do isys=1,nsys
            fchi2 = fchi2 + rsys(isys)**2
            if (ersys(isys).gt.0) then
               ersys(isys) = dsqrt(2. / ersys(isys))
            endif
         enddo

         

c.... added by JF 05/03/07
c....print out the correlated chi2
         if (iflag.eq.3) then
            fcorchi2=0.0d0
            do isys=1,nsys
               fcorchi2= fcorchi2 +rsys(isys)**2
            enddo
         endif
c...........................



      elseif (ICHI2.eq.3) then  ! full covariance matrix

         do ipoint=1,n0

            d_i = daten(ipoint)
            t_i = theo(ipoint)
            error_i = alpha(ipoint)


            do jpoint=1,n0
               d_j = daten(jpoint)
               t_j = theo(jpoint)
               error_j = alpha(jpoint)

               chisq = (d_i-t_i)*cov(ipoint,jpoint)*(d_j-t_j)
               fchi2 = fchi2 + chisq

               if (iflag.eq.3) then
                  h1iset = JSET(ipoint)
                  pchi2(h1iset) = pchi2(h1iset)+chisq
               endif

            enddo

         enddo


*     ---------------------------------------------------------

      elseif (ICHI2.eq.2) then  !  CTEQ-like chi2

         do ipoint=1,n0

            d_i = daten(ipoint)
            t_i = theo(ipoint)
            error_i = alpha(ipoint)

            chisq = (d_i-t_i)**2 / error_i**2
            fchi2 = fchi2 + chisq

            do jsys=1,nsys
               bsys(jsys) = bsys(jsys)
     +              + beta(jsys,ipoint)*(d_i-t_i)/error_i**2
            enddo

            if (iflag.eq.3) then
               h1iset = JSET(ipoint)
               pchi2(h1iset) = pchi2(h1iset)+chisq
            endif


         enddo

*     
*     - Now calculate the term to subtract from chi2 (cf CTEQ)

         sub = 0.d0

         do i=1,nsys
            do j=1,nsys
               sub = sub + bsys(i) * sysa(i,j) * bsys(j)
            enddo
         enddo

         fchi2 = fchi2 - sub


*     -- and get the systematic shifts :

         do i=1,nsys
            RSYS(i) = 0.
            do j=1,nsys
               rsys(i) = rsys(i) + sysa(i,j) * bsys(j)
            enddo
         enddo

c.... added by JF 05/03/07
c....print out the correlated chi2
         if (iflag.eq.3) then
            fcorchi2=0.0d0
            do isys=1,nsys
               fcorchi2= fcorchi2 +rsys(isys)**2
            enddo
         endif
c...........................

      endif

*     ---------------------------------------------------------


      if (iflag.eq.1) close(87)

      if (iflag.eq.3) then
         mpar = npar
         do i=1,mpar0
            pkeep(i) = p(i)
         enddo

c write out data points with fitted theory
         call writefittedpoints

      endif

C  June 27 2009, add pdf lenght term:
      if (ILenPdf.gt.0) then
         call  PDFLength(DeltaLength)
      
         print *,'Chi2 from PDF length:',DeltaLength,fchi2
      else
         DeltaLength = 0.
      endif

      f = fchi2 + DeltaLength

      icount = icount + 1
      if (ONLINE.and.lprint) then
         call cpu_time(time3)
         print '(''cpu_time'',3F10.2)', time1, time3, time3-time1 
         write(6,'(A20,i6,F12.2,i4,F12.2)') ' FitPDF f,ndf,f/ndf ',icount, f, ndf, f/ndf
*START DEBUG

C-1- 27/07/2010: added condition 'Itheory<>2' -----
       if (Itheory.ne.2) then
C-2- 27/10/2010: end of the addition --------------

            if (icount.eq.1.and.cdebug.eq.1) then
               do i = 1,kpr
                  x = xpr(i)
                  write(6,*) ' '
                  call Get_Partons(x,q0,quv,qdv,qus,qds,qst,qch,qbt,qgl)
                  write(6,600) x,q0,quv,qdv,qus,qds,qst,qch,qbt,qgl
                  write(6,*) 'gluon,singlet,uval,dval,sea'
     +                 ,gluon(x),singlet(x),Uminus(x),Dminus(x),sea(x)
               enddo
 600     format(1x, 10(F16.8,2x))
            endif

C-1- 27/07/2010: added condition 'Itheory<>2' -----
        endif   ! Itheory.ne.2
C-2- 27/10/2010: end of the addition --------------

      endif

      if (iflag.eq.1) then
         write(85,*) 'First iteration ',f,ndf,f/ndf
      endif
      if (iflag.eq.3) then
         write(85,*) 'After minimisation ',f,ndf,f/ndf
         write(85,*)
         write(6,*)
         write(6,*) 'After minimisation ',f,ndf,f/ndf
         write(6,*)
      endif


c      if (lprint.and.ONLINE) 
c     +     write(6,'(''RSYS'',4g8.2)') rsys(1),rsys(2),rsys(3),rsys(4)
      if (lNORMA.and.ONLINE) then
         if (lprint) 
     +        write(6,*) ' fitted norm ',
     +        p(15),p(16),p(17),p(18),p(19),p(20),
     +        p(21),p(22),p(23)
      endif


      if (iflag.eq.3) then
         if(itheory.eq.0) then      
         endif
         write(85,*) ' Partial chi2s '
         do h1iset=1,nset
            write(6,*) 'Dataset ',h1iset,pchi2(h1iset),npts(h1iset)
            write(85,*) 'Dataset ',h1iset,pchi2(h1iset),npts(h1iset)
         enddo
         write(85,*)

c.... added by JF 05/03/07
         write(6,*) 'Correlated Chi2 ', fcorchi2
         write(85,*) 'Correlated Chi2 ', fcorchi2
c.................................

       base_pdfname = 'output/pdfs_q2val_'

C-1- 22/07/2010: added the condition -------------
       if (ITheory.ne.2) then
          call Evolution
          call store_pdfs(base_pdfname)
       endif
C-2- 22/07/2010: end of the addition -------------

       write(85,*) 'Systematic shifts '
       open(unit=77,file='output/systematics_polar.txt')
       do jsys=1,nsys
          write(77,*)CompressIdx(jsys),' ', SYSTEM(CompressIdx(jsys)),rsys(jsys),' +/- ',ersys(jsys)
          write(85,*)CompressIdx(jsys),' ', SYSTEM(CompressIdx(jsys)),rsys(jsys),' +/- ',ersys(jsys)
       enddo
       close(77)

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



