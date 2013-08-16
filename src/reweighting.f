c----------------------------------------------------------------------
      subroutine  pdfreweighting
*     ---------------------------------------------------------
*     Main subroutine for PDF reweighting
*     ---------------------------------------------------------
      implicit none
      external fcn

      include 'steering.inc'
      include 'thresholds.inc'
      include 'couplings.inc'
      include 'for_debug.inc'
      include 'fcn.inc'
      include 'reweighting.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'alphas.inc'
      include 'theo.inc'

      character*72 minfile
      integer rwpdfsets
      integer replicas
      integer in
      integer ios
      character*180 rwpdfsteerfile
      character*180 outfilenames
      character*180 pdfreplicagridname

      character*180 fileName
      character*4  intfileName
      character*5  stringnreplicas
      character*1  points

      character delimiter

      logical isMCPDFMethod,isSymmetricErrors
      logical ex
      double precision chi2tot
      integer isSymmErrors
      integer iset
      integer idatapoint
      integer irow
      integer intfile
      
C Function:
      double precision chi2data_theory
      double precision theorypred
      double precision alphasPDF
      double precision corrmatrix(npoints,npoints)

C---------------------------------------------------------------
c     initialize number of PDF replicas
C---------------------------------------------------------------
      call GetPDFUncType(isMCPDFMethod,isSymmetricErrors)
      call numberPDF(rwpdfsets)
      print '(''isSymmetricErrors  '',L1)',
     $     isSymmetricErrors

      if ((isSymmetricErrors)) then
         isSymmErrors=1
      else 
         isSymmErrors=0
      endif

      print '(''isSymmetricErrors  '',I1)',
     $     isSymmErrors

*     ------------------------------------------------
*     initialize RWPDF steering file 
*     ------------------------------------------------ 

      RWPDFSET=RWPDFSET(1:index(RWPDFSET,'.LHgrid')-1)

      if ((RWMETHOD.eq.1)) then
         rwpdfsteerfile='input_steering/'//TRIM(RWPDFSET)
     $        //'_'//TRIM(RWDATA)//'_chi2.in'
      else if ((RWMETHOD.eq.2)) then
         rwpdfsteerfile='input_steering/'//TRIM(RWPDFSET)
     $     //'_'//TRIM(RWDATA)//'_data.in'
      endif

      if ((RWMETHOD.eq.1)) then
         outfilenames=''//TRIM(RWPDFSET)//'_'//TRIM(RWDATA)
     $        //'_'//'chi2'
      else if ((RWMETHOD.eq.2)) then
         outfilenames=''//TRIM(RWPDFSET)//'_'//
     $        TRIM(RWDATA)//'_'//'data'

      endif

      if (isMCPDFMethod) then
         print *,'PDF set is already NNPDF set 
     $     proceed directly to reweighting'
      else 
         replicas=RWREPLICAS
         rwpdfsets=RWREPLICAS
        call create_randompdfreplicas(TRIM(RWPDFSET)//CHAR(0),
     $        TRIM(outfilenames)//CHAR(0),replicas,
     $        isSymmErrors)
        print *,TRIM(RWPDFSET)

        write( stringnreplicas, '(i5)' )  rwpdfsets
        pdfreplicagridname='output/'//TRIM(outfilenames)//'/'//
     $       TRIM(RWPDFSET)//'_'//
     $       TRIM(adjustl(stringnreplicas))//
     $       'InputReplicas.LHgrid'
        print *,TRIM(pdfreplicagridname)

        call InitPDFset(pdfreplicagridname);
      endif

*     ------------------------------------------------
*     check whether file input is specified in hera fitter steering
*     ------------------------------------------------ 

      if ((npoints.eq.0)) then
         Call HF_errlog(12042302,
     $        'W:WARNING: No input files given in herafittersteering.'// 
     $        'Assume they are externally provided in NNPDF format.'// 
     $        'Proceed to reweighting without overwritting')

         if ((RWMETHOD.eq.1)) then
            fileName = 'NNPDF/data/'//
     $           TRIM(RWDATA)//'/'//TRIM(RWPDFSET)//'/chi2/'
     $           //TRIM(RWDATA)//'_'//TRIM(RWPDFSET)//'-chi2.res'

            inquire(file=TRIM(fileName),exist=ex) 
            
            if (ex) then
               Call HF_errlog(12042301,
     $              'W:WARNING: File '//TRIM(fileName)//
     $              ' exists. Proceed to reweighting without'//
     $              ' overwritting')
               goto 36
            endif

         else if ((RWMETHOD.eq.2)) then
            fileName = 'NNPDF/data/'//
     $           TRIM(RWDATA)//'/'
     $           //TRIM(RWDATA)//'_exp.res'
            
            inquire(file=TRIM(fileName),exist=ex) 
                      
            if (ex) then

               npoints = 0
               open(86,file=TRIM(fileName))
               do idatapoint=1,100001
                  read(86,*,IOSTAT=ios) points
                  npoints = npoints + 1
                  if (ios.ne.0) then
                     goto 26
                  endif
                  if (idatapoint == 100000) then
                     print '(''Error: Maximum number of data points'//
     $                    'for NNPDF reweighting exceeded.'//
     $                    'Fix in line 115 of src/nnpdfreweighting.f'')'
                     call HF_stop
                  endif
               enddo     
            endif
         endif
      endif

 26   close(86)

      print *,npoints

      print *,'line 152'

*     ------------------------------------------------
*     write out rweighting input steering file
*     ------------------------------------------------ 

      open(87,file=rwpdfsteerfile,status='unknown')

      write(87,'(A)') '% Prior & Data Settings'
      write(87,'(A,A)') 'NNPDFPATH: NNPDF'
      write(87,'(A,A)') 'PRIOR:  ',TRIM(RWPDFSET)
      write(87,'(A,A)') 'RWDATA: ',TRIM(RWDATA)
      if ((RWMETHOD.eq.1)) then
         write(87,'(A)') 'RWMODE: CHI2' 
      else if ((RWMETHOD.eq.2)) then
         write(87,'(A)') 'RWMODE: DATA' 
      else 
         print '(''Unknown reweighting method for NNPDF requested'')'
         call HF_stop
      endif
      write(87,'(A,I6)') 'NDATA:  ',npoints
      write(87,'(A)') '' 
      write(87,'(A)') '% Output Settings' 
      write(87,'(A,A)') 'OUTDIR: ',TRIM(outfilenames) 
      write(87,'(A)') 'OUTDESC: Bayesian reweighted PDF' 
      write(87,'(A)') 'PLOTFORMAT: ' 
      write(87,'(A)') '' 
      write(87,'(A)') '% LHGrid Settings' 
      write(87,'(A)') 'LHGRID: Y' 

      write( stringnreplicas, '(i5)' ) RWOUTREPLICAS 

      write(87,'(A,A,A,A,A)') 'LHGRIDFILE:  ',TRIM(outfilenames),
c     $     '.LHgrid' 
     $     '_nRep',TRIM(adjustl(stringnreplicas)),'.LHgrid' 
      write(87,'(A,I4)') 'LHGRIDNREP: ',RWOUTREPLICAS 
      write(87,'(A)') 'DESCRIPTION: ' 
      write(87,'(A)') 'ENDDESC' 
      close(87)

*     ------------------------------------------------
*     open output chi or data File - first create target directory             
*     ------------------------------------------------ 

      CALL system('mkdir -p NNPDF/data/'//
     $     TRIM(RWDATA))
      CALL system('mkdir -p NNPDF/data/'//
     $     TRIM(RWDATA)//'/'//TRIM(RWPDFSET))
      if ((RWMETHOD.eq.1)) then
      CALL system('mkdir -p NNPDF/data/'//
     $     TRIM(RWDATA)//'/'//TRIM(RWPDFSET)//'/chi2/')
      else if ((RWMETHOD.eq.2)) then
      CALL system('mkdir -p NNPDF/data/'//
     $     TRIM(RWDATA)//'/'//TRIM(RWPDFSET)//'/obs/')
      endif

*     ------------------------------------------------
*     if chi2 or data files already exist, 
*           skip calling minuit and theory calculations
*     ------------------------------------------------ 

      if ((RWMETHOD.eq.1)) then
         fileName = 'NNPDF/data/'//
     $        TRIM(RWDATA)//'/'//TRIM(RWPDFSET)//'/chi2/'
     $        //TRIM(RWDATA)//'_'//TRIM(RWPDFSET)//'-chi2.res'

      else if ((RWMETHOD.eq.2)) then
         fileName = 'NNPDF/data/'//
     $        TRIM(RWDATA)//'/'
     $        //TRIM(RWDATA)//'_exp.res'
      endif

      inquire(file=TRIM(fileName),exist=ex) 

      if (ex) then
         Call HF_errlog(12042301,
     $        'W:WARNING: File '//TRIM(fileName)//
     $        ' exists. Proceed to reweighting without overwritting')
         goto 36
      endif


      if ((RWMETHOD.eq.2)) then
*     ------------------------------------------------
* Write out data results - Format: bin number | bin border | value
*     ------------------------------------------------
         fileName = 'NNPDF/data/'//
     $        TRIM(RWDATA)//'/'
     $        //TRIM(RWDATA)//'_exp.res'

         open(86,file=TRIM(fileName))
         do idatapoint=1,npoints
            write(86,*) idatapoint, ' ',AbstractBins(1,idatapoint),'  ',
     $           DATEN(idatapoint)
         enddo
         close(86)
         print *,'(ios.ne.0)'
      
*     ------------------------------------------------
* Write out covariance matrix (assume for the moment no correlation)
* for it to have proper format (no newline), need to fill an array first,  
* then do implicit loop
*     ------------------------------------------------
         fileName = 'NNPDF/data/'//
     $        TRIM(RWDATA)//'/'
     $        //TRIM(RWDATA)//'_CovMat.res'

         open(86,file=TRIM(fileName))
 
         do idatapoint=1,npoints
            do irow=1,npoints
               if ((irow.eq.idatapoint)) then
                  corrmatrix(idatapoint,irow)=E_TOT(idatapoint)
               else 
                  corrmatrix(idatapoint,irow)=0
               endif
            enddo
         enddo

         do, idatapoint=1,npoints
            write(86,'(10000g10.2)') ( corrmatrix(idatapoint,irow), 
     $           irow=1,npoints )
         enddo

         close(86)

         fileName = 'NNPDF/data/'//
     $        TRIM(RWDATA)//'/'//TRIM(RWPDFSET)//'/obs/'
     $        //TRIM(RWDATA)//'_1001.res'

      endif

*
* Prepare output file:
*
      open(86,file=TRIM(fileName))
      
*     ------------------------------------------------
*     loop over all Reweighting sets to write out chi2
*     ------------------------------------------------ 
*     RWREPLICAS

      do iset=1, rwpdfsets
         call InitPDF(iset)

         alphas = alphasPDF(Mz)
         chi2tot = chi2data_theory(min(2,iset),RWMETHOD)  

         if ((RWMETHOD.eq.1)) then    
            print '(''Got MC set='',i5,'' chi2='',F10.1,'' ndf='',i5)',
     $           iset,chi2tot,ndfmini
            write(86,*) iset, ' ', chi2tot/ndfmini

         else if ((RWMETHOD.eq.2)) then

            intfile=(1000+iset)

            print *,intfile
            write(intfileName,'(I4)') intfile

            fileName = 'NNPDF/data/'//
     $           TRIM(RWDATA)//'/'//TRIM(RWPDFSET)//'/obs/'
     $           //TRIM(RWDATA)//'_'//TRIM(intfileName)//'.res'

            print *,TRIM(fileName)

            open(86,file=TRIM(fileName))

            print '(''Got MC set='',i5,
     $           ''  theory prediction(0)='',F10.5)',
     $           iset,THEO(1)

            do, idatapoint=1,npoints
               theorypred = THEO(idatapoint)
    
               write(86,*) idatapoint, ' ', 
     $              AbstractBins(1,idatapoint),' ', theorypred
            enddo            
            close(86)
         endif 
      enddo

      close(86)

 36   continue

      call nnnpdf_reweight(rwpdfsets,TRIM(rwpdfsteerfile)//CHAR(0))

      end
