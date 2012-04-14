c----------------------------------------------------------------------
      subroutine  nnpdfreweighting
*     ---------------------------------------------------------
*     Main subroutine for NNPDF reweighting
*     ---------------------------------------------------------
      implicit none
      external fcn

      include 'steering.inc'
      include 'thresholds.inc'
      include 'couplings.inc'
      include 'for_debug.inc'
      include 'fcn.inc'
      include 'nnpdf.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'alphas.inc'

      character*72 minfile
      integer npdfsets
      integer in
      character*180 nnpdfsteerfile
      character*180 outfilenames

      character*180 fileName

      character delimiter
      logical ex
      double precision chi2tot
      integer iset
      
C Function:
      double precision chi2data_theory
      double precision alphasPDF
C---------------------------------------------------------------
c     initialize number of PDF replicas
      call numberPDF(npdfsets)


*     ------------------------------------------------
*     initialize NNPDF steering file 
*     ------------------------------------------------ 

      NNPDFSET=NNPDFSET(1:index(NNPDFSET,'.LHgrid')-1)

      if ((NNPDFREWEIGHTMETHOD.eq.1)) then
         nnpdfsteerfile='input_steering/'//TRIM(NNPDFSET)
     $        //'_'//TRIM(NNPDFRWDATA)//'_chi2.in'
      else if ((NNPDFREWEIGHTMETHOD.eq.2)) then
         nnpdfsteerfile='input_steering/'//TRIM(NNPDFSET)
     $     //'_'//TRIM(NNPDFRWDATA)//'_data.in'
      endif



      if ((NNPDFREWEIGHTMETHOD.eq.1)) then
         outfilenames=''//TRIM(NNPDFSET)//'_'//TRIM(NNPDFRWDATA)
     $        //'_'//'chi2'
      else if ((NNPDFREWEIGHTMETHOD.eq.2)) then
         outfilenames=''//TRIM(NNPDFSET)//'_'//
     $        TRIM(NNPDFRWDATA)//'_'//'data'
      endif

      if ((NNPDFREWEIGHTMETHOD.eq.2)) then 
         print *,
     $ 'ERROR: reweighting based on data/theory not yet implemented'
         call hf_stop
      endif 
*     ------------------------------------------------
*     write out NNPDF input steering file
*     ------------------------------------------------ 

      open(87,file=nnpdfsteerfile,status='unknown')

      write(87,'(A)') '% Prior & Data Settings'
      write(87,'(A,A)') 'NNPDFPATH: NNPDF'
      write(87,'(A,A)') 'PRIOR:  ',TRIM(NNPDFSET)
      write(87,'(A,A)') 'RWDATA: ',TRIM(NNPDFRWDATA)
      if ((NNPDFREWEIGHTMETHOD.eq.1)) then
         write(87,'(A)') 'RWMODE: CHI2' 
      else if ((NNPDFREWEIGHTMETHOD.eq.2)) then
         write(87,'(A)') 'RWMODE: DATA' 
      else 
         print '(''Unknown reweighting method for NNPDF requested'')'
         call HF_stop
      endif
      write(87,'(A,I6)') 'NDATA:  ',npoints
      write(87,'(A)') '' 
      write(87,'(A)') '% Output Settings' 
      write(87,'(A,A)') 'OUTDIR: ',TRIM(outfilenames) 
      write(87,'(A)') 'OUTDESC: ' 
      write(87,'(A)') 'PLOTFORMAT: ' 
      write(87,'(A)') '' 
      write(87,'(A)') '% LHGrid Settings' 
      write(87,'(A)') 'LHGRID: Y' 
      write(87,'(A,A,A,I2,A)') 'LHGRIDFILE:  ',TRIM(outfilenames),
     $     '_',NNPDFOUTREPLICAS,'.LHgrid' 
      write(87,'(A,I2)') 'LHGRIDNREP: ',NNPDFOUTREPLICAS 
      write(87,'(A)') 'DESCRIPTION: ' 
      write(87,'(A)') 'ENDDESC' 
      close(87)

*     ------------------------------------------------
*     open output chi or data File - first create target directory             
*     ------------------------------------------------ 

      CALL system('mkdir -p NNPDF/data/'//
     $     TRIM(NNPDFRWDATA))
      CALL system('mkdir -p NNPDF/data/'//
     $     TRIM(NNPDFRWDATA)//'/'//TRIM(NNPDFSET))
      if ((NNPDFREWEIGHTMETHOD.eq.1)) then
      CALL system('mkdir -p NNPDF/data/'//
     $     TRIM(NNPDFRWDATA)//'/'//TRIM(NNPDFSET)//'/chi2/')
      else if ((NNPDFREWEIGHTMETHOD.eq.2)) then
      CALL system('mkdir -p NNPDF/data/'//
     $     TRIM(NNPDFRWDATA)//'/'//TRIM(NNPDFSET)//'/obs/')
      endif

*     ------------------------------------------------
*     if chi2 or data files already exist, 
*           skip calling minuit and theory calculations
*     ------------------------------------------------ 

      if ((NNPDFREWEIGHTMETHOD.eq.1)) then
         fileName = 'NNPDF/data/'//
     $        TRIM(NNPDFRWDATA)//'/'//TRIM(NNPDFSET)//'/chi2/'
     $        //TRIM(NNPDFRWDATA)//'_'//TRIM(NNPDFSET)//'-chi2.res'

      else if ((NNPDFREWEIGHTMETHOD.eq.2)) then
         fileName = 'NNPDF/data/'//
     $        TRIM(NNPDFRWDATA)//'/'//TRIM(NNPDFSET)//'/obs/'
     $        //TRIM(NNPDFRWDATA)//'_'//TRIM(NNPDFSET)//'-data.res'
      endif

      inquire(file=TRIM(fileName),exist=ex) 

      if (ex) then
         Call HF_errlog(12042301,
     $        'W:WARNING: File '//TRIM(fileName)//
     $        ' exists. Proceed to reweighting without overwritting')
         goto 36
      endif

*
* Prepare output file:
*
      open(86,file=TRIM(fileName))


      
*     ------------------------------------------------
*     loop over all NNPDF sets to write out chi2
*     ------------------------------------------------ 

      do iset=1, npdfsets
         call InitPDF(iset)
         alphas = alphasPDF(Mz)
         chi2tot = chi2data_theory(min(2,iset))        
         print '(''Got MC set='',i5,'' chi2='',F10.1,'' ndf='',i5)',
     $        iset,chi2tot,ndfmini
         write(86,*) iset, ' ', chi2tot/ndfmini
      enddo

      close(86)



 36   continue

      call nnnpdf_reweight(npdfsets,TRIM(nnpdfsteerfile)//CHAR(0))


      end
