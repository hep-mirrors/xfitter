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

      character*72 minfile
      integer npdfsets
      integer in
      character*163 nnpdfsteerfile
      character*163 outfilenames
      character delimiter
      logical ex
      double precision chi2tot
      
C Function:
      double precision chi2data_theory

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
         print *,'ERROR: reweighting based on data/theory not yet implemented'
         goto 37
      endif 
*     ------------------------------------------------
*     write out NNPDF input steering file
*     ------------------------------------------------ 

      open(87,file=nnpdfsteerfile)

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
         inquire(file='NNPDF/data/'//
     $        TRIM(NNPDFRWDATA)//'/'//TRIM(NNPDFSET)//'/chi2/'
     $        //TRIM(NNPDFRWDATA)//'_'//TRIM(NNPDFSET)//'-chi2.res', 
     $        exist=ex) 
      else if ((NNPDFREWEIGHTMETHOD.eq.2)) then
         inquire(file='NNPDF/data/'//
     $        TRIM(NNPDFRWDATA)//'/'//TRIM(NNPDFSET)//'/obs/'
     $        //TRIM(NNPDFRWDATA)//'_'//TRIM(NNPDFSET)//'-data.res',
     $        exist=ex) 
      endif

      if (ex) goto 36

      if ((NNPDFREWEIGHTMETHOD.eq.1)) then
         open(86,file='NNPDF/data/'//
     $        TRIM(NNPDFRWDATA)//'/'//TRIM(NNPDFSET)//'/chi2/'
     $        //TRIM(NNPDFRWDATA)//'_'//TRIM(NNPDFSET)//'-chi2.res')
      else if ((NNPDFREWEIGHTMETHOD.eq.2)) then
         open(86,file='NNPDF/data/'//
     $        TRIM(NNPDFRWDATA)//'/'//TRIM(NNPDFSET)//'/obs/'
     $        //TRIM(NNPDFRWDATA)//'_'//TRIM(NNPDFSET)//'-data.res')
      endif

      
*     ------------------------------------------------
*     loop over all NNPDF sets to write out chi2
*     ------------------------------------------------ 

      do FLAGNNPDF=1,npdfsets

c         open ( 125, file='input_steering/minuit.out.nnpdf.txt' )
c         minfile='input_steering/minuit.in.nnpdf.txt' 
c         write(6,*) ' read minuit input params from file ',minfile
c         call HF_errlog(12020504,
c     +     'I: read minuit input params from file '//minfile) 
c         open ( 124, file=minfile )
c         open (  17, file='input_steering/minuit.save.nnpdf.txt' )

c         call mintio(124,125,17)

*     ------------------------------------------------
*     initialize MINUIT
*     ------------------------------------------------ 

         call InitPDF(FLAGNNPDF)
C         call minuit(fcn,0)

         chi2tot = chi2data_theory(1)
         
c         close( 124) 
c         close( 125) 
c         close( 17) 

      enddo

      close(86)

36    call nnnpdf_reweight(npdfsets,nnpdfsteerfile)

      close(87)

 37   end
