      Subroutine GetNCxskt(IDataSet, XSecType)
C----------------------------------------------------------------
C
C  NC double differential reduced cross section calculation 
C  for dataset IDataSet. Fills global array THEO.
C  Needs 'Q2', 'x', 'y' columns in a data file and following CInfo fields:
C                'sqrt(S)','e charge', 'reduced', 'e polarity' 
C
C---------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'fcn.inc'
      include 'couplings.inc'
      include 'qcdnumhelper.inc'

      character*(*) XSecType
      integer IDataSet
      integer idxQ2, idxX, idxY, i,  idx
      
      double precision X(NPMaxDIS),Y(NPMaxDIS),Q2(NPMaxDIS),XSec(NPMaxDIS)
      double precision XSecqpm(NPMaxDIS),Xsecglu(NPMaxDIS)
      double precision Charge, polarity, alphaem_run, factor
      logical IsReduced
      
      double Precision xx,q2x,phi,phie,phiqpm,phieqpm

C Functions:
      integer GetBinIndex
      integer GetInfoIndex
      double precision AEMRUN
      
c cascade stuff
      Logical Firstgrid
      Data Firstgrid/.true./
      Logical ex
      Integer Iread,Err,Irr
      Data Iread/1/
      logical firsth
	double precision auh 
      common/f2fit/auh(50),firsth
      Integer Iseed
      Common/random/ISEED
      
      ISEED = 2313134

C---------------------------------------------------------

      if (NDATAPOINTS(IDataSet).gt.NPMaxDIS) then
         print *,'ERROR IN GetDisXsection'
         print *,'INCREASE NPMaxDIS to ',NDATAPOINTS(IDataSet)
         call HF_stop
      endif
      If(Firstgrid) then
         inquire(FILE='theoryfiles/updf/f2qpm-grid.dat',EXIST=ex)
c create grid
         if(ex) then
            open(4,FILE='theoryfiles/updf/f2qpm-grid.dat', FORM='formatted',
     +      STATUS='OLD',IOSTAT=IRR,ERR=220)
            write(6,*) ' theoryfiles/updf/f2qpm-grid.dat existing '
         else 
            write(6,*) '  creating f2qpm-grid.dat ' 
            open(4,FILE='theoryfiles/updf/f2qpm-grid.dat', FORM='formatted',STATUS='NEW',
     +      IOSTAT=IRR,ERR=220)
         Endif
         Firstgrid = .false.
       Endif

C
C Get indexes for Q2, x and y bins:
C
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX  = GetBinIndex(IDataSet,'x')
      idxY = GetBinIndex(IDataSet,'y')
      IsReduced = DATASETInfo( GetInfoIndex(IDataSet,'reduced'), IDataSet).gt.0


      if (idxQ2.eq.0 .or. idxX.eq.0 .or. idxY.eq.0) then
         Return
      endif
C prepare bins:
      do i=1,NDATAPOINTS(IDataSet)
C
C Reference from the dataset to a global data index:
C
         idx =  DATASETIDX(IDataSet,i)
C
C Local X,Y,Q2 arrays, used for QCDNUM SF caclulations:
C
         X(i)   = AbstractBins(idxX,idx)
         Y(i)   = AbstractBins(idxY,idx)
         Q2(i)  = AbstractBins(idxQ2,idx)
         xx = x(i)
         q2x = q2(i)
c         write(6,*) ' getncxskt ',i,Idataset,idx
c call kt facotrisation sigma_red         
         if(itheory.eq.101) then
            call sigcalc(xx,q2x,phi,phie)
            elseif(Itheory.eq.102) then 
            if(Xsecglu(idx).eq.0) then
              call sigcalc(xx,q2x,phi,phie)
              Xsecglu(idx) = phi
              write(6,*) 'getncxskt ',i,Idataset,xx,q2x,y(i), phi,phie
            endif
            
            else
c call kt facotrisation sigma_red         
            call siggrid(xx,q2x,phi,phie)
         Endif
         if(Xsecglu(idx).ne.0) then
           Xsec(i) = Xsecglu(idx)
           else
           Xsec(i) = phi
c         write(6,*) idx, Xsecqpm(idx)
         Endif

         if(Xsecqpm(idx).eq.0) then
            If(Itheory.ne.103) then
               if(auh(7).gt.0.001) call sigqpm(xx,q2x,phiqpm,phieqpm)
            else
               if(iread.eq.0) then
                  if(auh(7).gt.0.001) call sigqpm(xx,q2x,phiqpm,phieqpm)
c               write(6,*) ' itheory1 qpm :x,q2,data,theo: ',xx,q2x,phi,phiqpm,auh(7)
                  write(4,*) phiqpm  
               else
                  read(4,*) phiqpm
               endif
            endif
            Xsecqpm(idx) = phiqpm
c            write(6,*) ' qpm xsection taken from memory '
         endif 
c         write(6,*) ' getncxskt ',xx,q2x,phi
      enddo
c      write(6,*) ' getncxskt ',itheory,auh(1),auh(7)
      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)

         factor=1.D0
         If(Itheory.eq.101) then
            THEO(idx) =  factor*(XSec(i) + Xsecqpm(idx)*max(0.,auh(7)))
         Elseif(Itheory.eq.102) then
            THEO(idx) =  factor*(XSec(i)*auh(1) + Xsecqpm(idx)*max(0.,auh(7)))
         Elseif(Itheory.eq.103) then
            THEO(idx) =  factor*(XSec(i) + Xsecqpm(idx)*max(0.,auh(7)))
         Endif

      enddo
      Return
220   write(6,*) ' getncxskt: f2qpm_grid.dat  on unit 3 cannot be opened '          

      end


