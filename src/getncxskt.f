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
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "fcn.inc"
#include "couplings.inc"
#include "qcdnumhelper.inc"

      character*(*) XSecType
      integer IDataSet
      integer idxQ2, idxX, idxY, i,  idx
      
      double precision X(NPMaxDIS),Y(NPMaxDIS),Q2(NPMaxDIS),XSec(NPMaxDIS)
      double precision XSecqpm(NPMaxDIS),Xsecglu(NPMaxDIS)
      double precision Charge, polarity, alphaem_run, factor
      logical IsReduced
      
      double Precision xx,q2x,yx,epsilon,phit,phi,phie,phil,phiqpm,phieqpm
      double precision phitl,phill,phic,philc,phib,philb

C Functions:
      integer GetBinIndex
      integer GetInfoIndex
      double precision AEMRUN
      
c cascade stuff
      Integer IcasHF
      Common/CasHF/IcasHF
      Logical Firstgrid
      Data Firstgrid/.true./
      Logical ex
      Integer Iread,Err,Irr
      Data Iread/0/
      logical firsth
	double precision auh 
      common/f2fit/auh(50),firsth
      Integer Iseed
      Common/random/ISEED
      Integer Itheory_ca
      Common/theory/Itheory_ca
      Integer j
      Double precision q2ini,xtest,q2test,f2qpm
      Integer ngridmax
      Parameter (ngridmax = 1000)
      Double Precision q2grid(ngridmax),xgrid(ngridmax),f2qpmgrid(ngridmax)
      Double precision small 
      Parameter (small=1.d-3)
      Integer jmax,jj,ifound
      
      ISEED = 2313134
cc      ISEED = 1234567

C---------------------------------------------------------
      itheory_ca = itheory
      if (NDATAPOINTS(IDataSet).gt.NPMaxDIS) then
         print *,'ERROR IN GetDisXsection'
         print *,'INCREASE NPMaxDIS to ',NDATAPOINTS(IDataSet)
         call HF_stop
      endif
      If(Firstgrid) then
         inquire(FILE='Cascade/updfGrids/f2qpm-grid.dat',EXIST=ex)
c create grid
         if(ex) then
            open(4,FILE='Cascade/updfGrids/f2qpm-grid.dat', FORM='formatted',
     +      STATUS='OLD',IOSTAT=IRR,ERR=220)
            write(6,*) ' Cascade/updfGrids/f2qpm-grid.dat existing '
         else 
            write(6,*) '  creating f2qpm-grid.dat ' 
            open(4,FILE='Cascade/updfGrids/f2qpm-grid.dat', FORM='formatted',STATUS='NEW',
     +      IOSTAT=IRR,ERR=220)
            iread = 0
         Endif
         Firstgrid = .false.
         if(iread.eq.1) then
            q2ini = 0.
            jmax=0
            J=0
330         read(4,*,END=331,ERR=331) xtest,q2test,f2qpm
            J = J + 1
            if(j.gt.ngridmax) then
              write(6,*) ' getncxskt: j > ngridmax -> STOP ',j,ngridmax
              stop
            endif
            q2grid(j) = q2test
            xgrid(j) = xtest
            f2qpmgrid(j)=f2qpm
c            write(6,*) ' test: ',xgrid(j),q2grid(j),f2qpmgrid(j)
            goto 330
331         continue 
            rewind 4
            jmax = j
         endif
       Endif

C
C Get indexes for Q2, x and y bins:
C
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX  = GetBinIndex(IDataSet,'x')
      idxY = GetBinIndex(IDataSet,'y')
      IsReduced = DATASETInfo( GetInfoIndex(IDataSet,'reduced'), IDataSet).gt.0

      if (idxQ2.eq.0 .or. idxX.eq.0 ) then
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
         if(idxY.ne.0) Y(i)   = AbstractBins(idxY,idx)
         Q2(i)  = AbstractBins(idxQ2,idx)
         xx = x(i)
         q2x = q2(i)
         yx = y(i)
         epsilon = yx**2/2d0/(1d0-yx+yx**2/2d0)
c         epsilon =0
         ICasHF = 0
c something needs to be done for charm
         if(XSecType.eq.'CHARMDIS') then
            IcasHF = 4
         endif
c         write(6,*) ' getncxskt ',XsecType,i,Idataset,idx
c call kt facotrisation sigma_red       
         if(itheory.eq.101) then
c somthing needs to be done here for FL
            call sigcalc(xx,q2x,phi,phie)
            elseif(Itheory.eq.102.or.Itheory.eq.105) then 
            if(Xsecglu(idx).eq.0) then
c something needs to be done here for FL
              IcasHF=0
              call sigcalc(xx,q2x,phil,phie)
              IcasHF=4
              call sigcalc(xx,q2x,phic,phie)
              IcasHF=5
              call sigcalc(xx,q2x,phib,phie)
              
              Xsecglu(idx) = phil+phic+phib
            endif
            
            else
c call kt factorisation sigma_red   
            phitl =0 
            phill =0 
            phic = 0
            philc = 0
            phib = 0
            philb = 0      
            call siggrid(xx,q2x,phitl,phill,phic,philc,phib,philb)
            if(XSecType.eq.'CHARMDIS') then
               phit = phic
               phil = philc
               else
               phit = phitl + phic + phib
               phil = phill + philc + philb
            endif
            if(XSecType.eq.'F2') then
               phil =0
            endif
            phi = phit - epsilon*phil
            if (XSecType.eq.'FL') then
                phi = phil
            endif
            if(phi.ne.phi) then
              write(6,*) ' getncxskt: xx,q2x',xx,q2x,phit,phil,epsilon
            endif
         Endif
         if(Xsecglu(idx).ne.0) then
           Xsec(i) = Xsecglu(idx)
           else
           Xsec(i) = phi
         Endif
         if(Xsecqpm(idx).eq.0) then
            if(IcasHF.lt.4.and.XsecType.ne.'FL' ) then
               If(Itheory.eq.101) then
                  phiqpm = 0
                  if(auh(7).gt.0.001) call sigqpm(xx,q2x,phiqpm,phieqpm)
               else
                  if(iread.eq.0) then
                     phiqpm = 1E-44
                     if(auh(7).gt.0.001) call sigqpm(xx,q2x,phiqpm,phieqpm)
c                     write(6,*) ' qpm :x,q2,qpm,norm: ',xx,q2x,phi,phiqpm,auh(7)
                     write(4,*) xx,q2x,phiqpm  
                  else
                     ifound = 0
                     do jj=1,jmax
                      if(xx.ge.xgrid(jj)-small*xgrid(jj).and.xx.le.xgrid(jj)+small*xgrid(jj)) then 
                         if(q2x.ge.q2grid(jj)-small.and.q2x.le.q2grid(jj)+small) then
                            phiqpm=f2qpmgrid(jj)  
                            ifound = 1 
                            goto 334
                         endif                      
                      endif
                     enddo
                  
  334                if(ifound.eq.1) then
c                     write(6,*) ' x,q2 point found ',xx,xgrid(jj),jj,q2x,q2grid(jj),jj,f2qpmgrid(jj),phiqpm
                     else
                        write(6,*) 'no x point found ',jj
                     endif
                  endif
               endif
               if (XSecType.eq.'FL') then
                   phiqpm = 0
               endif
               Xsecqpm(idx) = phiqpm
            endif
         endif 
      enddo
      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)

         factor=1.D0

         If(Itheory.eq.101) then
            THEO(idx) =  factor*(XSec(i) + Xsecqpm(idx)*max(0.,auh(7)))
	      if(XSecType.eq.'CHARMDIS') then
	       THEO(idx) =  factor*(XSec(i))
	      endif
         Elseif(Itheory.eq.102.or.Itheory.eq.105) then
            THEO(idx) =  factor*(XSec(i)*auh(1) + Xsecqpm(idx)*max(0.,auh(7)))
	      if(XSecType.eq.'CHARMDIS') then
	        THEO(idx) =  factor*(XSec(i)*auh(1))
	      endif
         Elseif(Itheory.eq.103.or.Itheory.eq.104) then
            THEO(idx) =  factor*(XSec(i) + Xsecqpm(idx)*max(0.,auh(7)))
	      if(XSecType.eq.'CHARMDIS') then
	        THEO(idx) =  factor*(XSec(i))
	      endif
         Endif
      enddo
      Return
220   write(6,*) ' getncxskt: f2qpm_grid.dat  on unit 3 cannot be opened '          

      end


