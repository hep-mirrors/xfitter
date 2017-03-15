      Subroutine GetTopFastNLOXsectionNormalised(IDataSet
     >  ,UseNormalisation)
C---------------------------------------------------------------------------
C
C  Created 25/06/2014.  Calculate top-pair cross sections 
C  By M. Guzzi and K. Lipka
C
C---------------------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"

      logical UseNormalisation
      integer IDataSet
      integer GetInfoIndex
      double precision UseZMVFNS
      integer local_hfscheme, idx

      idx = GetInfoIndex(IDataSet,'UseZMVFNS')
      
      if (idx.gt.0) then
         UseZMVFNS=DATASETInfo(idx,IDataSet)
      else
         UseZMVFNS=0
      endif

      if(UseZMVFNS.eq.0.) then
         local_hfscheme = HFSCHEME
      else
         local_hfscheme = 0
      endif
      
      Call GetTopFastNLOXsection(IDataSet, UseNormalisation)
      
      end










      Subroutine GetTopFastNLOXsection(IDataSet, UseNormalisation)

C---------------------------------------------------------------------------
C
C  Created 25/06/2014. Calculate top-pair production cross sections 
C
C---------------------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"

      logical UseNormalisation
      integer IDataSet
      integer NFmax, NdataBins, Npth
      integer i, k, idx, idxVarlow, idxVarhigh
      integer NGauss

C Functions:
      integer GetBinIndex

      parameter(NdataBins=100)
      parameter(NFmax=8000)
      parameter(NGauss=25)

      double precision Varlow(NdataBins), Varhigh(NdataBins)
      double precision xmin,xmax,XsecV,TOTXsec,total
      
      double precision xgauss(NGauss), wgauss(NGauss)

      integer NFmax1
      double precision XSec(NFmax), XSec1(NFmax),ThBin(NFmax)
      common / integral / ThBin, XSec, NFmax1, Npth



      NFmax1=NFmax
C  Get cross section

c      print*,'+++++++++++++++++++++++++++++++++++++++++++++++'
c      print*,'Top pair prod Xsec is being computed by DiffTop'
c      print*,'+++++++++++++++++++++++++++++++++++++++++++++++'

      call fastnlocalctop(IDataSet,XSec,ThBin,TOTXsec,Npth)
      

      idxVarlow = GetBinIndex(IDataSet,'ptlow')
      idxVarhigh = GetBinIndex(IDataSet,'pthigh')

      xmin=ThBin(1)
      xmax=ThBin(Npth)

c      print*,'DiffTop tot Xsec: xmin, xmax, TOTXsec=',xmin,xmax,TOTXsec



        
!MK14
!       Initialization to zero of the diff. Xsec array 
!       before THE bin-by-bin integration 
!       to obtain the normalized Xsec
      do i=1,NDATAPOINTS(IDataSet)
        Xsec1(i)=0.0d0
      enddo
      
      
      if(NDATAPOINTS(IDataSet) .eq. 1) then
        
        idx =  DATASETIDX(IDataSet,NDATAPOINTS(IDataSet))
        
        Varlow(1) = AbstractBins(idxVarlow,idx)
        Varhigh(1) = AbstractBins(idxVarhigh,idx)
        
!MK14 hard wired: theory predictions are unsafe for PT<1 GeV           
        if (Varlow(1).eq.0.) then 
          Varlow(1)=1.0d0
        endif
        
        THEO(idx) = TOTXsec
     
      elseif(NDATAPOINTS(IDataSet) .gt. 1) then
        
         do i=1,NDATAPOINTS(IDataSet)
          
          idx =  DATASETIDX(IDataSet,i)
          
          Varlow(i) = AbstractBins(idxVarlow,idx)
          Varhigh(i) = AbstractBins(idxVarhigh,idx)
          
!MK14 hard wired: theory predictions are unsafe for PT<1 GeV           
          if (Varlow(1).eq.0.) then 
            Varlow(1)=1.0d0
          endif
          
         call gauleg(Varlow(i),Varhigh(i),xgauss,wgauss,NGauss)
          
          
          do 11 k=1,NGauss
            Xsec1(i) = Xsec1(i) + XsecV(xgauss(k))*wgauss(k)
            Xsec(i) = Xsec1(i)
 11       continue
          
          
          if (UseNormalisation) then 
            THEO(idx) = Xsec(i)/TOTXsec/(Varhigh(i)-Varlow(i))
C            print*,Varlow(i),Varhigh(i),DATEN(idx),THEO(idx)
          else
            THEO(idx) = Xsec(i)/(Varhigh(i)-Varlow(i))
C            print*,Varlow(i),Varhigh(i),DATEN(idx),THEO(idx)
          endif
          
        enddo
        
      endif    
      
      end
      








      double precision function XsecV(z)
      implicit none

      integer NFmax,i, Npth
      parameter(NFmax=8000)
      double precision z,res


      integer NFmax1
      double precision ThBin(NFmax), XSec(NFmax)
      common / integral / ThBin,XSec, NFmax1, Npth

      call interp(XSec,z,ThBin, Npth, res)

      XsecV=res

      end








      SUBROUTINE gauleg(x1,x2,x,w,n)
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.14159265359d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 



















