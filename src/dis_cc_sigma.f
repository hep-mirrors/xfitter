      Subroutine GetCCXsection(IDataSet)
C----------------------------------------------------------------
C
C  Created by SG, 26/05/2011
C
C----------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'couplings.inc'
      include 'qcdnumhelper.inc'
      include 'polarity.inc'


      integer IDataSet

      integer idxQ2, idxX, idxY, i,  idx
      integer NPmax
      parameter(NPmax=1000)

      double precision X(NPmax),Y(NPmax),Q2(NPmax)
      double precision yplus, yminus

      double precision F2(NPmax),xF3(NPmax),FL(NPmax)

      double precision XSec

      double precision Xv,Yv,Q2v

      double precision Charge, S, FactorCC, polarity 
      double precision err_pol_unc, shift_pol
      double precision err_pol_corL
      double precision err_pol_corT

      logical IsReduced


C Functions:
      integer GetBinIndex
      integer GetInfoIndex


C-------------------------------------------------------------

      polarity=0.d0
      err_pol_unc=0.d0
      err_pol_corL=0.d0
      err_pol_corT=0.d0
      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN GetReducedNCXsection'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)
         stop
      endif

C Get indexes for Q2, x and y:
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX  = GetBinIndex(IDataSet,'x')
      idxY = GetBinIndex(IDataSet,'y')

      if (idxQ2.eq.0 .or. idxX.eq.0 .or. idxY.eq.0) then
         Return
      endif
C prepare bins:
      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)

         X(i)   = AbstractBins(idxX,idx)
         Y(i)   = AbstractBins(idxY,idx)
         Q2(i)  = AbstractBins(idxQ2,idx)
      enddo
     
C QCDNUM, caclulate FL, F2 and xF3 for all bins:
      charge = DATASETInfo( GetInfoIndex(IDataSet,'e charge'), IDataSet)
      polarity = DATASETInfo( GetInfoIndex(IDataSet,'e polarity'), IDataSet)

      if (polarity.ne.0) then
         err_pol_unc = 
     $        DATASETInfo( GetInfoIndex(IDataSet,'pol err unc'), IDataSet)
         err_pol_corL = 
     $        DATASETInfo( GetInfoIndex(IDataSet,'pol err corLpol'), IDataSet)
         err_pol_corT = 
     $        DATASETInfo( GetInfoIndex(IDataSet,'pol err corTpol'), IDataSet)
      endif
      
      S = (DATASETInfo( GetInfoIndex(IDataSet,'sqrt(S)'), IDataSet))**2
      IsReduced = DATASETInfo( GetInfoIndex(IDataSet,'reduced'), IDataSet).gt.0

C----------------------------------------------------------

      if (charge.lt.0.) then
         if (polarity.gt.0) then
            shift_pol=shift_polRHm
         else
            shift_pol=shift_polLHm
         endif
      else
         if (polarity.gt.0) then
            shift_pol=shift_polRHp
         else
            shift_pol=shift_polLHp
         endif
      endif
      
      polarity=polarity*(1+err_pol_unc/100*shift_pol+
     $     err_pol_corL/100*shift_polL+
     $     err_pol_corT/100*shift_polT)
      
      
      
      if(polarity.ne.0.d0) then
         print '( ''charge:  '', F8.4, 
     $        '' pol: '', F16.4, 
     $        ''shift pol: '', F16.4 , 
     $        ''shift Lpol: '', F16.4 , 
     $        ''shift Tpol: '', F16.4 )', 
     $        charge, polarity,shift_pol,shift_polL,shift_polT
      endif
      



      if (charge.gt.0) then
         CALL ZMSTFUN(1,CCEP2F,X,Q2,FL,NDATAPOINTS(IDataSet),0)
         CALL ZMSTFUN(2,CCEP2F,X,Q2,F2,NDATAPOINTS(IDataSet),0)
         CALL ZMSTFUN(3,CCEP3F,X,Q2,XF3,NDATAPOINTS(IDataSet),0)    
      else
         CALL ZMSTFUN(1,CCEM2F,X,Q2,FL,NDATAPOINTS(IDataSet),0)
         CALL ZMSTFUN(2,CCEM2F,X,Q2,F2,NDATAPOINTS(IDataSet),0)      
         CALL ZMSTFUN(3,CCEM3F,X,Q2,XF3,NDATAPOINTS(IDataSet),0) 
      endif

      do i=1,NDATAPOINTS(IDataSet)
         yplus  = 1+(1-y(i))**2
         yminus = 1-(1-y(i))**2
         if (charge.gt.0) then
            XSec = 0.5*(yplus*F2(i) - yminus*xF3(i) - y(i)*y(i)*FL(i))
            Xsec = Xsec*(1+polarity)
         else
            XSec = 0.5*(yplus*F2(i) + yminus*xF3(i) - y(i)*y(i)*FL(i))
            Xsec = Xsec*(1-polarity)
         endif

         if (.not. IsReduced) then
            factorCC=(Mw**4/(Mw**2+q2(i))**2)*Gf**2/(2*pi*x(i))*convfac
            XSec = XSec*factorCC
         endif

         idx =  DATASETIDX(IDataSet,i)
         THEO(idx) =  XSec

      enddo

      end
