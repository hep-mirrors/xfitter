C-----------------------------------------------
C> \brief Get DIS NC cross section
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C-----------------------------------------------
      Subroutine GetNCXsection(IDataSet, local_hfscheme)
#include "steering.inc"

      if(itheory.ge.100) then
         call GetNCxskt(IDataSet, 'NCDIS')
      else
         call GetDisXsection(IDataSet, 'NCDIS', local_hfscheme)
      endif

      end
      
C-----------------------------------------------
C> \brief Get DIS CC cross section
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C-----------------------------------------------
      Subroutine GetCCXsection(IDataSet, local_hfscheme)
      call GetDisXsection(IDataSet, 'CCDIS', local_hfscheme)
      end

C----------------------------------------------------
C> \brief Get DIS NC charm production cross section
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C----------------------------------------------------
      Subroutine GetNCCharmXsection(IDataSet, local_hfscheme)
#include "steering.inc"
      if(itheory.ge.100) then
         call GetNCxskt(IDataSet, 'CHARMDIS')
      else
         call GetDisXsection(IDataSet, 'CHARMDIS', local_hfscheme)
      endif
      end
C----------------------------------------------------
C> \brief Get DIS NC beauty production cross section
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C----------------------------------------------------
      Subroutine GetNCBeautyXsection(IDataSet, local_hfscheme)
#include "steering.inc"
      if(itheory.ge.100) then
         call GetNCxskt(IDataSet, 'BEAUTYDIS')
      else
         call GetDisXsection(IDataSet, 'BEAUTYDIS', local_hfscheme)
      endif
      end
C----------------------------------------------------
C> \brief Get DIS NC FL
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C----------------------------------------------------
      Subroutine GetNCFL(IDataSet, local_hfscheme)
#include "steering.inc"
      if(itheory.ge.100) then
         call GetNCxskt(IDataSet, 'FL')
      else
         call GetDisXsection(IDataSet, 'FL', local_hfscheme)
      endif
      end

C----------------------------------------------------
C> \brief Get DIS NC F2
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C----------------------------------------------------
      Subroutine GetNCF2(IDataSet, local_hfscheme)
#include "steering.inc"
      if(itheory.ge.100) then
         call GetNCxskt(IDataSet, 'F2')
      else
         call GetDisXsection(IDataSet, 'F2', local_hfscheme)
      endif
      end

C----------------------------------------------------
C> \brief Get DIS NC cross section
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C----------------------------------------------------      
      Subroutine GetIntegratedNCXsection(IDataSet, local_hfscheme)
      call GetIntegratedDisXsection(IDataSet, 'NCDIS', local_hfscheme)
      end

C----------------------------------------------------
C> \brief Get DIS CC cross section
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C----------------------------------------------------       
      Subroutine GetIntegratedCCXsection(IDataSet, local_hfscheme)
C This is not fully tested
      call GetIntegratedDisXsection(IDataSet, 'CCDIS', local_hfscheme)
      end


C----------------------------------------------------------------
C
C> \brief Double differential integrated DIS cross section calculation
C> \details Fills global array THEO.
C> Searches for xmin, xmax, ymin, ymax, q2min, q2max column names 
C> in a data file as well as 'sqrt(S)' and 'lumi(e-)/lumi(tot)' in
C> additional information
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C> \param XSecType type of the DIS process
C
C> \author Krzysztof Nowak
C> \date 20/01/2012
C---------------------------------------------------------------
      Subroutine  GetIntegratedDisXsection(IDataSet, XSecType, local_hfscheme)

      implicit none
#include "ntot.inc"
#include "couplings.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "qcdnumhelper.inc"

      character*(*) XSecType
      integer IDataSet
      integer idxQ2min, idxQ2max, idxYmin, idxYmax, idxXmin, idxXmax
      integer i,  idx, iq2, ix, j, kkk
      
      integer local_hfscheme
      integer nq2split
      parameter(nq2split=25)
      integer nxsplit
      parameter(nxsplit=25)

      double precision X(NPMaxDIS),Y(NPMaxDIS),Q2(NPMaxDIS),XSecT(NPMaxDIS)
      double precision XSecP(NPMaxDIS),XSecM(NPMaxDIS)
      double precision q2_1,q2_2,x_1,x_2, alphaem_run, factor, XSec, YPlus
      double precision dq2(NPMaxDIS), dx(NPMaxDIS)
      double precision Charge, polarity, S
      double precision q2min,q2max,xmin,xmax,ymin,ymax
      double precision EmToEtotRatio, alm_mz
      logical CopyValue
      integer NSubBins


C Functions:
      integer GetBinIndex
      integer GetInfoIndex
      double precision AEMRUN

C---------------------------------------------------------

      if ((nq2split+1)*(nxsplit+1).gt.NPMaxDIS) then
         print *,'ERROR IN GetIntegratedNCXsection'
         print *,'INCREASE NPMax to ',(nq2split+1)*(nxsplit+1)
         call HF_stop
      endif

C
C Get indexes for Q2, x and y bins:
C
      idxQ2min = GetBinIndex(IDataSet,'q2min')
      idxYmin = GetBinIndex(IDataSet,'ymin')
      idxQ2max = GetBinIndex(IDataSet,'q2max')
      idxYmax = GetBinIndex(IDataSet,'ymax')
      S = (DATASETInfo( GetInfoIndex(IDataSet,'sqrt(S)'), IDataSet))**2
      EmToEtotRatio=(DATASETInfo(GetInfoIndex(IDataSet,
     $     'lumi(e-)/lumi(tot)'),IDataSet))

      if (idxQ2min.eq.0 .or. idxQ2max.eq.0) then
         call HF_errlog(12100800,
     + 'F: Q2 bins not well defined in a data file that expecting to calculate DIS cross sections')
      endif
      if (idxYmin.eq.0 .or. idxYmax.eq.0) then
         call HF_errlog(12100801,
     + 'F: y bins need to be defined in a data file to calculate DIS cross sections')
      endif

      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)

         CopyValue = .true.
         if(idx.gt.1) then      ! maybe I can just copy previous entry
            if((idxQ2min.gt.0).and.(idxQ2max.gt.0)) then
               if((AbstractBins(idxQ2min,idx).ne.AbstractBins(idxQ2min,idx-1)).or.
     +              (AbstractBins(idxQ2max,idx).ne.AbstractBins(idxQ2max,idx-1))) then
                  CopyValue = .false.
               endif
            endif
            if((idxYmin.gt.0).and.(idxYmax.gt.0)) then
               if((AbstractBins(idxYmin,idx).ne.AbstractBins(idxYmin,idx-1)).or.
     +              (AbstractBins(idxYmax,idx).ne.AbstractBins(idxYmax,idx-1))) then
                  CopyValue = .false.
               endif
            endif
            
            if(CopyValue) then
               THEO(idx) = THEO(idx-1)
               cycle
            endif
         endif

         q2min = AbstractBins(idxQ2min,idx)
         q2max = AbstractBins(idxQ2max,idx)
         ymin  = AbstractBins(idxYmin,idx)
         ymax  = AbstractBins(idxYmax,idx)
         xmin = q2min / (S * ymax)
         xmax = q2max / (S * ymin)
c add checks on x kinematics
         if(xmin.lt.1D-5) then
            xmin = 1D-5
         endif
         if(xmax.gt.0.98) then
           xmax = 0.98
         endif 


c do integration in log space for Q2 and x, if linear scale wanted, see
c example below for x (commented out)         
         j=0
         do iq2=0,nq2split
            q2_1 = q2_2
            q2_2 = exp( log(q2min) + (log(q2max) - log(q2min)) / nq2split*dble(iq2))            
            if(iq2.gt.0) then
               do ix=0, nxsplit
                  x_1 = x_2
                  x_2 = exp( log(xmin) + (log(xmax) - log(xmin)) / nxsplit*dble(ix))            
                  if(ix.gt.0) then
                     j=j+1
                     dq2(j) = q2_2 - q2_1
                     q2(j) = exp( log(q2_1) + 0.5*(log(q2_2) - log(q2_1)) ) 

                     dx(j) = x_2 - x_1
                     x(j) = exp( log(x_1) + 0.5*(log(x_2) - log(x_1)) ) 
  
                     y(j) = q2(j) / (S * x(j))
c check that calculated y values agree with limits given in data file                  
                     if(y(j).lt.AbstractBins(idxYmin,idx).or.y(j).gt.AbstractBins(idxYmax,idx)) then
                       y(j)  = 0.d0
                     endif 

                  endif
               enddo
            endif
         enddo                  ! loop over q2 subgrid
         NSubBins = j

c -------  same with linear integration in x:         
c         j=0
c         do iq2=0,nq2split
c            q2_1 = q2_2
c            q2_2 = exp( log(q2min) + (log(q2max) - log(q2min)) / nq2split*dble(iq2))            
c            if(iq2.gt.0) then
c               do ix=0, nxsplit-1
c                  j=j+1
c                  dq2(j) = q2_2 - q2_1
c                  q2(j) = exp( log(q2_1) + 0.5*(log(q2_2) - log(q2_1)) ) 
c
c                  
c                  if(LoopOverYBins) then
c                     xmax = q2(j) / (S * AbstractBins(idxYmin,idx))
c                     xmin = q2(j) / (S * AbstractBins(idxYmax,idx))
c                  else
c                     xmax = AbstractBins(idxXmax,idx)
c                     xmin = AbstractBins(idxXmin,idx)
c                  endif
c                  
c                  x(j) = xmin + (xmax-xmin)/dble(nxsplit)*(dble(ix)+0.5)
c                  dx(j) = (xmax-xmin) / dble(nxsplit)
c                  y(j) = q2(j) / (S * x(j))
c               enddo
c            endif
c         enddo                  ! loop over q2 subgrid
c         NSubBins = j


         call ReadPolarityAndCharge(idataset,charge,polarity)
         call CalcReducedXsectionForXYQ2(X,Y,Q2,NSubBins,charge,
     $        polarity,IDataSet,XSecType, local_hfscheme, XSecT)
c also get x-os sections for e+ and e- to calculate ratio later         
         call CalcReducedXsectionForXYQ2(X,Y,Q2,NSubBins, 1.D0,
     $        polarity,IDataSet,XSecType, local_hfscheme, XSecP)
         call CalcReducedXsectionForXYQ2(X,Y,Q2,NSubBins,-1.D0,
     $        polarity,IDataSet,XSecType, local_hfscheme, XSecM)

         XSec = 0.D0
         do j=1, NSubBins
            if(EmToEtotRatio.ne.0.D0) then
               XSecT(j) = EmToEtotRatio*XSecM(j) + (1.D0-EmToEtotRatio)*XSecP(j)
            endif

            Yplus  = 1. + (1.-y(j))**2

c do add cross section if y was outside of cuts
            if(y(j).eq.0.D0) then
               XSecT(j) = 0.D0
            endif

            factor=1.D0
            if (XSecType.eq.'CCDIS') then
               factor=(Mw**4/(Mw**2+q2(j))**2)*Gf**2/(2*pi*x(j))*convfac
            else if (XSecType.eq.'NCDIS') then
c               alm_mz = 1.d0 / 128.9d0
c               alphaem_run = alm_mz/(1. - alm_mz * 2/(3.*pi)*log(q2(j)/mz**2))
               alphaem_run=alphaem ! not a running alpha_em!
               if (DATASETREACTION(IDataSet)
     $              .eq.'FastNLO ep jets normalised') then
                  alphaem_run = aemrun(q2(j))
               endif

               factor=2*pi*alphaem_run**2*YPlus/(x(j)*q2(j)**2)*convfac
            else
               print *, 'GetIntegratedDisXsection, XSecType',XSecType,
     $              'not supported' 
               stop
            endif


            XSecT(j) = XSecT(j) * factor
            XSecT(j) = XSecT(j) * dq2(j)
            XSecT(j) = XSecT(j) * dx(j)

            XSec = XSec+XSecT(j)
         enddo
         
         THEO(idx) =  XSec
c temporary divide over dq2
c         THEO(idx) =  XSec / (q2max - q2min)
c         print *, idx, ':', THEO(idx)
c         call HF_stop


      enddo   ! loop over data points

      end



C----------------------------------------------------------------
C
C> \brief NC and CC double differential reduced cross section calculation 
C>  for dataset IDataSet. 
C> \details Fills global array THEO.
C>  Needs 'Q2', 'x', 'y' columns in a data file and following CInfo fields:
C>                'sqrt(S)','e charge', 'reduced', 'e polarity' 
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C> \param XSecType type of the DIS process
C
C  Created by SG, 25/05/2011
C  Start with zero mass implementation
C                 14/06/2011 : re-introduce RT code
C---------------------------------------------------------------
      Subroutine GetDisXsection(IDataSet, XSecType, local_hfscheme)

      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "fcn.inc"
#include "couplings.inc"
#include "qcdnumhelper.inc"
#include "for_debug.inc"

      character*(*) XSecType
      integer IDataSet, local_hfscheme
      integer idxQ2, idxX, idxY, i,  idx, idxS
      
      double precision X(NPMaxDIS),Y(NPMaxDIS),Q2(NPMaxDIS),XSec(NPMaxDIS)
      double precision Charge, polarity, alphaem_run, factor, S, YPlus
      logical IsReduced

C Functions:
      integer GetBinIndex
      integer GetInfoIndex
      double precision AEMRUN
! [--- KK 2015-08-30
      double precision dsigma_red_t4
! ---]

c H1qcdfunc
      integer ifirst
      data ifirst /1/
C---------------------------------------------------------
      if(debug) then
        print*,'GetDisXsection: XSEC TYPE = ', XSecType
      endif

      if (NDATAPOINTS(IDataSet).gt.NPMaxDIS) then
         print *,'ERROR IN GetDisXsection'
         print *,'INCREASE NPMaxDIS to ',NDATAPOINTS(IDataSet)
         call HF_stop
      endif

C
C Get indexes for Q2, x and y bins:
C
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX  = GetBinIndex(IDataSet,'x')
      idxY = GetBinIndex(IDataSet,'y')
      IsReduced = DATASETInfo( GetInfoIndex(IDataSet,'reduced'), IDataSet).gt.0



      if (idxY.eq.0) then
         idxS =  GetInfoIndex(IDataSet,'sqrt(S)')
         if (idxS .gt. 0) then
            S = (DATASETInfo( GetInfoIndex(IDataSet,'sqrt(S)')
     $           , IDataSet))**2
         else
            print *,
     $ 'ERROR: DIS sample, neigher S nor y are defined !'
            print *,' !!! STOP STOP STOP STOP !!!'
            call HF_stop
         endif
      endif

      if (idxQ2.eq.0 .or. idxX.eq.0) then
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
         Q2(i)  = AbstractBins(idxQ2,idx)
         if (idxY.eq.0) then
            Y(i)   = Q2(i) / ( X(i) * S )
         else
            Y(i)   = AbstractBins(idxY,idx)
         endif

      enddo

      call ReadPolarityAndCharge(idataset,charge,polarity)
      call CalcReducedXsectionForXYQ2(X,Y,Q2,NDATAPOINTS(IDataSet),
     $     charge,polarity,IDataSet,XSecType, local_hfscheme,XSec)


      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)

         Yplus  = 1. + (1.-Y(i))**2

         factor=1.D0
         if(.not. IsReduced) then
            if (XSecType.eq.'CCDIS') then
               factor=(Mw**4/(Mw**2+q2(i))**2)*Gf**2/(2*pi*x(i))*convfac
            else if (XSecType.eq.'NCDIS'.or.XSecType.eq.'CHARMDIS'.or.
     $              XSecType.eq.'BEAUTYDIS') then
!               alphaem_run = aemrun(q2(i))
               alphaem_run = alphaem
               factor=2*pi*alphaem_run**2*Yplus/(x(i)*q2(i)**2)*convfac
            else if (XSecType.eq.'FL') then
               factor=1.D0
            else if (XSecType.eq.'F2') then
               factor=1.D0
            else
               print *, 'GetDisXsection, XSecType',XSecType,
     $              'not supported'
               stop
            endif
         endif

         if (local_hfscheme==27) factor=1d0

         THEO(idx) =  XSec(i)*factor
         
! [--- KK 2015-08-30, WS 2015-10-10
C         print*,'Twist study: Prepare to implemented. doHiTwist = ', doHiTwist
         if(doHiTwist) then
           if (XSecType.eq.'NCDIS') then 
!              print*,'Twist study: Past doHiTwist check. HiTwistType = ',HiTwistType
             if(HiTwistType.eq.'Twist4') then
!               print*,'Twist study: Twist4 from LM implemented'
               THEO(idx) = THEO(idx) + dsigma_red_t4(Y(i),X(i),Q2(i))
             endif
           endif
         endif
! ---]

      enddo

      if ((iflagFCN.eq.3).and.(h1QCDFUNC).and.(XSecType.eq.'NCDIS')) then
         if (ifirst.eq.1) then
            print*,'getting output for the H1QCDFUNC'
        
            call GetH1qcdfuncOutput(charge, 0.36D0 ) ! polarity)
            ifirst=0
          
         endif
      endif
      end



C----------------------------------------------------------------
C
C> \brief Get polarity and charge from the data file
C> \param IDataSet index of data set
C> \param charge of the lepton beam
C> \param polarity of the lepton beam
C
C> \author Krzysztof Nowak
C> \date   18/01/2012
C---------------------------------------------------------------
      Subroutine ReadPolarityAndCharge(IDataSet,charge,polarity)

      implicit none
#include "ntot.inc"
#include "datasets.inc"
#include "polarity.inc"

C Input:
      integer IDataSet

C Output:
      double precision charge, polarity

      double precision err_pol_unc, shift_pol
      double precision err_pol_corL
      double precision err_pol_corT

C Functions:
      integer GetInfoIndex

      polarity=0.d0
      err_pol_unc=0.d0
      err_pol_corL=0.d0
      err_pol_corT=0.d0
      if(GetInfoIndex(IDataSet,'e polarity').ne.0) then
         polarity = DATASETInfo( GetInfoIndex(IDataSet,'e polarity'),
     $     IDataSet)
      endif
      if (polarity.ne.0) then
         err_pol_unc = 
     $        DATASETInfo( GetInfoIndex(IDataSet,'pol err unc'), IDataSet)
         err_pol_corL = 
     $        DATASETInfo( GetInfoIndex(IDataSet,'pol err corLpol'), IDataSet)
         err_pol_corT = 
     $        DATASETInfo( GetInfoIndex(IDataSet,'pol err corTpol'), IDataSet)
      endif
      if(GetInfoIndex(IDataSet,'e charge').ne.0) then
         charge = DATASETInfo( GetInfoIndex(IDataSet,'e charge'), IDataSet)
      endif

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

c
c      if(polarity.ne.0.d0) then
c         print '( ''charge:  '', F8.4, 
c     $        '' pol: '', F16.4, 
c     $        ''shift pol: '', F16.4 , 
c     $        ''shift Lpol: '', F16.4 , 
c     $        ''shift Tpol: '', F16.4 )', 
c     $        charge, polarity,shift_pol,shift_polL,shift_polT
c      endif
c
      end


C----------------------------------------------------------------
C> \brief Double differential reduced cross section calculation 
C>  for a table given by X, Y, Q2. 
C> \details Fills array XSec
C> \param[in] X, Y, Q2 kinematic bins
C> \param[in] npts number of data points
C> \param[in] charge of the lepton beam
C> \param[in] polarity of the lepton beam
C> \param[in] XSecType DIS process type
C> \param[in] IDataSet index of data set
C> \param[in] local_hfscheme heavy flavour scheme
C> \param[out] XSec calculated theoretical cross section
C
C> \author Krzysztof Nowak
C> \date 18/01/2012
C---------------------------------------------------------------
      Subroutine CalcReducedXsectionForXYQ2(X,Y,Q2,npts,charge,polarity,
     $     IDataSet,XSecType,local_hfscheme,XSec)

      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "qcdnumhelper.inc"

C Input:
      integer npts, IDataSet, local_hfscheme
      character*(*) XSecType
      double precision X(NPMaxDIS),Y(NPMaxDIS),Q2(NPMaxDIS)
      double precision Charge, polarity
C Output: 
      double precision XSec(NPMaxDIS),XSecs(NPMaxDIS)
      integer i, idx
      double precision yplus(NPMaxDIS), yminus(NPMaxDIS)
      double precision F2(NPMaxDIS),xF3(NPMaxDIS),FL(NPMaxDIS)
      double precision F2gamma(NPMaxDIS),FLgamma(NPMaxDIS)

      double precision F2in(NPMaxDIS),xF3in(NPMaxDIS),FLin(NPMaxDIS)
      double precision F2gammain(NPMaxDIS),FLgammain(NPMaxDIS)

      double precision F2c(NPMaxDIS),FLc(NPMaxDIS),F2b(NPMaxDIS),FLb(NPMaxDIS)


C---------------------------------------------------------


C Protect against overflow of internal arrays:
C
      if (npts.gt.NPMaxDIS) then
         print *,'ERROR IN CalculateReducedXsection'
         print *,'INCREASE NPMaxDIS (qcdnumhelper.inc) to ',npts
         call HF_stop
      endif

      if (itheory.lt.50) then

         call UseZmvnsScheme(F2, FL, xF3, F2gamma, FLgamma,
     $        q2, x, npts, polarity, charge, XSecType, local_hfscheme)

         F2in=F2
         FLin=FL
         xF3in=xF3
         F2gammain=F2gamma
         FLgammain=FLgamma

         if     (mod(local_hfscheme,10).eq.1) then
            
            call UseAcotScheme(F2, FL, XF3, F2c, FLc, F2b, FLb, 
     $           x, q2, npts, polarity, XSecType, F2in, FLin, XF3in, 
     $           charge, local_hfscheme, IDataSet)
            
         elseif (mod(local_hfscheme,10).eq.2) then
            
            call UseRtScheme(F2, FL, XF3, F2c, FLc, F2b, FLb, 
     $           x, q2, npts, XSecType, F2gammain, FLgammain,
     $        local_hfscheme, IDataSet)
            
         elseif (mod(local_hfscheme,10).eq.3) then 
            
            call UseHqstfScheme(F2, FL, XF3, F2c, FLc, F2b, FLb, 
     $           x, q2, npts, XSecType)
            
         elseif (mod(local_hfscheme,10).eq.4) then 
            
            call UseABKMFFScheme(F2, FL, XF3, F2c, FLc, F2b, FLb, 
     $           x, q2, npts, XSecType, charge, polarity, F2gamma,
     $           FLgamma, IDataSet)

         elseif (mod(local_hfscheme,10).eq.5) then

            call UseFONLLScheme(F2,FL,xF3,F2c,FLc,F2b,FLb, 
     1                          x,q2,npts,polarity,XSecType,
     2                          charge,local_hfscheme,IDataSet)

         elseif (mod(local_hfscheme,10).eq.6) then

            call UseMELAZmvnsScheme(F2,FL,xF3,F2c,FLc,F2b,FLb, 
     1                              x,q2,npts,XSecType,
     2                              charge,IDataSet)

         elseif (mod(local_hfscheme,10).eq.7) then

            call UseSAcotchiScheme(F2,FL,XF3,F2c,FLc,F2b,FLb,xsecs,
     $           x, q2, y, npts, polarity, XSecType, F2in, FLin, XF3in,
     $           charge, local_hfscheme, IDataSet)

         endif

      elseif (itheory.eq.50) then
         call UseFractalFit(F2, FL, XF3, x,q2,npts,XSecType, IDataSet)
         
      endif
      

C all the transformations below are array operations!
      yplus = 1+(1-y)**2
      yminus = 1-(1-y)**2

      if(XSecType.eq.'CCDIS') then
         if (charge.gt.0) then
            XSec = 0.5*(yplus*F2 - yminus*xF3 - y*y*FL)
            Xsec = Xsec*(1+polarity)
         else
            XSec = 0.5*(yplus*F2 + yminus*xF3 - y*y*FL)
            Xsec = Xsec*(1-polarity)
         endif
      else if(XSecType.eq.'NCDIS') then
         XSec = F2 + yminus/yplus*xF3 - y*y/yplus*FL


c         XSec = FL !hp
      else if(XSecType.eq.'CHARMDIS') then
         XSec = F2c - y*y/yplus*FLc
      else if(XSecType.eq.'BEAUTYDIS') then
         XSec = F2b - y*y/yplus*FLb
      else if(XSecType.eq.'FL') then
         XSec = FL
      else if(XSecType.eq.'F2') then
         XSec = F2
      else
         print *, 'CalcReducedXsectionForXYQ2, XSecType',
     $        XSecType,'not supported'
         stop
      endif
      if (local_hfscheme==27) XSec = XSecs
      
      end


C----------------------------------------------------------------
C> \brief Calculates F2, FL, xF3, F2gamma and FLgamma using ZMVFNS from QCDNUM
C> \param[out] f2, fl, xf3, f2gamma, flgamma structure functions
C> \param[in] q2, x kinematic bin
C> \param[in] npts total number of points
C> \param[in] polarity of the lepton beam
C> \param[in] charge of the lepton beam
C> \param[in] XSecType DIS process type
C
C  Created by Krzysztof Nowak, 31/01/2012
C---------------------------------------------------------------
      subroutine UseZmvnsScheme(f2, fl, xf3, f2gamma, flgamma,
     $     q2, x, npts, polarity, charge, XSecType, local_hfscheme)

      implicit none
#include "steering.inc"
#include "couplings.inc"
#include "qcdnumhelper.inc"


C Input:
      double precision X(NPMaxDIS),Q2(NPMaxDIS)
      double precision charge, polarity
      integer npts, local_hfscheme
      character*(*) XSecType

C Output: 
      double precision F2(NPMaxDIS), FL(NPMaxDIS), xF3(NPMaxDIS)
      double precision F2gamma(NPMaxDIS), FLgamma(NPMaxDIS)
      double precision F2c(NPMaxDIS),FLc(NPMaxDIS)
      double precision F2b(NPMaxDIS),FLb(NPMaxDIS)

C--------------------------------------------------------
C Temporary variables:
      double precision F2m(NPMaxDIS),xF3m(NPMaxDIS),FLm(NPMaxDIS)
      integer i
      double precision ve,ae,au,ad,vu,vd,A_u,A_d,B_u,B_d,pz


C EW param

      double precision sin2th_eff, xkappa, epsilon
      double precision deltar,sweff, sin2thw2
      double precision cau, cad, cvu, cvd




C QCDNUM ZMVFNS, caclulate FL, F2 and xF3 for d- and u- type quarks all bins:

      if(XSecType.eq.'CCDIS') then
         if (charge.gt.0) then
            CALL ZMSTFUN(1,CCEP2F,X,Q2,FL,npts,0)
            CALL ZMSTFUN(2,CCEP2F,X,Q2,F2,npts,0)
            CALL ZMSTFUN(3,CCEP3F,X,Q2,XF3,npts,0)    
         else
            CALL ZMSTFUN(1,CCEM2F,X,Q2,FL,npts,0)
            CALL ZMSTFUN(2,CCEM2F,X,Q2,F2,npts,0)      
            CALL ZMSTFUN(3,CCEM3F,X,Q2,XF3,npts,0) 
         endif
      elseif (XSecType.eq.'NCDIS'.or.XSecType.eq.'CHARMDIS'
     $        .or.XSecType.eq.'F2'
     $        .or.XSecType.eq.'FL'.or.XSecType.eq.'BEAUTYDIS') then
C     u-type ( u+c ) contributions 
         CALL ZMSTFUN(1,CNEP2F,X,Q2,FL,npts,0)
         CALL ZMSTFUN(2,CNEP2F,X,Q2,F2,npts,0)
         CALL ZMSTFUN(3,CNEP3F,X,Q2,XF3,npts,0)    
         
C     d-type (d + s + b) contributions
         CALL ZMSTFUN(1,CNEM2F,X,Q2,FLm,npts,0)
         CALL ZMSTFUN(2,CNEM2F,X,Q2,F2m,npts,0)
         CALL ZMSTFUN(3,CNEM3F,X,Q2,XF3m,npts,0) 
      else
         print *, 'UseZmvnsScheme, XSecType',XSecType,
     $        'not supported'
         stop
      endif

c     for NC needs to combine F2p with F2m etc.        

      if(XSecType.eq.'NCDIS'.or.XSecType.eq.'CHARMDIS'.or.
     $     XSecType.eq.'F2'.or.
     $     XSecType.eq.'FL'.or.XSecType.eq.'BEAUTYDIS') then

         do i=1, npts

            if(EWFIT.eq.0) then
C
C EW couplings of the electron
C
               ve = -0.5d0 + 2.*sin2thw
               ae = -0.5d0         
            
C
C and quarks
C         
         
               au = 0.5d0
               ad = -0.5d0
            
               vu = au - (4.d0/3.d0)*sin2thw
               vd = ad + (2.d0/3.d0)*sin2thw
               
               PZ = 4.d0 * sin2thw * cos2thw * (1.+Mz**2/Q2(i))
            else 

               call wrap_ew(q2(i),sweff,deltar,cau,cad,cvu,cvd,polarity,charge)
               sin2thw2 = 1.d0 - MW**2/MZ**2
               sin2th_eff = sweff
               xkappa = sin2th_eff/sin2thw2
               epsilon = xkappa -1.0
               ve = -0.5d0 + 2.*sin2th_eff
               ae = -0.5d0
               
               vu = cvu - (4.d0/3.d0)*epsilon*sin2thw2
               vd = cvd + (2.d0/3.d0)*epsilon*sin2thw2
               au = cau
               ad = cad
*     
*     Feed the EW parameters to APFEL 
*     
               if (mod(local_hfscheme,10).eq.5) then
                  call SetSin2ThetaW(sin2th_eff)
                  call SetPropagatorCorrection(deltar)
                  call SetEWCouplings(vd,vu,ad,au)
               endif
               
C     Propagator factor PZ
               PZ = 4.d0*sin2thw2*(1.d0 - sin2thw2)*(1.+Mz**2/Q2(i))
               PZ = PZ*(1.d0 - Deltar)
            endif               
            PZ = 1./Pz
C     EW couplings of u-type and d-type quarks at the scale Q2
               
            if (charge.gt.0) then
               A_u = e2u        ! gamma
     $              + (-ve-polarity*ae)*PZ*2.*euq*vu !gamma-Z
     $              + (ve**2 + ae**2+2*polarity*ve*ae)*PZ**2*(vu**2+au**2) !Z
               
               A_d = e2d 
     $              + (-ve-polarity*ae)*PZ*2.*edq*vd 
     $              + (ve**2 + ae**2+2*polarity*ve*ae)*PZ**2*(vd**2+ad**2)
               
               B_u = (ae+polarity*ve)*PZ*2.*euq*au !gamma-Z
     $              + (-2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vu*au !Z
               B_d = (ae+polarity*ve)*PZ*2.*edq*ad 
     $              + (-2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vd*ad
            else
               A_u = e2u        ! gamma
     $              + (-ve+polarity*ae)*PZ*2.*euq*vu !gamma-Z
     $              + (ve**2 + ae**2-2*polarity*ve*ae)*PZ**2*(vu**2+au**2) !Z
               
               A_d = e2d 
     $              + (-ve+polarity*ae)*PZ*2.*edq*vd 
     $              + (ve**2 + ae**2-2*polarity*ve*ae)*PZ**2*(vd**2+ad**2)
               
               B_u = (-ae+polarity*ve)*PZ*2.*euq*au !gamma-Z
     $              + (2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vu*au !Z
               B_d = (-ae+polarity*ve)*PZ*2.*edq*ad 
     $              + (2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vd*ad
               
            endif
            
cv for polarised case should reduce to:
cv         A_u = e2u - ve*PZ*2.*euq*vu +(ve**2 + ae**2)*PZ**2*(vu**2+au**2)
cv         A_d = e2d - ve*PZ*2.*edq*vd +(ve**2 + ae**2)*PZ**2*(vd**2+ad**2)
cv         B_u = -ae*PZ*2.*euq*au + 2.*ve*ae*(PZ**2)*2.*vu*au
cv         B_d = -ae*PZ*2.*edq*ad + 2.*ve*ae*(PZ**2)*2.*vd*ad

            
            F2Gamma(i) = 4.D0/9.D0 * F2(i)  + 1.D0/9.D0 * F2m(i)
            FLGamma(i) = 4.D0/9.D0 * FL(i)  + 1.D0/9.D0 * FLm(i)
            XF3(i)  = B_U*XF3(i)  + B_D*XF3m(i)
            F2(i)   = A_U*F2(i)   + A_D*F2m(i)
            FL(i)   = A_U*FL(i)   + A_D*FLm(i)
            
         enddo
      elseif(XSecType.eq.'CCDIS') then
         do i=1,npts
            F2Gamma(i) = 4.D0/9.D0 * F2(i)  + 1.D0/9.D0 * F2m(i)
            FLGamma(i) = 4.D0/9.D0 * FL(i)  + 1.D0/9.D0 * FLm(i)
         enddo
      else 
         print *, 'UseZmvnsScheme, XSecType',XSecType,
     $        'not supported'
         stop
      endif
      
      end

C----------------------------------------------------------------
C> Calculates F2, FL, XF3, F2c, FLc, F2b, FLb according to ACOT scheme
C> \param[out] F2, FL, xF3, F2c, FLc, F2b, FLb structure functions
C> \param[in] q2, x kinematic bin
C> \param[in] npts total number of points
C> \param[in] polarity of the lepton beam
C> \param[in] charge of the lepton beam
C> \param[in] XSecType DIS process type
C> \param[in] local_hfscheme heavy flavour scheme
C> \param[in] IDataSet data set index
C> \param[in] F2in, FLin, XF3in structure functions calculated by QCDNUM
C
C  Created by Krzysztof Nowak, 23/01/2012
C---------------------------------------------------------------
      subroutine UseAcotScheme(F2, FL, XF3, F2c, FLc, F2b, FLb, 
     $     x, q2, npts, polarity, XSecType,  F2in, FLin, XF3in, 
     $     charge, local_hfscheme, IDataSet)

      implicit none
#include "ntot.inc"
#include "datasets.inc"
#include "steering.inc"
#include "fcn.inc"
#include "qcdnumhelper.inc"

C Input:
      double precision x(NPMaxDIS), q2(NPMaxDIS)
      double precision charge, polarity
      integer npts,IDataSet
      character*(*) XSecType
      double precision F2in(NPMaxDIS), FLin(NPMaxDIS), xF3in(NPMaxDIS)
C Output:
      double precision F2(NPMaxDIS), FL(NPMaxDIS), xF3(NPMaxDIS)
      double precision F2c(NPMaxDIS),FLc(NPMaxDIS)
      double precision F2b(NPMaxDIS),FLb(NPMaxDIS)

C Additional variables:
      integer icharge, i, idx, local_hfscheme
      logical UseKFactors
      double precision f123l(4),f123lc(4),f123lb(4)
      
      if (mod(local_hfscheme,10).eq.1) then
         UseKFactors = .true.   !ACOT Full , ACOT Chi, ACOT ZM
!      else
!         UseKFactors = .false.  !ACOT ZM
      endif

c     icharge_in: 0 NC: photon exchange only
c     icharge_in: 4 NC: e+ gamma+gammaZ+Z 
c     icharge_in: 5 NC: e- gamma+gammaZ+Z 
c     icharge_in:-1 CC e-
c     icharge_in:+1 CC e+
      if(XSecType.eq.'CCDIS') then
         icharge = int(charge)
      else if (XSecType.eq.'NCDIS') then
         if (charge.gt.0) icharge = 4
         if (charge.lt.0) icharge = 5

      else if (XSecType.eq.'CHARMDIS'.or.XSecType.eq.'BEAUTYDIS') then
         icharge = 0
      else
         print *, 'UseAcotScheme, XSecType', XSecType,
     $        'not supported'
         stop
      endif

      do i=1,npts
         idx =  DATASETIDX(IDataSet,i)

C ---------------------------------------
C     F2, FL, XF3 are already computed by QCDNUM:  (FIO 15 Dec 2012)
c     Pass inside sf_acot_wrap for K-Factor method:
C     Important: This relies on the call to UseZmvnsScheme
         f123l(2)= f2in(i)
         f123l(3)=xf3in(i)
         f123l(4)= fLin(i)
         f123l(1)= (f2in(i)-fLin(i))/(2.0d0*x(i))
C ---------------------------------------
         
         call sf_acot_wrap(x(i),q2(i),
     $        f123l,f123lc,f123lb,
     $        local_hfscheme, icharge, 
     $        iFlagFCN, idx,
     $        UseKFactors, polarity)
      
         FL(i)  = F123L(4)
         F2(i)  = F123L(2)
         XF3(i) = x(i) * F123L(3)
c         if ((charge.gt.0).and.(XSecType.eq.'NCDIS')) then
c            XF3(i) = - XF3(i)
c         endif

         FLc(i)  = F123Lc(4)
         F2c(i)  = F123Lc(2)

         FLb(i)  = F123Lb(4)
         F2b(i)  = F123Lb(2)
         
      enddo

      end

C----------------------------------------------------------------      
C> Calculates F2, FL, XF3, F2c, FLc, F2b, FLb according to ACOT scheme
C> \param[out] F2, FL, xF3, F2c, FLc, F2b, FLb structure functions    
C> \param[in] q2, x kinematic bin                                     
C> \param[in] npts total number of points                             
C> \param[in] polarity of the lepton beam                             
C> \param[in] charge of the lepton beam                               
C> \param[in] XSecType DIS process type                              
C> \param[in] local_hfscheme heavy flavour scheme                     
C> \param[in] IDataSet data set index                                
C> \param[in] F2in, FLin, XF3in structure functions calculated by QCDNUM
C                                                                       
C  Created by Krzysztof Nowak, 23/01/2012                               
C---------------------------------------------------------------        
Cmarco2014 IN PROGRESS---------------------------------------------    

      subroutine UseSAcotchiScheme(F2,FL,XF3,F2c,FLc,F2b,FLb,XSec,
     $     x, q2, y, npts, polarity, XSecType,  F2in, FLin, XF3in,
     $     charge, local_hfscheme, IDataSet)

      implicit none
#include "ntot.inc"
#include "datasets.inc"
#include "steering.inc"
#include "fcn.inc"
#include "qcdnumhelper.inc"
C Input: 
      double precision x(NPMaxDIS), q2(NPMaxDIS),y(NPMaxDIS)
      double precision charge, polarity
      integer npts,IDataSet
      character*(*) XSecType
      double precision F2in(NPMaxDIS), FLin(NPMaxDIS), xF3in(NPMaxDIS)
C Output: 
      double precision F2(NPMaxDIS), FL(NPMaxDIS), xF3(NPMaxDIS)
      double precision F2c(NPMaxDIS),FLc(NPMaxDIS)
      double precision F2b(NPMaxDIS),FLb(NPMaxDIS)
      double precision XSec(NPMaxDIS)

C Additional variables:
      integer icharge, i, idx, local_hfscheme,ID_WZ
      logical UseKFactors
      double precision StrFn(0:3)   !f123l(4),f123lc(4),f123lb(4)

c communication with S-ACOT-chi module 
      integer lio               !flag to distinguish among the 4 H1ZEUS data sets 
      integer iLptn1,iLptn2     !labels for init partons
      double precision xsecs
c
      integer idis
c test 
      character t1*23,t2*23

      !call t2f(t1)   

      if(XSecType.eq.'CCDIS') then
        if(charge.eq.-1) then
          lio=21
        elseif(charge.eq.1) then
          lio=-21
        endif
      else if (XSecType.eq.'NCDIS') then
         if (charge.eq.-1) then
           lio=22
         elseif (charge.eq.1) then
           lio=-22
         endif
      else if (XSecType.eq.'CHARMDIS') then
         !if (charge.eq.-1) then
         !  lio=22 
         !elseif (charge.eq.1) then
         !  lio=-22
         !endif
         print*,'Not supported for the moment'
         stop
      else if (XSecType.eq.'FL') then
         print*,'FL is not supported by SACOT-chi for the moment'
         stop
      else if (XSecType.eq.'F2') then
         print*,'F2 is not supported by SACOT-chi for the moment'
         stop
      else
         print *, 'UseSACOTChischeme, XSecType', XSecType,
     $     'not supported'
         stop
      endif

      iLptn1 = lio/10
      iLptn2 = sign(mod(abs(lio),10),lio)

C for the moment
      if (mod(local_hfscheme,10).eq.7) then
        UseKFactors = .true.    !S-ACOT-chi implements several definitions  
                                !for the K-factors 
      else
        print*,'UNDER CONSTRUCTION! ======================='
        STOP
      endif

      ID_WZ=1
      Call SetEwk(ID_WZ)

      idis =mod(local_hfscheme/10,10)
          !=1, using structure as output 
          !=2, using reduce cross section as output    

      do i=1,npts
c        idx =  DATASETIDX(IDataSet,i) 
c        print*,'Idataset, idx = ', IdataSet, idx 

        call sf_sacotchi_wrap(idis,XSecType,x(i),q2(i),y(i)
     &  ,StrFn,iLptn1,iLptn2,xsecs)

        if (idis==1) then !using structure as output 
          F2(i) = StrFn(2)
          FL(i) = StrFn(0)
          xF3(i)= StrFn(3)*x(i)
        else if (idis==2) then !using reduce cross section as output 
          xsec(i)=xsecs
        end if
      enddo
                                                                  
      end
C----------------------------------------------------------------
C> \brief Calculates F2, FL, XF3, F2c, FLc, F2b, FLb 
C>  according to Robert Thorne scheme
C> \param[out] F2, FL, xF3, F2c, FLc, F2b, FLb structure functions
C> \param[in] q2, x kinematic bin
C> \param[in] npts total number of points
C> \param[in] XSecType DIS process type
C> \param[in] local_hfscheme heavy flavour scheme
C> \param[in] IDataSet data set index
C> \param[in] F2gamma, FLgamma structure functions calculated by mstwnc_wrap
C>
C  Created by Krzysztof Nowak, 23/01/2012
C---------------------------------------------------------------
      subroutine UseRtScheme(F2, FL, XF3, F2c, FLc, F2b, FLb, 
     $        x, q2, npts, XSecType, F2gamma, FLgamma, local_hfscheme, IDataSet)

      implicit none
#include "ntot.inc"
#include "datasets.inc"
#include "steering.inc"
#include "fcn.inc"
#include "qcdnumhelper.inc"

C Input:
      double precision x(NPMaxDIS), q2(NPMaxDIS)
      double precision F2gamma(NPMaxDIS), FLgamma(NPMaxDIS)
      integer npts,IDataSet, local_hfscheme
      character*(*) XSecType

C Output:
      double precision F2(NPMaxDIS), FL(NPMaxDIS), xF3(NPMaxDIS)
      double precision F2c(NPMaxDIS),FLc(NPMaxDIS)
      double precision F2b(NPMaxDIS),FLb(NPMaxDIS)

C Additional variables:
      integer i, idx
      logical UseKFactors
      Double precision f2RT,flRT
      Double precision f2cRT,flcRT,f2bRT,flbRT

C RT code good only for NC case      
      if (XSecType.eq.'CCDIS') return
      
      
      if (local_hfscheme.eq.22.or.local_hfscheme.eq.222) then 
         UseKFactors = .true.    ! RT (Standard OR OPT)Fast
      else
         UseKFactors = .false.   ! RT
      endif


C RT does not provide terms beyond gamma exchange. Since they occur at high Q2,
C use QCDNUM to take them into account as a "k"-factor 
C
C  F2_total^{RT} =  F2_{\gamma}^{RT}  *  (  F2_{total}^{QCDNUM}/F2_{\gamma}^{QCDNUM}   
C
!$OMP PARALLEL DO SCHEDULE(GUIDED)
!$OMP& SHARED(DATASETIDX,IDataSet,x,q2,F2Gamma,FLGamma,UseKFactors,F2,FL,F2c,FLc,F2b,FLb,npts) 
!$OMP& DEFAULT(PRIVATE) 
      do i=1,npts
         idx =  DATASETIDX(IDataSet,i)
            call  mstwnc_wrap(
     $        x(i),q2(i),1,
           ! Output:
     $        f2RT,f2cRT,f2bRT,flRT,flcRT,flbRT
           ! Input:
     $          ,iFlagFCN,idx    ! fcn flag, data point index
     $          ,F2Gamma(i),FLGamma(i)
     $          ,UseKFactors
     $          )


      
C     Replace F2,FL from QCDNUM by RT values
C     Keep xF3 from QCDNUM

         F2(i) = F2RT * (F2(i)/F2Gamma(i))
         if (I_Fit_Order.NE.1) then
            FL(i) = FLRT * (FL(i)/FLGamma(i))
         else
            FL(i) = FLRT
         endif

         F2c(i) = f2cRT
         FLc(i) = flcRT
         F2b(i) = f2bRT
         FLb(i) = flbRT

         
      enddo
      end




C----------------------------------------------------------------
C> \brief Calculates F2, FL, XF3, F2c, FLc, F2b and FLb using HQST scheme
C> \param[out] F2, FL, xF3, F2c, FLc, F2b, FLb structure functions
C> \param[in] q2, x kinematic bin
C> \param[in] npts total number of points
C> \param[in] XSecType DIS process type
C
C     Created by Krzysztof Nowak, 31/01/2012
C---------------------------------------------------------------
      subroutine UseHqstfScheme(F2, FL, XF3, F2c, FLc, F2b, FLb, 
     $     x, q2, npts, XSecType)

      implicit none
#include "qcdnumhelper.inc"
      
C Input:
      double precision X(NPMaxDIS),Q2(NPMaxDIS)
      integer npts
      character*(*) XSecType
      
C Output: 
      double precision F2(NPMaxDIS), FL(NPMaxDIS), xF3(NPMaxDIS)
      double precision F2c(NPMaxDIS),FLc(NPMaxDIS),F2b(NPMaxDIS),FLb(NPMaxDIS)

C Additional Variables:
      double precision NC2FHF(-6:6)
      integer i
  
C     HQSTF code good only for NC case      
      if (XSecType.eq.'CCDIS') return

      NC2FHF = 4.D0/9.D0 * CNEP2F  + 1.D0/9.D0 * CNEM2F

c      print*, 'voica is here', CNEP2F, NC2FHF
      CALL HQSTFUN(2,1,NC2FHF,X,Q2,F2c,npts,0)
      CALL HQSTFUN(1,1,NC2FHF,X,Q2,FLc,npts,0)
      CALL HQSTFUN(2,-2,NC2FHF,X,Q2,F2b,npts,0)
      CALL HQSTFUN(1,-2,NC2FHF,X,Q2,FLb,npts,0)
      
      do i=1,npts
         F2(i) = F2(i) + F2c(i) + F2b(i) 
         FL(i) = FL(i) + FLc(i) + FLb(i)
      enddo

      end




C----------------------------------------------------------------
C> \brief  Calculates F2, FL, XF3, F2c, FLc, F2b and FLb using ABKM FF scheme
C> \param[out] F2, FL, xF3, F2c, FLc, F2b, FLb structure functions
C> \param[in] q2, x kinematic bin
C> \param[in] npts total number of points
C> \param[in] XSecType DIS process type
C> \param[in] charge of the lepton beam
C> \param[in] polarity of the lepton beam
C> \param[in] IDataSet data set index
C> \param[in] F2gamma, FLgamma structure functions calculated by QCDNUM
C>
C---------------------------------------------------------------
      subroutine UseABKMFFScheme(F2, FL, XF3, F2c, FLc, F2b, FLb, 
     $     x, q2, npts, XSecType, charge, polarity,  
     $     F2gamma, FLgamma, IDataSet)

      implicit none
#include "ntot.inc"
#include "datasets.inc"
#include "couplings.inc"
#include "steering.inc"
#include "fcn.inc"
#include "qcdnumhelper.inc"
      
C Input:
      double precision X(NPMaxDIS),Q2(NPMaxDIS)
      double precision F2gamma(NPMaxDIS), FLgamma(NPMaxDIS)
      double precision charge, polarity
      integer npts,IDataSet
      character*(*) XSecType
      
C Output: 
      double precision F2(NPMaxDIS), FL(NPMaxDIS), xF3(NPMaxDIS)
      double precision F2c(NPMaxDIS),FLc(NPMaxDIS),F2b(NPMaxDIS),FLb(NPMaxDIS)

C ABKM 
      Double precision f2abkm,flabkm,f3abkm
      Double precision f2cabkm,flcabkm,f3cabkm
      Double precision f2babkm,flbabkm,f3babkm
      integer ncflag

C Additional Variables:
      double precision NC2FHF(-6:6)
      integer i, idx
  

      if (XSecType.eq.'CCDIS') then
        ncflag = 0 
      else
        ncflag = 1
      endif  

      do i=1,npts
        idx =  DATASETIDX(IDataSet,i)

        call sf_abkm_wrap(x(i),q2(i)
     $      ,f2abkm,flabkm,f3abkm,f2cabkm,flcabkm,f3cabkm
     $      ,f2babkm,flbabkm,f3babkm,ncflag,charge
     $      ,polarity,sin2thw,cos2thw,Mz)

        F2(i) = f2abkm + f2cabkm + f2babkm
        FL(i) = flabkm + flcabkm + flbabkm
        XF3(i) = x(i)*(f3abkm+f3cabkm)

        F2c(i) = f2cabkm
        FLc(i) = flcabkm
        F2b(i) = f2babkm
        FLb(i) = flbabkm

c         write(6,*) 'ABKM:x,q2,F2,FL,xF3,f2c,flc,xf3c', 
c     &   i,x(i),q2(i),f2abkm,flabkm,x(i)*f3abkm,f2cabkm,flcabkm,
c     &   x(i)*f3cabkm    


      enddo
      end


C----------------------------------------------------------------
C> Calculates F2, FL, XF3 in parametric form
C> \param[out] F2, FL, xF3  structure functions
C> \param[in] q2, x kinematic bin
C> \param[in] npts total number of points
C> \param[in] IDataSet data set index
C> \param[in] XSecType DIS process type
C
C  Created by Voica Radescu, 24/01/2014
C---------------------------------------------------------------
      subroutine UseFractalFit(F2, FL, XF3, 
     $     x, q2, npts, XSecType, IDataSet)
      

      implicit none
#include "ntot.inc"
#include "datasets.inc"
#include "steering.inc"
#include "fcn.inc"
#include "fractal.inc"
#include "qcdnumhelper.inc"

C     Input:
      double precision x(NPMaxDIS), q2(NPMaxDIS)
      integer npts,IDataSet
      character*(*) XSecType
C     Output:
      double precision F2(NPMaxDIS), FL(NPMaxDIS), xF3(NPMaxDIS)
C     Additional variables:
      integer i, idx
      
      if (XSecType.eq.'NCDIS') then
         do i=1,npts
            idx =  DATASETIDX(IDataSet,i)


!            print*,'frac:', f_D0, f_Q02, f_D3, f_D1, f_D2, f_R

            F2(i)=f_D0*f_Q02*((Q2(i)/(Q2(i)+f_Q02))**(f_D2-1))
     $           *(x(i)**(-f_D2+1))/(1+f_D3-f_D1*log(x(i)))
            F2(i)=F2(i)*(x(i)**(-f_D1*log(1+Q2(i)/f_Q02))*
     $           ((1+Q2(i)/f_Q02)**(f_D3+1))-1)

            FL(i)= F2(i)*f_R/(1+f_R)
            
            XF3(i) =0.d0
         enddo
         
      else
         print *, 'UseFractalFit, XSecType', XSecType,
     $        'not supported'
         stop
      endif
      

C     ---------------------------------------
      
      end

C----------------------------------------------------------------
C> Calculates F2, FL, XF3, F2c, FLc, F2b, FLb according to FONLL scheme
C> \param[out] F2, FL, xF3, F2c, FLc, F2b, FLb structure functions
C> \param[in] q2, x kinematic bin
C> \param[in] npts total number of points
C> \param[in] polarity of the lepton beam
C> \param[in] charge of the lepton beam
C> \param[in] XSecType DIS process type
C> \param[in] local_hfscheme heavy flavour scheme
C> \param[in] IDataSet data set index
C
C  Created by Valerio Bertone, 25/03/2015
C---------------------------------------------------------------
      subroutine UseFONLLScheme(F2,FL,xF3,F2c,FLc,F2b,FLb, 
     1                          x,q2,npts,polarity,XSecType,
     2                          charge,local_hfscheme,IDataSet)
*
      implicit none
*
#include "ntot.inc"
#include "datasets.inc"
#include "steering.inc"
#include "fcn.inc"
#include "qcdnumhelper.inc"
**
*     Input Varibales
*
      integer npts,IDataSet
      integer local_hfscheme
      double precision x(NPMaxDIS),q2(NPMaxDIS)
      double precision charge,polarity
      character*(*) XSecType
**
*     Internal Variables
*
      integer i,nh
      double precision muoQ
      double precision s2,m2,HeavyQuarkMass
      double precision dummy
**
*     Output Variables
*
      double precision F2(NPMaxDIS),FL(NPMaxDIS),xF3(NPMaxDIS)
      double precision F2c(NPMaxDIS),FLc(NPMaxDIS)
      double precision F2b(NPMaxDIS),FLb(NPMaxDIS)
*
      call SetPolarizationDIS(polarity)
*
      if(XSecType.eq.'CCDIS')then
         call SetProcessDIS("CC")
      elseif(XSecType.eq.'NCDIS')then
         call SetProcessDIS("NC")
      elseif(XSecType.eq.'CHARMDIS'.or.XSecType.eq.'BEAUTYDIS') then
         call SetProcessDIS("NC")
      else
         write(6,*) 'UseFONLLScheme, XSecType ',XSecType,
     1              ' not supported'
         call HF_stop
      endif
*
      if(charge.lt.0d0)then
         call SetProjectileDIS("electron")
      else
         call SetProjectileDIS("positron")
      endif
*
      call SetPDFSet("external1")
*
*     First compute structure functions with m2 / Q2 = scalea1
*
      muoQ = dsqrt(scalea1)
c      muoQ = 1d0
      call SetRenQRatio(muoQ)
      call SetFacQRatio(muoQ)
      do i=1,npts
         call sf_fonll_wrap(x(i),q2(i),muoQ,
     1                      F2(i),FL(i),xF3(i),
     2                      F2c(i),FLc(i),F2b(i),FLb(i))
      enddo
*
*     If scaleb1 is different from zero ...
*
      if(scaleb1.ne.0d0)then
c      if(scalea1.ne.1d0.or.scaleb1.ne.0d0)then
*     Subtract heavy quark contributions
         do i=1,npts
            F2(i) = F2(i) - F2c(i) - F2b(i)
            FL(i) = FL(i) - FLc(i) - FLb(i)
         enddo
*     number of heavy quarks in the final state (squared)
         nh = 4                           ! NC
         if(XSecType.eq.'CCDIS') nh = 1   ! CC
*     Recompute the charm contribution at the right scales and add it back
*     to the total structure functions
         do i=1,npts
            s2   = scalea1 * q2(i)
            m2   = scaleb1 * HF_MASS(1)**2 !HeavyQuarkMass(4,dsqrt(q2(i)))**2d0
            muoQ = dsqrt( ( s2 + nh * m2 ) / q2(i) )
            call SetRenQRatio(muoQ)
            call SetFacQRatio(muoQ)
            call sf_fonll_wrap(x(i),q2(i),muoQ,
     1                         dummy,dummy,dummy,
     2                         F2c(i),FLc(i),dummy,dummy)
            F2(i) = F2(i) + F2c(i)
            FL(i) = FL(i) + FLc(i)
         enddo
*     Recompute the bottom contribution at the right scales and add it back
*     to the total structure functions
         do i=1,npts
            s2   = scalea1 * q2(i)
            m2   = scaleb1 * HF_MASS(2)**2 !HeavyQuarkMass(5,dsqrt(q2(i)))**2d0
            muoQ = dsqrt( ( s2 + nh * m2 ) / q2(i) )
            call SetRenQRatio(muoQ)
            call SetFacQRatio(muoQ)
            call sf_fonll_wrap(x(i),q2(i),muoQ,
     1                         dummy,dummy,dummy,
     2                         dummy,dummy,F2b(i),FLb(i))
            F2(i) = F2(i) + F2b(i)
            FL(i) = FL(i) + FLb(i)
         enddo
      endif
*     Adjust sign of F3
      xF3 = - charge * xF3
*     Divide structure functions by 2 for the CC process
      if(XSecType.eq.'CCDIS')then
         F2  = 0.5d0 * F2
         FL  = 0.5d0 * FL
         xF3 = - 0.5d0 * charge * xF3
         F2c = 0.5d0 * F2c
         FLc = 0.5d0 * FLc
         F2b = 0.5d0 * F2b
         FLb = 0.5d0 * FLb
      endif
*
      return
      end

C----------------------------------------------------------------
C> Calculates F2, FL, XF3, F2c, FLc, F2b, FLb in the ZM-VFNS scheme using MELA
C> \param[out] F2, FL, xF3, F2c, FLc, F2b, FLb structure functions
C> \param[in] q2, x kinematic bin
C> \param[in] npts total number of points
C> \param[in] polarity of the lepton beam
C> \param[in] charge of the lepton beam
C> \param[in] XSecType DIS process type
C> \param[in] local_hfscheme heavy flavour scheme
C> \param[in] IDataSet data set index
C> \param[in] F2in, FLin, XF3in structure functions calculated by QCDNUM
C
C  Created by Valerio Bertone, 08/09/2015
C---------------------------------------------------------------
      subroutine UseMELAZmvnsScheme(F2,FL,xF3,F2c,FLc,F2b,FLb, 
     1                              x,q2,npts,XSecType,
     2                              charge,IDataSet)
*
      implicit none
*
#include "ntot.inc"
#include "datasets.inc"
#include "steering.inc"
#include "fcn.inc"
#include "qcdnumhelper.inc"
**
*     Input Varibales
*
      integer npts,IDataSet
      integer local_hfscheme
      double precision x(NPMaxDIS),q2(NPMaxDIS)
      double precision charge
      character*(*) XSecType
**
*     Internal Variables
*
      integer i
      integer nQ
      double complex Q(2)
      double complex SFx(3,0:6)
**
*     Output Variables
*
      double precision F2(NPMaxDIS),FL(NPMaxDIS),xF3(NPMaxDIS)
      double precision F2c(NPMaxDIS),FLc(NPMaxDIS)
      double precision F2b(NPMaxDIS),FLb(NPMaxDIS)
*
      if(XSecType.ne.'NCDIS'.and.XSecType.ne.'CHARMDIS'.and.
     1   XSecType.ne.'BEAUTYDIS') then
         write(6,*) 'UseMELAZmvnsScheme, XSecType ',XSecType,
     1              ' not supported'
         call HF_stop
      endif
*
      if(charge.gt.0d0)then
         write(6,*) 'UseMELAZmvnsScheme, charge ',charge,
     1              ' not supported'
         call HF_stop
      endif
*
      nQ   = 2
      Q(1) = sqrt(starting_scale)
      do i=1,npts
         Q(2) = dsqrt(q2(i))
         call xStructureFunctions(x(i),nQ,Q,SFx)
*
         F2(i)  = abs(SFx(1,0))
         FL(i)  = abs(SFx(2,0))
         xF3(i) = abs(SFx(3,0))

         F2c(i) = abs(SFx(1,4))
         FLc(i) = abs(SFx(2,4))

         F2b(i) = abs(SFx(1,5))
         FLb(i) = abs(SFx(2,5))
      enddo
*
      return
      end
