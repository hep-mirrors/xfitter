C-----------------------------------------------
C> \brief Get DIS NC cross section
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C-----------------------------------------------
      Subroutine GetNCXsection(IDataSet, local_hfscheme)
      include 'steering.inc'

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
      include 'steering.inc'
      if(itheory.ge.100) then
         call GetNCxskt(IDataSet, 'CHARMDIS')
      else
         call GetDisXsection(IDataSet, 'CHARMDIS', local_hfscheme)
      endif
      end
C----------------------------------------------------
C> \brief Get DIS NC FL
C> \param IDataSet index of data set
C> \param local_hfscheme heavy flavour scheme
C----------------------------------------------------
      Subroutine GetNCFL(IDataSet, local_hfscheme)
      include 'steering.inc'
      if(itheory.ge.100) then
         call GetNCxskt(IDataSet, 'FL')
      else
         call GetDisXsection(IDataSet, 'FL', local_hfscheme)
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
C> \autor Krzisztof Nowak
C> \date 20/01/2012
C---------------------------------------------------------------
      Subroutine  GetIntegratedDisXsection(IDataSet, XSecType, local_hfscheme)

      implicit none
      include 'ntot.inc'
      include 'couplings.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'qcdnumhelper.inc'

      character*(*) XSecType
      integer IDataSet
      integer idxQ2min, idxQ2max, idxYmin, idxYmax, idxXmin, idxXmax
      integer i,  idx, iq2, ix, j, kkk
      
      integer local_hfscheme
      integer nq2split
      parameter(nq2split=25)
      integer nxsplit
      parameter(nxsplit=25)

      double precision X(NPMaxDIS),Y(NPMaxDIS),Q2(NPMaxDIS),XSecP(NPMaxDIS),XSecN(NPMaxDIS)
      double precision q2_1, q2_2, alphaem_run, factor, XSec, YPlus
      double precision dq2(NPMaxDIS), dx(NPMaxDIS)
      double precision Charge, polarity, S
      double precision q2min,q2max,xmin,xmax
      double precision EmToEtotRatio, alm_mz
      logical LoopOverYBins, CopyValue
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
      idxXmin  = GetBinIndex(IDataSet,'xmin')
      idxYmin = GetBinIndex(IDataSet,'ymin')
      idxQ2max = GetBinIndex(IDataSet,'q2max')
      idxXmax  = GetBinIndex(IDataSet,'xmax')
      idxYmax = GetBinIndex(IDataSet,'ymax')
      LoopOverYBins = .false.
      S = (DATASETInfo( GetInfoIndex(IDataSet,'sqrt(S)'), IDataSet))**2
      EmToEtotRatio=(DATASETInfo(GetInfoIndex(IDataSet,
     $     'lumi(e-)/lumi(tot)'),IDataSet))

      if (idxQ2min.eq.0 .or. idxQ2max.eq.0) then
         call HF_errlog(12100800,
     + 'F: Q2 bins not well defined in a data file that expecting to calculate DIS cross sections')
      endif
      if (idxXmin.eq.0 .or. idxXmax.eq.0) then
         if (idxYmin.eq.0 .or. idxYmax.eq.0) then
         call HF_errlog(12100801,
     + 'F: x or y bins need to be defined in a data file to calculate DIS cross sections')
         endif
         LoopOverYBins = .true.
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
            if((idxXmin.gt.0).and.(idxXmax.gt.0)) then
               if((AbstractBins(idxXmin,idx).ne.AbstractBins(idxXmin,idx-1)).or.
     +              (AbstractBins(idxXmax,idx).ne.AbstractBins(idxXmax,idx-1))) then
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

         j=0
         do iq2=0,nq2split
            q2_1 = q2_2
            q2_2 = exp( log(q2min) + (log(q2max) - log(q2min)) / nq2split*dble(iq2))            
            if(iq2.gt.0) then
               do ix=0, nxsplit-1
                  j=j+1
                  dq2(j) = q2_2 - q2_1
                  q2(j) = exp( log(q2_1) + 0.5*(log(q2_2) - log(q2_1)) ) 

                  
                  if(LoopOverYBins) then
                     xmax = q2(j) / (S * AbstractBins(idxYmin,idx))
                     xmin = q2(j) / (S * AbstractBins(idxYmax,idx))
                  else
                     xmax = AbstractBins(idxXmax,idx)
                     xmin = AbstractBins(idxXmin,idx)
                  endif
                  
                  x(j) = xmin + (xmax-xmin)/dble(nxsplit)*(dble(ix)+0.5)
                  dx(j) = (xmax-xmin) / dble(nxsplit)
                  y(j) = q2(j) / (S * x(j))
               enddo
            endif
         enddo                  ! loop over q2 subgrid
         NSubBins = j

         polarity = 0.D0


         call CalcReducedXsectionForXYQ2(X,Y,Q2,NSubBins, 1.D0,
     $        polarity,IDataSet,XSecType, local_hfscheme, XSecP)
         call CalcReducedXsectionForXYQ2(X,Y,Q2,NSubBins,-1.D0,
     $        polarity,IDataSet,XSecType, local_hfscheme, XSecN)

         XSec = 0.D0
         do j=1, NSubBins
            XSecP(j) = EmToEtotRatio*XSecN(j) + (1.D0-EmToEtotRatio)*XSecP(j)

            Yplus  = 1. + (1.-y(j))**2
            XSecP(j) = XSecP(j) * YPlus

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

               factor=2*pi*alphaem_run**2/(x(j)*q2(j)**2)*convfac
            else
               print *, 'GetIntegratedDisXsection, XSecType',XSecType,
     $              'not supported' 
               stop
            endif


            XSecP(j) = XSecP(j) * factor
            XSecP(j) = XSecP(j) * dq2(j)
            XSecP(j) = XSecP(j) * dx(j)

            XSec = XSec+XSecP(j)
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
      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'fcn.inc'
      include 'couplings.inc'
      include 'qcdnumhelper.inc'
      include 'for_debug.inc'

      character*(*) XSecType
      integer IDataSet, local_hfscheme
      integer idxQ2, idxX, idxY, i,  idx, idxS
      
      double precision X(NPMaxDIS),Y(NPMaxDIS),Q2(NPMaxDIS),XSec(NPMaxDIS)
      double precision Charge, polarity, alphaem_run, factor, S
      logical IsReduced

C Functions:
      integer GetBinIndex
      integer GetInfoIndex
      double precision AEMRUN

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

         factor=1.D0
         if(.not. IsReduced) then
            if (XSecType.eq.'CCDIS') then
               factor=(Mw**4/(Mw**2+q2(i))**2)*Gf**2/(2*pi*x(i))*convfac
            else if (XSecType.eq.'NCDIS'.or.XSecType.eq.'CHARMDIS') then
!               alphaem_run = aemrun(q2(i))
               alphaem_run = alphaem
               factor=2*pi*alphaem_run**2/(x(i)*q2(i)**2)*convfac
            else if (XSecType.eq.'FL') then
               factor=1.D0
            else
               print *, 'GetDisXsection, XSecType',XSecType,
     $              'not supported'
               stop
            endif
         endif

         THEO(idx) =  XSec(i)*factor


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
      include 'ntot.inc'
      include 'datasets.inc'
      include 'polarity.inc'

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
      polarity = DATASETInfo( GetInfoIndex(IDataSet,'e polarity'),
     $     IDataSet)
      if (polarity.ne.0) then
         err_pol_unc = 
     $        DATASETInfo( GetInfoIndex(IDataSet,'pol err unc'), IDataSet)
         err_pol_corL = 
     $        DATASETInfo( GetInfoIndex(IDataSet,'pol err corLpol'), IDataSet)
         err_pol_corT = 
     $        DATASETInfo( GetInfoIndex(IDataSet,'pol err corTpol'), IDataSet)
      endif
      charge = DATASETInfo( GetInfoIndex(IDataSet,'e charge'), IDataSet)

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
      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'qcdnumhelper.inc'

C Input:
      integer npts, IDataSet, local_hfscheme
      character*(*) XSecType
      double precision X(NPMaxDIS),Y(NPMaxDIS),Q2(NPMaxDIS)
      double precision Charge, polarity
C Output: 
      double precision XSec(NPMaxDIS)
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
     $        q2, x, npts, polarity, charge, XSecType)

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
      else if(XSecType.eq.'FL') then
         XSec = FL
      else
         print *, 'CalcReducedXsectionForXYQ2, XSecType',
     $        XSecType,'not supported'
         stop
      endif
      
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
     $     q2, x, npts, polarity, charge, XSecType)

      implicit none
      include 'steering.inc'
      include 'couplings.inc'
      include 'qcdnumhelper.inc'

C Input:
      double precision X(NPMaxDIS),Q2(NPMaxDIS)
      double precision charge, polarity
      integer npts
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
     $        .or.XSecType.eq.'FL') then
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
     $     XSecType.eq.'FL') then
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
         else
            call wrap_ew(q2,sweff,deltar,cau,cad,cvu,cvd,polarity,charge)
            sin2thw2 = 1.d0 - MW**2/MZ**2
            sin2th_eff = 0.23134d0
            xkappa = sin2th_eff/sin2thw
            epsilon = xkappa -1.0
            ve = -0.5d0 + 2.*sin2th_eff
            ae = -0.5d0
               
            vu = cvu - (4.d0/3.d0)*epsilon*sin2thw2
            vd = cvd + (2.d0/3.d0)*epsilon*sin2thw2
            au = cau
            ad = cad
         endif
         
         do i=1,npts
            
C     Propagator factor PZ
            PZ = 4.d0 * sin2thw * cos2thw * (1.+Mz**2/Q2(i))
               
C     modify propagator for EW corrections -- only for EW fit
            if (EWfit.ne.0) PZ = PZ * (1.d0 - Deltar)
               
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
      include 'ntot.inc'
      include 'datasets.inc'
      include 'steering.inc'
      include 'fcn.inc'
      include 'qcdnumhelper.inc'

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

      else if (XSecType.eq.'CHARMDIS') then
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
      include 'ntot.inc'
      include 'datasets.inc'
      include 'steering.inc'
      include 'fcn.inc'
      include 'qcdnumhelper.inc'

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
      include 'qcdnumhelper.inc'
      
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
      include 'ntot.inc'
      include 'datasets.inc'
      include 'couplings.inc'
      include 'steering.inc'
      include 'fcn.inc'
      include 'qcdnumhelper.inc'
      
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
      include 'ntot.inc'
      include 'datasets.inc'
      include 'steering.inc'
      include 'fcn.inc'
      include 'fractal.inc'      
      include 'qcdnumhelper.inc'

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
