      Subroutine  GetIntegratedNCXsection(IDataSet)
C-----------------------------------------------------------------
C
C  Created by SG, 24/05/2011
C  Calculate x-integrated  cross section.
C
C-----------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'couplings.inc'
      include 'qcdnumhelper.inc'
      
      integer IDataSet

      integer idxQ2, idxX1, idxX2, i, j, idx

      integer NSplit
      parameter (NSplit=100)

      integer NQ2Max
      parameter (NQ2Max=60)

      integer Npt

      double precision charge,S,y,yplus,yminus,polarity
      double precision xmin,xmax,q2v,x(0:NSplit,NQ2Max)
     $     ,q2(0:NSplit,NQ2Max)
      double precision F2p(0:NSplit,NQ2Max),xF3p(0:NSplit,NQ2Max)
     $     ,FLp(0:NSplit,NQ2Max)
      double precision F2m(0:NSplit,NQ2Max),xF3m(0:NSplit,NQ2Max)
     $     ,FLm(0:NSplit,NQ2Max)
      double precision F2,xF3,FL
      double precision XSec(0:NSplit),Xsint
      

      double precision ve,ae,au,ad,vu,vd,A_u,A_d,B_u,B_d,pz
      double precision alpha_run


C Functions:
      integer GetBinIndex
      integer GetInfoIndex
      double precision DSIMPS
      double precision AEMRUN
C-----------------------------------------------------------------

C
C Electron couplings 
C
      ve = -0.5d0 + 2.*sin2thw
      ae = -0.5d0                  

C
C u-type and d-type couplings:
C
      au = 0.5d0
      ad = -0.5d0
      vu = au - (4.d0/3.d0)*sin2thw
      vd = ad + (2.d0/3.d0)*sin2thw

      if (NDATAPOINTS(IDataSet).gt.NQ2Max) then
         print *,'ERROR IN GetIntegratedNCXsection'
         print *,'INCREASE NQ2Max to ',NDATAPOINTS(IDataSet)
         stop
      endif

C Get indexes for Q2, x1 and x2:
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX1 = GetBinIndex(IDataSet,'xmin')
      idxX2 = GetBinIndex(IDataSet,'xmax')

      if (idxX1.eq.0 .or. idxX2.eq.0 .or. idxQ2.eq.0) then
         Return
      endif


C prepare bins:
      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)

         xmin = AbstractBins(idxX1,idx)
         xmax = AbstractBins(idxX2,idx)
         q2v  = AbstractBins(idxQ2,idx)

         do j=0,NSplit
            q2(j,i) = q2v
            x(j,i)  = (xmax-xmin)/NSplit * j + xmin
         enddo

      enddo

      Npt = (NSplit+1)*NDATAPOINTS(IDataSet)


C QCDNUM, caclulate FL, F2 and xF3 for all bins:
C 
C Initialise polarisation, in case info is not provided.
C
      polarity=0.d0

      polarity = DATASETInfo( GetInfoIndex(IDataSet,'e polarity'), IDataSet)
      charge = DATASETInfo( GetInfoIndex(IDataSet,'e charge'), IDataSet)
      S = (DATASETInfo( GetInfoIndex(IDataSet,'sqrt(S)'), IDataSet))**2


      CALL ZMSTFUN(1,CNEP2F,X,Q2,FLp,Npt,0)
      CALL ZMSTFUN(2,CNEP2F,X,Q2,F2p,Npt,0)
      CALL ZMSTFUN(3,CNEP3F,X,Q2,XF3p,Npt,0)    
      CALL ZMSTFUN(1,CNEM2F,X,Q2,FLm,Npt,0)
      CALL ZMSTFUN(2,CNEM2F,X,Q2,F2m,Npt,0)
      CALL ZMSTFUN(3,CNEM3F,X,Q2,XF3m,Npt,0) 


      do i=1,NDATAPOINTS(IDataSet)
         
         PZ = 4.d0 * sin2thw * cos2thw * (1.+Mz**2/Q2(0,i))
         PZ = 1./Pz


         if (charge.gt.0) then
            A_u = e2u           ! gamma
     $           + (-ve-polarity*ae)*PZ*2.*euq*vu !gamma-Z
     $           + (ve**2 + ae**2+2*polarity*ve*ae)*PZ**2*(vu**2+au**2) !Z
           
            A_d = e2d 
     $           + (-ve-polarity*ae)*PZ*2.*edq*vd 
     $           + (ve**2 + ae**2+2*polarity*ve*ae)*PZ**2*(vd**2+ad**2)
            
            B_u = (ae+polarity*ve)*PZ*2.*euq*au !gamma-Z
     $           + (-2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vu*au !Z
            B_d = (ae+polarity*ve)*PZ*2.*edq*ad 
     $           + (-2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vd*ad


         else
            A_u = e2u           ! gamma
     $           + (-ve+polarity*ae)*PZ*2.*euq*vu !gamma-Z
     $           + (ve**2 + ae**2-2*polarity*ve*ae)*PZ**2*(vu**2+au**2) !Z
           
            A_d = e2d 
     $           + (-ve+polarity*ae)*PZ*2.*edq*vd 
     $           + (ve**2 + ae**2-2*polarity*ve*ae)*PZ**2*(vd**2+ad**2)
            
            B_u = (-ae+polarity*ve)*PZ*2.*euq*au !gamma-Z
     $           + (2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vu*au !Z
            B_d = (-ae+polarity*ve)*PZ*2.*edq*ad 
     $           + (2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vd*ad
            
         endif


cv in case polarity=0 it reduces to: (eminus case)
cv         A_u = e2u - ve*PZ*2.*euq*vu +(ve**2 + ae**2)*PZ**2*(vu**2+au**2)
cv         A_d = e2d - ve*PZ*2.*edq*vd +(ve**2 + ae**2)*PZ**2*(vd**2+ad**2)
cv         B_u = -ae*PZ*2.*euq*au + 2.*ve*ae*(PZ**2)*2.*vu*au
cv         B_d = -ae*PZ*2.*edq*ad + 2.*ve*ae*(PZ**2)*2.*vd*ad


         alpha_run = AEMRUN(q2(0,i))

C Get x-sections:
         do j=0,NSplit
            y = q2(j,i)/(S*X(j,i))
            yplus  = 1+(1-y)**2
            yminus = 1-(1-y)**2


            XF3  = B_U*XF3p(j,i)  + B_D*XF3m(j,i)
            F2   = A_U*F2p(j,i)   + A_D*F2m(j,i)
            FL   = A_U*FLp(j,i)   + A_D*FLm(j,i)


C polarisation already taken into account in the coupling
            XSec(j) = yplus*F2 + yminus*xF3 - y*y*FL

cv to get back to the old unpolarised case then one has to uncomment this:
cv         if (charge.gt.0) then            
cv            XSec(j) = yplus*F2 - yminus*xF3 - y*y*FL
cv         else
cv            XSec(j) = yplus*F2 + yminus*xF3 - y*y*FL
cv         endif

            XSec(j) = XSec(j) * (2*pi*alpha_run**2)/ 
     $           (q2(j,i)**2*x(j,i))*convfac

         enddo


C Integrate:
         xsint = DSIMPS(XSec(0),x(0,i),x(NSplit,i),NSplit)


         idx =  DATASETIDX(IDataSet,i)
         THEO(idx) =  xsint

      enddo


C-----------------------------------------------------------------
      end


      Subroutine GetReducedNCXsection(IDataSet)
C----------------------------------------------------------------
C
C  NC Double differential reduced cross section calculation for dataset IDataSet
C  Fills global array THEO.
C
C  Created by SG, 25/05/2011
C  Start with zero mass implementation
C                 14/06/2011 : re-introduce RT code
C---------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'couplings.inc'
      include 'qcdnumhelper.inc'
      include 'fcn.inc'
      include 'polarity.inc'

      integer IDataSet

      double precision ve,ae,au,ad,vu,vd,A_u,A_d,B_u,B_d,pz

      integer idxQ2, idxX, idxY, i,  idx
      
      integer NPmax
      parameter(NPmax=1000)

      double precision X(NPmax),Y(NPmax),Q2(NPmax)
      double precision yplus, yminus

      double precision F2p(NPmax),xF3p(NPmax),FLp(NPmax)
      double precision F2m(NPmax),xF3m(NPmax),FLm(NPmax)
      
      double precision F2,xF3,FL,XSec

      double precision Xv,Yv,Q2v

      double precision Charge, S, polarity
      double precision err_pol_unc, shift_pol
      double precision err_pol_corL
      double precision err_pol_corT

      double precision FLGamma, F2Gamma
C RT:
      Double precision f2pRT,flpRT,f1pRT,rpRT,f2nRT,flnRT,f1nRT,rnRT,
     $     f2cRT,flcRT,f1cRT,f2bRT,flbRT,f1bRT, F2rt, FLrt
      logical UseKFactors

C HF:
      double precision NC2FHF(-6:6)
      double precision FLc, F2c, FLb, F2b
      dimension FLc(NPmax),F2c(NPmax),FLb(NPmax),F2b(NPmax)

C Functions:
      integer GetBinIndex
      integer GetInfoIndex

c EW param

      double precision sin2th_eff, xkappa, epsilon
      double precision deltar,sweff, sin2thw2
      double precision cau, cad, cvu, cvd


c ACOT 
      double precision f123l(4),f123lc(4),f123lb(4),f2nc
      integer icharge,ij
      
C---------------------------------------------------------

      if (IFlagFCN.eq.1) then
C
C Execute for the first iteration only.
C
         if (HFSCHEME.eq.22.or.HFSCHEME.eq.11.or.HFSCHEME.eq.1) then
            UseKFactors = .true.
         else
            UseKFactors = .false.
         endif
      endif

C
C Protect against overflow of internal arrays:
C

      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN GetReducedNCXsection'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)
         stop
      endif


C
C Get indexes for Q2, x and y bins:
C
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX  = GetBinIndex(IDataSet,'x')
      idxY = GetBinIndex(IDataSet,'y')

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
      enddo

C
C Get polarity, charge and CME information:
C     

C 
C Initialise polarisation, in case info is not provided.
C
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
      S = (DATASETInfo( GetInfoIndex(IDataSet,'sqrt(S)'), IDataSet))**2



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


         if (idataset.eq.2) then
            print '(''cau: cvu: cad: cvd'',5F16.6)',cau,cvu,cad,cvd,deltar
         endif


         sin2thw2 = 1.d0 - MW**2/MZ**2
         sin2th_eff = 0.23134d0
         xkappa = sin2th_eff/sin2thw
         epsilon = xkappa -1.0
c         print*, sin2thw, sweff, xkappa, epsilon
         ve = -0.5d0 + 2.*sin2th_eff
         ae = -0.5d0

         vu = cvu - (4.d0/3.d0)*epsilon*sin2thw2
         vd = cvd + (2.d0/3.d0)*epsilon*sin2thw2
         au = cau
         ad = cad

      endif


C QCDNUM ZMVFNS, caclulate FL, F2 and xF3 for d- and u- type quarks all bins:

C u-type ( u+c ) contributions 
      CALL ZMSTFUN(1,CNEP2F,X,Q2,FLp,NDATAPOINTS(IDataSet),0)
      CALL ZMSTFUN(2,CNEP2F,X,Q2,F2p,NDATAPOINTS(IDataSet),0)
      CALL ZMSTFUN(3,CNEP3F,X,Q2,XF3p,NDATAPOINTS(IDataSet),0)    

C d-type (d + s + b) contributions
      CALL ZMSTFUN(1,CNEM2F,X,Q2,FLm,NDATAPOINTS(IDataSet),0)
      CALL ZMSTFUN(2,CNEM2F,X,Q2,F2m,NDATAPOINTS(IDataSet),0)
      CALL ZMSTFUN(3,CNEM3F,X,Q2,XF3m,NDATAPOINTS(IDataSet),0) 

C heavy quark contribution (c and b)      
      if (mod(HFSCHEME,10).eq.3) then 

         NC2FHF = 4.D0/9.D0 * CNEP2F  + 1.D0/9.D0 * CNEM2F
c        write(6,*) 'NC2FHF:', NC2FHF

         CALL HQSTFUN(2,1,NC2FHF,X,Q2,F2c,NDATAPOINTS(IDataSet),0)
         CALL HQSTFUN(1,1,NC2FHF,X,Q2,FLc,NDATAPOINTS(IDataSet),0)
         CALL HQSTFUN(2,-2,NC2FHF,X,Q2,F2b,NDATAPOINTS(IDataSet),0)
         CALL HQSTFUN(1,-2,NC2FHF,X,Q2,FLb,NDATAPOINTS(IDataSet),0)
c         write(6,*) 'HQSTFUN: X,Q2,F2c,FLc', X(1),Q2(1),F2c(1),FLc(1)
c         write(6,*) 'HQSTFUN: X,Q2,F2b,FLb', X(1),Q2(1),F2b(1),FLb(1)

      endif
C
C Prepare theory prediction for chi2 calculation:
C
      do i=1,NDATAPOINTS(IDataSet)

C Get the index of the point in the global data table:
         idx =  DATASETIDX(IDataSet,i)

C Propagator factor PZ
         PZ = 4.d0 * sin2thw * cos2thw * (1.+Mz**2/Q2(i))

C modify propagator for EW corrections -- only for EW fit
         if (EWfit.ne.0) PZ = PZ * (1.d0 - Deltar)

         PZ = 1./Pz
C EW couplings of u-type and d-type quarks at the scale Q2

         if (charge.gt.0) then
            A_u = e2u           ! gamma
     $           + (-ve-polarity*ae)*PZ*2.*euq*vu !gamma-Z
     $           + (ve**2 + ae**2+2*polarity*ve*ae)*PZ**2*(vu**2+au**2) !Z
            
            A_d = e2d 
     $           + (-ve-polarity*ae)*PZ*2.*edq*vd 
     $           + (ve**2 + ae**2+2*polarity*ve*ae)*PZ**2*(vd**2+ad**2)
            
            B_u = (ae+polarity*ve)*PZ*2.*euq*au !gamma-Z
     $           + (-2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vu*au !Z
            B_d = (ae+polarity*ve)*PZ*2.*edq*ad 
     $           + (-2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vd*ad
         else
            A_u = e2u           ! gamma
     $           + (-ve+polarity*ae)*PZ*2.*euq*vu !gamma-Z
     $           + (ve**2 + ae**2-2*polarity*ve*ae)*PZ**2*(vu**2+au**2) !Z
            
            A_d = e2d 
     $           + (-ve+polarity*ae)*PZ*2.*edq*vd 
     $           + (ve**2 + ae**2-2*polarity*ve*ae)*PZ**2*(vd**2+ad**2)
            
            B_u = (-ae+polarity*ve)*PZ*2.*euq*au !gamma-Z
     $           + (2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vu*au !Z
            B_d = (-ae+polarity*ve)*PZ*2.*edq*ad 
     $           + (2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vd*ad



         endif

cv for polarised case should reduce to:
cv         A_u = e2u - ve*PZ*2.*euq*vu +(ve**2 + ae**2)*PZ**2*(vu**2+au**2)
cv         A_d = e2d - ve*PZ*2.*edq*vd +(ve**2 + ae**2)*PZ**2*(vd**2+ad**2)
cv         B_u = -ae*PZ*2.*euq*au + 2.*ve*ae*(PZ**2)*2.*vu*au
cv         B_d = -ae*PZ*2.*edq*ad + 2.*ve*ae*(PZ**2)*2.*vd*ad

C Get x-sections:
         yplus  = 1+(1-y(i))**2
         yminus = 1-(1-y(i))**2

C
C xF3, F2, FL from QCDNUM:
C
         XF3  = B_U*XF3p(i)  + B_D*XF3m(i)
         F2   = A_U*F2p(i)   + A_D*F2m(i)
         FL   = A_U*FLp(i)   + A_D*FLm(i)

C-----------------------------------------------------------------------
C  Extra heavy flavour schemes
C


c
C ACOT scheme 
C                     
         if (mod(HFSCHEME,10).eq.1) then
            call sf_acot_wrap(x(i),q2(i),
     $           f123l,f123lc,f123lb,
     $           hfscheme, 4, 
     $           iFlagFCN, idx,
     $           UseKFactors)
            
            FL  = F123L(4)
            F2  = F123L(2)
            
            if (charge.gt.0) then
               XF3 = - x(i)*F123L(3)
            else
               XF3 = x(i)*F123L(3)
            endif
            

C RT scheme 
C            
         elseif (mod(HFSCHEME,10).eq.2) then 


C RT does not provide terms beyond gamma exchange. Since they occur at high Q2,
C use QCDNUM to take them into account as a "k"-factor 
C
C  F2_total^{RT} =  F2_{\gamma}^{RT}  *  (  F2_{total}^{QCDNUM}/F2_{\gamma}^{QCDNUM}   
C

           F2Gamma = 4.D0/9.D0 * F2p(i)  + 1.D0/9.D0 * F2m(i)
           FLGamma = 4.D0/9.D0 * FLp(i)  + 1.D0/9.D0 * FLm(i)


           call sfun_wrap(x(i),q2(i)
     $          ,f2pRT,flpRT,f1pRT,
     +          rpRT,f2nRT,flnRT,
     +          f1nRT,rnRT,f2cRT,
     +          flcRT,f1cRT,f2bRT,
     +          flbRT,f1bRT
           ! Input:
     $          ,iFlagFCN,idx    ! fcn flag, data point index
     $          ,F2Gamma,FLGamma
     $          ,UseKFactors
     $          )
           



           F2rt = F2pRT * (F2/F2Gamma)
           FLrt = FLpRT * (FL/FLGamma)


C Replace F2,FL from QCDNUM by RT values
C Keep xF3 from QCDNUM

           F2 = F2rt
           FL = FLrt
        endif


C FFNS, heavy quark contribution (c and b) to F2 and FL   
        if (mod(HFSCHEME,10).eq.3) then 

c           write(6,*) 'FFNS before sum: F2', i,F2
           F2 = F2 + F2c(i) + F2b(i) 
           FL = FL + FLc(i) + FLb(i)
c           write(6,*) 'FFNS: F2,F2c,F2b', i,F2,F2c(i),F2b(i)

        endif
         

C polarisation already taken into account in the couplings (A_q,B_q)

        XSec = F2 + yminus/yplus*xF3 - y(i)*y(i)/yplus*FL

cv to get back to the old unpolarised case then one has to uncomment this:
cv        if (charge.gt.0) then
cv
cv           XSec = F2 - yminus/yplus*xF3 - y(i)*y(i)/yplus*FL
cv
cv       else
cv           XSec = F2 + yminus/yplus*xF3 - y(i)*y(i)/yplus*FL
cv        endif



C
C Store cross-section prediction in the global cross-sections table:
C
        THEO(idx) =  XSec
      enddo

      end
