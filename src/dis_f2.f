      Subroutine GetF2p(IDataSet)
C----------------------------------------------------------------
C
C  F2 structure function calculation for dataset IDataSet
C  Fills global array THEO.
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

      integer IDataSet

      double precision ve,ae,au,ad,vu,vd,A_u,A_d,B_u,B_d,pz

      integer idxQ2, idxX, idxY, i,  idx
      
      integer NPmax
      parameter(NPmax=1000)

      double precision X(NPmax),Y(NPmax),Q2(NPmax)
      double precision yplus, yminus

      double precision F2p(NPmax)
      double precision F2m(NPmax)
      
      double precision F2

      double precision Xv,Yv,Q2v

      double precision Charge, S, polarity


      double precision FLGamma, F2Gamma
C RT:
      Double precision f2pRT,flpRT,f1pRT,rpRT,f2nRT,flnRT,f1nRT,rnRT,
     $     f2cRT,flcRT,f1cRT,f2bRT,flbRT,f1bRT, F2rt, FLrt
      logical UseKFactors


C Functions:
      integer GetBinIndex
      integer GetInfoIndex

C---------------------------------------------------------

      if (IFlagFCN.eq.1) then
C
C Execute for the first iteration only.
C
         if (HFSCHEME.eq.22) then
            UseKFactors = .true.
         else
            UseKFactors = .false.
         endif
      endif

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
 
C
C Protect against overflow of internal arrays:
C

      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN GetF2p'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)
         call HF_stop
      endif


C
C Get indexes for Q2, x and y bins:
C
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX  = GetBinIndex(IDataSet,'x')
c      idxY = GetBinIndex(IDataSet,'y')

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
C         Y(i)   = AbstractBins(idxY,idx)
         Q2(i)  = AbstractBins(idxQ2,idx)
      enddo

C
C Get polarity, charge and CME information:
C     

C 
C Initialise polarisation, in case info is not provided.
C
      polarity=0.d0

      polarity = DATASETInfo( GetInfoIndex(IDataSet,'e polarity'),
     $     IDataSet)
      charge = DATASETInfo( GetInfoIndex(IDataSet,'e charge'), IDataSet)
c      S = (DATASETInfo( GetInfoIndex(IDataSet,'sqrt(S)'), IDataSet))**2

C QCDNUM ZMVFNS, caclulate FL, F2 and xF3 for d- and u- type quarks all bins:

C u-type ( u+c ) contributions 

      CALL ZMSTFUN(2,CNEP2F,X,Q2,F2p,NDATAPOINTS(IDataSet),0)
c      print *,'F2p', F2p

C d-type (d + s + b) contributions

      CALL ZMSTFUN(2,CNEM2F,X,Q2,F2m,NDATAPOINTS(IDataSet),0)


C
C Prepare theory prediction for chi2 calculation:
C
      do i=1,NDATAPOINTS(IDataSet)

C Get the index of the point in the global data table:
         idx =  DATASETIDX(IDataSet,i)

C Propagator factor PZ
         PZ = 4.d0 * sin2thw * cos2thw * (1.+Mz**2/Q2(i))
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

C
C  F2 from QCDNUM:C
        
         F2   = A_U*F2p(i)   + A_D*F2m(i)
ch         print *,'F2p', F2

C-----------------------------------------------------------------------
C  Extra heavy flavour schemes
C

c
C RT scheme 
C
        if (mod(HFSCHEME,10).eq.2) then 

C RT does not provide terms beyond gamma exchange. Since they occur at high Q2,
C use QCDNUM to take them into account as a "k"-factor 
C
C  F2_total^{RT} =  F2_{\gamma}^{RT}  *  (  F2_{total}^{QCDNUM}/F2_{\gamma}^{QCDNUM}   
C

           F2Gamma = 4.D0/9.D0 * F2p(i)  + 1.D0/9.D0 * F2m(i)


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
ch          print *,'F2p', F2pRT


C Replace F2,FL from QCDNUM by RT values
C Keep xF3 from QCDNUM

           F2 = F2rt
        endif

         


        THEO(idx) =  F2
ch        print *,'F2', F2
      enddo

      end
