      Subroutine PDF_param_iteration(p,iflag)
C-------------------------------------------------------
C
C Created 5 June 2011. Move PDF parameterisation setting from FCN 
C
C  Input:  p(*)  -- input minuit parameters
C          iflag -- minuit flag 
C
C--------------------------------------------------------
      implicit none


      double precision p(*)
      integer iflag
#include "pdfparam.inc"
#include "steering.inc"
#include "alphas.inc"
#include "thresholds.inc"
#include "extrapars.inc"
#include "polarity.inc"
#include "couplings.inc"
      integer i

      double precision fs,rs
      double precision fshermes
      double precision alphasPDF, StepAlphaS, StepFs
      double precision Aph, Bph, Cph, Dph, Eph
      double precision StepAph, StepBph, StepCph, StepDph, StepEph
      	
C-------------------------------------------------------
      logical LFirstTime
      data LFirstTime /.true./
      integer idxAlphaS, idxRS, idxFCharm, idxFS   !> indices for alphas and fs
      integer idphA, idphB, idphC, idphD, idphE !> indices for photon Pdf
      data idphA/0/
      data idphB/0/
      data idphC/0/
      data idphD/0/
      data idphE/0/

      integer GetParameterIndex  !> function to read parameter index
C-------------------------------------------------------
      integer idxShiftPolLHp   !> indices shiftpol
      integer idxShiftPolRHp   !> indices shiftpol
      integer idxShiftPolLHm   !> indices shiftpol
      integer idxShiftPolRHm    !> indices shiftpol
      integer idxShiftPolT   !> indices shiftpol
      integer idxShiftPolL   !> indices shiftpol
      
C-------------------------------------------------------
      !> Additional scaling parameter "temperature"
      integer idxTemperature
      data idxTemperature/0/

C-------------------------------------------------------
      integer idxAuEW, idxAdEW, idxVuEW, idxVdEW !> indices for EW param
      
      integer idxVcs

      logical LPolFits       !> Logical to init polarisation fits
      data LPolFits/.false./

C-------------------------------------------------------
      if(ITheory.ge.100) return 
      if(Itheory.eq.50) then
         call DecodeFractal(p) 
         return
      endif

      if (LFirstTime) then    
         LFirstTime = .false.
         if (ITheory.eq.11.or.ITheory.eq.35) then
            idphA = GetParameterIndex('Aph')
            idphB = GetParameterIndex('Bph')
            idphC = GetParameterIndex('Cph')
            idphD = GetParameterIndex('Dph')
            idphE = GetParameterIndex('Eph')

            if ((idphA.eq.0).and.(idphB.eq.0).and.(idphC.eq.0).and.
     $           (idphD.eq.0).and.(idphE.eq.0)) then
               if( PDFStyle.ne.'LHAPDF'.and.PDFStyle.ne.'LHAPDFQ0') then

                  print *,'Did not find photon parameters'
                  print *,'Add to ExtraParamters: Aph,Bph,Cph,Dph,Eph'
                  Call HF_errlog(15052700,
     $                 'W: Add to ExtraParamters: Aph,Bph,Cph,Dph,Eph')
               endif

            else
               idphA = iExtraParamMinuit(idphA)
               idphB = iExtraParamMinuit(idphB)
               idphC = iExtraParamMinuit(idphC)
               idphD = iExtraParamMinuit(idphD)
               idphE = iExtraParamMinuit(idphE)

               StepAph=ExtraParamStep(idphA)
               StepBph=ExtraParamStep(idphB)
               StepCph=ExtraParamStep(idphC)
               StepDph=ExtraParamStep(idphD)
               StepEph=ExtraParamStep(idphE)
            endif
            
         endif

         idxAlphaS = GetParameterIndex('alphas')     
         if (idxAlphaS.eq.0) then
            print *,'Did not find alpha_S parameter'
            print *,'Add to ExtraParamters with the name alphas'
            call HF_stop
         else
            StepAlphaS=ExtraParamStep(idxAlphaS)
            idxAlphaS = iExtraParamMinuit(idxAlphaS)
         endif

         idxFS = GetParameterIndex('fs')
         idxRS = GetParameterIndex('rs')

         if ((idxRS.eq.0).and.(idxFS.eq.0)) then
            print *,'Did not find fs nor rs parameter'
            print *,'Add to ExtraParamters with the name rs or fs'
            Call HF_errlog(13050800,
     $           'S: Add to ExtraParamters with the name rs or fs')
         elseif ((idxRS.ne.0).and.(idxFS.ne.0)) then
            print *,'Use either rs or fs, NOT both:
     $           Both rs and fs are defined' 
            Call HF_errlog(13050801,
     $           'S: Use either rs or fs, NOT both')
         elseif (idxFS.ne.0) then
c            print*,'idxFs', idxFs, ExtraParamStep(idxFS),ExtraParamValue(idxFS)
            idxFS = iExtraParamMinuit(idxFS)
!            StepFs = ExtraParamStep(idxFS)
         elseif (idxRS.ne.0) then
            idxRS = iExtraParamMinuit(idxRS)
!            StepFs = ExtraParamStep(idxRS)
         endif


         idxFCharm = GetParameterIndex('fcharm')
         if (idxFCharm.gt.0) then
            idxFCharm = iExtraParamMinuit(idxFCharm)
         endif

         idxTemperature = GetParameterIndex('Temperature')
         if (idxTemperature.gt.0) then
            idxTemperature = iExtraParamMinuit(idxTemperature)
         endif

         idxAuEW = GetParameterIndex('auEW')
         idxAdEW = GetParameterIndex('adEW')
         idxVuEW = GetParameterIndex('vuEW')
         idxVdEW = GetParameterIndex('vdEW')

         if (idxAuEW.eq.0.or.idxAdEW.eq.0.or.
     $        idxVuEW.eq.0.or.idxVdEW.eq.0) then
         else
            idxAuEW = iExtraParamMinuit(idxAuEW)
            idxAdEW = iExtraParamMinuit(idxAdEW)
            idxVuEW = iExtraParamMinuit(idxVuEW)
            idxVdEW = iExtraParamMinuit(idxVdEW)
         endif

         idxVcs = GetParameterIndex('Vcs')
         if (idxVcs.gt.0) then
            idxVcs = iExtraParamMinuit(idxVcs)
            call hf_errlog(20112013,'I:Float Vcs')
         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         idxShiftPolLHp = GetParameterIndex('shiftpolLHp')
         idxShiftPolRHp = GetParameterIndex('shiftpolRHp')
         idxShiftPolLHm = GetParameterIndex('shiftpolLHm')
         idxShiftPolRHm = GetParameterIndex('shiftpolRHm')
         idxShiftPolL = GetParameterIndex('shiftpolL')
         idxShiftPolT = GetParameterIndex('shiftpolT')
 
        if (idxShiftPoLRHm.eq.0.or.
     $        idxShiftPolRHm.eq.0.or.
     $        idxShiftPolLHm.eq.0.or.
     $        idxShiftPolL.eq.0.or.
     $        idxShiftPolT.eq.0.or.
     $        idxShiftPolLHp.eq.0) then
         else
            idxShiftPolRHm = iExtraParamMinuit(idxShiftPolRHm)
            idxShiftPolLHm = iExtraParamMinuit(idxShiftPolLHm)
            idxShiftPolRHp = iExtraParamMinuit(idxShiftPolRHp)
            idxShiftPolLHp = iExtraParamMinuit(idxShiftPolLHp)
            idxShiftPolT = iExtraParamMinuit(idxShiftPolT)
            idxShiftPolL = iExtraParamMinuit(idxShiftPolL)
            LPolFits = .true.
         endif



      endif

C Polarisation shifts extra param

      if (LPolFits) then
         shift_polRHp=p(idxShiftPolRHp)
         shift_polLHp=p(idxShiftPolLHp)
         shift_polLHm=p(idxShiftPolLHm)
         shift_polRHm=p(idxShiftPolRHm)
         shift_polL=p(idxShiftPolL)
         shift_polT=p(idxShiftPolT)
      else
         shift_polRHp=0.0
         shift_polLHp=0.0
         shift_polLHm=0.0
         shift_polRHm=0.0
         shift_polL  =0.0
         shift_polT  =0.0
      endif




C EW extra param
      IF (idxAuEw.gt.0) then
         cau_ew=p(idxAuEW)
         cad_ew=p(idxAdEW)
         cvu_ew=p(idxVuEW)
         cvd_ew=p(idxVdEW)
      else
         cau_ew = 0.
         cad_ew = 0.
         cvu_ew = 0.
         cvd_ew = 0.
      endif


C "Temperature"
      if (idxTemperature.gt.0) then
         Temperature = p(idxTemperature)
         print *,'Temperature=',Temperature
      endif


C Get from extra pars:
      alphas=p(idxAlphaS)      
      
      if (idxFS.ne.0) then
         fstrange=p(idxFS)
      elseif (idxRS.ne.0) then
         fstrange=p(idxRS)/(p(idxRS)+1)
      endif
C In case PDF and alphas needs to be read from LHAPDF (iparam=0, ipdfset=5)
C maybe instead warning message should be issued
      
      if( PDFStyle.eq.'LHAPDF'.or.PDFStyle.eq.'LHAPDFQ0') then
         if (StepAlphaS.eq.0.) then
            alphas=alphasPDF(Mz)
            Call HF_errlog(13051401,
     $              'W: alphas is fixed and taken from LHAPDF file') 
         else
            Call HF_errlog(13051402,
     $           'W: alphas is free and taken from steering.txt  ') 
         endif
      endif
         


C Hermes strange prepare:
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(0.D0)
      endif
     

      if (q0.ge.qc.and.idxfcharm.gt.0) then
         fcharm=p(idxFCharm)
      else
         fcharm=0.
      endif


      if (idxVcs.gt.0) then
         Vcs = p(idxVcs)
      endif



!!!!!!!!!!!!!!!!!!!!!!!

C 10 Aug 2011: Standard parametrisation:
      if ((idphA.ne.0).or.(idphB.ne.0).or.(idphC.ne.0).or. 
     $    (idphD.ne.0).or.(idphE.ne.0))then
         parphoton(1)= p(idphA)
         parphoton(2)= p(idphB)
         parphoton(3)= p(idphC)
         parphoton(4)= p(idphD)
         parphoton(5)= p(idphE)
      endif
      Call DecodePara(p)


C  25 Jan 2011: Poly params for valence:
      if (NPOLYVAL.gt.0) then
         Call StorePoly(p,iflag)
      endif

C  22 Apr 2011: CT parameterisation:
      if (PDFStyle.eq.'CTEQ'.or.PDFStyle.eq.'CTEQHERA') then
         Call DecodeCtPara(p)
         if ((idphA.ne.0).or.(idphB.ne.0).or.(idphC.ne.0).or.
     $    (idphD.ne.0).or.(idphE.ne.0))then
         ctphoton(1)= p(idphA)
         ctphoton(2)= p(idphB)
         ctphoton(3)= p(idphC)
         ctphoton(4)= p(idphD)
         ctphoton(5)= p(idphE)

         endif

      endif

C  22 Sep 2011: AS parameterisation:
        if ((PDFStyle.eq.'AS').or.(PDFStyle.eq.'BiLog')) then
         Call DecodeASPara(p)
      endif



C
C Chebyshev for the gluon:
C
      if (NCHEBGLU.gt.0) then
         do i=1,NCHEBGLU
            ChebPars(i) = p(30+i)
         enddo
      endif
C
C  Chebyshev param. for the sea:
C
      if (NCHEBSea.gt.0) then
         do i=1,NCHEBSea
C  Offset is now steering parameter (default = 70, params start from 41)
            ChebParsSea(i) = p(30+IOFFSETCHEBSEA+i)
         enddo
      endif

      if (NChebGlu.gt.0 .or. NChebSea.gt.0) then
         call ChebToPoly
      endif


C  22 Nov 2011: dipole model parameters
      if (DipoleModel.gt.0) then
         call DecodeDipolePar(p)
      endif

      end


* -------------------------------------------------------
      double precision function flav_number(q)
* -------------------------------------------------------

      implicit none
#include "thresholds.inc"

      double precision q
      
      flav_number = 3.d0
      if (q.ge.qc) flav_number = 4.d0
      if (q.ge.qb) flav_number = 5.d0

      return
      end


      subroutine DecodePara(pars)
C-------------------------------------------------------
C Created 22 Apr 11 by SG. Decode minuit input for CTEQ-like param.
C   pars(1-10)  - gluon
C   pars(11-20)  - Uv
C   pars(21-30)  - Dv
C   pars(31-40)  - Ubar, U
C   pars(41-50)  - Dbar, D
C   pars(51-60)  - sea, delta
C   pars(91-100)  - others
C------------------------------------------------------
      implicit none 
#include "pdfparam.inc"
#include "steering.inc"
#include "for_debug.inc"
      double precision pars(*)
      integer i,j
      logical lfirstt
      data lfirstt /.true./

      double precision fs
      double precision fshermes
C---------------------------------------------------------
      if (lfirstt) then
         lfirstt = .false.
         print *,'DecodePara INFO: First time call'        
      endif

C Hermes strange prepare:
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(0.D0)
      endif

C     simple copy first:
      do i=1,10
         parglue(i) = pars(i)
         paruval(i) = pars(10+i)
         pardval(i) = pars(20+i)
         parubar(i) = pars(30+i)
         pardbar(i) = pars(40+i)
         paru(i) = pars(50+i)
         pard(i) = pars(60+i)
         parsea(i) = pars(70+i)
         parstr(i) = pars(80+i)
         parother(i) = pars(90+i)
      enddo

      if (PDF_DECOMPOSITION.eq.'D_U_Dbar_Ubar') then     !  H1PDF2k like

         if (pard(2).eq.0)    pard(2)=paru(2)
         if (parubar(2).eq.0) parubar(2)=paru(2)
         if (pardbar(2).eq.0) pardbar(2)=paru(2)
         if (pardbar(1).eq.0) pardbar(1)=pard(1)
         if (paru(1).eq.0)    parU(1)=pard(1)*(1.D0-fs)/(1.D0-fcharm)
         if (parUbar(1).eq.0) parUbar(1)=parU(1)
cv        elseif (iparam.eq.2) then
cv           if (pardval(2).eq.0)   pardval(2)=paruval(2) !  Bud    = Buv 
cv           if (parubar(2).eq.0)   parubar(2)=pardbar(2) !  Bubar  = Bdbar
cv         pardval(2)=paruval(2)
cv         parubar(2)=pardbar(2)
cv         parUbar(1)=pardbar(1)*(1.D0-fs)/(1.D0-fcharm)


!> this style is common to HERAPDF, ATLASPDF:
      elseif (index(PDF_DECOMPOSITION,'Dv_Uv_Dbar_Ubar_Str').ne.0) then

         if (pardval(2).eq.0)   pardval(2)=paruval(2) !  Bud    = Buv 
         if (parubar(2).eq.0)   parubar(2)=pardbar(2)  !  Bubar  = Bdbar 
         if (parstr(1).eq.0.and.
     $        parstr(2).eq.0.and.
     $        parstr(3).eq.0) then
            
!> use coupled strange to Dbar
            FreeStrange=.false.  
            
         else
            FreeStrange=.true.
         endif
!> couple Bstr and Cstr to dbar when zero:
         if (FreeStrange) then
            if (parstr(2).eq.0.and.parstr(3).ne.0) parstr(2)=pardbar(2)
            if (parstr(3).eq.0.and.parstr(2).ne.0) parstr(3)=pardbar(3)
            
         endif
         
         if (fs.ne.-10000.and.(FreeStrange)) then
!> then use ubar and dbar (not Dbar and Ubar)
            parstr(1)=fs/(1.-fs)*pardbar(1)
            if (parubar(1).eq.0) parubar(1) = pardbar(1)
         else                  
!> then use Dbar and Ubar
            if (parubar(1).eq.0)   parubar(1)=pardbar(1)*(1.D0-fs) 
     $           /(1.D0-fcharm) !then use Ubar=Dbar
         endif



c      elseif (iparam.eq.3) then ! g,uval,dval,sea as in ZEUS-S 2002 fit
c         
c         paruval(2)=0.5
c         pardval(2)=0.5
c         parstr(2)=0.5
c         parstr(3)=parsea(3)+2.
         
      elseif (PDFStyle.eq.'CHEB'.or.PDFStyle.eq.'ZEUS Jet') then
         pardval(2)=paruval(2)
        

*  dbar-ubar (not Ubar - Dbar), Adel fixed to output of ZEUS-S fit   
         parstr(1)=0.27
         parstr(2)=0.5
         parstr(3)=parsea(3)+2.

c      elseif (iparam.eq.24) then ! g,uval,dval,sea as in ZEUS-JET fit
c         pardval(2)=paruval(2)
c         parstr(1)=0.27
c         parstr(2)=0.5
c         parstr(3)=parsea(3)+2.
      endif         

      if (debug) then
         print '(''1uv:'',11F10.4)',(paruval(i),i=1,10)
         print '(''1dv:'',11F10.4)',(pardval(i),i=1,10)
         print '(''1Ub:'',11F10.4)',(parubar(i),i=1,10)
         print '(''1Db:'',11F10.4)',(pardbar(i),i=1,10)
         print '(''1GL:'',11F10.4)',(parglue(i),i=1,10)
         print '(''1ST:'',11F10.4)',(parstr(i),i=1,10)
          if (ITheory.eq.11) then
             print '(''1PH:'',11F10.4)',(parphoton(i),i=1,10)
          endif
      endif

C---------------------------------------------------------
      end

C---------------------------------------------------------

      double precision function para(x,a)
C----------------------------------------------------
C
C standard-like parameterisation: 
C  AF = (a*x**b)*(1 - x)**c*(1 + d*x + e*x**2+f*x**3)-
C     - (ap*x**bp)*(1-x)**cp
C
C-----------------------------------------------------
      implicit none
      double precision x,a(1:10)
      double precision AF
      
      AF = a(1)*x**a(2)*(1 - x)**a(3)*(1 + a(4)*x
     $     + a(5)*x**2+a(6)*x**3+a(10)*x**0.5)-a(7)*x**a(8)*(1-x)**a(9)
      
      para = AF

      end



      subroutine DecodeCtPara(pars)
C-------------------------------------------------------
C Created 22 Apr 11 by SG. Decode minuit input for CTEQ-like param.
C   pars(1-6)  - gluon
C   pars(11-16)  - Uv
C   pars(21-26)  - Dv
C   pars(31-36)  - Ubar
C   pars(41-46)  - Dbar
C   pars(81-86)  - Str
C   pars(95-100)  - alphas, fstrange, fcharm
C------------------------------------------------------
      implicit none 
#include "pdfparam.inc"
#include "steering.inc"
      double precision pars(*)
      integer i
      logical lfirstt
      data lfirstt /.true./

      double precision fs
      double precision fshermes
C---------------------------------------------------------
      if (lfirstt) then
         lfirstt = .false.
         print *,'DecodeCtPara INFO: First time call'        
      endif

C Hermes strange prepare:
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(0.D0)
      endif

C simple copy first:
      do i=1,9
         ctglue(i) = pars(i)
         ctuval(i) = pars(10+i)
         ctdval(i) = pars(20+i)
         ctubar(i) = pars(30+i)
         ctdbar(i) = pars(40+i)
cv add str
         ctstr(i)  = pars(80+i)
         ctother(i)= pars(90+i)
      enddo

c      UF = a(1)*exp(a(4)*x)*(1 - x)**a(3)*x**(a(2))*(1 + exp(a(5))*x 
c     $     + exp(a(6))*x**2)


C Extra constrains:
      if (ctubar(2).eq.0) ctubar(2) = ctdbar(2)  ! Bubar = Bdbar
      if (ctuval(2).eq.0) ctuval(2) = ctdval(2)  ! Buv = Bdv

!> use ubar and dbar (not Dbar and Ubar)
      ctstr(1)=fs/(1.-fs)*ctdbar(1)
      if (ctubar(1).eq.0) ctubar(1) = ctdbar(1)

!> use coupled strange to Dbar ! 
      if (ctstr(2).eq.0) ctstr(2)=ctdbar(2)
      if (ctstr(3).eq.0) ctstr(3)=ctdbar(3)
      if (ctstr(4).eq.0) ctstr(4)=ctdbar(4)
            

C (other constraints from sum-rules)


C---------------------------------------------------------
      end

      double precision function ctpara(x,a)
C----------------------------------------------------
C
C cteq-like parameterisation: 
C  UF = a0*E**(a3*x)*(1 - x)**a2*x**(a1 + n)*(1 + E**a4*x + E**a5*x**2)
C
C-----------------------------------------------------
      implicit none
      double precision x,a(1:9)
      double precision UF
      UF = a(1)*exp(a(4)*x)*(1 - x)**a(3)*x**(a(2))*(1 + exp(a(5))*x 
     $     + exp(a(6))*x**2)-a(7)*x**a(8)*(1-x)**a(9)
      
      ctpara = UF

      end

      double precision function ctherapara(x,a)
C----------------------------------------------------
C
C hybrid cteq-hera parameterisation: 
c  UF = (A*x**B)*(1 - x)**C * exp(A4*x) * (1 + D*x + E*x**2) 
c     -AP*x**BP*(1-x)**CP
C
C-----------------------------------------------------
      implicit none
      double precision x,a(1:9)
      double precision UF
      UF = a(1)*(1 - x)**a(3)*x**(a(2))*exp(a(6)*x)*(1 + a(4)*x 
     $     + a(5)*x**2)-a(7)*x**a(8)*(1-x)**a(9)
      
      ctherapara = UF

      end

      subroutine DecodeASPara(pars)
C-------------------------------------------------------
C Created 20 Jul 11 by VR. Decode minuit input for AS param.
C   pars(21-25)  - Uv
C   pars(32-35)  - Dv
C   pars(43-45)  - Ubar
C   pars(51-55)  - Dbar
C   pars(1-5)  - gluon
C------------------------------------------------------
      implicit none 
#include "pdfparam.inc"
#include "steering.inc"
      double precision pars(*)
      integer i
      logical lfirstt
      data lfirstt /.true./
      double precision fs
      double precision fshermes
C---------------------------------------------------------
      if (lfirstt) then
         lfirstt = .false.
         print *,'DecodeASPara INFO: First time call'        
      endif


C Hermes strange prepare:
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(0.D0)
      endif

C simple copy first:
      do i=1,5
         asglue(i) = pars(i)
         asuval(i) = pars(10+i)
         asdval(i) = pars(20+i)
         asubar(i) = pars(30+i)
         asdbar(i) = pars(40+i)
         asother(i) = pars(94+i)
      enddo

C Extra constrains:
      if (pars(31).eq.0) then
         asubar(1) = asdbar(1) * (1.D0-fs)/(1.D0-fcharm) !
         asubar(2) = asdbar(2)  ! Bubar = Bdbar
         asubar(3) = asdbar(3)  ! Bubar = Bdbar
      endif

C Impose Buv = Bdv if parameter for Buv = 0.
c      if (pars(12).eq.0) then
c         asuval(2) = asdval(2)  ! Buv = Bdv
c      endif

c      print '(''2uv:'',5F10.4)',(asuval(i),i=1,5)
c      print '(''2dv:'',5F10.4)',(asdval(i),i=1,5)
c      print '(''2Ub:'',5F10.4)',(asubar(i),i=1,5)
c      print '(''2Db:'',5F10.4)',(asdbar(i),i=1,5)
c      print '(''2GL:'',5F10.4)',(asglue(i),i=1,5)



C---------------------------------------------------------
      end


      double precision function splogn(x,a)
C----------------------------------------------------
c     Special lognormal function 
c
c
c     A.Schoening, University Heidelberg, Physikalisches Institut
c     Creation: 12.6.2011
c   A1*x**(A2-A3*log(x))*(1-x)**(A4-A5*log(1-x))
C-----------------------------------------------------
      implicit none
      double precision x,a(1:5)
      double precision  splogn1
      
      splogn1=0.0d0
      if (x.gt.0.d0.and.x.lt.1.d0) then
         splogn1=A(1)*x**(A(2)-A(3)*log(x))*
     $        (1.d0-x)**(A(4)-A(5)*log(1.d0-x))
      endif


      if (abs(splogn1).lt.1d30 .and. abs(splogn1).gt.1d-30) then
c value in allowed range
         splogn=splogn1
      else
         splogn=0.0d0
      endif

      return
      end

      




* -------------------------------------------------------
      double precision function gluon(x)
* -------------------------------------------------------
* x *g(x,Q2)

      implicit none
#include "pdfparam.inc"
#include "steering.inc"
      double precision x
      integer i
C External function:
      double precision PolyParam,ctpara,ctherapara,para,splogn
C-------------------------------------------------


C    22 Apr 11, SG, Add CTEQ-like
      if (PDFStyle.eq.'CTEQ') then
         gluon = ctpara(x,ctglue)
         return
      endif

      if (PDFStyle.eq.'CTEQHERA') then
         gluon = ctherapara(x,ctglue)
         return
      endif

C    22 Sept 11, VR, Add AS
      if ((PDFStyle.eq.'AS').or.(PDFStyle.eq.'BiLog')) then
         gluon = splogn(x,asglue)
         return
      endif
      if (nchebglu.eq.0) then

!> HERAPDF style goes in here:
         gluon=para(x,parglue)
      else
C
C Use polynomial representation of cheb. 
C
         gluon = parglue(1) * PolyParam(x,nchebGlu,polyPars,chebxminlog)
         if (ichebtypeGlu.eq.0) then
C Do nothing
         else if (ichebtypeGlu.eq.1) then
            gluon = gluon * (1 - x) ! force PDFs=0 for x=1
         endif
         
      endif
      end

* -------------------------------------------------------
      double precision function photon(x)
* -------------------------------------------------------
* x *photon(x,Q2)

      implicit none
#include "pdfparam.inc"
#include "steering.inc"
      double precision x
      integer i
C External function:
      double precision ctpara,ctherapara,para
C-------------------------------------------------


C    22 Apr 11, SG, Add CTEQ-like
      if (PDFStyle.eq.'CTEQ') then
         photon = ctpara(x,ctphoton)
         return
      elseif (PDFStyle.eq.'CTEQHERA') then
         photon = ctherapara(x,ctphoton)
         return
      else
         photon=para(x,parphoton)
         
      endif
      
      end


C-------------------------------------------------
      Subroutine ChebToPoly()
C-------------------------------------------------
C
C Utility to convert chebyshev to standard polynomial expansion.
C
      implicit none
#include "pdfparam.inc"
#include "steering.inc"
      integer i
C---------------------------------
      if (nchebGlu.gt.0) then
         do i=1,nchebmax
            polyPars(i) = 0
         enddo
         call DCHPWS(nchebGlu,ChebPars,polyPars)
      endif
      if (nchebSea.gt.0) then
         do i=1,nchebmax
            polyParsSea(i) = 0
         enddo
         call DCHPWS(nchebSea,ChebParsSea,polyParsSea)
      endif
      end      


* -------------------------------------------------------
      double precision function H1U(x)
* -------------------------------------------------------
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision x,para



      H1U=para(x,paru)


      return
      end

* -------------------------------------------------------
      double precision function H1D(x)
* -------------------------------------------------------
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision x,para


      H1D=para(x,pard)


      return
      end

* -------------------------------------------------------
      double precision function Uval(x)
* -------------------------------------------------------
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision x,x23
      double precision PolyVal,ctpara,ctherapara,para,splogn
C---------------------------------------------------

C    22 Apr 11, SG, Add CTEQ-like
      if (PDFStyle.eq.'CTEQ') then
         UVal = ctpara(x,ctuval)
         return
      endif

      if (PDFStyle.eq.'CTEQHERA') then
         UVal = ctherapara(x,ctuval)
         return
      endif

C    22 Sep 11, VR, Add AS
      if ((PDFStyle.eq.'AS').or.(PDFStyle.eq.'BiLog')) then
         UVal = splogn(x,asuval)
         return
      endif

C
C 25 Jan 2011: add polynomial param 
C

      if (NPOLYVAL.eq.0) then
!> HERAPDF style goes in here
         Uval=para(x,paruval)
      else
C 
C PDFs are parameterised as a function of x23 = x^{2/3}
C
         x23 = x**(2.D0/3.D0)
         Uval = paruval(1) * PolyVal(x23,NPOLYVALINT,PolyUval)
      endif

      return
      end


* -------------------------------------------------------
      double precision function Dval(x)
* -------------------------------------------------------
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision x,x23
      double precision PolyVal,ctpara,ctherapara,para,splogn
C--------------------------------------------------------

C    22 Apr 11, SG, Add CTEQ-like
      if (PDFStyle.eq.'CTEQ') then
         DVal = ctpara(x,ctdval)
         return
      endif
      if (PDFStyle.eq.'CTEQHERA') then
         DVal = ctherapara(x,ctdval)
         return
      endif
C    22 Sep 11, VR, Add AS
      if ((PDFStyle.eq.'AS').or.(PDFStyle.eq.'BiLog')) then
         DVal = splogn(x,asdval)
         return
      endif


C
C 25 Jan 2011: add polynomial param 
C

      if (NPOLYVAL.eq.0) then
!> HERAPDF style goes in here
         Dval=para(x,pardval)
      else
C
C PDFs are parameterised as a function of x23 = x^{2/3}
C
         x23 = x**(2.D0/3.D0)
         Dval = pardval(1) * PolyVal(x23,NPOLYVALINT,PolyDval)
      endif

      return
      end

* -------------------------------------------------------
      double precision function PolyVal(x,NPOLY,Poly)
C-----------------------------------------------------
C 25 Jan 2011
C
C Evaluate polynomial sum fast
C
      implicit none
      integer NPOLY
      double precision x,Poly(NPOLY)
      integer i
      double precision sum
C------------------------------------------------
      sum = 0
      do i=NPoly,1,-1
         sum = x*(sum + Poly(i))
      enddo
      PolyVal = sum
      end


* -------------------------------------------------------
      double precision function sea(x)
* -------------------------------------------------------
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision x,Ubar,Dbar
C External function:
      double precision PolyParam,para
C--------------------------------------------------

      if (Index(PDF_DECOMPOSITION,'Dbar_Ubar').gt.0) then

         sea = Ubar(x) + Dbar(x)
      elseif (Index(PDF_DECOMPOSITION,'Sea').gt.0) then 
!         print*,'heeeere'
* warning for iparam = 3 or 4, the sea is 2 * sum (ubar +dbar + sbar + cbar)

         if (nchebSea.eq.0) then
            sea=para(x,parsea)
         else
            sea = PolyParam(x,nchebSea,polyParsSea,chebxminlog)
            if (ichebtypeSea.eq.0) then
C Do nothing
            else if (ichebtypeSea.eq.1) then
               sea = sea * (1 - x)   ! force PDFs=0 for x=1
            endif
         endif
      endif
 
      return
      end
* -------------------------------------------------------
      double precision Function PolyParam(x,ncheb,poly,xminlog)
* -------------------------------------------------------
C
C  SG: Use polynomial representation of cheb. 
C
      implicit none
      integer ncheb
      double precision x,poly(ncheb),xminlog
      double precision xx,sum
      integer i
C-------------------------------------------------
      xx = (2.D0*log(x)-xminlog)/(-xminlog)
      sum = poly(ncheb)
      do i=ncheb-1,1,-1
         sum = sum*xx + poly(i)
      enddo
      

      PolyParam = sum
      end


* -------------------------------------------------------
      double precision function dbmub(x)
* -------------------------------------------------------
* new jf , added to fit a la ZEUS, dbmub = dbar-ubar (not Dbar - Ubar)
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision x,para
      

      if (PDFStyle.eq.'CHEB'.or.PDFStyle.eq.'ZEUS Jet') then
         dbmub=para(x,parstr)
      endif
      return
      end

* -------------------------------------------------------
      double precision function qstrange (x)
* -------------------------------------------------------
* new jf , added to fit a la ZEUS, qstrange = 0.1 (i.e. fstrange *.5) * sea
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision x,sea,Dbar, para, ctpara,ctherapara
C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif
      
      if (PDFStyle.eq.'CTEQ') then
         qstrange = ctpara(x, ctstr)
         return
      endif

      if (PDFStyle.eq.'CTEQHERA') then
         qstrange = ctherapara(x, ctstr)
         return
      endif

      if (PDFStyle.eq.'CHEB'.or.PDFStyle.eq.'ZEUS Jet') then
         qstrange = 0.5 * fs * sea(x)
         return
         
      endif

      if (FreeStrange) then
         qstrange = para(x, parstr)
      else
         qstrange = fs * Dbar(x)
      endif


c      elseif (iparam.eq.222222.or.iparam.eq.222223) then
c         qstrange = fs * Dbar(x)/(1-fs)
cv      else
cv      endif

      return
      end

* -------------------------------------------------------
      double precision function cbar(x)
* -------------------------------------------------------
* new2 jf    
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
#include "thresholds.inc"
      double precision x,pdf,q2
      dimension pdf(-N_CHARGE_PDF:N_CHARGE_PDF+N_NEUTRAL_PDF)

      double precision sing,flav_number,QPDFXQ,vcplus,vicplus,cplus
      integer iflag,iq0,iqc,iqfromq

      if (PDFStyle.eq.'CHEB'.or.PDFStyle.eq.'ZEUS Jet') then
         if (x.eq.1) goto 999
         call hf_get_pdfs(x,q2,pdf)
* charm a la ZEUS 
         if (q0.lt.qc) then
            cbar = 0.
         else

            cbar=pdf(-4)
         endif

      endif
 999  continue
      return
      end

* -------------------------------------------------------
      double precision function Ubar(x)
* -------------------------------------------------------
* new2 jf    
* corrected for iparam=2 and iparam = 3 or 4
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
#include "thresholds.inc"
      double precision x,sea,dbmub,qstrange,cbar
      double precision sing,flav_number,QPDFXQ
      integer iflag,iq0,iqb,iqc,iqfromq,jtest
      double precision ctpara,ctherapara,para, splogn
C----------------------------------------------
* new2 jf SPECIAL TEST with dubar

C    22 Sep 11, VR, Add AS
      if ((PDFStyle.eq.'AS').or.(PDFStyle.eq.'BiLog')) then
         Ubar = splogn(x,asubar)
         return
      endif
      

C    22 Apr 11, SG, Add CTEQ-like
      if (PDFStyle.eq.'CTEQ') then
         Ubar=ctpara(x,ctubar)/(1-fcharm)
         return
      endif

      if (PDFStyle.eq.'CTEQHERA') then
         Ubar=ctherapara(x,ctubar)/(1-fcharm)
         return
      endif


      if (PDFStyle.eq.'CHEB'.or.PDFStyle.eq.'ZEUS Jet') then 
         Ubar = (0.5d0 * sea(x) - dbmub(x) - qstrange (x) + cbar(x))/2.d0
         return
      endif
     
cv      elseif (iparam.eq.222222.or.iparam.eq.222223
cv     $        .or.iparam.eq.2011) then
! Ubar=ubar/(1-fc) --> ubar=Ubar*(1-fc)         


      if (FreeStrange) then
         Ubar=para(x,parubar)
      else
         Ubar=para(x,parubar)/(1-fcharm)
      endif
      

 
      return
      end

* -------------------------------------------------------
      double precision function Dbar(x)
* -------------------------------------------------------
* new2 jf
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision x,sea,Ubar
      double precision ctpara,ctherapara,para,splogn
C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif
C------------------------------------------------------------

C    22 Sep 11, VR, Add AS

      if ((PDFStyle.eq.'AS').or.(PDFStyle.eq.'BiLog')) then
         Dbar = splogn(x,asdbar)
         return
      endif

C    22 Apr 11, SG, Add CTEQ-like
      if (PDFStyle.eq.'CTEQ') then
         Dbar=ctpara(x,ctdbar)+ctpara(x,ctstr)
         return
      endif

      if (PDFStyle.eq.'CTEQHERA') then
         Dbar=ctherapara(x,ctdbar)+ctherapara(x,ctstr)
         return
      endif

      if (PDFStyle.eq.'CHEB'.or.PDFStyle.eq.'ZEUS Jet') then
         Dbar = sea(x) * 0.5d0 - Ubar(x)
         return
      endif

!      elseif(iparam.eq.222222.or.iparam.eq.222223) then
!         Dbar=para(x,pardbar)/(1-fstrange)
      if (FreeStrange) then
         Dbar=para(x,pardbar)+para(x,parstr)
      else
         Dbar=para(x,pardbar)
      endif



      end


C---------------------------------------------------      
      double precision function fshermes(x)
C---------------------------------------------------
C
C X-dependent strange fraction, inspired by HERMES data
C Created 31 Oct 2009 by SG.
C
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision x
      double precision hermes_xcent,hermes_xrise
      parameter (hermes_xcent= 0.07)
      parameter (hermes_xrise=20.0)
      logical lfirstlocal
      data lfirstlocal/.true./
C---------------------------------------------------
      if (lfirstlocal) then
         lfirstlocal = .false.
         print *,'----------------------------------------------------'
         print *,' Called FSHERMES. Hermes-inspired strange density '
         print *,' Amp,Xmean,Xslope=',fstrange,hermes_xcent
     $        ,hermes_xrise
         print *,'----------------------------------------------------'
      endif

      if (x.lt.1.0D-8) then
         fshermes = fstrange*( 
     $        0.5D0*(1.D0+tanh(-(1.0D-8-hermes_xcent)*hermes_xrise)))
      else
         fshermes = fstrange*( 
     $        0.5D0*(1.D0+tanh(-(x-hermes_xcent)*hermes_xrise)))
      endif
C      print *,'DEBUG:',x,fshermes,fstrange
C---------------------------------------------------
      end




C---------------------------------------------------
      subroutine StorePoly(p,iflag)
C---------------------------------------------------
C
C Created 27 Jan 2011 by SG. Transfer parameters from MINUIT array p to
C internal arrays PolyUval and PolyDval
C
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision p(*)
      integer iflag,i,nn
      double precision ptmp(100)

C-----------------------------------      
      if (iflag.eq.1) then
         print *,
     $        'INFO: USING POLY parameterisation for valence quarks'
         print '(''Require PDFs at x=1 to vanish as (1-x)**'',I1)'
     $        ,IZPOPOLY
         if (IPOLYSQR.eq.1) then
            print *,
     $      'Square valence parameterisation, enforce positivity'
         else if (IPOLYSQR.eq.0) then
         else
            print *,'Invalid IPOLYSQR=',IPOLYSQR
            print *,'Can be 1 or 0, stop'
            call HF_stop
         endif
         
      endif
      
      if (IPOLYSQR.eq.0) then
         Call DecodePoly(p(21),PolyUVal,NPOLYVAL,IZPOPOLY)
         Call DecodePoly(p(31),PolyDVal,NPOLYVAL,IZPOPOLY)
         NPOLYVALINT = NPOLYVAL+IZPOPOLY
      else if  (IPOLYSQR.eq.1) then
         Call DecodePoly(p(21),Ptmp,NPOLYVAL,IZPOPOLY)
         Call SquarePoly(NPOLYVAL+IZPOPOLY,PTmp,NPOLYVALINT,PolyUVal)
         Call DecodePoly(p(31),Ptmp,NPOLYVAL,IZPOPOLY)
         Call SquarePoly(NPOLYVAL+IZPOPOLY,PTmp,NPOLYVALINT,PolyDVal)
      else
         print *,'Invalid IPOLYSQR=',IPOLYSQR
         print *,'Can be 1 or 0, stop'
         call HF_stop
      endif

      end      

C---------------------------------------------------
      subroutine DecodePoly(pars,poly,np,iz)
C---------------------------------------------------
C
C Created 29 Jan 2011. Decode input minuit parameters to internal PolyVal
C arrays for different order of (1-x) 
C     
      implicit none
      integer np,iz
      double precision pars(*),poly(*)
      integer i
C----------------------------------------------
      if (iz.eq.1) then
* \times (1-x)
         Poly(1) = pars(1)
         do i=2,np
            Poly(i) = pars(i)-pars(i-1)
         enddo
         Poly(np+1) = -pars(np)
      elseif (iz.eq.2) then
* \times (1-x)^2
         Poly(1) = pars(1)
         Poly(2) = pars(2)-2*pars(1)
         do i=3,np
            Poly(i) = pars(i)-2*pars(i-1)+pars(i-2)
         enddo
         Poly(np+1) = -2*pars(np)+pars(np-1)
         Poly(np+2) = pars(np)
      endif
      
      end
C---------------------------------------------------
      subroutine squarepoly(Npar,PolyIn,Npar2,PolyOut)
C-----------------------------------------
C
C     1 Feb 2011 by SG. Square a polynom. 
C
C-----------------------------------------
      implicit none
      integer Npar, Npar2
      double precision PolyIn(Npar), PolyOut(*)
      integer j,k,i
C--------------------------------------------------
      NPar2 = NPar*2-1

      do k=1,NPar2
         PolyOut(k) = 0.
         do j=1,k
            i = k-j+1
            if (i.le.NPar.and.j.le.NPar) then
               PolyOut(k) = PolyOut(k) + PolyIn(i)*PolyIn(j)
            endif
         enddo
      enddo

      end




C-----------------------------------------------------
      Subroutine PDFLength(DeltaChi2)
C---------------------------------------------------
C
C Created 27 June 2009 by SG
C Add extra constraint for the PDF "length" = int_wmin^wmax \sqrt{1+pdf'(W)**2} dw
C
      implicit none
      
#include "steering.inc"
#include "pdfparam.inc"
#include "pdflength.inc"
      double precision DeltaChi2
      double precision pdflen(5)
      double precision zero
      integer i
      integer ngrid
      parameter(ngrid=500)
      double precision grid(ngrid+1),d0,dst,val,wmin,wmax
C External functions
      double precision glulen,sealen
      double precision ssdint
      external glulen,sealen
      logical LFirstIn
      data LFirstIn /.true./
C----------------------------------------------------
C
C
      if (LFirstIn) then
        Wmin = WMNLEN
        Wmax = WMXLEN
        LFirstIn = .false.
        print '(''PDFLENGTH INITIALIZATION: Set Wmin,Wmax to '',2F8.2)',
     $        Wmin,Wmax
      endif

      DeltaChi2 = 0

      if (pdfLenWeight(1).gt.0) then
         pdflen(1) = ssdint(Wmin,glulen,Wmax)
     $     -(Wmax-Wmin)
         
         DeltaChi2 = DeltaChi2 + pdflen(1)*pdfLenWeight(1)
      else
         pdflen(1) = 0.
      endif

      if (pdfLenWeight(2).gt.0) then
         pdflen(2) = ssdint(Wmin,sealen,Wmax)
     $     -(Wmax-Wmin)
         
         DeltaChi2 = DeltaChi2 + pdflen(2)*pdfLenWeight(2)
      else
         pdflen(2) = 0.
      endif
 

      print *,'Gluon length=',pdflen(1)
      print *,'Sea length=',pdflen(2)


C---------------------------------------------------------
      end


C---------------------------------------------------
      double precision function powerLen(W,a,b,c,d,e,f,ap,bp,cp)
C---------------------------------------------------
C
C Utility to calculate pdf length element in W for power parameterization.
C
      implicit none
#include "steering.inc"
      double precision W,a,b,c,d,e,f,q2,x,der,derw,p,ap,bp,cp
C----------------------------------------------------
C Assume Q2=4
!      Q2 = 4.D0
      Q2 = starting_scale
      X = Q2/(Q2 + W*W)

      p   = (1.D0+d*x+e*x*x+f*x*x*x)
      der = a*x**b*(1.D0-x)**c*p*
     $     (b/x 
     $     - c/(1.0D0-x) 
     $     + (d+2.D0*e*x+3.D0*f*x*x)/p)
     $     -ap*x**bp*(1.D0-x)**cp*
     $     (bp/x
     $     -cp/(1.0D0-x))

C W derrivative:
      derw = - der* (2*W*Q2)/((W*W+Q2)*(W*W+Q2))

      PowerLen = sqrt(1.D0+derw*derw)      
C----------------------------------------------------
      end

C---------------------------------------------------      
      double precision function ChebLen(W,ncheb,poly,a,xminlog,iType)
C---------------------------------------------------
C
C Utility to calculate pdf length element in W for chebyshev parameterization.
C
      implicit none
#include "steering.inc"
      integer ncheb,iType
      double precision W,poly(ncheb),a,xminlog
      double precision Q2,X,XX,Sum,der,derw,sum2
      integer i
      logical LFirst
      data LFirst /.true./
C------------------------------------------------------
      if (LFirst) then
         print *,'First time in ChebLen. IType=',itype
         LFirst = .false.
      endif
C Assume Q2=4
      Q2 = starting_scale
      X = Q2/(Q2 + W*W)
C
C get derrivative:
C
      xx =  (2.D0*log(x)-xminlog)/(-xminlog)
      sum = poly(ncheb)*(ncheb-1)
         
      do i=ncheb-1,2,-1
         sum = sum*xx + poly(i)*(i-1.D0)
      enddo

C SG: Fix for (1-x) dumping:
      if (iType.eq.0) then
C do nothing
      else if (iType.eq.1) then        
C subract term corresponding to -x:
         sum2 = poly(ncheb)*ncheb
         do i=ncheb-1,1,-1
            sum2 = sum2*xx + poly(i)*i
         enddo
         sum = sum - sum2
      endif

      der = sum * a * 2.D0/(-xminlog) / x
C W derrivative:
      derw = - der* (2*W*Q2)/((W*W+Q2)*(W*W+Q2))

      ChebLen = sqrt(1.D0 + derw*derw)
C-------------------------------------------------------
      end
      

C---------------------------------------------------
      double precision function glulen(W)
C---------------------------------------------------
      implicit none
      double precision W
#include "pdfparam.inc"
#include "steering.inc"
C
      double precision PowerLen,ChebLen
C----------------------------------------------------


      if (nchebglu.eq.0) then
         glulen = powerlen(W,parglue(1)
     $        ,parglue(2),parglue(3),parglue(4),parglue(5),parglue(6),
     $        parglue(7),parglue(8),parglue(9))
      else
         glulen = cheblen(W,nchebGlu,polyPars,parglue(1),chebxminlog,
     $        ichebtypeGlu)
      endif
      end

C---------------------------------------------------
      double precision function Sealen(W)
C---------------------------------------------------
      implicit none
      double precision W
#include "pdfparam.inc"
#include "steering.inc"
C
      double precision PowerLen,ChebLen
C----------------------------------------------------


      if (nchebsea.eq.0) then
         Sealen = powerlen(W,parsea(1),parsea(2)
     $        ,parsea(3),parsea(4),parsea(5),parsea(6),
     $        parsea(7),parsea(8),parsea(9))
      else
         Sealen = cheblen(W,nchebSea,polyParsSea,1.D0,chebxminlog,
     $        ichebtypeSea)

      endif
      end
C---------------------------------------------------
      subroutine SaveRestorePdfs(imode)
C-------------------------------------------------------
C Leave only the contribution of the valence quarks
C for DGLAP+Dipole model fits.
C This subroutine ius called by subroutine fcn
C-------------------------------------------------------
      implicit none

#include "steering.inc"
#include "pdfparam.inc"
      integer imode

      double precision savglue(10)
      double precision savsea(10)
      double precision savdv(10)
      double precision savuv(10)
      integer i
C------------------------------

      if (imode.eq.0) then
C save all
         do i=1,10
            savglue(i) = parglue(i)
            savsea(i)  = parsea(i)
            savdv(i)   = pardval(i)
            savuv(i)   = paruval(i)
         enddo
      elseif (imode.eq.1) then
         do i=1,10
C Restore all:
            parglue(i) = savglue(i) 
            parsea(i)  = savsea(i)  
            pardval(i) = savdv(i) 
            paruval(i) = savuv(i)  

C For dglap fit, reset gluon and sea, keep valence:
            parglue(i) = 0.
            parsea(i)  = 0.
         enddo
      elseif (imode.eq.2) then
         do i=1,10
C restore
            parglue(i) = savglue(i) 
            parsea(i)  = savsea(i)  
            pardval(i) = savdv(i) 
            paruval(i) = savuv(i)  

C For dipole, reset valence:
            pardval(i)  = 0.0
            paruval(i)  = 0.0
         enddo         

      elseif (imode.eq.3) then
         do i=1,10
C restore
            parglue(i) = savglue(i) 
            parsea(i)  = savsea(i)  
            pardval(i) = savdv(i) 
            paruval(i) = savuv(i)  

         enddo         
      endif

      print *,'in save restore',imode,parglue(1),parglue(2),parglue(3)

      end

C !> Read PDF from a text file.

      double precision function pdf_from_text(x,id)
      implicit none
#include "steering.inc"
#include "ntot.inc"
#include "pdfparam.inc"
      double precision x
      integer id
      logical lfirst
      save lfirst
      data lfirst/.true./
      
      integer NXgrid
      double precision Q20
      namelist/XGrid/NXgrid, Q20

      double precision xx(NxgridMax)
      double precision xuv(NxgridMax)
      double precision xdv(NxgridMax)
      double precision xUbar(NxgridMax)
      double precision xDbar(NxgridMax)
      double precision xg(NxgridMax)
      save xx, xuv, xdv, xUbar, xDbar, xg
      integer i,ix

      integer ixfrmx
      logical XXATIX
C------------------------------------------------------
      if (lfirst) then
         lfirst = .false.
         open (51,file=LHAPDFSET,status='old',err=91)
         read (51,nml=XGrid,err=92,end=93) 
C Read the data
         read (51,*,err=94,end=95) (xx(i),i=1,NXgrid)  ! x 
C .... Add check ....
         do i=1,NXgrid
            ix = ixfrmx(xx(i))
            if (.not. xxatix(xx(i),ix)) then
               call hf_errlog(6,'F:Mis-match of the QCDNUM and text
     $ file grid in file '//Trim(LHAPDFSET))
            endif
         enddo
C ... read the tables ...
         read (51,*,err=94,end=95) (xuv(i),i=1,NXgrid)  
         read (51,*,err=94,end=95) (xdv(i),i=1,NXgrid)  
         read (51,*,err=94,end=95) (xubar(i),i=1,NXgrid)  
         read (51,*,err=94,end=95) (xdbar(i),i=1,NXgrid)  
         read (51,*,err=94,end=95) (xg(i),i=1,NXgrid)  
         
c         print '(4F12.6)',(xx(i),xg(i),xdbar(i),xubar(i),i=1,nxgrid)
c         stop

         call hf_errlog(301213,'I:Read PDF data from '
     $        //Trim(LHAPDFSET))

         close (51)
      endif
C------------------------------------------------------
      pdf_from_text = 0.0D0

C Get grid point:      
      ix = ixfrmx(x)
      if (id.eq.0) then
         pdf_from_text = xg(ix)
      elseif (id.eq.1) then
         pdf_from_text = xdv(ix)
      elseif (id.eq.2) then
         pdf_from_text = xuv(ix)
      elseif (id.eq.3) then
         pdf_from_text = 2*xdbar(ix) * fstrange ! /(1-fstrange)
      elseif (id.eq.4) then
         pdf_from_text = xubar(ix)
      elseif (id.eq.5) then
         pdf_from_text = xdbar(ix) ! * 1/(1-fstrange)
      elseif (id.eq.6) then
         pdf_from_text = 0.d0
      endif

C      print *,ix,xx(ix),x,fstrange
C      stop
      

      return
 91   call hf_errlog(1,'F:pdf_from_text: Can not open file '
     $     //Trim(LHAPDFSET))
 92   call hf_errlog(2,
     $     'F:pdf_from_text: Error reading namelist XGrid in '
     $     //Trim(LHAPDFSET))
 93   call hf_errlog(3,
     $     'F:pdf_from_text: Can not find namelist XGrid in '
     $     //Trim(LHAPDFSET))
 94   call hf_errlog(4,'F:pdf_from_text: Error reading PDF data in '
     $     //Trim(LHAPDFSET))
 95   call hf_errlog(5,'F:pdf_from_text: End of file for PDF data in '
     $     //Trim(LHAPDFSET))
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine DecodeFractal(pars)

#include "fractal.inc"
#include "extrapars.inc"
            
      
      integer GetParameterIndex
      double precision pars(*)
      integer idpx

      write(*,*) ' read parameters for fractal model >>'
      
!      idpx = GetParameterIndex('frac_1')
!      idpx = iExtraParamMinuit(idpx)
      f_D0 = pars(1)
      f_Q02 = pars(2)
      f_D2 = pars(3)
      f_D3 = pars(4)
      f_D1 = pars(5)
      f_R  = pars(6)
      end
      
      
