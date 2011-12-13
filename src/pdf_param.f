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
      include 'pdfparam.inc'
      include 'steering.inc'
      include 'alphas.inc'
      include 'thresholds.inc'
      include 'extrapars.inc'
      include 'polarity.inc'
      integer i

      double precision fs
      double precision fshermes

C-------------------------------------------------------
      logical LFirstTime
      data LFirstTime /.true./
      integer idxAlphaS, idxFs,   !> indices for alphas and fs
     $     idxFCharm
      integer GetParameterIndex  !> function to read parameter index
C-------------------------------------------------------
      integer idxShiftPolLHp   !> indices shiftpol
      integer idxShiftPolRHp   !> indices shiftpol
      integer idxShiftPolLHm   !> indices shiftpol
      integer idxShiftPolRHm   !> indices shiftpol
      integer idxShiftPolT   !> indices shiftpol
      integer idxShiftPolL   !> indices shiftpol
C-------------------------------------------------------
      integer idxAuEW, idxAdEW, idxVuEW, idxVdEW !> indices for EW param

C-------------------------------------------------------
      if (LFirstTime) then
         LFirstTime = .false.
         idxAlphaS = GetParameterIndex('alphas')
         if (idxAlphaS.eq.0) then
            print *,'Did not find alpha_S parameter'
            print *,'Add to ExtraParamters with the name alphas'
            print *,'Stop'
            stop
         else
            idxAlphaS = iExtraParamMinuit(idxAlphaS)
         endif
         idxFS = GetParameterIndex('fs')
         if (idxFS.eq.0) then
            print *,'Did not find fs parameter'
            print *,'Add to ExtraParamters with the name fs'
            print *,'Stop'
            stop
         else
            idxFS = iExtraParamMinuit(idxFS)
         endif
         
         idxFCharm = GetParameterIndex('fcharm')
         if (idxFCharm.gt.0) then
            idxFCharm = iExtraParamMinuit(idxFCharm)
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
         endif



      endif

C Polarisation shifts extra param


      shift_polRHp=p(idxShiftPolRHp)
      shift_polLHp=p(idxShiftPolLHp)
      shift_polLHm=p(idxShiftPolLHm)
      shift_polRHm=p(idxShiftPolRHm)
      shift_polL=p(idxShiftPolL)
      shift_polT=p(idxShiftPolT)




C EW extra param
      cau_ew=p(idxAuEW)
      cad_ew=p(idxAdEW)
      cvu_ew=p(idxVuEW)
      cvd_ew=p(idxVdEW)


C PDF param    

      Ag=p(1)
      Bg=p(2)
      Cg=p(3)
      Dg=p(4)
      Eg=p(5)
      Fg=p(6)
  
      Apg=p(7)
      Bpg=p(8)
      Cpg=p(9)

      Auv=p(11)
      Buv=p(12)
      Cuv=p(13)
      Duv=p(14)
      Euv=p(15)
      Fuv=p(16)

      Adv=p(21)
      Bdv=p(22)
      Cdv=p(23)
      Ddv=p(24)
      Edv=p(25)
      Fdv=p(26)

      Aubar=p(31)
      Bubar=p(32)
      Cubar=p(33)
      Dubar=p(34)


      AU=p(51)
      BU=p(52)
      CU=p(53)
      DU=p(54)

      Adbar=p(41)
      Bdbar=p(42)
      Cdbar=p(43)
      Ddbar=p(44)

      AD=p(61)
      BD=p(62)
      CD=p(63)
      DD=p(64)


      Asea=p(71)
      Bsea=p(72)
      Csea=p(73)
      Dsea=p(74)

      Adel=p(81)
      Bdel=p(82)
      Cdel=p(83)
      Ddel=p(84)

      Astr=p(81)
      Bstr=p(82)

      if (Bstr.eq.0) then
         Bstr = Bdbar
      endif

      Cstr=p(83)      
      if (Cstr.eq.0) then
         Cstr = Cdbar
      endif


C Get from extra pars:
      alphas=p(idxAlphaS)      
      fstrange=p(idxFS)



C Hermes strange prepare:
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(0.)
      endif
     

      if (q0.ge.qc.and.idxfcharm.gt.0) then
         fcharm=p(idxFCharm)
      else
         fcharm=0.
      endif


!!!!!!!!!!!!!!!!!!!!!!!

C 10 Aug 2011: Standard parametrisation:
      Call DecodePara(p)


C  25 Jan 2011: Poly params for valence:
      if (NPOLYVAL.gt.0) then
         Call StorePoly(p,iflag)
      endif

C  22 Apr 2011: CT parameterisation:
      if (IPARAM.eq.171717) then
         Call DecodeCtPara(p)
      endif


C different cases of parametrisation type
      if ((iparam.eq.1).or.(iparam.eq.21)) then     !  H1PDF2k like

         Bd = Bu
         Bubar = Bu
         Bdbar = Bu
         Adbar = Ad
         aU = aD * (1.-fs)/(1.-fcharm)
         aUbar = aU

C     Chebyshev param. for the gluon:
C
         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(30+i)
            enddo
            call ChebToPoly
         endif
         
      elseif (iparam.eq.2) then


         Bdv = Buv
         Aubar = Adbar * (1.-fs)/(1.-fcharm)
         Bubar = Bdbar
         

C
C  Chebyshev param. for the gluon:
C
         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(30+i)
            enddo
            call ChebToPoly
         endif

      elseif (iparam.eq.222222) then
 
         Aubar=Adbar

         Bubar=Bdbar
         
         Bdv=Buv


      elseif (iparam.eq.222223) then

 
         Aubar=Adbar
         Bubar=Bdbar         
         cpg=25.      

      elseif (iparam.eq.2011) then

         Aubar=Adbar
         Bubar=Bdbar        
C
C Fix strange to dbar if zero:
C
         if (Bstr.eq.0) then
            Bstr=Bdbar
         endif
         Astr=fs/(1-fs)*Adbar
         cpg=25.

      elseif ((iparam.eq.22).or.(iparam.eq.221).or.(iparam.eq.229)
     $        .or.(iparam.eq.222)) then


         if (iparam.ne.221.and.iparam.ne.229) then               
            Bdv = Buv
         endif
               
         Aubar = Adbar * (1.-fs)/(1.-fcharm)         
         Bubar = Bdbar
         if (iparam.eq.222.or.iparam.eq.229) then
            Cpg=25.
         endif


C
C  Chebyshev param. for the gluon:
C     
         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(30+i)
            enddo
            call ChebToPoly
         endif


      elseif (iparam.eq.225) then

C
C  Chebyshev param. for the gluon:
C
         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(30+i)
            enddo
            call ChebToPoly
         endif



      elseif (iparam.eq.3) then ! g,uval,dval,sea as in ZEUS-S 2002 fit
         

         Buv = 0.5
         Bdv = 0.5
         Bdel = 0.5
         Cdel = Csea +2.
         
      elseif ((iparam.eq.4).or.(iparam.eq.24)) then ! g,uval,dval,sea as in ZEUS-JET fit
         
         Bdv = Buv

*  dbar-ubar (not Ubar - Dbar), Adel fixed to output of ZEUS-S fit   
 
         Adel = 0.27          
         Bdel = 0.5
         Cdel = Csea +2.

C  Chebyshev param. for the gluon:
C
         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(30+i)
            enddo
            call ChebToPoly
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
      include 'thresholds.inc'

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
      include 'pdfparam.inc'
      include 'steering.inc'
      include 'for_debug.inc'
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
         fs = fshermes(0.)
      endif

C     simple copy first:
      do i=1,10
         parglue(i) = pars(i)
         paruval(i) = pars(10+i)
         pardval(i) = pars(20+i)
         parubar(i) = pars(30+i)
         pardbar(i) = pars(40+i)
         parsea(i) = pars(70+i)
         paru(i) = pars(50+i)
         pard(i) = pars(60+i)
         pardel(i) = pars(80+i)
         parother(i) = pars(90+i)
      enddo


      if ((iparam.eq.1).or.(iparam.eq.21)) then     !  H1PDF2k like

         pard(2)=paru(2)
         parubar(2)=paru(2)
         pardbar(2)=paru(2)
         pardbar(1)=pard(1)
         parU(1)=pard(1)*(1.D0-fs)/(1.D0-fcharm)
         parUbar(1)=parU(1)
         
      elseif (iparam.eq.2) then

         pardval(2)=paruval(2)
         parubar(2)=pardbar(2)
         parUbar(1)=pardbar(1)*(1.D0-fs)/(1.D0-fcharm)

      elseif ((iparam.eq.22).or.(iparam.eq.221).or.
     $        (iparam.eq.222)) then
         
         if (iparam.ne.221) then               
            pardval(2)=paruval(2)
         endif
         parubar(2)=pardbar(2)
         parUbar(1)=pardbar(1)*(1.D0-fs)/(1.D0-fcharm)

      elseif (iparam.eq.222222) then
         pardval(2)=paruval(2)
         parubar(2)=pardbar(2)
         parubar(1)=pardbar(1)

      elseif (iparam.eq.222223) then
         parubar(2)=pardbar(2)
         parubar(1)=pardbar(1)


      elseif (iparam.eq.2011) then
         parubar(2)=pardbar(2)
         parubar(1)=pardbar(1)

         if (pardel(2).eq.0) then
            pardel(2)=pardbar(2)
         endif

         if (pardel(3).eq.0) then
            pardel(3)=pardbar(3)
         endif

         if (fs.ne.-10000) then
            pardel(1)=fs/(1.-fs)*pardbar(1)
         endif

      elseif (iparam.eq.229) then
cv         parglue(9)=25.
         parUbar(2)=parDbar(2)
         parUbar(1)=parDbar(1)*(1.D0-fs)/(1.D0-fcharm)

      elseif (iparam.eq.3) then ! g,uval,dval,sea as in ZEUS-S 2002 fit
         
         paruval(2)=0.5
         pardval(2)=0.5
         pardel(2)=0.5
         pardel(3)=parsea(3)+2.
         
      elseif (iparam.eq.4) then ! g,uval,dval,sea as in ZEUS-JET fit
         pardval(2)=paruval(2)
        

*  dbar-ubar (not Ubar - Dbar), Adel fixed to output of ZEUS-S fit   
         pardel(1)=0.27
         pardel(2)=0.5
         pardel(3)=parsea(3)+2.

      elseif (iparam.eq.24) then ! g,uval,dval,sea as in ZEUS-JET fit
         pardval(2)=paruval(2)
         pardel(1)=0.27
         pardel(2)=0.5
         pardel(3)=parsea(3)+2.
         
      endif         

      if (debug) then
         print '(''1uv:'',11F10.4)',(paruval(i),i=1,10)
         print '(''1dv:'',11F10.4)',(pardval(i),i=1,10)
         print '(''1Ub:'',11F10.4)',(parubar(i),i=1,10)
         print '(''1Db:'',11F10.4)',(pardbar(i),i=1,10)
         print '(''1GL:'',11F10.4)',(parglue(i),i=1,10)
         print '(''1ST:'',11F10.4)',(pardel(i),i=1,10)
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
C   pars(95-100)  - alphas, fstrange, fcharm
C------------------------------------------------------
      implicit none 
      include 'pdfparam.inc'
      include 'steering.inc'
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
         fs = fshermes(0.)
      endif

C simple copy first:
      do i=1,6
         ctglue(i) = pars(i)
         ctuval(i) = pars(10+i)
         ctdval(i) = pars(20+i)
         ctubar(i) = pars(30+i)
         ctdbar(i) = pars(40+i)
         ctother(i)= pars(94+i)
      enddo


C Extra constrains:
      if (pars(31).eq.0) then
         ctubar(1) = ctdbar(1) * (1.D0-fs)/(1.D0-fcharm) ! normalization ubar = dbar 
         ctubar(2) = ctdbar(2)  ! Bubar = Bdbar
      endif

C Impose Buv = Bdv if parameter for Buv = 0.
      if (pars(12).eq.0) then
         ctuval(2) = ctdval(2)  ! Buv = Bdv
      endif
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
      double precision x,a(1:6)
      double precision UF
      UF = a(1)*exp(a(4)*x)*(1 - x)**a(3)*x**(a(2))*(1 + exp(a(5))*x 
     $     + exp(a(6))*x**2)
      
      ctpara = UF

      end




* -------------------------------------------------------
      double precision function gluon(x)
* -------------------------------------------------------
* x *g(x,Q2)

      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision x
      integer i
C External function:
      double precision PolyParam,ctpara,para
C-------------------------------------------------


C    22 Apr 11, SG, Add CTEQ-like
      if (iparam.eq.171717) then
         gluon = ctpara(x,ctglue)
         return
      endif


      if (nchebglu.eq.0) then

         gluon=para(x,parglue)
         
      else
C
C SG: Use polynomial representation of cheb. 
C
         gluon = ag * PolyParam(x,nchebGlu,polyPars,chebxminlog)
         if (ichebtypeGlu.eq.0) then
C Do nothing
         else if (ichebtypeGlu.eq.1) then
            gluon = gluon * (1 - x) ! force PDFs=0 for x=1
         endif
         
      endif
      end

      Subroutine ChebToPoly()
C
C Utility to convert chebyshev to standard polynomial expansion.
C
      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
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
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,para


      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.21) then
         H1U=para(x,paru)
cv       H1U = au * x**bu * (1.-x)**cu * (1. + du*x +  eu *x**2+ fu*x**3)
      endif

      return
      end

* -------------------------------------------------------
      double precision function H1D(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,para

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.21) then
         H1D=para(x,pard)
cv       H1D = ad * x**bd * (1.-x)**cd * (1. + dd*x)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Uval(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,x23
      double precision PolyVal,ctpara,para
C---------------------------------------------------

C    22 Apr 11, SG, Add CTEQ-like
      if (iparam.eq.171717) then
         UVal = ctpara(x,ctuval)
         return
      endif

      if (iparam.eq.2.or.iparam.eq.3.or.iparam.eq.4.or.iparam.eq.222222
     $     .or.iparam.eq.22.or.iparam.eq.24.or.iparam.eq.225.
     $     .or.iparam.eq.222223.or.iparam.eq.2011
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then

C
C 25 Jan 2011: add polynomial param 
C
         if (NPOLYVAL.eq.0) then
            Uval=para(x,paruval)

cv            Uval = aUv * x**bUv * (1.-x)**cUv
cv     +           * (1. + dUv*x + eUv *x**2+  fUV *x**3)
         else
C 
C PDFs are parameterised as a function of x23 = x^{2/3}
C
            x23 = x**(2.D0/3.D0)
            Uval = aUv * PolyVal(x23,NPOLYVALINT,PolyUval)
         endif
      endif

      return
      end


* -------------------------------------------------------
      double precision function Dval(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,x23
      double precision PolyVal,ctpara,para
C--------------------------------------------------------

C    22 Apr 11, SG, Add CTEQ-like
      if (iparam.eq.171717) then
         DVal = ctpara(x,ctdval)
         return
      endif


      if (iparam.eq.2.or.iparam.eq.3.or.iparam.eq.4
     $     .or.iparam.eq.22.or.iparam.eq.24.or.iparam.eq.225
     $     .or.iparam.eq.222222.or.iparam.eq.222223.or.iparam.eq.2011
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
C
C 25 Jan 2011: add polynomial param 
C
         if (NPOLYVAL.eq.0) then
            Dval=para(x,pardval)

cv            Dval = aDv * x**bDv * (1.-x)**cDv
cv     +           * (1. + dDv*x + fDv * x**3)
         else
C
C PDFs are parameterised as a function of x23 = x^{2/3}
C
            x23 = x**(2.D0/3.D0)
            Dval = aDv * PolyVal(x23,NPOLYVALINT,PolyDval)
         endif
      endif

      return
      end


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
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,Ubar,Dbar
C External function:
      double precision PolyParam,para
C--------------------------------------------------

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.2
     $     .or.iparam.eq.21.or.iparam.eq.22.or.iparam.eq.225
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
         sea = Ubar(x) + Dbar(x)
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then 
* warning for iparam = 3 or 4, the sea is 2 * sum (ubar +dbar + sbar + cbar)

         if (nchebSea.eq.0) then
            sea=para(x,parsea)

cv            sea = Asea * x**Bsea * (1.-x)**Csea
cv     +           * (1. + Dsea*x)
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

      double precision Function PolyParam(x,ncheb,poly,xminlog)
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
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,para
      
      if (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then 
         dbmub=para(x,pardel)
cv         dbmub = Adel * x**Bdel * (1.-x)**Cdel
      endif
 
      return
      end
* -------------------------------------------------------
      double precision function qstrange (x)
* -------------------------------------------------------
* new jf , added to fit a la ZEUS, qstrange = 0.1 (i.e. fstrange *.5) * sea
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,sea,Dbar

C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif
      
      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.2
     $     .or.iparam.eq.21.or.iparam.eq.22.or.iparam.eq.225
     $     .or.iparam.eq.171717
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then

         qstrange = fs * Dbar(x)

      elseif (iparam.eq.222222
     $     .or.iparam.eq.222223) then

         qstrange = fs * Dbar(x)/(1-fs)

      elseif (iparam.eq.2011) then
         qstrange = Astr*x**Bstr*(1-x)**Cstr
cv         print*,'voica is here', Astr, bstr, cstr


      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then 
         qstrange = 0.5 * fs * sea(x)
      endif
 

      return
      end
* -------------------------------------------------------
      double precision function cbar(x)
* -------------------------------------------------------
* new2 jf    
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,pdf,q2
      dimension pdf(-6:6)

      double precision sing,flav_number,QPDFXQ,vcplus,vicplus,cplus
      integer iflag,iq0,iqc,iqfromq

      if (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         if (x.eq.1) goto 999
         call fpdfxq(IPDFSET,x,q2,pdf,0)
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
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,sea,dbmub,qstrange,cbar
      double precision sing,flav_number,QPDFXQ
      integer iflag,iq0,iqb,iqc,iqfromq,jtest
      double precision ctpara,para
C----------------------------------------------
* new2 jf SPECIAL TEST with dubar


C    22 Apr 11, SG, Add CTEQ-like
      if (iparam.eq.171717) then
         Ubar = ctpara(x,ctubar)
         return

      elseif (iparam.eq.222222.or.iparam.eq.222223
     $        .or.iparam.eq.2011) then
! Ubar=ubar/(1-fc) --> ubar=Ubar*(1-fc)         
         Ubar=para(x,parubar)/(1-fcharm)
      endif


      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.2
     $     .or.iparam.eq.21.or.iparam.eq.22.or.iparam.eq.225
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
cv         Ubar = aubar * x**bubar * (1.-x)**cubar * (1. + dubar *x)
         Ubar=para(x,parubar)

      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
 
         Ubar = (0.5d0 * sea(x) - dbmub(x) - qstrange (x) + cbar(x))/2.d0

      endif
 
      return
      end

* -------------------------------------------------------
      double precision function Dbar(x)
* -------------------------------------------------------
* new2 jf
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,sea,Ubar
      double precision ctpara,para
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

C    22 Apr 11, SG, Add CTEQ-like
      if (iparam.eq.171717) then
         Dbar = ctpara(x,ctdbar)
         return
      elseif(iparam.eq.222222.or.iparam.eq.222223) then

         
         Dbar=para(x,pardbar)/(1-fstrange)
      elseif (iparam.eq.2011) then
         Dbar=para(x,pardbar)+para(x,pardel)

      endif


* SPECIAL TEST with ddbar      
      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.2
     $     .or.iparam.eq.21.or.iparam.eq.22.or.iparam.eq.225
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
         
cv         Dbar = adbar * x**bdbar * (1.-x)**cdbar * (1. + ddbar *x)
         Dbar=para(x,pardbar)

      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         Dbar = sea(x) * 0.5d0 - Ubar(x)
      endif

      return
      end

* -------------------------------------------------------
      double precision function singlet(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,H1U,H1D,Ubar,Dbar,sea,Uval,Dval

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.21) then
         singlet = H1U(x)+H1D(x)+Ubar(x)+Dbar(x)

      elseif (iparam.eq.2.or.iparam.eq.22.or.iparam.eq.225
     $        .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
         singlet = 2.d0 * sea(x) + Uval(x) + Dval(x)
* new jf
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         singlet = sea(x) + Uval(x) + Dval(x)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Uminus(x)
* -------------------------------------------------------
      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision x,H1U,Ubar,Uval

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.21) then
       Uminus = H1U(x) - Ubar(x)
* new jf
      elseif (iparam.eq.2.or.iparam.eq.3.or.iparam.eq.4
     $      .or.iparam.eq.22.or.iparam.eq.225.or.iparam.eq.24
     $      .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
       Uminus = Uval(x)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Dminus(x)
* -------------------------------------------------------
      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision x,H1D,Dbar,Dval

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.21) then
       Dminus = H1D(x) - Dbar(x)
* new jf
      elseif (iparam.eq.2.or.iparam.eq.3.or.iparam.eq.4
     $      .or.iparam.eq.22.or.iparam.eq.225.or.iparam.eq.24
     $      .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
       Dminus = Dval(x)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Deltaud(x)
* -------------------------------------------------------
      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision x,H1U,Ubar,H1D,Dbar
      double precision Uval,Dval,qstrange,cbar
C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif

      if (iparam.eq.1.or.iparam.eq.21) then
        Deltaud = H1U(x) - fcharm * Ubar(x)
     +          + Ubar(x) * (1. - fcharm)
     +          - (H1D(x) - fs * Dbar(x))
     +          - (Dbar(x) * (1. - fs))
      elseif (iparam.eq.2.or.iparam.eq.22.or.iparam.eq.225
     $       .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
* doubious factor two removed -- REVERT CHANGE (Sasha Glazov)
         Deltaud =  Uval(x) + 2.*(1.-fcharm)*Ubar(x)
     +        - (Dval(x) + 2.*(1.-fs)*Dbar(x))
CSG        Deltaud =  Uval(x) + (1.-fcharm)*Ubar(x)
CSG     +          - (Dval(x) + (1.-fs)*Dbar(x))
* new jf
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
*  doubious factor two removed  -- revert change
         Deltaud =  Uval(x) + 2.*(Ubar(x)-cbar(x))
     +        - (Dval(x) + 2.*(Dbar(x)-qstrange(x)))
CSG         Deltaud =  Uval(x) + (Ubar(x)-cbar(x))
CSG     +          - (Dval(x) + (Dbar(x)-qstrange(x)))
      elseif (iparam.eq.11) then
        Deltaud = H1U(x) + Ubar(x)
     +          - (H1D(x) - fs * Dbar(x))
     +          - (Dbar(x) * (1. - fs))
      endif

      return
      end

* -------------------------------------------------------
      double precision function Splus(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,Dbar,singlet,flav_number,qstrange
C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif


      Splus = 2.d0 * fs * Dbar(x) 
     +      - singlet(x)/flav_number(q0)
*new jf
      If (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         Splus = 2.d0 * qstrange(x) 
     +      - singlet(x)/flav_number(q0)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Uplus(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,Ubar,H1U,singlet,flav_number
      double precision Uval,cbar


* to be modified for iparam = 2 or 1 when q0.lt.qc

      if (iparam.eq.1.or.iparam.eq.21) then
        Uplus = Ubar(x)*(1.-fcharm)
     +          + H1U(x) - fcharm*Ubar(x) 
     +          - singlet(x)/flav_number(q0)

* new jf
      elseif (iparam.eq.2.or.iparam.eq.22.or.iparam.eq.225
     $       .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
        Uplus = Uval(x) + 2.*(1.-fcharm)*Ubar(x) 
     +        - singlet(x)/flav_number(q0)
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         Uplus = Uval(x) + 2.d0 * (Ubar(x)-cbar(x))
     +    - singlet(x)/flav_number(q0)
      elseif (iparam.eq.11) then
        Uplus = Ubar(x)
     +          + H1U(x) 
     +          - singlet(x)/flav_number(q0)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Dplus(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,Dbar,H1D,singlet,flav_number,Dval,qstrange
C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif
*mis print corrected 30/04/2008 (jf)
      if (iparam.eq.1.or.iparam.eq.21) then
        Dplus = Dbar(x)*(1.-fs)
     +          + H1D(x) - fs*Dbar(x) 
     +          - singlet(x)/flav_number(q0)
      elseif (iparam.eq.2.or.iparam.eq.22.or.iparam.eq.225
     $       .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
        Dplus = Dval(x) + 2.*(1.-fs)*Dbar(x) 
     +        - singlet(x)/flav_number(q0)
* new jf
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
        Dplus = Dval(x) + 2.*(Dbar(x)-qstrange(x)) 
     +        - singlet(x)/flav_number(q0)
      endif

      return
      end


* -------------------------------------------------------
      double precision function Cplus(x)
* -------------------------------------------------------
* given at starting scale if q0 > qc (allows to write
* that c+cbar is a fraction fcharm of Ubar for charm H1 style);
* else defined at qc (new2)

      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,singlet,flav_number,Ubar,factor
      integer iflag,jtest,iqfrmq,iq0,iqc,iqb

* Charm H1 style 

      if (iparam.eq.1.or.iparam.eq.2.or.iparam.eq.21.or.iparam.eq.22.or.iparam.eq.225
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
         if (q0.ge.qc) then
            Cplus = 2.d0 * fcharm * Ubar(x) 
     +       - singlet(x)/flav_number(q0)
         else

            Cplus = - singlet(x)/flav_number(qc)
           
         endif
      endif

* Charm ZEUS style

      if (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.11.or.iparam.eq.24) then

         Cplus = - singlet(x)/flav_number(qc)
         if(q0.ge.qc) then
            IQ0 = iqfrmq(q0)
            IQC = iqfrmq(qc)
            IQB = iqfrmq(qb)
            factor = -1. /flav_number(qc)
         endif
      endif

      return
      end

* -------------------------------------------------------
      double precision function Bplus(x)
* -------------------------------------------------------
* given at scale = b threshold
      implicit none

      include 'thresholds.inc'
      double precision x,singlet,flav_number
      integer iflag



      Bplus = - singlet(x)/flav_number(qb)
      return
      end

* -------------------------------------------------------
      double precision function Deltacs(x)
* -------------------------------------------------------
      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision x,Ubar,Dbar,qstrange,cbar
C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif

      Deltacs = 2.*fcharm*Ubar(x) 
     +        - 2.*fs*Dbar(x)
*new jf
      If (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         Deltacs = 2.*cbar(x) 
     +        - 2.*qstrange(x)
      endif

      return
      end

      
      double precision function fshermes(x)
C
C X-dependent strange fraction, inspired by HERMES data
C Created 31 Oct 2009 by SG.
C
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
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

      subroutine StorePoly(p,iflag)
C
C Created 27 Jan 2011 by SG. Transfer parameters from MINUIT array p to
C internal arrays PolyUval and PolyDval
C
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
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
            stop
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
         stop
      endif

      end      


      subroutine DecodePoly(pars,poly,np,iz)
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
C
C Created 27 June 2009 by SG
C Add extra constraint for the PDF "length" = int_wmin^wmax \sqrt{1+pdf'(W)**2} dw
C
      implicit none
      
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'pdflength.inc'
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
 

      print *,'Gluon length=',pdflen(1),ag,bg,cg,dg,fg
      print *,'Sea length=',pdflen(2)


C---------------------------------------------------------
      end

      double precision function powerLen(W,a,b,c,d,f)
C
C Utility to calculate pdf length element in W for power parameterization.
C
      implicit none
      double precision W,a,b,c,d,f,q2,x,der,derw,p
C----------------------------------------------------
C Assume Q2=4
      Q2 = 4.D0
      X = Q2/(Q2 + W*W)

      p   = (1.D0+d*x+f*x*x*x)
      der = a*x**b*(1.D0-x)**c*p*
     $     (b/x 
     $     - c/(1.0D0-x) 
     $     + (d+3.D0*f*x*x)/p)
C W derrivative:
      derw = - der* (2*W*Q2)/((W*W+Q2)*(W*W+Q2))

      PowerLen = sqrt(1.D0+derw*derw)      
C----------------------------------------------------
      end

      
      double precision function ChebLen(W,ncheb,poly,a,xminlog,iType)
C
C Utility to calculate pdf length element in W for chebyshev parameterization.
C
      implicit none
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
      Q2 = 4.D0
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
      


      double precision function glulen(W)
      implicit none
      double precision W
      include 'pdfparam.inc'
      include 'steering.inc'
C
      double precision PowerLen,ChebLen
C----------------------------------------------------


      if (nchebglu.eq.0) then
         glulen = powerlen(W,ag,bg,cg,dg,fg)
      else
         glulen = cheblen(W,nchebGlu,polyPars,ag,chebxminlog,
     $        ichebtypeGlu)
      endif
      end


      double precision function Sealen(W)
      implicit none
      double precision W
      include 'pdfparam.inc'
      include 'steering.inc'
C
      double precision PowerLen,ChebLen
C----------------------------------------------------


      if (nchebsea.eq.0) then
         Sealen = powerlen(W,aSea,bSea,cSea,dSea,0.0D0)
      else
         Sealen = cheblen(W,nchebSea,polyParsSea,1.D0,chebxminlog,
     $        ichebtypeSea)

      endif
      end
