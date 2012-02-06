      subroutine sfunct_b(xbj,Q2,ft_l,fl_l,ft_c,fl_c)
c...
c...ft and fl sf computation
c...
      implicit double precision(a-h,o-z)

      logical DEBUG_PRINT
      data DEBUG_PRINT/.false./
      
      common/chsteer/ cm,icharm
      common/pass1/ sig0,xlam,x0,xm

c...light quark contribution
        call masscorr(xnew,xbj,xm,Q2)       !mass correction


C==== 26/07/2010 ADDED OUTPUT =====================

        if (DEBUG_PRINT) then
           print *,' xnew,xbj,xm,Q2 ',xnew,xbj,xm,Q2
           WRITE(*,*) 'cm     = ', cm
           WRITE(*,*) 'icharm = ', icharm
           WRITE(*,*) 'sig0   = ', sig0
           WRITE(*,*) 'xlam   = ', xlam
           WRITE(*,*) 'x0     = ', x0
           WRITE(*,*) 'xm     = ', xm
        endif
C==================================================

        chr2 =  2.d0/3.d0                    !(2/3)^2+(1/3)^2+(1/3)^2
        ft_l =  chr2 * F_t(xnew,Q2) 
        fl_l =  chr2 * F_l(xnew,Q2) 

C        print *,'hi',ft_l,fl_l

        ft_c =  0.d0
        fl_c =  0.d0

c...charm and bottom quark contribution
       if (icharm.eq.1) then
            tmp = xm
             xm = cm
         call masscorr(xnew,xbj,xm,Q2)     !mass correction
           chr2 = 4.d0/9.d0                 !(2/3)^2
           ft_c = chr2 * F_T(xnew,Q2)
           fl_c = chr2 * F_L(xnew,Q2)
             xm = tmp 
       endif

      return
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine masscorr(xnew,xold,xm,Q2)
c
c...mass correction
c
      implicit double precision(a-h,o-z)

      if (Q2.ne.0.0d0) then
        xnew = xold  * (1.d0+4*xm**2/Q2)
      else
        xnew = xold
      endif 

      return
      end
C-----------------------------------------------------------
       FUNCTION F_T(XX,Q2)
c
c... transverse sf
c
       IMPLICIT REAL*8(A-H,O-Z)
       PARAMETER(PI2=9.869588,ALFINV=137.036)

       F_T = Q2*ALFINV/(4*PI2) * CROSS_T(XX,Q2) 

       RETURN
       END
C-----------------------------------------------------------
       FUNCTION F_L(XX,Q2)
c
c... longitudinal sf
c
       IMPLICIT REAL*8(A-H,O-Z)
       PARAMETER(PI2=9.869588,ALFINV=137.036)

       F_L = Q2*ALFINV/(4*PI2) * CROSS_L(XX,Q2)

       RETURN
       END

C---------------------------------------------------------- 
       FUNCTION WAVDIP_T_r(r)
       IMPLICIT REAL*8(A-H,O-Z)
       EXTERNAL RINT_T

         rmi = 1.d-3
         rma = r
       CALL SIMPS(rmi,rma,1.D-3,1.D-3,1.D-4,RINT_T,DUM,WYNIK,DUM1,DUM2)

       WAVDIP_T_r = WYNIK

       RETURN
       END  

C---------------------------------------------------------- 
       FUNCTION WAVDIP_L_r(r)
       IMPLICIT REAL*8(A-H,O-Z)
       EXTERNAL RINT_L

         rmi = 1.d-3
         rma = r
       CALL SIMPS(rmi,rma,1.D-3,1.D-3,1.D-4,RINT_L,DUM,WYNIK,DUM1,DUM2)

       WAVDIP_L_r = WYNIK

       RETURN
       END  

c-----------------------------------------------------------
       FUNCTION CROSS_T(XX,Q2)
c
c...without charge square e^2 ; 
c...dimension in GeV^-2 ; multiply by 0.389379 to get in mbarns
c 
       IMPLICIT REAL*8(A-H,O-Z)
       PARAMETER(PI=3.14159D0,PI2=9.869588,ALFINV=137.036)

           FAC = 3.d0 / (2*PI2*ALFINV) * (2*PI)
       CROSS_T = FAC * WAVDIP_T(XX,Q2)
 
       RETURN
       END  

c-----------------------------------------------------------
       FUNCTION CROSS_L(XX,Q2)
c
c...without charge square e^2 ; 
c...dimension in GeV^-2 ; multiply by 0.389379 to get in mbarns
c
       IMPLICIT REAL*8(A-H,O-Z)
       PARAMETER(PI=3.14159D0,PI2=9.869588,ALFINV=137.036)

           FAC = 3.d0 / (2*PI2*ALFINV) * (2*PI)
       CROSS_L = FAC * WAVDIP_L(XX,Q2)

       RETURN
       END   

C---------------------------------------------------------- 
       FUNCTION WAVDIP_T(XX,Q2)
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /PASSZR/ XPASS,Q22,RR
       EXTERNAL RINT_T

       XPASS = XX
         Q22 = Q2

         rmi = 1.d-3
         rma = 1.d4
       CALL SIMPS(rmi,rma,1.D-3,1.D-3,1.D-4,RINT_T,DUM,WYNIK,DUM1,DUM2)

       WAVDIP_T = WYNIK

       RETURN
       END  
C----------------------------------------------------------- 
       FUNCTION WAVDIP_L(XX,Q2)
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /PASSZR/ XPASS,Q22,RR
       EXTERNAL RINT_L

       XPASS = XX
         Q22 = Q2

         rmi = 1.d-3  
         rma = 1.d4
       CALL SIMPS(rmi,rma,1.D-3,1.D-3,1.D-4,RINT_L,DUM,WYNIK,DUM1,DUM2)

       WAVDIP_L = WYNIK
       
       RETURN
       END 
c------------------------------------------------------------
       FUNCTION RINT_T(r)
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /PASSZR/ XPASS,Q22,RR
       EXTERNAL ZINT_T


       RR = r

         zmi = 0.0d0
         zma = 1.0d0
       WYNIK = DGQUAD(ZINT_T,zmi,zma,96)

       RINT_T = WYNIK * r * DIP_CS(RR,XPASS)

       RETURN
       END        
c------------------------------------------------------------
       FUNCTION RINT_L(r)
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /PASSZR/ XPASS,Q22,RR
       EXTERNAL ZINT_L

       RR = r

         zmi = 0.0d0
         zma = 1.0d0
       WYNIK = DGQUAD(ZINT_L,zmi,zma,96)

       RINT_L = WYNIK * r * DIP_CS(RR,XPASS)

       RETURN
       END 

c------------------------------------------------------------
       FUNCTION ZINT_T(z)
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /PASSZR/ XPASS,Q22,RR

       ZINT_T = WAVE_T(z,RR,Q22) 

       RETURN
       END        
c------------------------------------------------------------
       FUNCTION ZINT_L(z)
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /PASSZR/ XPASS,Q22,RR

       ZINT_L = WAVE_L(z,RR,Q22) 

       RETURN
       END 

C------------------------------------------------------------
       FUNCTION WAVE_T(Z,R,Q2)
       IMPLICIT REAL*8(A-H,O-Z)
       common/pass1/ sig0,xlam,x0,xm

C       print *,sig0,xlam,x0,xm
C       stop

        QBAR2 = Z*(1.-Z)*Q2 + XM*XM
         QBAR = DSQRT(QBAR2)
         BES0 = DBESK0(QBAR*R)
         BES1 = DBESK1(QBAR*R)

       WAVE_T = (Z**2+(1.d0-Z)**2) * QBAR2 * BES1**2 
     $        +  XM*XM * BES0**2

       RETURN
       END 
C------------------------------------------------------------ 
       FUNCTION WAVE_L(Z,R,Q2)
       IMPLICIT REAL*8(A-H,O-Z)
       common/pass1/ sig0,xlam,x0,xm

        QBAR2 = Z*(1.d0-Z)*Q2 + XM*XM
         QBAR = DSQRT(QBAR2)
         BES0 = DBESK0(QBAR*R)
       
       WAVE_L = 4.d0 * Q2 * (Z*(1.d0-Z))**2 * BES0**2

       RETURN
       END 

C------------------------------------------------------------
      FUNCTION DIP_CS(r,xx)
c
c...dimension in GeV^-2 : multiply by 0.389379 to get it in mbarns
c
      IMPLICIT REAL*8(A-H,O-Z)
c      common/pass1/ sig0,xlam,x0,xm
      logical lfirst
      data lfirst/.true./
      include 'dipole.inc'
      include 'steering.inc'

      if (lfirst) then
         lfirst = .false.
         print *,'DIP_CS: Dipole model =',idipole
      endif

      SIGMA0 = sig0 * 13.15d0   ![GeV^-2]

      if (DipoleModel.eq.1.or.DipoleModel.eq.3) then
c...saturation model with evolution

c...GBW saturation scale
         r02 = R_02_sat(r,xx)
         DIP_CS = SIGMA0 * (1.0D0 - dexp(-r02))
      elseif (DipoleModel.eq.2.or.DipoleModel.eq.4) then
         Y = log(1/xx)
         DIP_CS = sigma_iim(r,y,
C Fit parameters:
     $     xlam, x0, sig0  )                 ! [GeV^-2]
C         r02 = R_02_sat(r,xx)
C         DIP_CS2 = SIGMA0 * (1.0D0 - dexp(-r02))
C         print *,'1',dip_cs,dip_cs2,xx,r,y
      elseif (DipoleModel.eq.5) then 
         Y = log(1/x)
         DIP_CS = sigma_proton(r,y) * SIGMA0
      else
         print *,'Unknown dipole x-section=',idipole
         print *,'Abort in  DIP_CS'
         call HF_errlog(301,'F: DIP_CS - Unknown dipole x-section')
      endif

      RETURN
      END

C--------------------------------------------------------------
      Function dip_csiim(r,xx)
C
C Dipole cross section from Iancu, Itakura, Munier paper
C
      IMPLICIT REAL*8(A-H,O-Z)
      common/pass1/ sig0,xlam,x0,xm
C-------------------------------------------------------------
      Qs  = (x0/x)**(xlam/2)
      rQs = r*Qs
      
      if (rQs.le.2) then
         
      else
         aN = 1 - exp(aa)
      endif

      end

C--------------------------------------------------------------



c---------------------------------------------------------------
      FUNCTION R_02_sat(r,xx)
c
c...saturation model without evolution
c
      IMPLICIT REAL*8(A-H,O-Z) 
      common/pass1/ sig0,xlam,x0,xm

           tmp = (xx/x0)**xlam 
      R_02_sat = r**2 / (4.d0*tmp)

      RETURN
      END

