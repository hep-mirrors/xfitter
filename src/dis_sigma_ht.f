! [--- KK, LM 2015-08-30, WS 2015-10-10   Twist4 corrections
!> @brief Twist-4 contribution to sigma reduced from the GBW dipole model without charm.
!>
!> @param[in] y Inelasticity,
!> @param[in] x Bjorken x,
!> @param[in] QQ Q^2 in GeV^2,
      double precision function dsigma_red_t4(y, x, QQ)
        implicit none
#include "steering.inc"
#include "couplings.inc"
                ! for alphaem, convfac, pi
        double precision y,x,QQ 
        double precision A0,x0,lambda,dmp  
        double precision rho,d4FT,d4FL
        
        ! --- cE = 1/15 + gamma_Euler
        double precision cE
        parameter(cE=1.d0/15+gammaEuler)
        
        ! --- \sum_q e_q^2
        double precision eq2sum
        parameter(eq2sum=2*e2d+e2u) ! sum e_q^2 for q=d,u,s.
        
        ! --- mbGeV2 = (\hbar c)^2 = 0.389... GeV^2 mbarn
        double precision mbGeV2  
        
        double precision tmp  
        
        ! Functions:
        double precision XParValueByName

        mbGeV2 = convfac*1.d-9

! <Twist params>
        x0   = XParValueByName('HT_x0')
        A0 = XParValueByName('HT_sig0')
        lambda = XParValueByName('HT_lambda')
        dmp = XParValueByName('HT_xbpow')
          
        rho = (x0/x)**lambda / QQ    

C   twist 4 corrections to FT, FL

C   generic GBW
        if(HiTwistSubType.eq.'lam-sig-x0') then
          ! --- Here A0 = sigma0 [mb]
          tmp = eq2sum/(20.*pi**3) * A0/mbGeV2 *rho**2*QQ
          d4FT = 3.*tmp
          d4FL = 4.*tmp *(dlog(rho) - cE)
        elseif(HiTwistSubType.eq.'lam-A0-x0') then
          ! --- Here A0 = eq2sum/(5*pi^3) x0^(2*lambda) * sigma0[mb]/mbGeV2
          tmp = A0/(QQ*x**(2.*lambda))
          d4FT = 3.d0/4. * tmp
          ! d4FL = tmp *(lambda*dlog(x0/x) -dlog(QQ) - cE)
          d4FL = tmp *(dlog(rho) - cE)
        elseif(HiTwistSubType.eq.'free') then
          ! --- Here A0 = eq2sum/(5*pi^3) x0^(2*lambda) * sigma0[mb]/mbGeV2
          tmp = -lambda*dlog(x) -dlog(QQ)
          ! d4FT = (cT0 + cTlog*tmp)/(QQ*x**(2.*lambda))
          ! d4FL = (cL0 + cLlog*tmp)/(QQ*x**(2.*lambda))
          A0 = QQ * x**(2.*lambda)
          d4FT = (XParValueByName('HT_cT0') + XParValueByName('HT_cTlog')*tmp)/A0
          d4FL = (XParValueByName('HT_cL0') + XParValueByName('HT_cLlog')*tmp)/A0
        else
          print *,'Illegal HiTwistSubType: "',HiTwistSubType,'"'
          call HF_stop
        endif

C  twist 4 correction to "sigma reduced"
C  --------------------------------------

        y = 1-y
        dsigma_red_t4 = d4FT + d4FL*2*y/(1+y**2)
        if(dmp.ne.0.0) dsigma_red_t4 = (1-x)**dmp * dsigma_red_t4

        return 
        end
! ---]
