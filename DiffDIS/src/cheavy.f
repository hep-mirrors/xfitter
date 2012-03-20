c=======================================================
      real*8 function cheavyzo(i,z,eps)
      implicit real*8(a-h,o-z)

c     this function returns the values of C_g(z,Q^2) 
c     and the deriv. wrt log Q^2. Here eps=m^2/Q^2.
c     If i=1  C_g for F2.  If i=2 deriv of C_g for F2 
c     If i=3  C_g for FL.  If i=4 (1-m2/Q2)*beta*Clq(massless)
c     If i=5  (1-m2/Q2)*beta*Clq(1)(massless)
      if(i.gt.5) stop
      z1=1.-z
      z2=z*z
      eps2=eps*eps
      beta2=1.-4.*eps*z/z1
cws
      if(beta2.lt.-1.d-8) go to 10
      if(beta2.lt.0.) then
        cheavyzo=0
        return
      endif

      beta=dsqrt(beta2)
      a=z2+z1*z1
      b=4.*z*(1.-3.*z)
      c=-8.*z2
      aa=8.*z*z1-1.
      bb=-4.*z*z1
      arg=(1.+beta)/(1.-beta)
      fac=dlog(arg)
      cf=4./3.
      ca=3.
      enf=4.
      ZETA2=1.64493406684823
      ZETA3=1.20205690315959
      go to (1,2,3,4,5) i
    1 cheavyzo=(a+b*eps+c*eps2)*fac+(aa+bb*eps)*beta
      return
    2 cheavyzo=(-b*eps-2.*c*eps2)*fac+(a+b*eps+c*eps2)/beta
     .      +(-bb*eps)*beta +(aa+bb*eps)*2.*z*eps/z1/beta
      return
    3 cheavyzo=-bb*beta+c*eps*fac
      return
    4 cheavyzo=(1.-eps)*beta*4.*cf*z
      return
    5 y=z
      Y1=1.-Y
      Y2=Y*Y
      Y3=Y*Y2
      Y4=Y3*Y
      Y5=Y4*Y
      DL=DLOG(Y)
      DL2=DL*DL
      DLM1=DLOG(Y1)
      YP1=Y+1.
      DLP1=DLOG(YP1)
      DLM2=DLOG(1.-Y2)
      DL3=DLOG(Y/Y1)
      DL4=DLOG(Y2/Y1)
      ALI2 =FASTLI2(Y)
      ALI2M=FASTLI2(-Y)
      ALI21=FASTLI2(Y1)
      ALI3 =FASTLI3(Y)
      ALI3M=FASTLI3(-Y)
      S12M =FASTS12(-Y)
      FACT1=4.*DL*(6.-3.*Y+47.*Y2-9.*Y3)/(15.*Y2)
     X-4.*ALI2M*(DL-2.*DLP1)-8.*ZETA3-2.*DL2*DLM2
     X+4.*DL*DLP1*DLP1-4.*DL*ALI2+0.4*(5.-3.*Y2)*DL2
     X-4.*(2.+10.*Y2+5.*Y3-3.*Y5)*(ALI2M+DL*DLP1)/(5.*Y3)
     X+4.*ZETA2*(DLM2-0.2*(5.-3.*Y2))+8.*S12M+4.*ALI3+4.*ALI3M
     X-23.*DLM1/3.-(144.+294.*Y-1729.*Y2+216.*Y3)/(90.*Y2)
      FACT2=ALI2+DL3*DL3-3.*ZETA2-(3.-22.*Y)*DL/(3.*Y)
     X+(6.-25.*Y)*DLM1/(6.*Y)-(78.-355.*Y)/(36.*Y)
      FACT3=DL4-(6.-25.*Y)/(6.*Y)
      FNS2LQ=4.*CF*(CA-2.*CF)*Y*FACT1
     X+8.*CF*CF*Y*FACT2-(8./3.)*CF*ENF*Y*FACT3
      cheavyzo=(1.-eps)*beta*FNS2LQ
      return

   10 print*,'T-R cheavy error: x > x0 for ',i,z,eps
   99 format(1x,'T-R cheavy: x > x0')
      stop
      end
