*
* all scale dependent hpij functions to be multiplied by powers of 
* dlog(((mq2+m2)/muF2))
*
! ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 1-loop CC g coefficient function h10_g,1

      real*8 function h1g10(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*

      H1g10=
     & +Llam*(1.D0-2.D0*x+2.D0*x**2)
      H1g10=H1g10+LQm*(1.D0-2.D0*x+2.D0*x**2)
      H1g10=H1g10-2.D0+8.D0*x-8.D0*x**2-2.D0*Hr1(0)+4.D0*Hr1(0)*x-4.D0*
     &    Hr1(0)*x**2-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*x**2

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 1-loop CC g coefficient function h10_g,2

      real*8 function h2g10(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*
      H2g10=
     & +Llam*(1.D0-2.D0*x+2.D0*x**2)
      H2g10=H2g10+LQm*(1.D0-2.D0*x+2.D0*x**2)
      H2g10=H2g10-2.D0+16.D0*x-16.D0*x**2-2.D0*Hr1(0)+4.D0*Hr1(0)*x-4.D0
     &    *Hr1(0)*x**2-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*x**2

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 1-loop CC g coefficient function h10_g,3

      real*8 function h3g10(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*
      H3g10=
     & +Llam*(1.D0-2.D0*x+2.D0*x**2)
      H3g10=H3g10+LQm*(-1.D0+2.D0*x-2.D0*x**2)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 1-loop CC g coefficient function h11_g

      real*8 function hg11(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*
      Hg11=
     & +1.D0-2.D0*x+2.D0*x**2

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC g coefficient function h20_g,1

      real*8 function h1g20(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*
      H1g20=
     & +Llam*(-4.D0/3.D0+16.D0/3.D0*x-16.D0/3.D0*x**2-4.D0/3.D0*Hr1(0)+
     &    8.D0/3.D0*Hr1(0)*x-8.D0/3.D0*Hr1(0)*x**2-4.D0/3.D0*Hr1(1)+8.D0
     &  /3.D0*Hr1(1)*x-8.D0/3.D0*Hr1(1)*x**2)
      H1g20=H1g20+Llam**2*(2.D0/3.D0-4.D0/3.D0*x+4.D0/3.D0*x**2)
      H1g20=H1g20+LQm*(4.D0/3.D0-16.D0/3.D0*x+16.D0/3.D0*x**2+4.D0/3.D0
     &    *Hr1(0)-8.D0/3.D0*Hr1(0)*x+8.D0/3.D0*Hr1(0)*x**2+4.D0/3.D0*
     &    Hr1(1)-8.D0/3.D0*Hr1(1)*x+8.D0/3.D0*Hr1(1)*x**2)
      H1g20=H1g20+LQm**2*(-2.D0/3.D0+4.D0/3.D0*x-4.D0/3.D0*x**2)
      H1g20=H1g20+ca*(188.D0/9.D0-215.D0/9.D0*x-361.D0/27.D0*x**2+280.D0
     &  /27.D0*dx+112.D0/3.D0*Hr1(0)+157.D0/3.D0*Hr1(0)*x-206.D0/3.D0*
     &    Hr1(0)*x**2+11.D0/3.D0*Hr1(1)+10.D0/3.D0*Hr1(1)*x-142.D0/9.D0
     &    *Hr1(1)*x**2-56.D0/9.D0*Hr1(1)*dx-24.D0*Hr2(-1,0)-28.D0*Hr2(
     &    -1,0)*x-4.D0/3.D0*Hr2(-1,0)*x**2-16.D0/3.D0*Hr2(-1,0)*dx-Hr2(
     &    0,0)+84.D0*Hr2(0,0)*x-365.D0/3.D0*Hr2(0,0)*x**2-8.D0*Hr2(0,1)
     &    +48.D0*Hr2(0,1)*x-115.D0*Hr2(0,1)*x**2-Hr2(1,0)+64.D0*Hr2(1,0
     &    )*x-79.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,0)*dx-5.D0*Hr2(1,1)+36.D0*
     &    Hr2(1,1)*x-133.D0/3.D0*Hr2(1,1)*x**2+16.D0/3.D0*Hr2(1,1)*dx+4.
     &  D0*Hr3(-1,-1,0)+8.D0*Hr3(-1,-1,0)*x-8.D0*Hr3(-1,-1,0)*x**2+10.D0
     &    *Hr3(-1,0,0)+20.D0*Hr3(-1,0,0)*x+28.D0*Hr3(-1,0,0)*x**2+8.D0*
     &    Hr3(-1,0,1)+16.D0*Hr3(-1,0,1)*x+16.D0*Hr3(-1,0,1)*x**2+16.D0*
     &    Hr3(0,-1,0)*x**2+18.D0*Hr3(0,0,0)+52.D0*Hr3(0,0,0)*x+8.D0*
     &    Hr3(0,0,1)+64.D0*Hr3(0,0,1)*x-16.D0*Hr3(0,0,1)*x**2+4.D0*Hr3(
     &    0,1,0)+64.D0*Hr3(0,1,0)*x-16.D0*Hr3(0,1,0)*x**2+48.D0*Hr3(0,1
     &    ,1)*x)
      H1g20=H1g20+ca*(-16.D0*Hr3(0,1,1)*x**2-12.D0*Hr3(1,0,0)+24.D0*
     &    Hr3(1,0,0)*x-16.D0*Hr3(1,0,0)*x**2-6.D0*Hr3(1,0,1)+12.D0*Hr3(
     &    1,0,1)*x-12.D0*Hr3(1,0,1)*x**2-14.D0*Hr3(1,1,0)+28.D0*Hr3(1,1
     &    ,0)*x-28.D0*Hr3(1,1,0)*x**2-2.D0*Hr3(1,1,1)+4.D0*Hr3(1,1,1)*x
     &    -4.D0*Hr3(1,1,1)*x**2)
      H1g20=H1g20+ca*Llam*(-50.D0/3.D0-94.D0/3.D0*x+434.D0/9.D0*x**2+16.
     &  D0/9.D0*dx-2.D0*Hr1(0)-48.D0*Hr1(0)*x+212.D0/3.D0*Hr1(0)*x**2+4.
     &  D0*Hr1(1)-56.D0*Hr1(1)*x+196.D0/3.D0*Hr1(1)*x**2-16.D0/3.D0*
     &    Hr1(1)*dx-4.D0*Hr2(-1,0)-8.D0*Hr2(-1,0)*x-8.D0*Hr2(-1,0)*x**2
     &    -12.D0*Hr2(0,0)-40.D0*Hr2(0,0)*x-48.D0*Hr2(0,1)*x+16.D0*Hr2(0
     &    ,1)*x**2+8.D0*Hr2(1,0)-16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+
     &    12.D0*Hr2(1,1)-24.D0*Hr2(1,1)*x+24.D0*Hr2(1,1)*x**2)
      H1g20=H1g20+ca*Llam**2*(1.D0+8.D0*x-31.D0/3.D0*x**2+4.D0/3.D0*dx+
     &    2.D0*Hr1(0)+8.D0*Hr1(0)*x-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(
     &    1)*x**2)
      H1g20=H1g20+ca*LQm*(-4.D0+50.D0*x-436.D0/9.D0*x**2+40.D0/9.D0*dx+
     &    2.D0*Hr1(0)+16.D0*Hr1(0)*x+88.D0/3.D0*Hr1(0)*x**2+8.D0*Hr1(1)
     &    *x-8.D0*Hr1(1)*x**2-4.D0*Hr2(-1,0)-8.D0*Hr2(-1,0)*x-8.D0*Hr2(
     &    -1,0)*x**2-4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x-4.D0*Hr2(1,1)+8.D0*
     &    Hr2(1,1)*x-8.D0*Hr2(1,1)*x**2)
      H1g20=H1g20+ca*LQm*Llam*(2.D0+16.D0*x-62.D0/3.D0*x**2+8.D0/3.D0*
     &    dx+4.D0*Hr1(0)+16.D0*Hr1(0)*x-4.D0*Hr1(1)+8.D0*Hr1(1)*x-8.D0*
     &    Hr1(1)*x**2)
      H1g20=H1g20+ca*LQm**2*(-1.D0-8.D0*x+31.D0/3.D0*x**2-4.D0/3.D0*dx-
     &    2.D0*Hr1(0)-8.D0*Hr1(0)*x+2.D0*Hr1(1)-4.D0*Hr1(1)*x+4.D0*Hr1(
     &    1)*x**2)
      H1g20=H1g20+cf*(-281.D0/10.D0+881.D0/10.D0*x-272.D0/5.D0*x**2-8.D0
     &  /5.D0*dx-64.D0/5.D0*Hr1(0)+597.D0/10.D0*Hr1(0)*x-372.D0/5.D0*
     &    Hr1(0)*x**2+8.D0/5.D0*Hr1(0)*dx-6.D0*Hr1(1)+77.D0*Hr1(1)*x-68.
     &  D0*Hr1(1)*x**2+48.D0*Hr2(-1,0)+32.D0*Hr2(-1,0)*x+32.D0/5.D0*
     &    Hr2(-1,0)*x**3-8.D0/5.D0*Hr2(-1,0)*dx**2-7.D0/2.D0*Hr2(0,0)+
     &    30.D0*Hr2(0,0)*x-62.D0*Hr2(0,0)*x**2-32.D0/5.D0*Hr2(0,0)*x**3
     &    -14.D0*Hr2(0,1)+84.D0*Hr2(0,1)*x-78.D0*Hr2(0,1)*x**2-25.D0*
     &    Hr2(1,0)+68.D0*Hr2(1,0)*x-62.D0*Hr2(1,0)*x**2-24.D0*Hr2(1,1)+
     &    84.D0*Hr2(1,1)*x-78.D0*Hr2(1,1)*x**2-32.D0*Hr3(-1,-1,0)-64.D0
     &    *Hr3(-1,-1,0)*x-32.D0*Hr3(-1,-1,0)*x**2+16.D0*Hr3(-1,0,0)+32.D
     &  0*Hr3(-1,0,0)*x+16.D0*Hr3(-1,0,0)*x**2+32.D0*Hr3(0,-1,0)+32.D0*
     &    Hr3(0,-1,0)*x**2-9.D0*Hr3(0,0,0)+18.D0*Hr3(0,0,0)*x-36.D0*
     &    Hr3(0,0,0)*x**2-16.D0*Hr3(0,0,1)+32.D0*Hr3(0,0,1)*x-48.D0*
     &    Hr3(0,0,1)*x**2-14.D0*Hr3(0,1,0)+28.D0*Hr3(0,1,0)*x-32.D0*
     &    Hr3(0,1,0)*x**2-18.D0*Hr3(0,1,1)+36.D0*Hr3(0,1,1)*x-44.D0*
     &    Hr3(0,1,1)*x**2)
      H1g20=H1g20+cf*(-2.D0*Hr3(1,0,0)+4.D0*Hr3(1,0,0)*x-20.D0*Hr3(1,0,
     &    0)*x**2-24.D0*Hr3(1,0,1)+48.D0*Hr3(1,0,1)*x-48.D0*Hr3(1,0,1)*
     &    x**2-16.D0*Hr3(1,1,0)+32.D0*Hr3(1,1,0)*x-32.D0*Hr3(1,1,0)*
     &    x**2-22.D0*Hr3(1,1,1)+44.D0*Hr3(1,1,1)*x-44.D0*Hr3(1,1,1)*
     &    x**2)
      H1g20=H1g20+cf*Llam*(5.D0-21.D0*x+12.D0*x**2+2.D0*Hr1(0)-20.D0*
     &    Hr1(0)*x+20.D0*Hr1(0)*x**2+7.D0*Hr1(1)-24.D0*Hr1(1)*x+20.D0*
     &    Hr1(1)*x**2+4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x+16.D0*Hr2(0,0)*x**2
     &    +6.D0*Hr2(0,1)-12.D0*Hr2(0,1)*x+16.D0*Hr2(0,1)*x**2+8.D0*Hr2(
     &    1,0)-16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,1)-16.D0
     &    *Hr2(1,1)*x+16.D0*Hr2(1,1)*x**2)
      H1g20=H1g20+cf*Llam**2*(-1.D0/2.D0+2.D0*x-Hr1(0)+2.D0*Hr1(0)*x-4.D
     &  0*Hr1(0)*x**2-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*x**2)
      H1g20=H1g20+cf*LQm*(5.D0-21.D0*x+12.D0*x**2+2.D0*Hr1(0)-20.D0*
     &    Hr1(0)*x+20.D0*Hr1(0)*x**2+7.D0*Hr1(1)-24.D0*Hr1(1)*x+20.D0*
     &    Hr1(1)*x**2+4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x+16.D0*Hr2(0,0)*x**2
     &    +6.D0*Hr2(0,1)-12.D0*Hr2(0,1)*x+16.D0*Hr2(0,1)*x**2+8.D0*Hr2(
     &    1,0)-16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,1)-16.D0
     &    *Hr2(1,1)*x+16.D0*Hr2(1,1)*x**2)
      H1g20=H1g20+cf*LQm**2*(-1.D0/2.D0+2.D0*x-Hr1(0)+2.D0*Hr1(0)*x-4.D0
     &    *Hr1(0)*x**2-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*x**2)
      H1g20=H1g20+z3*ca*(12.D0-20.D0*x+24.D0*x**2)
      H1g20=H1g20+z3*cf*(30.D0+4.D0*x+76.D0*x**2)
      H1g20=H1g20+z2*ca*(8.D0-76.D0*x+115.D0*x**2-16.D0/3.D0*dx-6.D0*
     &    Hr1(-1)-12.D0*Hr1(-1)*x-20.D0*Hr1(-1)*x**2-8.D0*Hr1(0)-64.D0*
     &    Hr1(0)*x+16.D0*Hr1(0)*x**2+8.D0*Hr1(1)-16.D0*Hr1(1)*x+8.D0*
     &    Hr1(1)*x**2)
      H1g20=H1g20+z2*ca*Llam*(40.D0*x-16.D0*x**2)
      H1g20=H1g20+z2*ca*LQm*(-8.D0*x)
      H1g20=H1g20+z2*cf*(14.D0-52.D0*x+78.D0*x**2+32.D0/5.D0*x**3-16.D0
     &    *Hr1(-1)-32.D0*Hr1(-1)*x-16.D0*Hr1(-1)*x**2+16.D0*Hr1(0)-32.D0
     &    *Hr1(0)*x+48.D0*Hr1(0)*x**2+8.D0*Hr1(1)-16.D0*Hr1(1)*x+32.D0*
     &    Hr1(1)*x**2)
      H1g20=H1g20+z2*cf*Llam*(-6.D0+12.D0*x-16.D0*x**2)
      H1g20=H1g20+z2*cf*LQm*(-6.D0+12.D0*x-16.D0*x**2)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC g coefficient function h20_g,2

      real*8 function h2g20(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*
      H2g20=
     & +Llam*(-4.D0/3.D0+32.D0/3.D0*x-32.D0/3.D0*x**2-4.D0/3.D0*Hr1(0)+
     &    8.D0/3.D0*Hr1(0)*x-8.D0/3.D0*Hr1(0)*x**2-4.D0/3.D0*Hr1(1)+8.D0
     &  /3.D0*Hr1(1)*x-8.D0/3.D0*Hr1(1)*x**2)
      H2g20=H2g20+Llam**2*(2.D0/3.D0-4.D0/3.D0*x+4.D0/3.D0*x**2)
      H2g20=H2g20+LQm*(4.D0/3.D0-32.D0/3.D0*x+32.D0/3.D0*x**2+4.D0/3.D0
     &    *Hr1(0)-8.D0/3.D0*Hr1(0)*x+8.D0/3.D0*Hr1(0)*x**2+4.D0/3.D0*
     &    Hr1(1)-8.D0/3.D0*Hr1(1)*x+8.D0/3.D0*Hr1(1)*x**2)
      H2g20=H2g20+LQm**2*(-2.D0/3.D0+4.D0/3.D0*x-4.D0/3.D0*x**2)
      H2g20=H2g20+ca*(236.D0/9.D0+601.D0/9.D0*x-2905.D0/27.D0*x**2+232.D
     &  0/27.D0*dx+160.D0/3.D0*Hr1(0)+541.D0/3.D0*Hr1(0)*x-830.D0/3.D0*
     &    Hr1(0)*x**2+59.D0/3.D0*Hr1(1)+442.D0/3.D0*Hr1(1)*x-1534.D0/9.D
     &  0*Hr1(1)*x**2-104.D0/9.D0*Hr1(1)*dx-24.D0*Hr2(-1,0)+4.D0*Hr2(-1
     &    ,0)*x+92.D0/3.D0*Hr2(-1,0)*x**2-16.D0/3.D0*Hr2(-1,0)*dx-Hr2(0
     &    ,0)+180.D0*Hr2(0,0)*x-365.D0/3.D0*Hr2(0,0)*x**2-8.D0*Hr2(0,1)
     &    +144.D0*Hr2(0,1)*x-147.D0*Hr2(0,1)*x**2-Hr2(1,0)+96.D0*Hr2(1,
     &    0)*x-111.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,0)*dx-5.D0*Hr2(1,1)+68.D0
     &    *Hr2(1,1)*x-229.D0/3.D0*Hr2(1,1)*x**2+16.D0/3.D0*Hr2(1,1)*dx+
     &    4.D0*Hr3(-1,-1,0)+8.D0*Hr3(-1,-1,0)*x-8.D0*Hr3(-1,-1,0)*x**2+
     &    10.D0*Hr3(-1,0,0)+20.D0*Hr3(-1,0,0)*x+28.D0*Hr3(-1,0,0)*x**2+
     &    8.D0*Hr3(-1,0,1)+16.D0*Hr3(-1,0,1)*x+16.D0*Hr3(-1,0,1)*x**2+
     &    16.D0*Hr3(0,-1,0)*x**2+18.D0*Hr3(0,0,0)+52.D0*Hr3(0,0,0)*x+8.D
     &  0*Hr3(0,0,1)+64.D0*Hr3(0,0,1)*x-16.D0*Hr3(0,0,1)*x**2+4.D0*Hr3(
     &    0,1,0)+64.D0*Hr3(0,1,0)*x-16.D0*Hr3(0,1,0)*x**2+48.D0*Hr3(0,1
     &    ,1)*x)
      H2g20=H2g20+ca*(-16.D0*Hr3(0,1,1)*x**2-12.D0*Hr3(1,0,0)+24.D0*
     &    Hr3(1,0,0)*x-16.D0*Hr3(1,0,0)*x**2-6.D0*Hr3(1,0,1)+12.D0*Hr3(
     &    1,0,1)*x-12.D0*Hr3(1,0,1)*x**2-14.D0*Hr3(1,1,0)+28.D0*Hr3(1,1
     &    ,0)*x-28.D0*Hr3(1,1,0)*x**2-2.D0*Hr3(1,1,1)+4.D0*Hr3(1,1,1)*x
     &    -4.D0*Hr3(1,1,1)*x**2)
      H2g20=H2g20+ca*Llam*(-98.D0/3.D0-334.D0/3.D0*x+1250.D0/9.D0*x**2+
     &    64.D0/9.D0*dx-2.D0*Hr1(0)-112.D0*Hr1(0)*x+212.D0/3.D0*Hr1(0)*
     &    x**2+4.D0*Hr1(1)-88.D0*Hr1(1)*x+292.D0/3.D0*Hr1(1)*x**2-16.D0/
     &    3.D0*Hr1(1)*dx-4.D0*Hr2(-1,0)-8.D0*Hr2(-1,0)*x-8.D0*Hr2(-1,0)
     &    *x**2-12.D0*Hr2(0,0)-40.D0*Hr2(0,0)*x-48.D0*Hr2(0,1)*x+16.D0*
     &    Hr2(0,1)*x**2+8.D0*Hr2(1,0)-16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*
     &    x**2+12.D0*Hr2(1,1)-24.D0*Hr2(1,1)*x+24.D0*Hr2(1,1)*x**2)
      H2g20=H2g20+ca*Llam**2*(1.D0+8.D0*x-31.D0/3.D0*x**2+4.D0/3.D0*dx+
     &    2.D0*Hr1(0)+8.D0*Hr1(0)*x-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(
     &    1)*x**2)
      H2g20=H2g20+ca*LQm*(-4.D0+50.D0*x-436.D0/9.D0*x**2+40.D0/9.D0*dx+
     &    2.D0*Hr1(0)+16.D0*Hr1(0)*x+88.D0/3.D0*Hr1(0)*x**2+8.D0*Hr1(1)
     &    *x-8.D0*Hr1(1)*x**2-4.D0*Hr2(-1,0)-8.D0*Hr2(-1,0)*x-8.D0*Hr2(
     &    -1,0)*x**2-4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x-4.D0*Hr2(1,1)+8.D0*
     &    Hr2(1,1)*x-8.D0*Hr2(1,1)*x**2)
      H2g20=H2g20+ca*LQm*Llam*(2.D0+16.D0*x-62.D0/3.D0*x**2+8.D0/3.D0*
     &    dx+4.D0*Hr1(0)+16.D0*Hr1(0)*x-4.D0*Hr1(1)+8.D0*Hr1(1)*x-8.D0*
     &    Hr1(1)*x**2)
      H2g20=H2g20+ca*LQm**2*(-1.D0-8.D0*x+31.D0/3.D0*x**2-4.D0/3.D0*dx-
     &    2.D0*Hr1(0)-8.D0*Hr1(0)*x+2.D0*Hr1(1)-4.D0*Hr1(1)*x+4.D0*Hr1(
     &    1)*x**2)
      H2g20=H2g20+cf*(-1099.D0/30.D0+273.D0/10.D0*x+64.D0/5.D0*x**2+8.D0
     &  /15.D0*dx-296.D0/15.D0*Hr1(0)+181.D0/10.D0*Hr1(0)*x-276.D0/5.D0
     &    *Hr1(0)*x**2-8.D0/15.D0*Hr1(0)*dx-14.D0*Hr1(1)+53.D0*Hr1(1)*x
     &    -36.D0*Hr1(1)*x**2+48.D0*Hr2(-1,0)+64.D0/3.D0*Hr2(-1,0)*x+96.D
     &  0/5.D0*Hr2(-1,0)*x**3+8.D0/15.D0*Hr2(-1,0)*dx**2-7.D0/2.D0*Hr2(
     &    0,0)+26.D0/3.D0*Hr2(0,0)*x-62.D0*Hr2(0,0)*x**2-96.D0/5.D0*
     &    Hr2(0,0)*x**3-14.D0*Hr2(0,1)+68.D0*Hr2(0,1)*x-78.D0*Hr2(0,1)*
     &    x**2-25.D0*Hr2(1,0)+68.D0*Hr2(1,0)*x-62.D0*Hr2(1,0)*x**2-24.D0
     &    *Hr2(1,1)+84.D0*Hr2(1,1)*x-78.D0*Hr2(1,1)*x**2-32.D0*Hr3(-1,
     &    -1,0)-64.D0*Hr3(-1,-1,0)*x-32.D0*Hr3(-1,-1,0)*x**2+16.D0*Hr3(
     &    -1,0,0)+32.D0*Hr3(-1,0,0)*x+16.D0*Hr3(-1,0,0)*x**2+32.D0*Hr3(
     &    0,-1,0)+32.D0*Hr3(0,-1,0)*x**2-9.D0*Hr3(0,0,0)+18.D0*Hr3(0,0,
     &    0)*x-36.D0*Hr3(0,0,0)*x**2-16.D0*Hr3(0,0,1)+32.D0*Hr3(0,0,1)*
     &    x-48.D0*Hr3(0,0,1)*x**2-14.D0*Hr3(0,1,0)+28.D0*Hr3(0,1,0)*x-
     &    32.D0*Hr3(0,1,0)*x**2-18.D0*Hr3(0,1,1)+36.D0*Hr3(0,1,1)*x-44.D
     &  0*Hr3(0,1,1)*x**2)
      H2g20=H2g20+cf*(-2.D0*Hr3(1,0,0)+4.D0*Hr3(1,0,0)*x-20.D0*Hr3(1,0,
     &    0)*x**2-24.D0*Hr3(1,0,1)+48.D0*Hr3(1,0,1)*x-48.D0*Hr3(1,0,1)*
     &    x**2-16.D0*Hr3(1,1,0)+32.D0*Hr3(1,1,0)*x-32.D0*Hr3(1,1,0)*
     &    x**2-22.D0*Hr3(1,1,1)+44.D0*Hr3(1,1,1)*x-44.D0*Hr3(1,1,1)*
     &    x**2)
      H2g20=H2g20+cf*Llam*(9.D0-17.D0*x+4.D0*x**2+2.D0*Hr1(0)-12.D0*
     &    Hr1(0)*x+20.D0*Hr1(0)*x**2+7.D0*Hr1(1)-24.D0*Hr1(1)*x+20.D0*
     &    Hr1(1)*x**2+4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x+16.D0*Hr2(0,0)*x**2
     &    +6.D0*Hr2(0,1)-12.D0*Hr2(0,1)*x+16.D0*Hr2(0,1)*x**2+8.D0*Hr2(
     &    1,0)-16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,1)-16.D0
     &    *Hr2(1,1)*x+16.D0*Hr2(1,1)*x**2)
      H2g20=H2g20+cf*Llam**2*(-1.D0/2.D0+2.D0*x-Hr1(0)+2.D0*Hr1(0)*x-4.D
     &  0*Hr1(0)*x**2-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*x**2)
      H2g20=H2g20+cf*LQm*(9.D0-17.D0*x+4.D0*x**2+2.D0*Hr1(0)-12.D0*Hr1(
     &    0)*x+20.D0*Hr1(0)*x**2+7.D0*Hr1(1)-24.D0*Hr1(1)*x+20.D0*Hr1(1
     &    )*x**2+4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x+16.D0*Hr2(0,0)*x**2+6.D0
     &    *Hr2(0,1)-12.D0*Hr2(0,1)*x+16.D0*Hr2(0,1)*x**2+8.D0*Hr2(1,0)-
     &    16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,1)-16.D0*Hr2(
     &    1,1)*x+16.D0*Hr2(1,1)*x**2)
      H2g20=H2g20+cf*LQm**2*(-1.D0/2.D0+2.D0*x-Hr1(0)+2.D0*Hr1(0)*x-4.D0
     &    *Hr1(0)*x**2-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*x**2)
      H2g20=H2g20+z3*ca*(12.D0-20.D0*x+24.D0*x**2)
      H2g20=H2g20+z3*cf*(30.D0+4.D0*x+76.D0*x**2)
      H2g20=H2g20+z2*ca*(8.D0-140.D0*x+147.D0*x**2-16.D0/3.D0*dx-6.D0*
     &    Hr1(-1)-12.D0*Hr1(-1)*x-20.D0*Hr1(-1)*x**2-8.D0*Hr1(0)-64.D0*
     &    Hr1(0)*x+16.D0*Hr1(0)*x**2+8.D0*Hr1(1)-16.D0*Hr1(1)*x+8.D0*
     &    Hr1(1)*x**2)
      H2g20=H2g20+z2*ca*Llam*(40.D0*x-16.D0*x**2)
      H2g20=H2g20+z2*ca*LQm*(-8.D0*x)
      H2g20=H2g20+z2*cf*(14.D0-140.D0/3.D0*x+78.D0*x**2+96.D0/5.D0*x**3
     &    -16.D0*Hr1(-1)-32.D0*Hr1(-1)*x-16.D0*Hr1(-1)*x**2+16.D0*Hr1(0
     &    )-32.D0*Hr1(0)*x+48.D0*Hr1(0)*x**2+8.D0*Hr1(1)-16.D0*Hr1(1)*x
     &    +32.D0*Hr1(1)*x**2)
      H2g20=H2g20+z2*cf*Llam*(-6.D0+12.D0*x-16.D0*x**2)
      H2g20=H2g20+z2*cf*LQm*(-6.D0+12.D0*x-16.D0*x**2)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC g coefficient function h20_g,3

      real*8 function h3g20(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*
      H3g20=
     & +Llam**2*(2.D0/3.D0-4.D0/3.D0*x+4.D0/3.D0*x**2)
      H3g20=H3g20+LQm*Llam*(-4.D0/3.D0+8.D0/3.D0*x-8.D0/3.D0*x**2)
      H3g20=H3g20+LQm**2*(2.D0/3.D0-4.D0/3.D0*x+4.D0/3.D0*x**2)
      H3g20=H3g20+ca*(1.D0/3.D0+157.D0/3.D0*x-1588.D0/27.D0*x**2+112.D0/
     &    27.D0*dx+14.D0/3.D0*Hr1(0)+43.D0/3.D0*Hr1(0)*x+400.D0/9.D0*
     &    Hr1(0)*x**2+Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*x**2+4.D0*Hr2(-1
     &    ,0)*x+4.D0*Hr2(-1,0)*x**2-Hr2(0,0)-4.D0*Hr2(0,0)*x-23.D0/3.D0
     &    *Hr2(0,0)*x**2-Hr2(0,1)*x**2-3.D0*Hr2(1,0)-16.D0*Hr2(1,0)*x+
     &    65.D0/3.D0*Hr2(1,0)*x**2-8.D0/3.D0*Hr2(1,0)*dx+Hr2(1,1)+4.D0*
     &    Hr2(1,1)*x-5.D0*Hr2(1,1)*x**2-4.D0*Hr3(-1,-1,0)-8.D0*Hr3(-1,
     &    -1,0)*x-8.D0*Hr3(-1,-1,0)*x**2+2.D0*Hr3(-1,0,0)+4.D0*Hr3(-1,0
     &    ,0)*x+4.D0*Hr3(-1,0,0)*x**2+2.D0*Hr3(0,0,0)+4.D0*Hr3(0,0,0)*x
     &    -4.D0*Hr3(0,1,0)-16.D0*Hr3(0,1,0)*x+2.D0*Hr3(1,0,1)-4.D0*Hr3(
     &    1,0,1)*x+4.D0*Hr3(1,0,1)*x**2+2.D0*Hr3(1,1,0)-4.D0*Hr3(1,1,0)
     &    *x+4.D0*Hr3(1,1,0)*x**2-2.D0*Hr3(1,1,1)+4.D0*Hr3(1,1,1)*x-4.D0
     &    *Hr3(1,1,1)*x**2)
      H3g20=H3g20+ca*Llam*(-4.D0+50.D0*x-436.D0/9.D0*x**2+40.D0/9.D0*dx
     &    +2.D0*Hr1(0)+16.D0*Hr1(0)*x+88.D0/3.D0*Hr1(0)*x**2+8.D0*Hr1(1
     &    )*x-8.D0*Hr1(1)*x**2+4.D0*Hr2(-1,0)+8.D0*Hr2(-1,0)*x+8.D0*
     &    Hr2(-1,0)*x**2-4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x-4.D0*Hr2(1,1)+8.D
     &  0*Hr2(1,1)*x-8.D0*Hr2(1,1)*x**2)
      H3g20=H3g20+ca*Llam**2*(1.D0+8.D0*x-31.D0/3.D0*x**2+4.D0/3.D0*dx+
     &    2.D0*Hr1(0)+8.D0*Hr1(0)*x-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(
     &    1)*x**2)
      H3g20=H3g20+ca*LQm*(4.D0-50.D0*x+436.D0/9.D0*x**2-40.D0/9.D0*dx-2.
     &  D0*Hr1(0)-16.D0*Hr1(0)*x-88.D0/3.D0*Hr1(0)*x**2-8.D0*Hr1(1)*x+8.
     &  D0*Hr1(1)*x**2-4.D0*Hr2(-1,0)-8.D0*Hr2(-1,0)*x-8.D0*Hr2(-1,0)*
     &    x**2+4.D0*Hr2(0,0)+8.D0*Hr2(0,0)*x+4.D0*Hr2(1,1)-8.D0*Hr2(1,1
     &    )*x+8.D0*Hr2(1,1)*x**2)
      H3g20=H3g20+ca*LQm*Llam*(-2.D0-16.D0*x+62.D0/3.D0*x**2-8.D0/3.D0*
     &    dx-4.D0*Hr1(0)-16.D0*Hr1(0)*x+4.D0*Hr1(1)-8.D0*Hr1(1)*x+8.D0*
     &    Hr1(1)*x**2)
      H3g20=H3g20+ca*LQm**2*(1.D0+8.D0*x-31.D0/3.D0*x**2+4.D0/3.D0*dx+2.
     &  D0*Hr1(0)+8.D0*Hr1(0)*x-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*
     &    x**2)
      H3g20=H3g20+cf*(-13.D0/2.D0+41.D0/2.D0*x-20.D0*x**2+4.D0*Hr1(0)+9.
     &  D0/2.D0*Hr1(0)*x+12.D0*Hr1(0)*x**2-13.D0*Hr1(1)*x+12.D0*Hr1(1)*
     &    x**2+1.D0/2.D0*Hr2(0,0)+6.D0*Hr2(0,0)*x-10.D0*Hr2(0,0)*x**2-2.
     &  D0*Hr2(0,1)-12.D0*Hr2(0,1)*x+6.D0*Hr2(0,1)*x**2-Hr2(1,0)+12.D0*
     &    Hr2(1,0)*x-10.D0*Hr2(1,0)*x**2-2.D0*Hr2(1,1)-4.D0*Hr2(1,1)*x+
     &    6.D0*Hr2(1,1)*x**2-Hr3(0,0,0)+2.D0*Hr3(0,0,0)*x-4.D0*Hr3(0,0,
     &    0)*x**2+2.D0*Hr3(0,1,0)-4.D0*Hr3(0,1,0)*x+2.D0*Hr3(0,1,1)-4.D0
     &    *Hr3(0,1,1)*x+4.D0*Hr3(0,1,1)*x**2-2.D0*Hr3(1,0,0)+4.D0*Hr3(1
     &    ,0,0)*x-4.D0*Hr3(1,0,0)*x**2+2.D0*Hr3(1,1,1)-4.D0*Hr3(1,1,1)*
     &    x+4.D0*Hr3(1,1,1)*x**2)
      H3g20=H3g20+cf*Llam*(9.D0-23.D0*x+10.D0*x**2+4.D0*Hr1(0)-16.D0*
     &    Hr1(0)*x+20.D0*Hr1(0)*x**2+7.D0*Hr1(1)-24.D0*Hr1(1)*x+20.D0*
     &    Hr1(1)*x**2+4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x+16.D0*Hr2(0,0)*x**2
     &    +6.D0*Hr2(0,1)-12.D0*Hr2(0,1)*x+16.D0*Hr2(0,1)*x**2+8.D0*Hr2(
     &    1,0)-16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,1)-16.D0
     &    *Hr2(1,1)*x+16.D0*Hr2(1,1)*x**2)
      H3g20=H3g20+cf*Llam**2*(-1.D0/2.D0+2.D0*x-Hr1(0)+2.D0*Hr1(0)*x-4.D
     &  0*Hr1(0)*x**2-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*x**2)
      H3g20=H3g20+cf*LQm*(-9.D0+23.D0*x-10.D0*x**2-4.D0*Hr1(0)+16.D0*
     &    Hr1(0)*x-20.D0*Hr1(0)*x**2-7.D0*Hr1(1)+24.D0*Hr1(1)*x-20.D0*
     &    Hr1(1)*x**2-4.D0*Hr2(0,0)+8.D0*Hr2(0,0)*x-16.D0*Hr2(0,0)*x**2
     &    -6.D0*Hr2(0,1)+12.D0*Hr2(0,1)*x-16.D0*Hr2(0,1)*x**2-8.D0*Hr2(
     &    1,0)+16.D0*Hr2(1,0)*x-16.D0*Hr2(1,0)*x**2-8.D0*Hr2(1,1)+16.D0
     &    *Hr2(1,1)*x-16.D0*Hr2(1,1)*x**2)
      H3g20=H3g20+cf*LQm**2*(1.D0/2.D0-2.D0*x+Hr1(0)-2.D0*Hr1(0)*x+4.D0
     &    *Hr1(0)*x**2+2.D0*Hr1(1)-4.D0*Hr1(1)*x+4.D0*Hr1(1)*x**2)
      H3g20=H3g20+z3*ca*(-10.D0-32.D0*x-4.D0*x**2)
      H3g20=H3g20+z3*cf*(2.D0-4.D0*x-4.D0*x**2)
      H3g20=H3g20+z2*ca*(5.D0*x**2-2.D0*Hr1(-1)-4.D0*Hr1(-1)*x-4.D0*
     &    Hr1(-1)*x**2)
      H3g20=H3g20+z2*ca*Llam*(4.D0+8.D0*x**2)
      H3g20=H3g20+z2*ca*LQm*(-4.D0-8.D0*x**2)
      H3g20=H3g20+z2*cf*(2.D0+12.D0*x-6.D0*x**2)
      H3g20=H3g20+z2*cf*Llam*(-6.D0+12.D0*x-16.D0*x**2)
      H3g20=H3g20+z2*cf*LQm*(6.D0-12.D0*x+16.D0*x**2)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC g coefficient function h21_g,1

      real*8 function h1g21(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*
      H1g21=
     & +Llam*(4.D0/3.D0-8.D0/3.D0*x+8.D0/3.D0*x**2)
      H1g21=H1g21+ca*(-50.D0/3.D0-94.D0/3.D0*x+434.D0/9.D0*x**2+16.D0/9.
     &  D0*dx-2.D0*Hr1(0)-48.D0*Hr1(0)*x+212.D0/3.D0*Hr1(0)*x**2+4.D0*
     &    Hr1(1)-56.D0*Hr1(1)*x+196.D0/3.D0*Hr1(1)*x**2-16.D0/3.D0*Hr1(
     &    1)*dx-4.D0*Hr2(-1,0)-8.D0*Hr2(-1,0)*x-8.D0*Hr2(-1,0)*x**2-12.D
     &  0*Hr2(0,0)-40.D0*Hr2(0,0)*x-48.D0*Hr2(0,1)*x+16.D0*Hr2(0,1)*
     &    x**2+8.D0*Hr2(1,0)-16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+12.D0
     &    *Hr2(1,1)-24.D0*Hr2(1,1)*x+24.D0*Hr2(1,1)*x**2)
      H1g21=H1g21+ca*Llam*(2.D0+16.D0*x-62.D0/3.D0*x**2+8.D0/3.D0*dx+4.D
     &  0*Hr1(0)+16.D0*Hr1(0)*x-4.D0*Hr1(1)+8.D0*Hr1(1)*x-8.D0*Hr1(1)*
     &    x**2)
      H1g21=H1g21+ca*LQm*(2.D0+16.D0*x-62.D0/3.D0*x**2+8.D0/3.D0*dx+4.D0
     &    *Hr1(0)+16.D0*Hr1(0)*x-4.D0*Hr1(1)+8.D0*Hr1(1)*x-8.D0*Hr1(1)*
     &    x**2)
      H1g21=H1g21+cf*(5.D0-21.D0*x+12.D0*x**2+2.D0*Hr1(0)-20.D0*Hr1(0)*
     &    x+20.D0*Hr1(0)*x**2+7.D0*Hr1(1)-24.D0*Hr1(1)*x+20.D0*Hr1(1)*
     &    x**2+4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x+16.D0*Hr2(0,0)*x**2+6.D0*
     &    Hr2(0,1)-12.D0*Hr2(0,1)*x+16.D0*Hr2(0,1)*x**2+8.D0*Hr2(1,0)-
     &    16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,1)-16.D0*Hr2(
     &    1,1)*x+16.D0*Hr2(1,1)*x**2)
      H1g21=H1g21+cf*Llam*(-1.D0+4.D0*x-2.D0*Hr1(0)+4.D0*Hr1(0)*x-8.D0*
     &    Hr1(0)*x**2-4.D0*Hr1(1)+8.D0*Hr1(1)*x-8.D0*Hr1(1)*x**2)
      H1g21=H1g21+z2*ca*(40.D0*x-16.D0*x**2)
      H1g21=H1g21+z2*cf*(-6.D0+12.D0*x-16.D0*x**2)
      H1g21=H1g21-4.D0/3.D0+16.D0/3.D0*x-16.D0/3.D0*x**2-4.D0/3.D0*Hr1(
     &    0)+8.D0/3.D0*Hr1(0)*x-8.D0/3.D0*Hr1(0)*x**2-4.D0/3.D0*Hr1(1)+
     &    8.D0/3.D0*Hr1(1)*x-8.D0/3.D0*Hr1(1)*x**2

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC g coefficient function h21_g,2

      real*8 function h2g21(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*
      H2g21=
     & +Llam*(4.D0/3.D0-8.D0/3.D0*x+8.D0/3.D0*x**2)
      H2g21=H2g21+ca*(-98.D0/3.D0-334.D0/3.D0*x+1250.D0/9.D0*x**2+64.D0/
     &    9.D0*dx-2.D0*Hr1(0)-112.D0*Hr1(0)*x+212.D0/3.D0*Hr1(0)*x**2+4.
     &  D0*Hr1(1)-88.D0*Hr1(1)*x+292.D0/3.D0*Hr1(1)*x**2-16.D0/3.D0*
     &    Hr1(1)*dx-4.D0*Hr2(-1,0)-8.D0*Hr2(-1,0)*x-8.D0*Hr2(-1,0)*x**2
     &    -12.D0*Hr2(0,0)-40.D0*Hr2(0,0)*x-48.D0*Hr2(0,1)*x+16.D0*Hr2(0
     &    ,1)*x**2+8.D0*Hr2(1,0)-16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+
     &    12.D0*Hr2(1,1)-24.D0*Hr2(1,1)*x+24.D0*Hr2(1,1)*x**2)
      H2g21=H2g21+ca*Llam*(2.D0+16.D0*x-62.D0/3.D0*x**2+8.D0/3.D0*dx+4.D
     &  0*Hr1(0)+16.D0*Hr1(0)*x-4.D0*Hr1(1)+8.D0*Hr1(1)*x-8.D0*Hr1(1)*
     &    x**2)
      H2g21=H2g21+ca*LQm*(2.D0+16.D0*x-62.D0/3.D0*x**2+8.D0/3.D0*dx+4.D0
     &    *Hr1(0)+16.D0*Hr1(0)*x-4.D0*Hr1(1)+8.D0*Hr1(1)*x-8.D0*Hr1(1)*
     &    x**2)
      H2g21=H2g21+cf*(9.D0-17.D0*x+4.D0*x**2+2.D0*Hr1(0)-12.D0*Hr1(0)*x
     &    +20.D0*Hr1(0)*x**2+7.D0*Hr1(1)-24.D0*Hr1(1)*x+20.D0*Hr1(1)*
     &    x**2+4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x+16.D0*Hr2(0,0)*x**2+6.D0*
     &    Hr2(0,1)-12.D0*Hr2(0,1)*x+16.D0*Hr2(0,1)*x**2+8.D0*Hr2(1,0)-
     &    16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,1)-16.D0*Hr2(
     &    1,1)*x+16.D0*Hr2(1,1)*x**2)
      H2g21=H2g21+cf*Llam*(-1.D0+4.D0*x-2.D0*Hr1(0)+4.D0*Hr1(0)*x-8.D0*
     &    Hr1(0)*x**2-4.D0*Hr1(1)+8.D0*Hr1(1)*x-8.D0*Hr1(1)*x**2)
      H2g21=H2g21+z2*ca*(40.D0*x-16.D0*x**2)
      H2g21=H2g21+z2*cf*(-6.D0+12.D0*x-16.D0*x**2)
      H2g21=H2g21-4.D0/3.D0+32.D0/3.D0*x-32.D0/3.D0*x**2-4.D0/3.D0*Hr1(
     &    0)+8.D0/3.D0*Hr1(0)*x-8.D0/3.D0*Hr1(0)*x**2-4.D0/3.D0*Hr1(1)+
     &    8.D0/3.D0*Hr1(1)*x-8.D0/3.D0*Hr1(1)*x**2

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC g coefficient function h21_g,3

      real*8 function h3g21(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*
      H3g21=
     & +Llam*(4.D0/3.D0-8.D0/3.D0*x+8.D0/3.D0*x**2)
      H3g21=H3g21+LQm*(-4.D0/3.D0+8.D0/3.D0*x-8.D0/3.D0*x**2)
      H3g21=H3g21+ca*(-4.D0+50.D0*x-436.D0/9.D0*x**2+40.D0/9.D0*dx+2.D0
     &    *Hr1(0)+16.D0*Hr1(0)*x+88.D0/3.D0*Hr1(0)*x**2+8.D0*Hr1(1)*x-8.
     &  D0*Hr1(1)*x**2+4.D0*Hr2(-1,0)+8.D0*Hr2(-1,0)*x+8.D0*Hr2(-1,0)*
     &    x**2-4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x-4.D0*Hr2(1,1)+8.D0*Hr2(1,1
     &    )*x-8.D0*Hr2(1,1)*x**2)
      H3g21=H3g21+ca*Llam*(2.D0+16.D0*x-62.D0/3.D0*x**2+8.D0/3.D0*dx+4.D
     &  0*Hr1(0)+16.D0*Hr1(0)*x-4.D0*Hr1(1)+8.D0*Hr1(1)*x-8.D0*Hr1(1)*
     &    x**2)
      H3g21=H3g21+ca*LQm*(-2.D0-16.D0*x+62.D0/3.D0*x**2-8.D0/3.D0*dx-4.D
     &  0*Hr1(0)-16.D0*Hr1(0)*x+4.D0*Hr1(1)-8.D0*Hr1(1)*x+8.D0*Hr1(1)*
     &    x**2)
      H3g21=H3g21+cf*(9.D0-23.D0*x+10.D0*x**2+4.D0*Hr1(0)-16.D0*Hr1(0)*
     &    x+20.D0*Hr1(0)*x**2+7.D0*Hr1(1)-24.D0*Hr1(1)*x+20.D0*Hr1(1)*
     &    x**2+4.D0*Hr2(0,0)-8.D0*Hr2(0,0)*x+16.D0*Hr2(0,0)*x**2+6.D0*
     &    Hr2(0,1)-12.D0*Hr2(0,1)*x+16.D0*Hr2(0,1)*x**2+8.D0*Hr2(1,0)-
     &    16.D0*Hr2(1,0)*x+16.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,1)-16.D0*Hr2(
     &    1,1)*x+16.D0*Hr2(1,1)*x**2)
      H3g21=H3g21+cf*Llam*(-1.D0+4.D0*x-2.D0*Hr1(0)+4.D0*Hr1(0)*x-8.D0*
     &    Hr1(0)*x**2-4.D0*Hr1(1)+8.D0*Hr1(1)*x-8.D0*Hr1(1)*x**2)
      H3g21=H3g21+z2*ca*(4.D0+8.D0*x**2)
      H3g21=H3g21+z2*cf*(-6.D0+12.D0*x-16.D0*x**2)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC g coefficient function h22_g

      real*8 function hg22(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*

      Hg22=
     & +ca*(1.D0+8.D0*x-31.D0/3.D0*x**2+4.D0/3.D0*dx+2.D0*Hr1(0)+8.D0*
     &    Hr1(0)*x-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*x**2)
      Hg22=Hg22+cf*(-1.D0/2.D0+2.D0*x-Hr1(0)+2.D0*Hr1(0)*x-4.D0*Hr1(0)*
     &    x**2-2.D0*Hr1(1)+4.D0*Hr1(1)*x-4.D0*Hr1(1)*x**2)
      Hg22=Hg22+2.D0/3.D0-4.D0/3.D0*x+4.D0/3.D0*x**2

*
       return
       end
*
* ---------------------------------------------------------------------



* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC ps coefficient function h20_ps,1

      real*8 function h1qps20(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*

      H1qPS20=
     & +cf*(107.D0/9.D0-323.D0/9.D0*x+368.D0/27.D0*x**2+280.D0/27.D0*dx
     &    +106.D0/3.D0*Hr1(0)-62.D0/3.D0*Hr1(0)*x+16.D0/3.D0*Hr1(0)*
     &    x**2+56.D0/3.D0*Hr1(1)-80.D0/3.D0*Hr1(1)*x+128.D0/9.D0*Hr1(1)
     &    *x**2-56.D0/9.D0*Hr1(1)*dx-16.D0*Hr2(-1,0)-16.D0*Hr2(-1,0)*x-
     &    16.D0/3.D0*Hr2(-1,0)*x**2-16.D0/3.D0*Hr2(-1,0)*dx-Hr2(0,0)+3.D
     &  0*Hr2(0,0)*x-56.D0/3.D0*Hr2(0,0)*x**2-16.D0*Hr2(0,1)*x-16.D0*
     &    Hr2(0,1)*x**2+6.D0*Hr2(1,0)-6.D0*Hr2(1,0)*x-8.D0*Hr2(1,0)*
     &    x**2+8.D0*Hr2(1,0)*dx+4.D0*Hr2(1,1)-4.D0*Hr2(1,1)*x-16.D0/3.D0
     &    *Hr2(1,1)*x**2+16.D0/3.D0*Hr2(1,1)*dx+18.D0*Hr3(0,0,0)+18.D0*
     &    Hr3(0,0,0)*x+16.D0*Hr3(0,0,1)+16.D0*Hr3(0,0,1)*x+12.D0*Hr3(0,
     &    1,0)+12.D0*Hr3(0,1,0)*x+8.D0*Hr3(0,1,1)+8.D0*Hr3(0,1,1)*x)
      H1qPS20=H1qPS20+cf*Llam*(-44.D0/3.D0+44.D0/3.D0*x-16.D0/9.D0*x**2
     &    +16.D0/9.D0*dx-2.D0*Hr1(0)+6.D0*Hr1(0)*x+32.D0/3.D0*Hr1(0)*
     &    x**2-4.D0*Hr1(1)+4.D0*Hr1(1)*x+16.D0/3.D0*Hr1(1)*x**2-16.D0/3.
     &  D0*Hr1(1)*dx-12.D0*Hr2(0,0)-12.D0*Hr2(0,0)*x-8.D0*Hr2(0,1)-8.D0
     &    *Hr2(0,1)*x)
      H1qPS20=H1qPS20+cf*Llam**2*(1.D0-x-4.D0/3.D0*x**2+4.D0/3.D0*dx+2.D
     &  0*Hr1(0)+2.D0*Hr1(0)*x)
      H1qPS20=H1qPS20+cf*LQm*(-4.D0+12.D0*x-112.D0/9.D0*x**2+40.D0/9.D0
     &    *dx+2.D0*Hr1(0)+10.D0*Hr1(0)*x+16.D0/3.D0*Hr1(0)*x**2-4.D0*
     &    Hr2(0,0)-4.D0*Hr2(0,0)*x)
      H1qPS20=H1qPS20+cf*LQm*Llam*(2.D0-2.D0*x-8.D0/3.D0*x**2+8.D0/3.D0
     &    *dx+4.D0*Hr1(0)+4.D0*Hr1(0)*x)
      H1qPS20=H1qPS20+cf*LQm**2*(-1.D0+x+4.D0/3.D0*x**2-4.D0/3.D0*dx-2.D
     &  0*Hr1(0)-2.D0*Hr1(0)*x)
      H1qPS20=H1qPS20+z2*cf*(16.D0*x**2-16.D0/3.D0*dx-16.D0*Hr1(0)-16.D0
     &    *Hr1(0)*x)
      H1qPS20=H1qPS20+z2*cf*Llam*(8.D0+8.D0*x)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC ps coefficient function h20_ps,2

      real*8 function h2qps20(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*

      H2qPS20=
     & +cf*(155.D0/9.D0-515.D0/9.D0*x+848.D0/27.D0*x**2+232.D0/27.D0*dx
     &    +154.D0/3.D0*Hr1(0)-110.D0/3.D0*Hr1(0)*x-80.D0/3.D0*Hr1(0)*
     &    x**2+104.D0/3.D0*Hr1(1)-80.D0/3.D0*Hr1(1)*x+32.D0/9.D0*Hr1(1)
     &    *x**2-104.D0/9.D0*Hr1(1)*dx-16.D0*Hr2(-1,0)-16.D0*Hr2(-1,0)*x
     &    -16.D0/3.D0*Hr2(-1,0)*x**2-16.D0/3.D0*Hr2(-1,0)*dx-Hr2(0,0)+
     &    35.D0*Hr2(0,0)*x-56.D0/3.D0*Hr2(0,0)*x**2-16.D0*Hr2(0,1)*x**2
     &    +6.D0*Hr2(1,0)-6.D0*Hr2(1,0)*x-8.D0*Hr2(1,0)*x**2+8.D0*Hr2(1,
     &    0)*dx+4.D0*Hr2(1,1)-4.D0*Hr2(1,1)*x-16.D0/3.D0*Hr2(1,1)*x**2+
     &    16.D0/3.D0*Hr2(1,1)*dx+18.D0*Hr3(0,0,0)+18.D0*Hr3(0,0,0)*x+16.
     &  D0*Hr3(0,0,1)+16.D0*Hr3(0,0,1)*x+12.D0*Hr3(0,1,0)+12.D0*Hr3(0,1
     &    ,0)*x+8.D0*Hr3(0,1,1)+8.D0*Hr3(0,1,1)*x)
      H2qPS20=H2qPS20+cf*Llam*(-92.D0/3.D0+44.D0/3.D0*x+80.D0/9.D0*x**2
     &    +64.D0/9.D0*dx-2.D0*Hr1(0)-10.D0*Hr1(0)*x+32.D0/3.D0*Hr1(0)*
     &    x**2-4.D0*Hr1(1)+4.D0*Hr1(1)*x+16.D0/3.D0*Hr1(1)*x**2-16.D0/3.
     &  D0*Hr1(1)*dx-12.D0*Hr2(0,0)-12.D0*Hr2(0,0)*x-8.D0*Hr2(0,1)-8.D0
     &    *Hr2(0,1)*x)
      H2qPS20=H2qPS20+cf*Llam**2*(1.D0-x-4.D0/3.D0*x**2+4.D0/3.D0*dx+2.D
     &  0*Hr1(0)+2.D0*Hr1(0)*x)
      H2qPS20=H2qPS20+cf*LQm*(-4.D0+12.D0*x-112.D0/9.D0*x**2+40.D0/9.D0
     &    *dx+2.D0*Hr1(0)+10.D0*Hr1(0)*x+16.D0/3.D0*Hr1(0)*x**2-4.D0*
     &    Hr2(0,0)-4.D0*Hr2(0,0)*x)
      H2qPS20=H2qPS20+cf*LQm*Llam*(2.D0-2.D0*x-8.D0/3.D0*x**2+8.D0/3.D0
     &    *dx+4.D0*Hr1(0)+4.D0*Hr1(0)*x)
      H2qPS20=H2qPS20+cf*LQm**2*(-1.D0+x+4.D0/3.D0*x**2-4.D0/3.D0*dx-2.D
     &  0*Hr1(0)-2.D0*Hr1(0)*x)
      H2qPS20=H2qPS20+z2*cf*(-16.D0*x+16.D0*x**2-16.D0/3.D0*dx-16.D0*
     &    Hr1(0)-16.D0*Hr1(0)*x)
      H2qPS20=H2qPS20+z2*cf*Llam*(8.D0+8.D0*x)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC ps coefficient function h20_ps,3

      real*8 function h3qps20(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*

      H3qPS20=
     & +cf*(1.D0/3.D0+31.D0/3.D0*x-400.D0/27.D0*x**2+112.D0/27.D0*dx+14.
     &  D0/3.D0*Hr1(0)+22.D0/3.D0*Hr1(0)*x+112.D0/9.D0*Hr1(0)*x**2-Hr2(
     &    0,0)-5.D0*Hr2(0,0)*x-8.D0/3.D0*Hr2(0,0)*x**2-2.D0*Hr2(1,0)+2.D
     &  0*Hr2(1,0)*x+8.D0/3.D0*Hr2(1,0)*x**2-8.D0/3.D0*Hr2(1,0)*dx+2.D0
     &    *Hr3(0,0,0)+2.D0*Hr3(0,0,0)*x-4.D0*Hr3(0,1,0)-4.D0*Hr3(0,1,0)
     &    *x)
      H3qPS20=H3qPS20+cf*Llam*(-4.D0+12.D0*x-112.D0/9.D0*x**2+40.D0/9.D0
     &    *dx+2.D0*Hr1(0)+10.D0*Hr1(0)*x+16.D0/3.D0*Hr1(0)*x**2-4.D0*
     &    Hr2(0,0)-4.D0*Hr2(0,0)*x)
      H3qPS20=H3qPS20+cf*Llam**2*(1.D0-x-4.D0/3.D0*x**2+4.D0/3.D0*dx+2.D
     &  0*Hr1(0)+2.D0*Hr1(0)*x)
      H3qPS20=H3qPS20+cf*LQm*(4.D0-12.D0*x+112.D0/9.D0*x**2-40.D0/9.D0*
     &    dx-2.D0*Hr1(0)-10.D0*Hr1(0)*x-16.D0/3.D0*Hr1(0)*x**2+4.D0*
     &    Hr2(0,0)+4.D0*Hr2(0,0)*x)
      H3qPS20=H3qPS20+cf*LQm*Llam*(-2.D0+2.D0*x+8.D0/3.D0*x**2-8.D0/3.D0
     &    *dx-4.D0*Hr1(0)-4.D0*Hr1(0)*x)
      H3qPS20=H3qPS20+cf*LQm**2*(1.D0-x-4.D0/3.D0*x**2+4.D0/3.D0*dx+2.D0
     &    *Hr1(0)+2.D0*Hr1(0)*x)
      H3qPS20=H3qPS20+z3*cf*(-8.D0-8.D0*x)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC ps coefficient function h21_ps,1

      real*8 function h1qps21(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*

      H1qPS21=
     & +cf*(-44.D0/3.D0+44.D0/3.D0*x-16.D0/9.D0*x**2+16.D0/9.D0*dx-2.D0
     &    *Hr1(0)+6.D0*Hr1(0)*x+32.D0/3.D0*Hr1(0)*x**2-4.D0*Hr1(1)+4.D0
     &    *Hr1(1)*x+16.D0/3.D0*Hr1(1)*x**2-16.D0/3.D0*Hr1(1)*dx-12.D0*
     &    Hr2(0,0)-12.D0*Hr2(0,0)*x-8.D0*Hr2(0,1)-8.D0*Hr2(0,1)*x)
      H1qPS21=H1qPS21+cf*Llam*(2.D0-2.D0*x-8.D0/3.D0*x**2+8.D0/3.D0*dx+
     &    4.D0*Hr1(0)+4.D0*Hr1(0)*x)
      H1qPS21=H1qPS21+cf*LQm*(2.D0-2.D0*x-8.D0/3.D0*x**2+8.D0/3.D0*dx+4.
     &  D0*Hr1(0)+4.D0*Hr1(0)*x)
      H1qPS21=H1qPS21+z2*cf*(8.D0+8.D0*x)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC ps coefficient function h21_ps,2

      real*8 function h2qps21(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*

      H2qPS21=
     & +cf*(-92.D0/3.D0+44.D0/3.D0*x+80.D0/9.D0*x**2+64.D0/9.D0*dx-2.D0
     &    *Hr1(0)-10.D0*Hr1(0)*x+32.D0/3.D0*Hr1(0)*x**2-4.D0*Hr1(1)+4.D0
     &    *Hr1(1)*x+16.D0/3.D0*Hr1(1)*x**2-16.D0/3.D0*Hr1(1)*dx-12.D0*
     &    Hr2(0,0)-12.D0*Hr2(0,0)*x-8.D0*Hr2(0,1)-8.D0*Hr2(0,1)*x)
      H2qPS21=H2qPS21+cf*Llam*(2.D0-2.D0*x-8.D0/3.D0*x**2+8.D0/3.D0*dx+
     &    4.D0*Hr1(0)+4.D0*Hr1(0)*x)
      H2qPS21=H2qPS21+cf*LQm*(2.D0-2.D0*x-8.D0/3.D0*x**2+8.D0/3.D0*dx+4.
     &  D0*Hr1(0)+4.D0*Hr1(0)*x)
      H2qPS21=H2qPS21+z2*cf*(8.D0+8.D0*x)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC ps coefficient function h21_ps,3

      real*8 function h3qps21(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*

      H3qPS21=
     & +cf*(-4.D0+12.D0*x-112.D0/9.D0*x**2+40.D0/9.D0*dx+2.D0*Hr1(0)+10.
     &  D0*Hr1(0)*x+16.D0/3.D0*Hr1(0)*x**2-4.D0*Hr2(0,0)-4.D0*Hr2(0,0)*
     &    x)
      H3qPS21=H3qPS21+cf*Llam*(2.D0-2.D0*x-8.D0/3.D0*x**2+8.D0/3.D0*dx+
     &    4.D0*Hr1(0)+4.D0*Hr1(0)*x)
      H3qPS21=H3qPS21+cf*LQm*(-2.D0+2.D0*x+8.D0/3.D0*x**2-8.D0/3.D0*dx-
     &    4.D0*Hr1(0)-4.D0*Hr1(0)*x)

*
       return
       end
*
* ---------------------------------------------------------------------
* .. The asymptotics of the regular piece of the 
* .. 2-loop CC ps coefficient function h22_ps

      real*8 function hqps22(x,lambda)
      implicit none
      real*8 x, dx, lambda, Llam, LQm
      real*8 z2, z3, z4, ca, cf
*
      complex*16 hc1, hc2, hc3, hc4
      real*8 hr1, hr2, hr3, hr4, hi1, hi2, hi3, hi4
      integer n1, n2, nw, i1, i2, i3
      parameter ( n1 = -1, n2 = 1, nw = 4 )
      dimension hc1(n1:n2),hc2(n1:n2,n1:n2),hc3(n1:n2,n1:n2,n1:n2),
     ,          hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hr1(n1:n2),hr2(n1:n2,n1:n2),hr3(n1:n2,n1:n2,n1:n2),
     ,          hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension hi1(n1:n2),hi2(n1:n2,n1:n2),hi3(n1:n2,n1:n2,n1:n2),
     ,          hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      parameter ( z2 = 1.6449 34066 84822 64365 d0,
     ,     z3 = 1.2020 56903 15959 42854 d0,
     ,     z4 = 1.0823 23233 71113 81916 d0 )

*
*  ...Colour factors and abbreviation
*
*     lambda=Q2/(Q2+m2), LQm = ln(Q2/m2)
      cf = 4./3.d0
      ca = 3.d0
      dx = 1.D0/x
      Llam = dlog(lambda)
      LQm = dlog(lambda/(1d0-lambda))
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       call hplog (x, nw, hc1,hc2,hc3,hc4, hr1,hr2,hr3,hr4,
     ,            hi1,hi2,hi3,hi4, n1, n2)
*
* ...The coefficient function in terms of the harmonic polylogs
*

      HqPS22=
     & +cf*(1.D0-x-4.D0/3.D0*x**2+4.D0/3.D0*dx+2.D0*Hr1(0)+2.D0*Hr1(0)*
     &    x)

*
       return
       end
*
* ---------------------------------------------------------------------
