c=======================================================
      real*8 function coeff3(y,eps)
      implicit real*8(a-h,o-z)
      dimension z1(17),z2(17)

      data z1/-0.857E+00,0.217E+02,-0.199E+02,0.162E+02,
     .0.204E+00,-0.142E+03,0.101E+03,-0.575E+02,0.832E+00,
     .0.175E+03,-0.139E+03,0.604E+02,0.344E-01,-0.648E+02,
     .0.538E+02,-0.160E+02,0.000E+00/
      data z2/-0.334E-03,0.677E-02,-0.108E-02,0.825E-03,
     .-0.425E+04,0.883E+05,-0.562E+05,0.601E+05,0.562E+07,
     .-0.122E+09,0.912E+08,-0.928E+08,-0.236E+10,0.499E+11,
     .-0.449E+11,0.412E+11,0.000E+00/

      epsl=eps
      eps4=4.d0*eps
      x0=1.d0/(1.d0+eps)
      yx0=y/x0
      yx01=1.d0-yx0
      epsl2=epsl*epsl
      epsl3=epsl2*epsl
      if(y.lt.0.05) go to 100    
      a0=z1(1)+z1(2)*epsl+z1(3)*epsl2+z1(4)*epsl3
      a1=z1(5)+z1(6)*epsl+z1(7)*epsl2+z1(8)*epsl3
      a2=z1(9)+z1(10)*epsl+z1(11)*epsl2+z1(12)*epsl3
      a3=z1(13)+z1(14)*epsl+z1(15)*epsl2+z1(16)*epsl3
      fac=a0+a1*yx0+a2*yx0*yx0+a3*yx0*yx0*yx0
      coeff3=fac*yx01**z1(17)
      return
  100 continue
      a0=z2(1)+z2(2)*epsl+z2(3)*epsl2+z2(4)*epsl3
      a1=z2(5)+z2(6)*epsl+z2(7)*epsl2+z2(8)*epsl3
      a2=z2(9)+z2(10)*epsl+z2(11)*epsl2+z2(12)*epsl3
      a3=z2(13)+z2(14)*epsl+z2(15)*epsl2+z2(16)*epsl3
      fac=a0+a1*yx0+a2*yx0*yx0+a3*yx0*yx0*yx0
      coeff3=fac*yx01**z2(17)
      return
      end
