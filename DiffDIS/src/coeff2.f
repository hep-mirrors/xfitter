c=======================================================
      real*8 function coeff2(y,eps)
      implicit real*8(a-h,o-z)
      dimension z1(17),z2(17)
  
      data z1/-0.183E+01,0.400E+01,0.159E+02,-0.357E+02,
     .-0.186E+00,-0.988E+02,0.712E+02,0.631E+02,-0.136E+01,
     .0.175E+03,-0.158E+03,-0.433E+02,0.375E+01,-0.913E+02,
     .0.842E+02,0.107E+02,0.629E+00/
      data z2/-0.204E+01,-0.127E+01,0.117E+01,-0.208E+00,
     .-0.228E+02,0.261E+03,-0.340E+03,0.153E+03,-0.591E+02,
     .-0.199E+04,-0.113E+04,0.299E+04,-0.408E+04,-0.446E+04,
     . 0.307E+05,-0.282E+05,0.144E+02/

      epsl=eps
      eps4=4.d0*eps
      x0=1.d0/(1.d0+eps4)
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
      coeff2=fac*yx01**z1(17)
      return
  100 continue
      a0=z2(1)+z2(2)*epsl+z2(3)*epsl2+z2(4)*epsl3
      a1=z2(5)+z2(6)*epsl+z2(7)*epsl2+z2(8)*epsl3
      a2=z2(9)+z2(10)*epsl+z2(11)*epsl2+z2(12)*epsl3
      a3=z2(13)+z2(14)*epsl+z2(15)*epsl2+z2(16)*epsl3
      fac=a0+a1*yx0+a2*yx0*yx0+a3*yx0*yx0*yx0
      coeff2=fac*yx01**z2(17)
      return
      end
