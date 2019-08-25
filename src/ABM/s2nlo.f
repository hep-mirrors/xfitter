      real*8 function s2nlo(x)
      implicit real*8 (a-h,o-z)

      a=x/(1.+x)
      b=1./(1.+x)
      s2nlo=ddilog(a)-ddilog(b)
     ++0.5d0*log(x)*log(x/(1.d0+x)**2)

      return
      end
 
