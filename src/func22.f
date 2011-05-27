      double precision function func22(id,x)

      implicit double precision (a-h,o-z)
      include 'pdfparam.inc'
      include 'steering.inc'

      func22 = 0.D0
      if (id.eq.0) func22=gluon(x)
      if (id.eq.1) func22=dval(x)
      if (id.eq.2) func22=uval(x)
      if (id.eq.3) func22=2*qstrange(x)
      if (id.eq.4) func22=Ubar(x)
      if (id.eq.5) func22=Dbar(x)
      if (id.eq.6) func22=0.d0

      return
      end
