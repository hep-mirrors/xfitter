      double precision function func1(id,x)

      implicit double precision (a-h,o-z)
      include 'pdfparam.inc'

      if (id.eq.0) func1=gluon(x)
      if (id.eq.1) func1=H1D(x)
      if (id.eq.2) func1=H1U(x)
      if (id.eq.3) func1=Ubar(x)
      if (id.eq.4) func1=Dbar(x)
      if (id.eq.6) func1=2*qstrange(x)
      if (id.eq.5) func1=0.d0


      return
      end
