      double precision function func24(id,x)

      implicit double precision (a-h,o-z)
      include 'pdfparam.inc'


      if (id.eq.0) func24=gluon(x)
      if (id.eq.1) func24=dval(x)
      if (id.eq.2) func24=uval(x)
      if (id.eq.3) func24=sea(x)
      if (id.eq.4) func24=dbmub(x)
      if (id.eq.5) func24=0.d0
      if (id.eq.6) func24=0.d0

      return
      end
