C     evolution parameters
      data mxord/4/                                   !maximum order of evolution
      data iordqcd/2/, iordqed/1/, nfin/0/, iosp/3/   !order, VFNS, spline
c      data q2c/2.001D0/,q2b/20.25D0/,q2t/30625.D0/,q0/2.0D0/
      data q2c/2.0449D0/, q2b/18.49D0/,q2t/32400.D0/, q0/2.0D0/     !thresh and mu20
      data amu2/1.0D0/, bmu2/0.0D0/                   !renor scale
      data as0/0.1205/, r20/8315.17839D0/             !alphas/(2*pi)
c      data aem0/0.00119306554042/,rem20/3.157729D0/   !alphaem/(2*pi)
      data aem0/0.00116171491308/,rem20/3.157729D0/   !alphaem/(2*pi)

C     x-grid parameters
      dimension xmin(5), iwt(5)
      data xmin/1.D-5,0.2D0,0.4D0,0.6D0,0.75D0/
      data iwt/1,2,4,8,16/
      data nxin/1000/

C     mu2-grid parameters
      dimension qq(4),wt(4)
c      data qq/2.D0,2.0449D0,18.49D0,32400.D0,1.D6/
      data qq/2.D0,2.0449D0,18.49D0,1.D4/
      data wt/1.D0,1D0,1D0,1D0/
c      data wt/10.D0,0.1D0,0.1D0,0.1D0,0.1D0/
      data nqin/1000/

C     user setup name for evolution
      character setup*(*), scheme*(*)
      parameter (setup = 'mrst')
      parameter (scheme = 'VFNS')

      common /aem/ aem0,rem20,q2b,q2t
