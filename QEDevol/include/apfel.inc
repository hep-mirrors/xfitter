C     evolution parameters
      data mxord/4/                                          !maximum order of evolution
      data iordqcd/3/,iordqed/1/,nfin/0/,iosp/3/             !order, VFNS, spline
      data q2c/2.001D0/,q2b/20.25D0/,q2t/30625.D0/,q0/2.0D0/ !thresh and mu20
      data amu2/1.0D0/,bmu2/0.0D0/                           !renor scale
      data as0/0.35/,r20/2.D0/                               !alphas/(2*pi)
      data aem0/0.00119306554042/,rem20/3.157729D0/          !alphaem/(2*pi)

C     x-grid parameters
      dimension xmin(5), iwt(5)
      data xmin/1.D-5,0.2D0,0.4D0,0.6D0,0.75D0/
      data iwt/1,2,4,8,16/
      data nxin/1000/

C     mu2-grid parameters
      dimension qq(5),wt(5)
      data qq/2.D0,2.0101D0,20.25D0,30625.D0,1.D6/
      data wt/10.D0,0.1D0,0.1D0,0.1D0,0.1D0/
      data nqin/1000/

C     user setup name for evolution
      character setup*(*), scheme*(*)
      parameter (setup = 'apfel')
      parameter (scheme = 'VFNS')

      common /aem/ aem0,rem20,q2b,q2t
