c---  ws
c---  for easy call from C
c=================================================
      subroutine FixQuark(iq)
      implicit none
      integer iq
      real*8 scale, alphas, nlf, eh, eh2
      character*10 q
      common/coupling/scale, alphas, nlf, eh2
      common/quark/q
      if (iq .eq. 4) then
        q = 'charm'
        eh = 2d0/3d0
      elseif (iq .eq. 5) then
        q = 'bottom'
        eh = -1d0/3d0
      elseif (iq .eq. 6) then
        q = 'top'
        eh = 2d0/3d0
      else
        PRINT*,'Illegal quark id',iq
        STOP
      endif
      eh2 = eh*eh
      END

c  calculates structure functions for heavy quark production
c  gamma* parton -> Q + Qbar + X folded with the parton densities of
c  the proton.  These routines can only be used with MSbar parton
c  densities of the proton.
c  NPB392 (1993) 162 - 229, PL 347B 143 - 151.
c  Last modified 23.7.96.  If you have questions, e-mail
c  riemersm@hades.ifh.de
c      program struct

c=================================================
      subroutine js(xbjfix,mq2fix,fhad)
      implicit none
c declare all variables for easier debugging
      integer ifch, pts, its, ihist
      integer nborn, ngcorr, nqcorr
      integer iloop, iloop1, nscale, ipdf, ichoose
      integer kount, max, nxch
      real*8 xbj, mq2, m2, tf, scale, alphas, mq2l(4)
      real*8 nlf, eh, eh2, c, n
      real*8 fhad, xfhad
      real*8 xbjfix, mq2fix
      real*8 pi, ca, cf
      real*8 qfac, mfac
      real*8 daind1, eps, zmin, zmax, est
      real*8 xchbin(11),fchfac
      data xchbin/0.00005d0,0.00013d0,0.0003d0,0.0005d0,0.0008d0,
     > 0.0012d0,0.002d0,0.0032d0,0.004d0,0.008d0,0.020d0/
      character*10 q
      parameter (pi = 3.1415926535897932384626433832795D0)
c group structure constants ca = N of SU(N) = 3, cf = (N^2 - 1)/(2 N) =
c 4/3 and tf = normalization for the color matrices = 1/2 in SU(N) theories
      parameter (ca = 3d0, cf = 4d0/3d0, tf = 0.5d0)
c integrating the the function ff gives the structure function at points 
c in xbj (Bjorken x) and mq2 (Q^2).
      external ffm
      common/include/nborn, ngcorr, nqcorr, ifch
      common/mass/m2, mq2, xbj
      common/coupling/scale, alphas, nlf, eh2
      common/compdf/ipdf
      common/quark/q
      COMMON/JSPARAM/qfac,mfac
c      ,ichoose

c call the input file struct_had.dat which sets the parameters for the
c calculation, and struct_had.plo which is a file that tells the
c integration routine VEGAS to give differential distributions in
c the variables chosen.  struct_had.out is the file to which the results 
c of the calculation are written.
cws      open (unit = 11, file ='/userdisk/zeusqcd/develop/struct.dat',
cws     + status = 'old')
cws      open (unit = 12, file ='struct.out', status = 'unknown')
c      open (unit = 13, file ='struct.plo', status = 'unknown')
cws      rewind 11
c      rewind 12
c      rewind 13
c      rewind 14
c      rewind 15
c      rewind 16
c      rewind 17
c      rewind 18
c
c      print*,' choice of distribution function: ipdf = ',ipdf
c      print*,' choice of f_l (0) f_t (1) f_2 (2): ifch = ',ifch
c      print*,' include the born terms (1) not (0): nborn = ',nborn
c      print*,' include the gluon corrections (1) not (0): ngcorr = '
c     #     ,ngcorr
c      print*,' include the quark corrections (1) not (0): nqcorr = '
c     #     ,nqcorr
c      print*,' factor multiplying Q^2 to give MF scale = ',qfac
c      print*,' factor multiplying m^2 to give MF scale = ',mfac
c      print*,' the mass of the produced heavy quark: mq = ',mq
c      print*,' the virtuality of the photon: mq2 = ', mq2fix
c      print*,' the number of light flavors: nlf = ', nlf
c      print*,' ichoose = (0) loop (1) one pt (2) mf dep = ',ichoose
c      print*,' the bjorken x fixed value: xbjfix = ', xbjfix
c      print*,' the produced quark is ',q
c mass, charge parameters
cws      if (q .eq. 'charm') then
cws         eh = 2d0/3d0
cws      elseif (q .eq. 'bottom') then
cws         eh = -1d0/3d0
cws      elseif (q .eq. 'top') then
cws         eh = 2d0/3d0
cws      endif
cws      eh2 = eh*eh
cws      m2 = mq*mq

cws -----------------------
      ichoose = 1

      if (ichoose .eq. 0) then
c     loop over chosen values of xbj and mq2
         mq2l(1) = 10d0
         mq2l(2) = 50d0
         mq2l(3) = 100d0
         mq2l(4) = 500d0
         do 20 iloop1 = 1,4
            mq2 = mq2l(iloop1)
            do 10 iloop = 1,36
               if (iloop .le. 9) then
                  xbj = dble(iloop)*10**(-4d0)
               elseif (iloop .ge. 10 .and. iloop .le. 18) then
                  xbj = dble(iloop-9)*10**(-3d0)
               elseif (iloop .ge. 19 .and. iloop .le. 27) then
                  xbj = dble(iloop-18)*10**(-2d0)
               elseif (iloop .ge. 28 .and. iloop .le. 36) then
                  xbj = dble(iloop-27)*10**(-1d0)
               endif
c     choose the mass factorization scale
               scale = dsqrt(qfac*mq2 + mfac*m2)
c     print*,'xbj = ',xbj
               if (mq2*(1.d0 - xbj)/xbj .le. 4.d0*m2) then 
                  print*,'variables out of range'
                  goto 10
               endif
c     call daind1
               zmax = mq2/(mq2 + 4d0*m2)
               zmin = xbj
               eps = 1d-6
               max = 10000
               xfhad = daind1(zmin,zmax,ffm,eps,2,max,kount,est)
               fhad = alphas*mq2/m2/4.d0/pi/pi*xfhad
cws               write(12,100) xbj, mq2, scale, fhad
cws               write(13,95) xbj, fhad
c        print*,xbj, mq2, scale, fhad
 10         continue
 20      continue
      elseif (ichoose .eq. 1) then
         xbj = xbjfix
         mq2 = mq2fix
         scale = dsqrt(qfac*mq2 + mfac*m2)
c     print*,'xbj = ',xbj
         if (mq2*(1.d0 - xbj)/xbj .le. 4.d0*m2) then 
            print*,'jacksmith: variables out of range'
            print*,mq2,xbj,m2,(mq2*(1.d0 - xbj)/xbj - 4.d0*m2)
            fhad = 0.
            goto 999
         endif
c     here we call the integration routine daind1 to integrate the function
c     F.  the first two arguments are the limits of the integration, then
c     the function to be integrated, the relative accuracy eps, because of
c     the two in the next argument, the maximum number of points to be taken
c     for the integral, the number adn the ESTimated error.
         zmax = mq2/(mq2 + 4d0*m2)
         zmin = xbj
         eps = 1d-6
         max = 10000
         xfhad = daind1(zmin,zmax,ffm,eps,2,max,kount,est)
         fhad = alphas*mq2/m2/4.d0/pi/pi*xfhad
cws         write(12,100) xbj, mq2, scale, fhad
c      print*,xbj, mq2, scale, fhad
      elseif (ichoose .eq. 2) then
         mq2 = mq2fix
         xbj = xbjfix
         do 30 iloop = 1,10
c     if (iloop .le. 10) qfac = dble(iloop)/10d0
c     if (iloop .gt. 10) qfac = dble(iloop - 9)
c     scale = dsqrt(qfac*mq2)
c     if (iloop .le. 10)
c     #           scale = dsqrt(dble(iloop)*mq2/10.d0)
c     if (iloop .ge. 11) 
c     #           scale = dsqrt(dble(iloop - 9)*mq2)
            if (mq2*(1.d0 - xbj)/xbj .le. 4.d0*m2) then 
               print*,'variables out of range'
               goto 999
            endif
            if (iloop .eq. 1) scale = 3d0
            if (iloop .eq. 2) scale = 5d0
            if (iloop .eq. 3) scale = 7d0
            if (iloop .eq. 4) scale = 9d0
            if (iloop .eq. 5) scale = 10d0
            if (iloop .eq. 6) scale = 15d0
            if (iloop .eq. 7) scale = 20d0
            if (iloop .eq. 8) scale = 25d0
            if (iloop .eq. 9) scale = 30d0
            if (iloop .eq. 10) scale = 100d0
c     here we call the integration routine daind1 to integrate the function
c     F.  the first two arguments are the limits of the integration, then
c     teh function to be integrated, the relative accuracy eps, because of
c     teh two in the next argument, the maximum number of points to be taken
c     for the integral, the number adn the ESTimated error.
            zmax = mq2/(mq2 + 4d0*m2)
            zmin = xbj
            eps = 1d-6
            max = 10000
            xfhad = daind1(zmin,zmax,ffm,eps,2,max,kount,est)
            fhad = alphas*mq2/m2/4.d0/pi/pi*xfhad
cws            write(12,100) xbj, mq2, scale*scale, fhad
 30      continue
c     print*,xbj, mq2, scale*scale, fhad
      elseif (ichoose .eq. 3) then
c     loop over chosen values of mq2 with fixed x
         xbj = xbjfix

ccc ******************************************************
ccc This is the part which  does the loop in x values and Q^2 values

         do 40 nxch=1,11
         xbj=xchbin(nxch)
         mq2=1.25d0/1.11d0
         do 40 iloop = 1,6
         mq2=mq2*1.11d0
         if(iloop.eq.6) mq2=1.43*1.43

c     choose the mass factorization scale
            scale = dsqrt(qfac*mq2 + mfac*m2)
c     print*,'xbj = ',xbj
            if (mq2*(1.d0 - xbj)/xbj .le. 4.d0*m2) then 
               print*,'variables out of range'
               goto 45
            endif
c     here we call the integration routine daind1 to integrate the function
c     F.  the first two arguments are the limits of the integration, then
c     teh function to be integrated, the relative accuracy eps, because of
c     teh two in the next argument, the maximum number of points to be taken
c     for the integral, the number and the ESTimated error.
            zmax = mq2/(mq2 + 4d0*m2)
            zmin = xbj
            eps = 1d-6
            max = 10000
            xfhad = daind1(zmin,zmax,ffm,eps,2,max,kount,est)
            fhad = alphas*mq2/m2/4.d0/pi/pi*xfhad
cws            write(12,100) xbj, mq2, fhad
c            write(13,95) mq2, fhad
c     print*,xbj, mq2, scale, fhad
 45         continue
 40      continue
      endif
c 100  format (1x,4(d10.4,3x))
 100  format(f10.6,f14.5,1x,e10.4)  
 95   format (1x,d10.4,3x,d10.4)
 999  return
      end

c  this routine gives the value for alpha_s and the parton densities to be
c  convoluted with the partonic cross sections to give the structure
c  functions or cross sections.
c  parton densities must be in MS-bar scheme
c  the densities used in the subsequent routines are number densities
c  be sure that is what you use.
c  also make sure lambda_qcd is set correctly if you choose to change
c  to a different parton density.
c  also there is a cut-off built in this routine if the mass factorization
c  scale is too low.  be sure that is consistent with the density
c  you use
c  last modified 10.9.96
c  more general calls to lambda_qcd, initialized from pdfs. (hold)
c=================================================
      subroutine fnew(xbj,z)
      implicit none
      real*8 Q2
      real*8 z, xbj, scale, alphas, xz, mu2, nlf
      real*8 eh2, pi
      real*8 gluexz, sumq3, sumelq3
      real*8 sumq4, sumelq4, sumq5, sumelq5
      real*8 dum, alam4, alam5, alam6, amas4, amas5, amas6
      real*8 xmin, qini, qmax
      real*8 qcdl, qcdl2
      real*8 uxz, dxz, ubxz, dbxz, csxz, ssxz, bsxz
      real*8 alambda,flavor,qsct,qsdt,talpha
      real*8 upv,dnv,usea,dsea,str,chm,bot,glu
      real*8 ALPHAzo
c      real*8 QPDFXQ
      character*10 q
      integer iord, isch, mxflv, ncount, irt
      integer ifch, nborn, nqcorr, ngcorr, ipdf
      integer mode
      integer IFLAG
      parameter (pi = 3.1415926535897932384626433832795D0)
      external ALPHAzo
c      external QPDFXQ
      common/include/nborn, ngcorr, nqcorr, ifch
      common/coupling/scale, alphas, nlf, eh2
      common/glue/gluexz
      common/lq3/sumq3,sumelq3
      common/lq4/sumq4,sumelq4
      common/lq5/sumq5,sumelq5
      common/compdf/ipdf
      common/quark/q
      common/TRSFPAR/alambda,flavor,qsct,qsdt,iord
      data ncount/0/
      save ncount
c  initialize input from pdfs.
c      if (ncount .eq. 0) then
c         do ilam  = 4,6
c         alam(ilam) = wlamd3(ipdf,iord,idint(xnf
c     mrst input

cws --- these set by TRinit
cws      alambda=0.3
cws      flavor=3.
cws      iord=1
cws      qsdt=7.29
cws      qsct=74.
      
      talpha=2.d0*dlog(scale/alambda)
      alphas=ALPHAzo(talpha)

c      print*,'mu2 = ',mu2
c      print*,'scale = ',scale
      xz = xbj/z
c  determining parton densities (these are number densities rather than
c  momentum densities, after division by xz)
      mode=1
c      call mrschm(xz,scale,1,upv,dnv,usea,dsea,str,chm,bot,glu)
      Q2=scale*scale
      call xpdfvs(xz,Q2,upv,dnv,usea,dsea,str,chm,bot,glu)
c      print*,'xz = ',xz
c      print*,'scale = ',scale
      gluexz = glu/xz
      if (nqcorr .ne. 0) then
         uxz = (upv+usea)/xz
         dxz = (dnv+dsea)/xz
         ubxz = usea/xz
         dbxz = dsea/xz
         ssxz = str/xz
         csxz = chm/xz
         bsxz = bot/xz
c  depending upon the number of light flavors (nlf) we choose the light
c  quark densities for charm nlf = 4, bottom nlf = 5, top nlf = 6.
         if (q .eq. 'charm') then
            sumq3 = uxz + ubxz + dxz + dbxz + 2.d0*ssxz
            sumelq3 = 4.0d0/9.0d0*(uxz + ubxz) + 1.0d0/9.0d0*(dxz
     #           + dbxz) + 1.0d0/9.0d0*(2.d0*ssxz)
         elseif (q .eq. 'bottom') then
            sumq4 = uxz + ubxz + dxz + dbxz + 2.d0*ssxz + 2.d0*csxz
            sumelq4 = 4.0d0/9.0d0*(uxz + ubxz) + 1.0d0/9.0d0*(dxz
     #      + dbxz) + 1.0d0/9.0d0*(2.d0*ssxz) + 4.d0/9.d0*(2.d0*csxz)
         elseif (q .eq. 'top') then
            sumq5 = uxz + ubxz + dxz + dbxz + 2.d0*ssxz
     #           + 2.d0*csxz + 2.d0*bsxz
            sumelq5 = 4.0d0/9.0d0*(uxz + ubxz) + 1.0d0/9.0d0*(dxz
     #      + dbxz) + 1.0d0/9.0d0*(2.d0*ssxz) + 4.d0/9.d0*(2.d0*csxz)
     #           + 1.d0/9.d0*(2.d0*bsxz)
         endif
      endif
      ncount = 1
      return
      end

