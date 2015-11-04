C     ===============================================
      double precision function AsVal1(iq,nf,ithresh)
C     ===============================================

C     Get alpha_s/2pi

      implicit double precision (a-h,o-z)
      common /qcdqedord/ iordqcd,iordqed

      call setord(iordqcd)

      mf     = nf             !avoid compiler warning
      AsVal1 = 0.D0
      if(ithresh .ge. 0) then
        AsVal1 = altabn(0,iq,1,ierr)
      elseif(ithresh .eq. -1) then
        AsVal1 = altabn(0,-iq,1,ierr)
      else
        stop 'AsVal1: wrong ithresh'
      endif

      return
      end

C     ===============================================
      double precision function AsVal2(iq,nf,ithresh)
C     ===============================================

C     Get (alpha_s/2pi)^2

      implicit double precision (a-h,o-z)
      common /qcdqedord/ iordqcd,iordqed

      call setord(iordqcd)

      mf     = nf             !avoid compiler warning
      AsVal2 = 0.D0
      if(ithresh .ge. 0) then
        AsVal2 = altabn(0,iq,2,ierr)
      elseif(ithresh .eq. -1) then
        AsVal2 = altabn(0,-iq,2,ierr)
      else
        stop 'AsVal2: wrong ithresh'
      endif

      return
      end

C     ===============================================
      double precision function AsVal3(iq,nf,ithresh)
C     ===============================================

C     Get (alpha_s/2pi)^3

      implicit double precision (a-h,o-z)
      common /qcdqedord/ iordqcd,iordqed

      call setord(iordqcd)

      mf     = nf             !avoid compiler warning
      AsVal3 = 0.D0
      if(ithresh .ge. 0) then
        AsVal3 = altabn(0,iq,3,ierr)
      elseif(ithresh .eq. -1) then
        AsVal3 = altabn(0,-iq,3,ierr)
      else
        stop 'AsVal3: wrong ithresh'
      endif

      return
      end

C     ===============================================
      double precision function AemVal1(iq,nf,ithresh)
C     ===============================================

C     Get alpha_em/2pi

      implicit double precision (a-h,o-z)
      common /aem/ aem0,rem20,q2b,q2t

      mthresh = ithresh          !avoid compiler warning

c      AemVal1 = 0d0
c      return

      aemb = aem0
     $  /(1.D0-4.D0/9.D0*(4.D0*int(4/2)+int((4+1)/2))
     $  *aem0*dlog(q2b/rem20))
      aemt = aemb
     $  /(1.D0-4.D0/9.D0*(4.D0*int(5/2)+int((5+1)/2))
     $  *aemb*dlog(q2t/q2b))
      if (nf.le.4) then
        AemVal1 = aem0
     $  /(1.D0-4.D0/9.D0*(4.D0*int(nf/2)+int((nf+1)/2))
     $  *aem0*dlog(qfrmiq(iq)/rem20))
      else if (nf.le.5) then
        AemVal1 = aemb
     $  /(1.D0-4.D0/9.D0*(4.D0*int(nf/2)+int((nf+1)/2))
     $  *aemb*dlog(qfrmiq(iq)/q2b))
      else if (nf.le.6) then
        AemVal1 = aemt
     $  /(1.D0-4.D0/9.D0*(4.D0*int(nf/2)+int((nf+1)/2))
     $  *aemt*dlog(qfrmiq(iq)/q2t))
      endif
      return
      end
