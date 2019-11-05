C--------------------------------------------------------------------------
C OZ 7.10.17 call user alpha_S and PDF routines (provided in sf_abkm_wrap.F)
#define USER
C-----------------
      SUBROUTINE PDFFILLGRID
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      include 'CONSTCOM.'
      INCLUDE 'PDFCOM.'

      real*8 fsp(nxtot),bs(nxtot),cs(nxtot),ds(nxtot),xx(nxtot)
      real*8 pdfs(0:8),dpdfs(0:8,27)

      character name*80
      real*8 f(-6:6)

!  Setup of the PDF grid

      smax=log(log(q2max/0.04)/log(q2ini(0)/0.04))
      smin=-log(log(q2min/0.04)/log(q2ini(0)/0.04))
      dels1(kschemepdf)=smax/nspgrid
      dels2(kschemepdf)=smin/nsmgrid
      xlog0=log(xbmin)
      delx2=-(xlog0-xlog1)/nxmgrid

      delx3=(log(1-xbmax)-log(1-x1))/(nxpgrid-1)

! Mapping with log(1-x) at large x
      do i=0,nxpgrid-1
        xgrid(i)=1-exp(log(1-xbmax)+delx3*(i-nxpgrid+1))
      end do
      xgrid(nxpgrid)=1.
! Mapping with log(x) at small x
      do i=-nxmgrid,0
        xgrid(i)=exp(xlog1+delx2*i)
      end do

      DO I=0,nspgrid
        SGRID(I,kschemepdf)=DELS1(kschemepdf)*I
      END DO

      DO I=-nsmgrid,0
        SGRID(I,kschemepdf)=DELS2(kschemepdf)*I
      END DO

!  Here we fill the PDFs grid to be used elsewhere

!  Select the PDF set from the LHAPDF library
      name=' '
#ifdef ABMP16
      if (kschemepdf.eq.2) then 
        if (kordpdf.eq.2) name='ABMP16_5_nnlo.LHgrid'   
      end if
      if (kschemepdf.eq.1) then 
        if (kordpdf.eq.2) name='ABMP16_4_nnlo.LHgrid'   
      end if
      if (kschemepdf.eq.0) then 
        if (kordpdf.eq.2) name='ABMP16_3_nnlo.LHgrid'   
      end if
#endif 
#ifdef ABM12
      if (kschemepdf.eq.2) then 
        if (kordpdf.eq.2) name='abm12lhc_5_nnlo.LHgrid'    
        if (kordpdf.eq.1) name='abm12lhc_5_nlo.LHgrid'    
      end if
      if (kschemepdf.eq.1) then 
        if (kordpdf.eq.2) name='abm12lhc_4_nnlo.LHgrid'    
      end if
      if (kschemepdf.eq.0) then 
        if (kordpdf.eq.2) name='abm12lhc_3_nnlo.LHgrid'    
      end if
#endif 
#ifdef ABM11
      if (kschemepdf.eq.2) then 
        if (kordpdf.eq.2) name='abm11_5n_nnlo.LHgrid'    
        if (kordpdf.eq.1) name='abm11_5n_nlo.LHgrid'    
      end if
      if (kschemepdf.eq.1) then 
        if (kordpdf.eq.2) name='abm11_4n_nnlo.LHgrid'    
        if (kordpdf.eq.1) name='abm11_4n_nlo.LHgrid'    
      end if
      if (kschemepdf.eq.0) then 
        if (kordpdf.eq.2) name='abm11_3n_nnlo.LHgrid'    
        if (kordpdf.eq.1) name='abm11_3n_nlo.LHgrid'    
      end if
#endif 
#ifdef ABKM09
      if (kschemepdf.eq.2) then 
        if (kordpdf.eq.2) name='abkm09_5_nnlo.LHgrid'    
        if (kordpdf.eq.1) name='abkm09_5_nlo.LHgrid'    
      end if
      if (kschemepdf.eq.1) then 
        if (kordpdf.eq.2) name='abkm09_4_nnlo.LHgrid'    
        if (kordpdf.eq.1) name='abkm09_4_nlo.LHgrid'    
      end if
      if (kschemepdf.eq.0) then 
        if (kordpdf.eq.2) name='abkm09_3_nnlo.LHgrid'    
        if (kordpdf.eq.1) name='abkm09_3_nlo.LHgrid'    
      end if
#endif 
#ifdef JR09
      if (kordpdf.eq.2) then
        if (kschemepdf.eq.2) name='JR09VFnnloE.LHgrid'    
        if (kschemepdf.eq.1) name='JR09VFnnloE.LHgrid'    
        if (kschemepdf.eq.0) name='JR09FFnnloE.LHgrid'    
      end if
#endif
#ifdef CT14
      if (kordpdf.eq.1) then 
        if (kschemepdf.eq.2) name='CT14nlo.LHgrid'    
       end if
      if (kordpdf.eq.2) then 
        if (kschemepdf.eq.2) name='CT14nnlo.LHgrid'    
      end if
#endif
#ifdef CT10
      if (kordpdf.eq.1) then 
        if (kschemepdf.eq.2) name='CT10nlo.LHgrid'    
        if (kschemepdf.eq.1) name='CT10nlo_nf4.LHgrid'   
        if (kschemepdf.eq.0) name='CT10nlo_nf3.LHgrid'    
       end if
      if (kordpdf.eq.2) then 
        if (kschemepdf.eq.2) name='CT10nnlo.LHgrid'    
      end if
#endif
#ifdef NN30
      if (kordpdf.eq.1) then
        if (kschemepdf.eq.2) name='NNPDF30_nnlo_as_0118.LHgrid'
      end if
      if (kordpdf.eq.2) then
        if (kschemepdf.eq.2) name='NNPDF30_nnlo_as_0118.LHgrid'
      end if
#endif
#ifdef NN23
      if (kordpdf.eq.1) then
        if (kschemepdf.eq.2) name='NNPDF23_nlo_FFN_NF5_as_0119.LHgrid'
        if (kschemepdf.eq.1) name='NNPDF23_nlo_FFN_NF4_as_0119.LHgrid'
        if (kschemepdf.eq.0) name='NNPDF23_nlo_as_0119.LHgrid'
      end if
      if (kordpdf.eq.2) then
        if (kschemepdf.eq.2) name='NNPDF23_nnlo_FFN_NF5_as_0119.LHgrid'
        if (kschemepdf.eq.1) name='NNPDF23_nnlo_FFN_NF4_as_0119.LHgrid'
        if (kschemepdf.eq.0) name='NNPDF23_nnlo_as_0119.LHgrid'
      end if
#endif
#ifdef NN21
      if (kordpdf.eq.1) then
        if (kschemepdf.eq.2) name='NNPDF21_FFN_NF5_100.LHgrid'
        if (kschemepdf.eq.1) name='NNPDF21_FFN_NF4_100.LHgrid'
        if (kschemepdf.eq.0) name='NNPDF21_FFN_NF3_100.LHgrid'
      end if
      if (kordpdf.eq.2) then
        name='NNPDF21_nnlo_100.LHgrid'
      end if
#endif
#ifdef MSTW08
      if (kschemepdf.eq.2) then 
        if (kordpdf.eq.2) name='MSTW2008nnlo68cl.LHgrid'    
        if (kordpdf.eq.1) name='MSTW2008nlo68cl.LHgrid'    
      end if
      if (kschemepdf.eq.1) then 
        if (kordpdf.eq.2) name='MSTW2008nnlo68cl_nf4.LHgrid'    
        if (kordpdf.eq.1) name='MSTW2008nlo68cl_nf4.LHgrid'    
      end if
      if (kschemepdf.eq.0) then 
        if (kordpdf.eq.2) name='MSTW2008nnlo68cl_nf3.LHgrid'    
        if (kordpdf.eq.1) name='MSTW2008nlo68cl_nf3.LHgrid'    
      end if
#endif
#ifdef MMHT14
      if (kschemepdf.eq.2) then 
        if (kordpdf.eq.1) name='MMHT2014nlo68cl.LHgrid'!//char(0)    
        if (kordpdf.eq.2) name='MMHT2014nnlo68cl.LHgrid'!//char(0)     
      end if
#endif
#ifdef CJ15
      if (kschemepdf.eq.2) then 
        if (kordpdf.eq.1) name='CJ15nlo.LHgrid'!//char(0)    
      end if
#endif
#ifdef HERAPDF1
      if (kschemepdf.eq.2) then 
        if (kordpdf.eq.1) name='HERAPDF10_EIG.LHgrid'
      end if
#endif
!  The NNLO toy LH PDFs 
#ifdef TOYLH
      CALL GETENV( 'GRIDS', name ) 
      name=name(:LNBLNK(name))//'pdfs/toylh/toyLH_NNLO_nf3.grid'
      call InitNNPDFwrap(name,NREP)
      call InitPDF(kpdfset)
#endif

! Set the PDF set for the case of LHAPDF interface  
#if defined(MSTW08)||defined(MMHT14)||defined(JR09)||defined(HERAPDF1)||defined(ABKM09)||defined(ABM11)||defined(ABM12)||defined(ABMP16)||defined(CT10)||defined(CT14)||defined(CJ15)||defined(NN21)||defined(NN23)||defined(NN30)  
      if (name.eq.' ') then 
        print *,'WARNING: the PDF set requested is not available'
        stop
      else
        call initPDFSetByName(name)
        call InitPDF(kpdfset)
c        do ii=1,10
c        kpdfset=ii
c        call initlha6(name,kpdfset)
c        end do
c        stop
      end if
#endif

! ...then fill the PDF grid....

      DO is=-nsmgrid,nspgrid
        q2=exp(exp(sgrid(is,kschemepdf))
     *        *log(q2ini(kschemepdf)/0.04))*0.04
        DO IX=-nxmgrid,nxpgrid
          xb=xgrid(ix)
#if defined(MSTW08)||defined(MMHT14)||defined(JR09)||defined(HERAPDF1)||defined(ABKM09)||defined(ABM11)||defined(ABM12)||defined(ABMP16)||defined(CT10)||defined(CT14)||defined(CJ15)||defined(NN21)||defined(NN23)||defined(NN30)||defined(TOYLH)  
! LHAPDF interface
!        call NumberPDF(kkk)
!        print *,kkk
#ifdef TOYLH
! A special case for the toy LH PDFs
           call NNPDFINTevolveLHA(xb,sqrt(q2),f,1)
           q20alphas=rmass(8)**2
           alphas0=0.35d0
           an=alphas_ffn(q2)
#else
           an=alphasPDF(sqrt(q2))
c           qq=sqrt(q2)
c           an=alpslha6(qq)      !    C++ interface to LHAPDF
#endif
c  \alpha_s  
           Y(kschemepdf,0,IX,is)=an
c  PDFs
           call evolvePDF(xb,sqrt(q2),f)
c           do ipdf=-5,5
c             f(ipdf) = xqglha6(ipdf, xb, qq)  !    C++ interface to LHAPDF
c           end do
c  gluon
           Y(kschemepdf,1,IX,is)=f(0)
c  d
           Y(kschemepdf,2,IX,is)=f(1)
c  dbar
           Y(kschemepdf,3,IX,is)=f(-1)
c  u
           Y(kschemepdf,4,IX,is)=f(2)
c  ubar
           Y(kschemepdf,5,IX,is)=f(-2)
c  s,sbar
           Y(kschemepdf,6,IX,is)=f(3)
           Y(kschemepdf,7,IX,is)=f(-3)
           Y(kschemepdf,8,IX,is)=0.
           Y(kschemepdf,9,IX,is)=0.
           Y(kschemepdf,10,IX,is)=0.
           Y(kschemepdf,11,IX,is)=0.
c  c,cbar
           if (kschemepdf.ge.1) then 
             Y(kschemepdf,8,IX,is)=f(-4)
             Y(kschemepdf,9,IX,is)=f(-4)
             Y(kschemepdf,10,IX,is)=0.
             Y(kschemepdf,11,IX,is)=0.
           end if
c  b,bbar
           if (kschemepdf.ge.2) then 
             Y(kschemepdf,10,IX,is)=f(-5)
             Y(kschemepdf,11,IX,is)=f(-5)
           end if
#endif
#ifdef USER
! A user PDF input
!  \alpha_s  
           Y(kschemepdf,0,IX,is)=useralphas(q2
     ,           ,kschemepdf,kordpdf,kpdfset)
!  gluon
           Y(kschemepdf,1,IX,is)=userpdfs(xb,q2,0
     ,           ,kschemepdf,kordpdf,kpdfset)
!  d
           Y(kschemepdf,2,IX,is)=userpdfs(xb,q2,1
     ,           ,kschemepdf,kordpdf,kpdfset)
!  dbar
           Y(kschemepdf,3,IX,is)=userpdfs(xb,q2,-1
     ,           ,kschemepdf,kordpdf,kpdfset)
!  u
           Y(kschemepdf,4,IX,is)=userpdfs(xb,q2,2
     ,           ,kschemepdf,kordpdf,kpdfset)
!  ubar
           Y(kschemepdf,5,IX,is)=userpdfs(xb,q2,-2
     ,           ,kschemepdf,kordpdf,kpdfset)
!  s,sbar
           Y(kschemepdf,6,IX,is)=userpdfs(xb,q2,3
     ,           ,kschemepdf,kordpdf,kpdfset)
           Y(kschemepdf,7,IX,is)=userpdfs(xb,q2,-3
     ,           ,kschemepdf,kordpdf,kpdfset)
           Y(kschemepdf,8,IX,is)=0.
           Y(kschemepdf,9,IX,is)=0.
           Y(kschemepdf,10,IX,is)=0.
           Y(kschemepdf,11,IX,is)=0.
!  c,cbar
           if (kschemepdf.ge.1) then 
             Y(kschemepdf,8,IX,is)=userpdfs(xb,q2,-4
     ,           ,kschemepdf,kordpdf,kpdfset)
             Y(kschemepdf,9,IX,is)=userpdfs(xb,q2,-4
     ,           ,kschemepdf,kordpdf,kpdfset)
             Y(kschemepdf,10,IX,is)=0.
             Y(kschemepdf,11,IX,is)=0.
           end if
!  b,bbar
           if (kschemepdf.ge.2) then 
             Y(kschemepdf,10,IX,is)=userpdfs(xb,q2,-5
     ,           ,kschemepdf,kordpdf,kpdfset)
             Y(kschemepdf,11,IX,is)=userpdfs(xb,q2,-5
     ,          ,kschemepdf,kordpdf,kpdfset)
           end if
#endif
        end do

! ...and finally, construct the spline interpolation coefficients for the grid.
        do iq=0,7+kschemepdf*2
          do ix=-nxmgrid,nxpgrid
            xx(ix+nxmgrid+1)=xgrid(ix)
            fsp(ix+nxmgrid+1)=y(kschemepdf,iq,ix,is)
          end do
          call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
          do ix=-nxmgrid,nxpgrid
              bcoeff(kschemepdf,iq,ix,is)=bs(ix+nxmgrid+1)
              ccoeff(kschemepdf,iq,ix,is)=cs(ix+nxmgrid+1)
              dcoeff(kschemepdf,iq,ix,is)=ds(ix+nxmgrid+1)
          end do
        end do
      end do

      RETURN
      END
C------------------
      REAL*8 FUNCTION XQG(IQ,xx,Q2,kp)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'

      DATA Q2SAV/-1./
      data xbsav/-1./
      save xbsav,q2sav,is,q,ix,p

      xb=min(xx,xbmax)

      if(q2.gt.q2max.or.q2.lt.q2min)
     -            print*,'q2 out of range in xqg:',q2,q2min,q2max
      s=log(log(q2/0.04)/log(q2ini(kp)/0.04))
      if (s.ge.0.) then
        is=int(s/dels1(kp))
        dss=(s-sgrid(is,kp))/dels1(kp)
      else
        is=int(s/dels2(kp))-1
        dss=(s-sgrid(is,kp))/dels2(kp)
      end if
      q=(s-sgrid(is,kp))/(sgrid(is+1,kp)-sgrid(is,kp))

      if(xb.ne.xbsav) then 
        if (xb.lt.xbmin)
     -              print *,'xb out of range in xqg:',xb

        do ix=-nxmgrid,nxpgrid-1
          if (xb.lt.xgrid(ix+1)) goto 300
        end do
 300  xbsav=xb
      end if

      aa=xb-xgrid(ix)

      if(is.ne.-nsmgrid.and.ix.ne.-nxmgrid) then 
        f0=y(kp,iq,ix,is)+aa*bcoeff(kp,iq,ix,is)
     +    +aa**2*ccoeff(kp,iq,ix,is)+aa**3*dcoeff(kp,iq,ix,is)
        fp=y(kp,iq,ix,is+1)+aa*bcoeff(kp,iq,ix,is+1)
     +    +aa**2*ccoeff(kp,iq,ix,is+1)+aa**3*dcoeff(kp,iq,ix,is+1)
        fm=y(kp,iq,ix,is-1)+aa*bcoeff(kp,iq,ix,is-1)
     +    +aa**2*ccoeff(kp,iq,ix,is-1)+aa**3*dcoeff(kp,iq,ix,is-1)
        xqg=fm*q*(q-1d0)/2d0 + f0*(1d0-q**2) + fp*q*(q+1d0)/2d0
      else
        f0=y(kp,iq,ix,is)+aa*bcoeff(kp,iq,ix,is)
     +    +aa**2*ccoeff(kp,iq,ix,is)+aa**3*dcoeff(kp,iq,ix,is)
        fp=y(kp,iq,ix,is+1)+aa*bcoeff(kp,iq,ix,is+1)
     +    +aa**2*ccoeff(kp,iq,ix,is+1)+aa**3*dcoeff(kp,iq,ix,is+1)
        xqg=f0*(1d0-q) + fp*q
      end if

      RETURN
      END
!------------------
      real*8 function xqgn(iq,xb,q2,kp)
      implicit double precision (a-h,o-z)

!  The neutron PDFs

      k=iq
      if (iq.eq.2) k=4
      if (iq.eq.3) k=5
      if (iq.eq.4) k=2
      if (iq.eq.5) k=3

      xqgn=xqg(k,xb,q2,kp)

      RETURN
      END
