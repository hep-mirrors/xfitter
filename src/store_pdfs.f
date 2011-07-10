
* --------------------------------------------------------------
      subroutine store_pdfs(base)
* --------------------------------------------------------------

      implicit none

      include 'steering.inc'
      include 'pdfparam.inc'

      integer i,ix,idx,iq2,iflag
      double precision q2,x,gval,sing,umin,dmin
*new jf
      double precision QPDFXQ,cplus,splus,bplus,uplus,dplus,U,D,sea,DbmUb
      double precision d_Ubar,d_Dbar,u_sea,d_sea,str,chm,bot
      double precision totstr,totDbar,totcha,totUbar,afs,afc,xbelow,delx
      double precision totusea,totdsea,afs_ud

*     some variables needed for getting structure functions
*     from QCDNUM+HF schemes
      double precision QSTFXQ, f2, fl, xf3, f2em, flem
      double precision f2qcdnum,flqcdnum,xf3qcdnum,f2emqcdnum,flemqcdnum
      integer mode
      double precision rt_f2p,rt_flp,rt_f1p,
     +     rt_rp,rt_f2n,rt_fln,
     +     rt_f1n,rt_rn,rt_f2c,
     +     rt_flc,rt_f1c,rt_f2b,
     +     rt_flb,rt_f1b

      double precision xnu, xrho, Qsimple      
      double precision F123(3)

      character*25 name
      character*25 h1name
      character*25 base
      character*25 fsfc
      character*25 namefsfc

* base = pdfs_q2val_
      character tag(NBANDS)*2
      data (tag(i),i=1,NBANDS) /'01','02','03','04','05',
     $     '06','07','08','09','10',
     $     '11','12','13','14','15','16','17','18','19','20',
     $     '21','22','23','24','25','26','27','28','29','30',
     $     '31','32','33','34','35','36','37','38','39','40' /

      include 'thresholds.inc'

      real*4 glun00h100,ubar00h100,uval00h100,dbar00h100
      real*4 dval00h100,cbar00h100,sbar00h100,bbar00h100

      double precision nf

      double precision pdf
      dimension pdf(-6:6)


      double precision pdfl,asfunc,alphas(0:160),q2valpdf(0:160)
     $     ,xvalpdf(0:160)
      dimension pdfl(-6:6)
      integer iq,jx,j,ierr,nfv

      write(6,*) '--------- in store-pdfs -------'

      idx = index(base,' ')-1
      namefsfc = 'output/fsfc'//'.txt'
      print *,namefsfc
      open (83,file=namefsfc)
      write(83,*) 'iparam,q2,strange,Dbar,fs,charm,Ubar,fc'
      do 999 iq2=1,NBANDS
cv         print*,'voicaaaaa', iq2, q2val(iq2), nbands
      q2 = Q2VAL(iq2)
      if (q2.lt.0) goto 999
      nf = 5.
      if (q2.lt.qb) nf=4.
      if (q2.lt.qc) nf=3.
      
      if (idx.gt.0) then
        name =base(1:idx)//tag(iq2)//'.txt'
        h1name = base(1:idx)//tag(iq2)//'.txt'
      else
        name =base//tag(iq2)//'.txt'
        h1name = base//tag(iq2)//'.txt'
      endif
        open(81,file=name)
c        open(82,file=h1name)

        totstr=  0.d0
        totDbar= 0.d0
        totcha = 0.d0
        totUbar = 0.d0
        totusea = 0.d0
        totdsea = 0.d0
        afs_ud = 0.d0
        afs = 0.d0
        afc = 0.d0

      delx = 0.d0
      x = 1.d0

      do ix=1,outnx
         xbelow = x
         x = log(outxrange(1)) + (log(outxrange(2))-log(outxrange(1)))
     $        * (ix-1.d0)/(outnx-1.d0)
         x = dexp(x)
* for integral calculation
         if(ix.gt.1) then
            delx = x - xbelow
         endif
         call fpdfxq(1,x,q2,pdf ,0)


         gval=pdf(0)


         U=pdf(2)+pdf(-4)
         if (q2.gt.qb) then
            D=pdf(1)+pdf(-3)+pdf(-5)
         else
            D=pdf(1)+pdf(-3)
         endif
         
         umin=pdf(2)-pdf(-2)
         dmin=pdf(1)-pdf(-1)
         d_Ubar=pdf(-2)+pdf(-4)
     
         if (q2.gt.qb) then
            d_Dbar=pdf(-1)+pdf(-3)+pdf(-5)
         else
            d_Dbar=pdf(-1)+pdf(-3)
         endif

         sea=d_Ubar+d_Dbar

*      DbmUb=d_Dbar-d_Ubar
         u_sea=pdf(-2)
         d_sea=pdf(-1)
         str=pdf(-3)

         chm = 0.0d0
         if (q2.gt.qc) then
            chm=pdf(-4)
         endif
      
         bot = 0.d0
         if (q2.gt.qb) then
            bot=pdf(-5)
         endif



* integral calculation to estimate fs and fc
      totstr = totstr + str*delx
      totDbar = totDbar + d_Dbar*delx
      totcha = totcha + chm*delx
      totUbar = totUbar +d_Ubar*delx
      totusea = totusea + u_sea*delx
      totdsea = totdsea + d_sea*delx
      if(q2.le.6.) then
      if(ix.eq.1) then
       write(83,*) 'ix,x,q2,totstr,totDbar,totcha,totUbar,
     +totusea,totdsea'
      endif
      write(83,*) ix,x,q2,totstr,totDbar,totcha,totUbar
     +,totusea,totdsea
      endif

      if (outform.eq.0) then
*     old "standard" writing - no structure functions!
        write(81,810)
     +     x,gval,U,D,d_Ubar,d_Dbar,umin,dmin,sea,u_sea,d_sea,str,chm,bot
      endif
 810     format(14(2x,G12.6))
 811     format(I3,2x,23(2x,G12.6))


      enddo
* special
      if(q2.le.6.) then
      afs = totstr/totDbar
      afc = totcha/totUbar
      afs_ud = totstr/(totusea+totdsea)
      write(83,*) iparam,q2,totstr,totDbar,afs,totcha,totUbar,afc
      write(83,*) totusea,totdsea,afs_ud
      endif
* end special
      close(81)
c      close(82)
      
999   continue


      close(83)


cv store for LHAPDF
cv  PDFs are Glue Uval Dval Ubar Dbar Str Chrm  Bot 

      DO Iq=0,160
         Q2=10**(8.30103/160D0*Iq )
         q2valpdf(iq) = q2
         alphas(iq)=ASFUNC(q2,nfv,ierr)
      enddo


      DO jx=0,160
         IF(Jx.LE.80)THEN
            X=10**(6D0/120D0*Jx-6D0)
         ELSE
            X=10**(2D0/80D0*(Jx-80)-201D-2)
         ENDIF
         xvalpdf(jx) = x
      enddo

C Prepare LHAPDF output

      open (76,file='output/lhapdf.block.txt',status='unknown')

      do iq2=1,23
         write (76,'(7E12.4)') (q2valpdf((iq2-1)*7+j),j=0,6)
      enddo

      do jx=1,23
         write (76,'(7E12.4)') (xvalpdf((jx-1)*7+j),j=0,6)
      enddo

      do iq2=1,23
         write (76,'(7E12.4)') (alphas((iq2-1)*7+j),j=0,6)
      enddo


      DO Iq=0,160
         Q2=10**(8.30103/160D0*Iq )
cv         grid(162+iq)=q2

         DO jx=0,160
            IF(Jx.LE.80)THEN
               X=10**(6D0/120D0*Jx-6D0)
            ELSE
               X=10**(2D0/80D0*(Jx-80)-201D-2)
            ENDIF
cv       grid(1+jx)=x
            call fpdfxq(1,x,q2,pdfl ,0)
            write(76,666) PDFl(0), PDFl(2)-PDFl(-2), PDFl(1)-PDFl(-1),
     $           PDFl(-2), PDFl(-1),
     $           PDFl(3), PDFl(4), PDFl(5)


         enddo
      ENDDO

      close (76)

 666   format(8(2x,G12.6))




      return
      end

