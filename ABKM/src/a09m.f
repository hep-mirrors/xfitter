      subroutine a09m(xb,q2,PDFS,DPDFS,NFLV,IPAR,IORD)
c
c     This is a code for the 3-, 4-, and 5-flavour NNLO nucleon parton 
c     distributions generated from the fit of [PRD 81, 014032 (2010)]
c     using the matching conditions of [Eur.Phys.J.C1:301-320,1998], 
c     with account of the PDF uncertainties. 
c
c     The Q**2 range is 0.8d0 < Q**2 < 2d8, the x range is 1d-7 < x < 1d0. 
c
c  Input parameters:
c    XB is the parton momentum fraction 
c    Q2 is the factorization scale
c    NFLV selects the PDFs sets with different number of flavours (3, 4, or 5)
c    IPAR controls output of the PDFs and \alpha_s uncertaintis (see 
c         description of the output parameters below)
c    IORD defines the order of PDFs (2: NNLO, 1: NLO)

c  Output parameters:
c     The array PDFS contains fitted values of the strong coupling constant 
c     and the parton distributions at given x and Q:
c        PDFS(0) -- \alpha_s
c        PDFS(1) -- valence u-quarks 
c        PDFS(2) -- valence d-quarks
c        PDFS(3) -- gluons 
c        PDFS(4) -- sea u-quarks 
c        PDFS(5) -- s-quarks 
c        PDFS(6) -- sea d-quarks 
c        PDFS(7) -- c-quarks
c        PDFS(8) -- b-quarks
c     Output array DPDFS(0:8,NVAR) contains derivatives of \alpha_s and
c     the PDFs on the parameters corresponding to the independent 
c     sources of the uncertainties, NVAR is the number of these sources, 
c     equal to 27 in the current version. The input parameter IPAR is used to 
c     optimize performance fo the code. If IPAR=0, no 
c     uncertainties are returned in DPDFS; if 0<IPAR<=NVAR, only 
c     the uncertainty due to the IPAR-th source is returned; if IPAR<0
c     all uncertainties for the sources from 1 to NVAR are returned. 
c     Using derivatives returned in DPDFS one can take into account the 
c     correlations between different PDFs and between PDFs and \alpha_s. 
c     All derivatives are transformed to the orthonormal basis of eigenvectors 
c     of the parameters error matrix therefore variation of the PDFs by 
c     the values of DPDFS is performed independently. For example, 
c     after the call of A09 with IPAR=-1 the dispersion of the i-th PDF can 
c     be stored in DELPDF using the code 
c
c-----------------
c          DELPDF=0.
c          do k=1,nvar
c            DELPDF=DELPDF+dpdfs(i,k)**2
c          end do
c-----------------
c     and its random value can be stored in RPDF using the code 
c-----------------
c          RPDF=pdfs(i)          
c          do k=1,nvar
c            s=0.
c            do l=1,96
c              s=s+(2*rndm(xxx)-1)/sqrt(32.)
c            end do
c            RPDF=RPDF+s*dpdfs(i,k)
c          end do
c-----------------
c         Comments: Sergey.Alekhin@ihep.ru                      
c                                                               
c     Initial version: Jul 2009      
c     Modification of Jan 2010: The NLO PDFs added

      implicit none

      integer nxb,nq,np,nvar
      parameter(nxb=99,nq=20,np=8,nvar=27)

      integer k,i,n,m,kx,nxbb
      integer NPDF,NFLV,IORD

      real*8 f(nxb,nq+1,0:np),xx(nxb)
      real*8 fsp(nxb),bs(nxb),cs(nxb),ds(nxb)
      real*8 bsp(nxb,nq+1,0:np),csp(nxb,nq+1,0:np),dsp(nxb,nq+1,0:np)
      real*8 bspd(nvar,nxb,nq+1,0:np),cspd(nvar,nxb,nq+1,0:np)
     ,      ,dspd(nvar,nxb,nq+1,0:np)
      real*8 pdfs(0:np),dpdfs(0:np,nvar)
      real*8 df(nvar,0:np,nxb,nq+1)
      real*8 x,qsq,dels,delx,x1,delx1,xlog1,xd,b,aa,ss,f0,fp,fm
      real*8 xb,q2,df0,dfp,dfm

      character pdford*1,suffix*4
      dimension pdford(3),suffix(2)

      real*8 xmin,xmax,qsqmin,qsqmax
      integer nflvs,npar1,npar2,ipar,iords
      integer lnblnk

c I/O channel to read the data
      integer nport
      character locdir0*128,locdir*128
      data nport/1/

      data pdford/'3','4','5'/
      data suffix/'_nlo','    '/

      data xmin,xmax,qsqmin,qsqmax/1d-7,1d0,0.8d0,2d8/
      data NFLVS, IORDS /-1, -1/

      save nflvs,f,df,dels,delx,x1,delx1,xlog1,nxbb,xx,iords

c put in your local address of the PDF files in LOCDIR0 
c      data locdir0 /'/usr/users/alekhin/grids/'/
c or in the system variable GRIDS
      CALL GETENV( 'GRIDS', locdir0 ) 
      locdir=locdir0(:LNBLNK(locdir0))//'pdfs/abm11/'

      npdf=nflv+3

      if (ipar.gt.nvar) print *,'Wrong call of the PDF uncertainties' 

      if (ipar.eq.0) then 
        npar1=0
        npar2=0
      end if
      if (ipar.ge.0) then 
        npar1=ipar
        npar2=ipar
      end if
      if (ipar.lt.0) then 
        npar1=1
        npar2=nvar
      end if

      if (nflvs.ne.nflv.or.iords.ne.iord) then 
        nflvs=nflv
        iords=iord
      else 
        goto 10
      end if

      dels=(dlog(dlog(qsqmax/0.04d0))-
     +      dlog(dlog(qsqmin/0.04d0)))/dble(nq-1)

      nxbb=nxb/2
      x1=0.3d0
      xlog1=dlog(x1)
      delx=(dlog(x1)-dlog(xmin))/dble(nxbb-1)
      DELX1=(1.d0-x1)**2/dble(nxbb+1)

*...X GRID
      do kx=1,nxbb
        xx(kx)=dexp(dlog(xmin)+delx*dble(kx-1))
      end do
      do kx=nxbb+1,nxb-1
        xx(kx)=1.d0-dsqrt(dabs((1.d0-x1)**2-delx1*dble(kx-nxbb)))
      end do
      xx(nxb)=1d0

*...Read input tables
      print *,'***** Reading PDFs from tables *****'
      open(unit=nport,status='old'
     ,    ,file=locdir(:LNBLNK(locdir))
     /     //'/abm11.dpdfs_'//pdford(nflv-2)//suffix(iord))
      do n=1,nxb-1
        do m=1,nq
          do i=0,npdf 
            read (nport,*) (df(k,i,n,m),k=1,nvar)
          end do
        end do
      end do
      close(unit=nport)

      do k=1,nvar
        do i=0,npdf
          do m=1,nq
            if (i.ne.0) then 
              df(k,i,nxb,m)=0d0
            else 
              df(k,i,nxb,m)=df(k,i,nxb-1,m)
            end if
            do n=1,nxb
              fsp(n)=df(k,i,n,m)
            end do
            call spline (nxb,xx,fsp,bs,cs,ds)
            do n=1,nxb
              bspd(k,n,m,i)=bs(n)
              cspd(k,n,m,i)=cs(n)
              dspd(k,n,m,i)=ds(n)
            end do
          end do
        end do
      end do

      open(unit=nport,status='old',err=199
     ,    ,file=locdir(:LNBLNK(locdir))
     /     //'/abm11.pdfs_'//pdford(nflv-2)//suffix(iord))
      do n=1,nxb-1
        do m=1,nq
          read(nport,*) (f(n,m,i),i=0,npdf)
        end do
      end do
      do i=0,npdf
        do m=1,nq
          if (i.ne.0) then 
            f(nxb,m,i)=0d0
          else 
            f(nxb,m,i)=f(nxb-1,m,i)
          end if
          do n=1,nxb-1
            f(n,m,i)=f(n,m,i)
          end do
          do n=1,nxb
            fsp(n)=f(n,m,i)
          end do
          call spline (nxb,xx,fsp,bs,cs,ds)
          do n=1,nxb
            bsp(n,m,i)=bs(n)
            csp(n,m,i)=cs(n)
            dsp(n,m,i)=ds(n)
          end do
        end do
      end do
      close(unit=nport)

  10  continue

      if((q2.lt.qsqmin).or.(q2.gt.qsqmax)) then
         print 99,q2,qsqmin,qsqmax
         return
      end if
      if((xb.lt.xmin).or.(xb.gt.xmax)) then
         print 98,xb,xmin,xmax
         return
      end if
  99  format('  AGRIDS WARNING:  Q^2 VALUE IS OUT OF RANGE   ',3g12.3)
  98  format('  AGRIDS WARNING:   X  VALUE IS OUT OF RANGE   ',3g12.3)

      x=max(xb,xmin)
      x=min(xb,xmax)
      qsq=max(q2,qsqmin)
      qsq=min(q2,qsqmax)

      if (x.gt.x1) then
        xd=(1d0-x1)**2-(1d0-x)**2
        n=int(xd/delx1)+nxbb
      else
        xd=dlog(x)-xlog1
        n=nxbb+int(xd/DELX)-1
      end if
      aa=x-xx(n)

      ss=dlog(dlog(qsq/0.04d0))-dlog(dlog(qsqmin/0.04d0))
      m=int(ss/dels)+1
      b=ss/dels-dble(m)+1.d0

      do i=0,npdf
        f0=f(n,m,i) + aa*bsp(n,m,i) + aa**2*csp(n,m,i) 
     +              + aa**3*dsp(n,m,i)
        fp=f(n,m+1,i) + aa*bsp(n,m+1,i) + aa**2*csp(n,m+1,i)
     +                + aa**3*dsp(n,m+1,i)
        if (m.ge.2) then 
          fm=f(n,m-1,i) + aa*bsp(n,m-1,i) + aa**2*csp(n,m-1,i)
     +                   +aa**3*dsp(n,m-1,i)
          pdfs(i)=fm*b*(b-1d0)/2d0 + f0*(1d0-b**2) + fp*b*(b+1d0)/2d0
        else 
          pdfs(i)= f0*(1d0-b) + fp*b
        end if
        if (npar1.gt.0) then 
          do k=npar1,npar2
            df0=df(k,i,n,m) + aa*bspd(k,n,m,i) + aa**2*cspd(k,n,m,i) 
     +                      + aa**3*dspd(k,n,m,i)
            dfp=df(k,i,n,m+1)+aa*bspd(k,n,m+1,i)+aa**2*cspd(k,n,m+1,i)
     +                        + aa**3*dspd(k,n,m+1,i)
            if (m.ge.2) then 
              dfm=df(k,i,n,m-1)+aa*bspd(k,n,m-1,i)+aa**2*cspd(k,n,m-1,i)
     +                          + aa**3*dspd(k,n,m-1,i)
              dpdfs(i,k)=dfm*b*(b-1d0)/2d0 
     +                  + df0*(1d0-b**2) +dfp*b*(b+1d0)/2d0
            else 
              dpdfs(i,k) = df0*(1d0-b) + dfp*b
            end if
          end do
        end if
      end do

      return

 199  print *,'The PDFs set is inavailable'

      return
      end
