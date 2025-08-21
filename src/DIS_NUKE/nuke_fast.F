c************************************************
c   Nuclear correction to structure functions
c   according to the calculation by
c 
c   S.A. Kulagin and R. Petti,
c       Nucl.Phys.A 765 (2006) 126-187 [hep-ph/0412425]
c 
c   Please refer to the above paper for results 
c   based upon the present code. 
c  
c   The routine interpolates tables for the ratio:
c      R_sf = (Nuclear SF) / (nucleon SF)
c   and therefore the result should be multiplied
c   by the corresponding value of the nucleon structure
c   function for isoscalar target, including target mass 
c   corrections. Note the isovector correction is already 
c   taken into account. The target mass corrections used 
c   for the present calculation are based on a 1/Q^2 expansion 
c   of the Georgi-Politzer formula, as described in the paper 
c   cited above (NPA 765, 126-187).
c
c   xb   = Bjorken x
c   q2   = Q**2 in GeV**2
c
c   nsf  = structure function
c        = 1 2*x*F1
c        = 2 F2
c        = 3 x*F3 
c            Note: for antineutrinos (kint<0) the sum of xF3(nu)+xF3(nubar) 
c                  is provided instead to avoid numerical problems due to 
c                  negative SF at small x values  
c
c   ityp = nucleon type
c        = 0  actual average in nucleus
c        = 1  proton in nucleus
c        = -1 neutron in nucleus
c
c   kint = interaction
c        = 1 charged lepton
c        = 2 neutrino (-2 antineutrino) CC
c        = 3 neutrino (-3 antineutrino) NC
c        = 4 neutrino CC charm production (-4 antineutrino)
c        = 5 neutrino CC not-charm (-5 antineutrino)
c
c   kord = order of QCD calculation
c        =  1 LO
c        =  2 NLO
c        =  3 NNLO
c
c  ftyp  = target nucleus 
c        = 1 Carbon (A=12, Z=6)   
c        = 2 Iron (A=56, Z=26)
c        = 3 Lead (A=207, Z=82)
c        = 4 Emulsion (mixture of nuclei) 
c    
c  syst  = systematic uncertainties
c        = 0.0 no systematics (default)  
c  
c************************************************

      real*8 function nuke_fast(xb,q2,nsf,ityp,kint1,kord,ftyp,syst)
      implicit none 
*  
      integer nsf,ityp,kint,kord,kint0,kint1,ftyp,Nsys,ifirst
      integer nxb,nq,np,ntarg,i,n,m,kk,tt,nxbb,nqb,kx,kq,Debug  
      parameter(Debug=0) 
      parameter(nxb=99,nq=20,np=81,ntarg=4)
      data ifirst /0/

      real syst,Nsigmas,rvec
      real*8 xb,q2,x,qsq,a,b,aa,f0,fp,fm
      real*8 qq(nq+1),xx(nxb),dnuke_fast  
      real*8 x1,delx,delx1,ss,xd,xlog1,dels 
* 
      real*8 xmin,xmax,qsqmin,qsqmax
      data xmin,xmax,qsqmin,qsqmax/1.d-4,1.0d0,0.50d0,1d3/

      external ranlux 

      real*8 fsp(nxb),bs(nxb),cs(nxb),ds(nxb)
      real*8 f_nuke(nxb,nq+1,0:np,ntarg),b_nuke(nxb,nq+1,0:np,ntarg) 
      real*8 c_nuke(nxb,nq+1,0:np,ntarg),d_nuke(nxb,nq+1,0:np,ntarg)

      real*8 f_dnuke(nxb,nq+1,0:np,ntarg),b_dnuke(nxb,nq+1,0:np,ntarg)
      real*8 c_dnuke(nxb,nq+1,0:np,ntarg),d_dnuke(nxb,nq+1,0:np,ntarg)

c define local address of the PDF files and SF tables
      character defdir*36
c      data defdir /'/usr/users/alekhin/grids/nuc/petti13'/
c define environment variable PDFTABLES for table location  
      CHARACTER(len=256)   ENVDATA, locfile
      integer dfhf,dfhf2

      integer nport
c I/O channel to read the data
      data nport/23/

      integer ityps,kints,kords
      data ityps,kints,kords /-999,-999,-999/
      save ityps,kints,kords,xx,qq,f_nuke
      save b_nuke,c_nuke,d_nuke,b_dnuke,c_dnuke,d_dnuke
      integer flaginit
      data flaginit/0/
      save flaginit
* 
      nuke_fast=0.d0
      if ((nsf.lt.1).or.(nsf.gt.3)) return 
      if (q2.gt.2.d8) return 
* 
*...Define path to tables
*
      envdata='CHAR(0)'
      call getenv('GRIDS',envdata)

      if (envdata.eq.' ') then
         envdata=defdir
      endif
      dfhf = index(envdata,' ') - 1
* 
*...Check input flags 
* 
      if (ifirst.eq.0) then 
        ifirst = 1
        if ((ftyp.le.0).or.(ftyp.gt.4)) then 
          write(*,*) 'Type of target out of range: ftyp = ',ftyp
          return
        end if 
      end if
*  
      kint = kint1
      if(IABS(kint1).EQ.4.OR.IABS(kint1).EQ.5) kint=2*kint1/IABS(kint1)
* 
      nxbb=nxb/2
      x1=0.3d0
      xlog1=dlog(x1)
      delx=(dlog(x1)-dlog(xmin))/dble(nxbb-1)
      DELX1=(1.d0-x1)**2/dble(nxbb+1)
* 
      dels=(dlog(dlog(qsqmax/0.04d0))-
     +      dlog(dlog(qsqmin/0.04d0)))/dble(nq-1)
*  
*...Load tables and grid 
* 
C      if (kords.eq.kord) goto 10 
      if (flaginit.eq.1) goto 10 
      !print*,ityp,kint1,kord,ftyp,kords
      !call flush(6)
      ityps=ityp
      kints=kint
      kords=kord
      flaginit=1
*  
*...X GRID 
      do kx=1,nxbb
        xx(kx)=dexp(dlog(xmin)+delx*dble(kx-1))
      end do
      do kx=nxbb+1,nxb-1
        xx(kx)=1.d0-dsqrt(dabs((1.d0-x1)**2-delx1*dble(kx-nxbb)))
      end do
      xx(nxb)=1.d0
* 
*...Q2 GRID 
      do kq=1,nq
        qq(kq)=0.04d0*dexp(dexp(dlog(dlog(qsqmin/0.04d0))
     +        +dels*dble(kq-1)))
      end do
*  
      do tt=1,ntarg 
        if (tt.eq.1) then 
          print *,'***** Reading C tables *****'
*          locfile =  envdata(1:dfhf)//'/nuc/petti14/nuke_C.dat'
          locfile =  '/home/zenaiev/soft/xfitter-abmp/src/DIS_NUKE/petti14/nuke_C.dat'
        else if (tt.eq.2) then 
          print *,'***** Reading Fe tables *****'
*          locfile =  envdata(1:dfhf)//'/nuc/petti14/nuke_C.dat'
          locfile =  '/home/zenaiev/soft/xfitter-abmp/src/DIS_NUKE/petti14/nuke_Fe.dat'
        else if (tt.eq.3) then 
          print *,'***** Reading Pb tables *****'
*          locfile =  envdata(1:dfhf)//'/nuc/petti14/nuke_Pb.dat'
          locfile =  '/home/zenaiev/soft/xfitter-abmp/src/DIS_NUKE/petti14/nuke_Pb.dat'
        else if (tt.eq.4) then 
          print *,'***** Reading emulsion tables *****'
*          locfile =  envdata(1:dfhf)//'/nuc/petti14/nuke_Emulsion.dat'
          locfile =  '/home/zenaiev/soft/xfitter-abmp/src/DIS_NUKE/petti14/nuke_Emulsion.dat'
        end if
        dfhf2 = index(locfile,' ') - 1
        print *,'Opening the file ',locfile(1:dfhf2)
        open(unit=nport,status='old',err=200,file=locfile)
* 
        do n=1,nxb-1
          do m=1,nq
            read(nport,*) (f_nuke(n,m,i,tt),i=1,np)
          end do
        end do
        close(unit=nport)
* 
        do i=0,np
          do m=1,nq
            if (i.ne.0) then
              f_nuke(nxb,m,i,tt)=0d0
            else
              f_nuke(nxb,m,i,tt)=f_nuke(nxb-1,m,i,tt)
            end if
            do n=1,nxb
              fsp(n)=f_nuke(n,m,i,tt)
            end do
            call spline (nxb,xx,fsp,bs,cs,ds)
            do n=1,nxb
              b_nuke(n,m,i,tt)=bs(n)
              c_nuke(n,m,i,tt)=cs(n)
              d_nuke(n,m,i,tt)=ds(n)
            end do
          end do
        end do
* 
        if(syst.ne.0.0) then
          Nsys=INT(abs(syst)/100.)
          Nsigmas=syst/abs(syst)*(abs(syst)-float(Nsys)*100.)
          if (Nsys.eq.2) then 
            if (ftyp.eq.1) then 
              locfile =  envdata(1:dfhf)//'/nuke_dC.dat'
            else if (ftyp.eq.2) then 
              locfile =  envdata(1:dfhf)//'/nuke_dFe.dat'
            end if
            dfhf2 = index(locfile,' ') - 1
            print *,'Opening the file ',locfile(1:dfhf2)
            open(unit=nport,status='old',err=200,file=locfile)
          end if
          do n=1,nxb-1
            do m=1,nq
              read(nport,*) (f_dnuke(n,m,i,tt),i=1,np)
            end do
          end do
          close(unit=nport)
* 
          do i=0,np
            do m=1,nq
              if (i.ne.0) then
                f_dnuke(nxb,m,i,tt)=0d0
              else
                f_dnuke(nxb,m,i,tt)=f_dnuke(nxb-1,m,i,tt)
              end if
              do n=1,nxb
                fsp(n)=f_dnuke(n,m,i,tt)
              end do
              call spline (nxb,xx,fsp,bs,cs,ds)
              do n=1,nxb
                b_dnuke(n,m,i,tt)=bs(n)
                c_dnuke(n,m,i,tt)=cs(n)
                d_dnuke(n,m,i,tt)=ds(n)
              end do
            end do
          end do

        end if
      end do 
  111 format (27f12.6)
* 
  10  continue
* 
      if((q2.lt.qsqmin).or.(q2.gt.qsqmax)) then
         if (Debug.gt.0) print 99,q2,qsqmin,qsqmax
         return
      end if
      if((xb.lt.xmin).or.(xb.gt.xmax)) then
         if (Debug.gt.0) print 98,xb,xmin,xmax
         return
      end if
  99  format(' NUKE_FAST WARNING:  Q^2 VALUE IS OUT OF RANGE   ',3g12.3)
  98  format(' NUKE_FAST WARNING:   X  VALUE IS OUT OF RANGE   ',3g12.3)
* 
*...Now actual interpolation
* 
      x=max(xb,xmin)
      x=min(xb,xmax)
      qsq=max(q2,qsqmin)
      qsq=min(q2,qsqmax)
* 
      if (x.gt.x1) then
        xd=(1d0-x1)**2-(1d0-x)**2
        n=int(xd/delx1)+nxbb
      else
        xd=dlog(x)-xlog1
        n=nxbb+int(xd/DELX)-1
      end if
      aa=x-xx(n)
* 
      ss=dlog(dlog(qsq/0.04d0))-dlog(dlog(qsqmin/0.04d0))
      m=int(ss/dels)+1
      b=ss/dels-dble(m)+1.d0
* 
      kint0=kint
      if(kint.eq.-2) kint0=6
      if(kint.eq.-3) kint0=7
      if(kint.eq.-4) kint0=8
      if(kint.eq.-5) kint0=9

      kk = 27*(nsf-1) + 9*(ityp+1) + kint0

      if ((f_nuke(n,m,kk,ftyp).ne.0d0).or.
     +    (f_nuke(n+1,m,kk,ftyp).ne.0d0)) then
        f0 = f_nuke(n,m,kk,ftyp) + aa*b_nuke(n,m,kk,ftyp) 
     +     + aa**2*c_nuke(n,m,kk,ftyp) + aa**3*d_nuke(n,m,kk,ftyp)
        fp = f_nuke(n,m+1,kk,ftyp) + aa*b_nuke(n,m+1,kk,ftyp)
     _     + aa**2*c_nuke(n,m+1,kk,ftyp) + aa**3*d_nuke(n,m+1,kk,ftyp)
        if (m.ge.2) then
          fm = f_nuke(n,m-1,kk,ftyp) + aa*b_nuke(n,m-1,kk,ftyp) 
     +       + aa**2*c_nuke(n,m-1,kk,ftyp) + aa**3*d_nuke(n,m-1,kk,ftyp)
          nuke_fast = fm*b*(b-1d0)/2d0 +f0*(1d0-b**2) +fp*b*(b+1d0)/2d0
        else
          nuke_fast = f0*(1d0-b) + fp*b
        end if
      else 
        nuke_fast = 0.d0
      end if
*. 
      if((syst.ne.0.0).and.(Nsys.eq.1))then 
        if ((f_dnuke(n,m,kk,ftyp).ne.0d0).or.
     +      (f_dnuke(n+1,m,kk,ftyp).ne.0d0)) then
          f0 = f_dnuke(n,m,kk,ftyp) + aa*b_dnuke(n,m,kk,ftyp) 
     +       + aa**2*c_dnuke(n,m,kk,ftyp) + aa**3*d_dnuke(n,m,kk,ftyp)
          fp = f_dnuke(n,m+1,kk,ftyp) + aa*b_dnuke(n,m+1,kk,ftyp) 
     +       + aa**2*c_dnuke(n,m+1,kk,ftyp)+aa**3*d_dnuke(n,m+1,kk,ftyp)
          if (m.ge.2) then
            fm = f_dnuke(n,m-1,kk,ftyp) + aa*b_dnuke(n,m-1,kk,ftyp)
     +        +aa**2*c_dnuke(n,m-1,kk,ftyp)+aa**3*d_dnuke(n,m-1,kk,ftyp)
            dnuke_fast = fm*b*(b-1d0)/2d0 + f0*(1d0-b**2) 
     +                                   + fp*b*(b+1d0)/2d0
          else
            dnuke_fast = f0*(1d0-b) + fp*b
          end if
        else 
          dnuke_fast = 0.d0
        end if

        if (Nsigmas.eq.0.0) then
          call RNORMX(rvec,1,ranlux)
        else 
          rvec = Nsigmas
        end if
        nuke_fast = nuke_fast +dble(rvec)*dnuke_fast   

      end if
* 
      if (nuke_fast.lt.0.d0) nuke_fast = 0.d0
* 
      return 
 200  print *,'The Nuke set is inavailable (FILE: ',locfile,')'

      end 


