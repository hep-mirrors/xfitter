      double precision function pbtmdsubr(ipdf, x, qmu2, mc, mb, nameC, len)
C-------------------------------------------------------
C
C External PDF reading for QCDNUM
C
C--------------------------------------------------------
      implicit none
      character nameC*512   
      Integer len
      integer ipdf
c      logical first
      double precision x,qmu2
      double precision xf(-6:11)
      double precision mc,mb,mc2,mb2


      double precision xkt
      double precision xTMD(-6:11)
      
      Integer i,icall,icall_old,icall_new,ic

      Double Precision xold, q2old
      data icall/0/, icall_new/0/,icall_old/0/
      character*128  LHAPDFSET
      Logical lhapdf
      data lhapdf/.false./
c      data lhapdf/.true./
      
      Logical First
      Data First/.true./

      Data xold/-1./
      data q2old /-1./


      character outdir*132
      character*132 mystring
      Common/ myoutdir/ outdir

      outdir = nameC(1:len)
      PBTMDsubr = 0.0d0
      icall = icall + 1
      
c      write(6,*) ' in pbtmdsubr ipdf = ',ipdf
c      write(6,*) ' in pbtmdsubr lhapdf = ',lhapdf
      
      if (lhapdf) then
c        write(6,*) ' before evolve ',x,sqrt(qmu2)
        if(First) then
           LHAPDFSET ='HERAPDF20_LO_EIG'
           call InitPDFsetByName(LHAPDFSET)
           First=.False.
        endif
        call evolvePDF(x, sqrt(qmu2), xf)      
c        write(6,*) ' after evolve xf = ',(xf(i),i=-4,4)
        PBTMDsubr = xf(ipdf)
        return
      endif
      
c      write(6,*) ' in PBTMDsubr ', trim(outdir),' test '
c      write(6,*) ' in PBTMDsubr ', x, qmu2,ipdf
      if (x.eq.xold.and.qmu2.eq.q2old) then
          if(xf(ipdf).ne.xf(ipdf)) then
             write(6,*) ' problem 1 in PBTMDsubr_old :',x,qmu2,xf(ipdf)
          endif
          if(xf(ipdf).ge.1E20) then
             write(6,*) ' problem 2 in PBTMDsubr_old :',x,qmu2,xf(ipdf)
          endif
c          if(xf(ipdf).le.1E-20.and.xf(ipdf).ge.1E-300) then
c             write(6,*) ' problem 3 in PBTMDsubr_old :',x,qmu2,xf(ipdf)
c          endif
c          write(6,*) ' PBTMDsubr old  ',x,qmu2,ipdf,xf(ipdf)
          PBTMDsubr = xf(ipdf) 
          icall_old =icall_old + 1
      else
         icall_new =icall_new + 1
c        write(6,*) ' PBTMDsubr new ',x,qmu2
c        write(6,*) ' PBTMDsubr xold, x ', xold ,x, ' q2old ,qmu2 ', q2old ,qmu2
        do i=-6,11
          xf(i)= 0.d0
          xTMD(i) = 0.d0
        end do
        if (x.ne.xold) then
          xold = x
        endif
        if(qmu2.ne.q2old) then
          q2old = qmu2
c          write(6,*) ' PBTMDsubr ',qmu2
        endif
        xkt = -9999.
        mc2 = mc**2
        mb2 = mb**2
c      write(6,*) ' PBTMDsubr x = ',x,' mu2 = ',qmu2
c      write(6,*) ' PBTMDsubr ', ixold,iq2old
        call TMDconv(x, xkt , qmu2, xTMD)
        do i=-6,11
c        do i=-6,6
           xf(i)= max(0.d0,xTMD(i))
          if(qmu2.le.mc2.and.iabs(i).eq.4)  xf(i) = 0. 
          if(qmu2.le.mb2.and.iabs(i).eq.5)  xf(i) = 0.
c check          
c          if(iabs(i).eq.4)  xf(i) = 0. 
          if(xf(i).ne.xf(i)) then
             write(6,*) ' problem 1 in PBTMDsubr :',x,qmu2,i,xf(i)
          endif
          if(xf(i).ge.1E10) then
             write(6,*) ' problem 2 in PBTMDsubr :',x,qmu2,i,xf(i)
          endif
c          write(6,*) '  pdf ',i,xf(i),ipdf
        end do
        icall_new = icall_new + 1
        ic=100000 
        PBTMDsubr = xf(ipdf)
      endif 
c      if(ipdf.eq.8) then
c           write(6,*) ' pdf ',ipdf,xf(ipdf)
c      endif
c      write(6,*) ' PBTMBsubr icall = ', icall, icall_new,icall_old, ipdf
c      write(6,*) ' PBTMBsubr icall = ',  ipdf,xf(ipdf)
      return 
      end
      
      
      subroutine TMDconv(x,xkt,q2,xTMD)
      implicit None
      Double Precision x,xkt,q2,xTMD(-6:11)
      character outdir*132
      Common/ myoutdir/ outdir
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      Integer ipart
      Double precision xx,q2x,xxkt
      Common /myvalues/xx,q2x,xxkt,ipart
      Double Precision iTMDg,iTMDq
      External iTMDg,iTMDq
      Double Precision A0,B0,eps,sum,eps_set 
      Integer nseg
      Double Precision abstol,err
      Double Precision DGAUSS
      External DGAUSS
      Real Gauss
      External Gauss 
      Integer Idebug
      Common/mychecks/Idebug
      Integer IudGridfile /0/
      character testfiles*132
      
      Logical test
c      data test /.true./
      data test /.false./
      
      Logical first,ex
      data first/.true./
c      IudGridfile = 0
c check whether we have separated u- and d- quark grid files
      if(first) then 
         testfiles=trim(trim(outdir)//'/tmd-grid-u-quark_int.dat')
         inquire(FILE=testfiles,EXIST=ex)
         if(ex) then
           IudGridfile = 1 
         endif
         write(6,*) ' check for separated u- and d- quark grid files ', IudGridfile
         first = .false. 
      endif

c      idebug =1

      do ipart=-5,11
          xTMD(ipart) = 0.
      end do 
      eps_set = 0.0001 ! default
c      eps_set = 0.00001 
c      eps_set = 0.000001 
c	eps_set = 0.001 

c      Write (6,*) ' in TMDconv ',x,sqrt(q2)
        nseg = 10 ! needed for released grids 
c        nseg = 30  ! Sara
        xx = x
        q2x = q2
        xxkt = xkt
        if(xkt.lt.0) then 
c           mygridfiles=trim(trim(outdir)//'/tmd-grid-gluon_int.dat'//char(0))
           mygridfiles_d=trim(trim(outdir)//'/tmd-grid-gluon_int.dat')
        else
           mygridfiles_d=trim(trim(outdir)//'/tmd-grid-gluon.dat')
        endif 
        if(idebug.eq.1) write(6,*) ' idebug: TMDconv mygridfiles = ',mygridfiles_d, ' q2x = ',q2x

      if(q2x.le.2.0.or.test) then
        do ipart=-5,5
          call TMDstart(xx,q2x,xTMD)
        end do
      else
c         write(6,*) ' TMDconv: beginn  x,q2,xTMD ',x,sqrt(q2),xTMD(0),xTMD(7),xTMD(0)-xTMD(7)
        
        do ipart=-5,11
c        do ipart=-5,5 
c        do ipart=-4,4
          A0 = 0.0
          B0 = 1.
          eps = eps_set
          abstol = 0
          sum = 0
          xx = x
          if (xx.ge.1) then
c            write(6,*) ' TMDconv x> 1 ', x
            xx = 0.99999999
          endif
          q2x = q2
c          call dgadap(A0,B0,iTMDg,eps,sum)
          CALL dADAPT(iTMDg,A0,B0,NSEG,eps,abstol,sum,ERR)
c          SUM = DGAUSS(iTMDg,A0,B0,EPS)
c          sum=sum*1.01
          xTMD(ipart) = sum  
c          write(6,*) ' TMDconv gluons ',ipart, sum
        end do 
c      write(6,*) ' TMD glu: x,q2,xTMD ',x,sqrt(q2),(xTMD(i),i=-3,3)  
        if(xkt.lt.0) then 
           mygridfiles_d=trim(trim(outdir)//'/tmd-grid-quark_int.dat')
           mygridfiles_u=' '
           if(IudGridfile.eq.1) then 
              mygridfiles_d=trim(trim(outdir)//'/tmd-grid-d-quark_int.dat')
              mygridfiles_u=trim(trim(outdir)//'/tmd-grid-u-quark_int.dat')
              mygridfiles_ubar=trim(trim(outdir)//'/tmd-grid-ubar-quark_int.dat')
              mygridfiles_dbar=trim(trim(outdir)//'/tmd-grid-dbar-quark_int.dat')
              inquire(FILE=mygridfiles_ubar,EXIST=ex)
              if(.not.ex) then
                mygridfiles_ubar=' '
                mygridfiles_dbar=' '
              endif

           endif 
        else
           mygridfiles_d=trim(trim(outdir)//'/tmd-grid-quark.dat')
           mygridfiles_u=' '
c           write(6,*) ' check for u-files ', IudGridfile
           if(IudGridfile.eq.1) then 
              mygridfiles_d=trim(trim(outdir)//'/tmd-grid-d-quark.dat')
              mygridfiles_u=trim(trim(outdir)//'/tmd-grid-u-quark.dat')
           endif
        endif 
        if(idebug.eq.1) write(6,*) ' idebug: TMDconv mygridfiles = ',trim(mygridfiles_u),' test'
c        write(6,*) ' TMDconv: glu  x,q2,xTMD ',x,sqrt(q2),xTMD(0),xTMD(7),xTMD(0)-xTMD(7)
        do ipart=-5,11
c        do ipart=-5,5
c        do ipart=-5,5 
          A0 = 0.0
          B0 = 1.
          eps = eps_set
          abstol = 0
          sum = 0
          xx = x 
          if (xx.ge.1) then
c            write(6,*) ' TMDconv x> 1 ', x
            xx = 0.99999999
          endif
          q2x = q2
c          call dgadap(A0,B0,iTMDq,eps,sum)
           CALL dADAPT(iTMDq,A0,B0,NSEG,eps,abstol,sum,ERR)
c           SUM = DGAUSS(iTMDq,A0,B0,EPS)
c         write(6,*) ' x,p, sum',x,sqrt(q2),sum,ipart
c          write(6,*) ' TMDconv quarks ',ipart, sum
          xTMD(ipart) = sum + xTMD(ipart)
        end do 
      endif
c      write(6,*) ' TMDconv: x,q2,xTMD ',x,sqrt(q2),xTMD(0),xTMD(7),xTMD(0)-xTMD(7)

      Return
      End


     
      Function iTMDg(xt)
      Implicit None
      Double Precision iTMDg
      Double Precision xt
      double precision xr,pr,xpqr(-6:11)
      Double Precision xstart,qstart,xpq0(-6:11)
      Double Precision x0,xmin,xmax,wt,test
      Integer iparton
      Double precision x,q2,xkt
      Common /myvalues/x,q2,xkt,iparton
      Double Precision func22
      External func22
      Logical First_dum,Fccfm1,Fccfm2
      Common/ myfirst/First_dum,Fccfm1,Fccfm2
c      write(6,*) ' iTMD pdf-composition ',PDF_DECOMPOSITION
      pr = sqrt(q2)
      xmin=x
      xmax=0.9999999999
c    generate x0 with 1/x
      wt = log(xmax/xmin)
      x0= xmin*(xmax/xmin)**xt
      xr = x/x0
      if(xkt.lt.0) then  
c        write(6,*) ' before iTMDgridg ',xr,pr
        call iTMDgridg(xr,pr,xpqr) 
c        write(6,*) ' after iTMDgridg ',xr,pr,xpqr(0)
        else
        call TMDgridg(xr,xkt,pr,xpqr) 
      endif

      xstart=x0 
      call TMDstart(xstart,qstart,xpq0)
      test=0
c average over parton - antiparton for better stability
      if(iabs(iparton).le.6) then
         test = (xpqr(iparton)+xpqr(-iparton))/2.*xpq0(0)/x0
      elseif(iparton.ge.7.and.iparton.le.11) then 
         test = xpqr(iparton)*xpq0(0)/x0
      else
         write(6,*) ' iTMDg: wrong iparton ',iparton
      endif
c      write(6,*) ' iTMDg ',xpqr(iparton),xpqr(-iparton),(xpqr(iparton)+xpqr(-iparton))/2.
c       if(iparton.eq.7) write(6,*) ' in iTMDg ',test-(xpqr(0)+xpqr(-0))/2.*xpq0(0)/x0

c      write(6,*) ' iTMD ',iparton,iker,x,q2,xr,pr**2,xpqr(iker),x0,xpq0(0)
      iTMDg = test * wt * x0
c      if(iabs(iparton).ge.8.and.xpqr(8).ne.0) then 
c         write(6,*) ' iTMDg ',xpqr(iparton),xpq0(0),x0
c      endif 
     
      return
      end
      Function iTMDq(xt)
      Implicit None
      Double Precision iTMDq
      Double Precision xt
      double precision xr,pr,xpqr(-6:11),xpqru(-6:11)
      double precision xpqrtest(-6:11)
      double precision xpqrdbar(-6:11),xpqrubar(-6:11)
      Double Precision x0,xmin,xmax,wt,test
      Double Precision xstart,qstart,xpq0(-6:11)
      Integer iparton
      Double precision x,q2,xkt
      Common /myvalues/x,q2,xkt,iparton
c      Double Precision func22
c      External func22
c      Logical First_dum,Fccfm1,Fccfm2
c      Common/ myfirst/First_dum,Fccfm1,Fccfm2

      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      
      Logical first/.true./
      
      Double Precision Kaa,Kab,Kac,Kad,Kae,Kaf,Kaabar,Kabbar,Kacbar,Kadbar,Kaebar,Kafbar
      Double Precision Kaa_u,Kab_u,Kac_u,Kad_u,Kae_u,Kaf_u,Kaabar_u,Kabbar_u,Kacbar_u,Kadbar_u,Kaebar_u,Kafbar_u
      Integer i
c      write(6,*) ' iTMD pdf-composition ',PDF_DECOMPOSITION
      pr = sqrt(q2)
      xmin=x
c	xmax=0.99999
      xmax=0.9999999999
c    generate x0 with 1/x
      wt = log(xmax/xmin)
      x0= xmin*(xmax/xmin)**xt

      test = 0.
      xr = x/x0
      if(xkt.lt.0) then  
c        write(6,*) ' iTMDq: ',First_dum,Fccfm1,Fccfm2
c        write(6,*) ' before iTMDgridq ',xr,pr
c        write(6,*) ' before iTMDgridq ', mygridfiles 
c        call iTMDgridq(xr,pr,xpqr)
        if(mygridfiles_u.ne.' ') then
c          write(6,*) ' we have u- and d- type quark grid files', mygridfiles, mygridfiles_u
c          mygridfiles_u=mygridfiles
          call iTMDgridq_d(xr,pr,xpqr)
          call iTMDgridq_u(xr,pr,xpqru)
          if(mygridfiles_ubar.ne.' ') then
            call iTMDgridq_dbar(xr,pr,xpqrdbar)
            call iTMDgridq_ubar(xr,pr,xpqrubar)
          else
            do i=-6,11
              xpqrdbar(i) = xpqr(i)
              xpqrubar(i) = xpqru(i)
            end do
c rearrange parton-antiparton in case of only parton pdfs            
            xpqrdbar(-1) = xpqr(1)
            xpqrdbar(1) = xpqr(-1)
            xpqrubar(-2) = xpqru(2)
            xpqrubar(2) = xpqru(-2)           
          endif 
        else
          call iTMDgridq(xr,pr,xpqrtest)
          do i=-6,11
c              xpqr(i) = xpqrtest(i)
              xpqr(i) = 0
              xpqru(i)= 0
              xpqrdbar(i) = 0
              xpqrubar(i) = 0
          end do
          
          xpqr(0) = xpqrtest(0)
          
          xpqr(1)=xpqrtest(-3)
          xpqr(2)=(xpqrtest(-2)+xpqrtest(-1))/2.
          xpqr(3)=(xpqrtest(-2)+xpqrtest(-1))/2.
          xpqr(-1)=xpqrtest(3)
          xpqr(-2)=(xpqrtest(2)+xpqrtest(1))/2.
          xpqr(-3)=(xpqrtest(2)+xpqrtest(1))/2.
          xpqr(4)=xpqrtest(-4)
          xpqr(5)=xpqrtest(-5)
          xpqr(-4)=xpqrtest(4)
          xpqr(-5)=xpqrtest(5)

          xpqru(1) = (xpqrtest(-2)+xpqrtest(-1))/2.
          xpqru(2) = xpqrtest(-3)
          xpqru(-1) = (xpqrtest(2)+xpqrtest(1))/2.
          xpqru(-2) = xpqrtest(3)
c          xpqru(-3) = (xpqrtest(2)+xpqrtest(1))/2. ! test
          xpqru(4) = xpqrtest(-4)
          xpqru(5) = xpqrtest(-5)
          xpqru(-4) = xpqrtest(4)
          xpqru(-5) = xpqrtest(5)
          
          xpqrdbar(1) = xpqrtest(3)
          xpqrdbar(-1) = xpqrtest(-3)
          xpqrdbar(2) = (xpqrtest(2)+xpqrtest(1))/2.
          xpqrdbar(-2) = (xpqrtest(-2)+xpqrtest(-1))/2.
          xpqrdbar(3) = (xpqrtest(2)+xpqrtest(1))/2.
          xpqrdbar(-3) = (xpqrtest(-2)+xpqrtest(-1))/2.
          xpqrdbar(4) = xpqrtest(4)
          xpqrdbar(5) = xpqrtest(5)
          xpqrdbar(-4) = xpqrtest(-4)
          xpqrdbar(-5) = xpqrtest(-5)


          xpqrubar(1) = (xpqrtest(2)+xpqrtest(1))/2. 
          xpqrubar(2) = xpqrtest(3)
          xpqrubar(-1) = (xpqrtest(-2)+xpqrtest(-1))/2.
          xpqrubar(-2) = xpqrtest(-3)
          xpqrubar(4) = xpqrtest(4)
          xpqrubar(5) = xpqrtest(5)
          xpqrubar(-4) = xpqrtest(-4)
          xpqrubar(-5) = xpqrtest(-5)
          
        endif
c        write(6,*) ' after iTMDgridq ',xr,pr,xpqr(0)
        else
c        call TMDgridq(xr,xkt,pr,xpqr)
        if(mygridfiles_u.ne.' ') then
c          write(6,*) ' we have u- and d- type quark TMD grid files' 
          call TMDgridq(xr,xkt,pr,xpqr)
          call TMDgridq_u(xr,xkt,pr,xpqru)
          xpqrdbar(-1) = xpqr(1)
          xpqrdbar(1) = xpqr(-1)
          xpqrubar(-2) = xpqru(2)
          xpqrubar(2) = xpqru(-2)           
        endif
      endif

      xstart=x0 
      call TMDstart(xstart,qstart,xpq0)

C electroweak: W+ (iparton=9 ):  -1, -3, -5, 2, 4, 6     -> u, c, t          u-quark
C electroweak: W+ (iparton=9 ):  -1, -3, -5, 2, 4, 6     -> dbar, sbar, bbar u-quark
C electroweak: W- (iparton=10):  -2, -4, -6, 1, 3, 5     -> d, s, b          d-quark
C electroweak: W- (iparton=10):  -2, -4, -6, 1, 3, 5     -> ubar, cbar, tbar ubar-quark

      if(iabs(iparton).le.11) then
         if(iparton.eq.0) then 
            test = xpqr(0)*(xpq0(1)+xpq0(2)+xpq0(3)+xpq0(-1)+xpq0(-2)+xpq0(-3))/x0
         elseif(iparton.eq.1) then 
            test = (xpqr(1)*xpq0(1) + xpqru(1)*xpq0(2)+ xpqr(3)*xpq0(3)
     +            + xpqrdbar(1)*xpq0(-1)+ xpqrubar(1)*xpq0(-2) + xpqrdbar(3)*xpq0(-3))/x0
         elseif(iparton.eq.-1) then 
            test = (xpqr(-1)*xpq0(1) + xpqru(-1)*xpq0(2)+ xpqr(-3)*xpq0(3)
     +            + xpqrdbar(-1)*xpq0(-1)+ xpqrubar(-1)*xpq0(-2) + xpqrdbar(-3)*xpq0(-3))/x0
         elseif(iparton.eq.2) then 
            test = (xpqr(2)*xpq0(1) + xpqru(2)*xpq0(2)+ xpqr(3)*xpq0(3)
     +            + xpqrdbar(2)*xpq0(-1)+ xpqrubar(2)*xpq0(-2) + xpqrdbar(3)*xpq0(-3))/x0
         elseif(iparton.eq.-2) then 
            test = (xpqr(-2)*xpq0(1) + xpqru(-2)*xpq0(2)+ xpqr(-3)*xpq0(3)
     +            + xpqrdbar(-2)*xpq0(-1)+ xpqrubar(-2)*xpq0(-2) + xpqrdbar(-3)*xpq0(-3))/x0
         elseif(iparton.eq.3) then 
            test = (xpqr(3)*xpq0(1) + xpqru(1)*xpq0(2)+ xpqr(1)*xpq0(3)
     +            + xpqrdbar(3)*xpq0(-1)+ xpqrubar(1)*xpq0(-2) + xpqrdbar(1)*xpq0(-3))/x0
c            test = ((xpqr(-3)+xpqr(-2))/2.*xpq0(1) + (xpqr(-3)+xpqr(-2))/2.*xpq0(2)+ xpqr(1)*xpq0(3)
c     +            + (xpqrdbar(-3)+xpqrdbar(-2))/2.*xpq0(-1)+ (xpqrdbar(-3)+xpqrdbar(-2))/2.*xpq0(-2) + xpqrdbar(1)*xpq0(-3))/x0
         elseif(iparton.eq.-3) then 
            test = (xpqr(-3)*xpq0(1) + xpqru(-1)*xpq0(2)+ xpqr(-1)*xpq0(3)
     +            + xpqrdbar(-3)*xpq0(-1)+ xpqrubar(-1)*xpq0(-2) + xpqrdbar(-1)*xpq0(-3))/x0
c            test = ((xpqr(-3)+xpqr(-2))/2.*xpq0(1) + (xpqr(-3)+xpqr(-2))/2.*xpq0(2)+ xpqr(-1)*xpq0(3)
c     +            + (xpqrdbar(-3)+xpqrdbar(-2))/2.*xpq0(-1)+ (xpqrdbar(-3)+xpqrdbar(-2))/2.*xpq0(-2) + xpqrdbar(-1)*xpq0(-3))/x0
         elseif(iabs(iparton).ge.4.and.iabs(iparton).le.6) then 
            test = (xpqr(iparton)*xpq0(1) + xpqru(iparton)*xpq0(2)+ xpqr(iparton)*xpq0(3)
     +            + xpqrdbar(iparton)*xpq0(-1)+ xpqrubar(iparton)*xpq0(-2) + xpqrdbar(iparton)*xpq0(-3))/x0
         elseif(iparton.ge.7.and.iparton.le.8) then
            test = (xpqr(iparton)*(xpq0(1)+xpq0(3)) + xpqru(iparton)*(xpq0(2)+xpq0(4)) 
     +            + xpqrdbar(iparton)*(xpq0(-1)+xpq0(-3)) + xpqrubar(iparton)*(xpq0(-2)+xpq0(-4)) )/x0   ! gamma, Z0
         elseif(iparton.eq.9) then  
            test = (xpqru(iparton)*xpq0(2)+xpqrdbar(iparton)*xpq0(-1)+xpqrdbar(iparton)*xpq0(-3))/x0 ! Wplus
         elseif(iparton.eq.10) then 
            test = (xpqr(iparton)*xpq0(1)+xpqr(iparton)*xpq0(3)+xpqrubar(iparton)*xpq0(-2))/x0 ! Wminus         
         elseif(iparton.eq.11) then
            test = xpqr(iparton)*(xpq0(1)+xpq0(3)+xpq0(-1)+xpq0(-3))/x0 +
     +          xpqru(iparton)*(xpq0(2)+xpq0(-2))/x0 
         else
            write(6,*) ' iTMDq; wrong iparton ', iparton 
         endif
       endif
234    continue       
c      write(6,*) ' iTMD ',iparton,iker,x,q2,xr,pr**2,xpqr(iker),x0
       iTMDq = test * wt * x0
c      if(iabs(iparton).ge.8.and.xpqr(8).ne.0 ) then 
c         write(6,*) ' iTMDq ',xpqr(iparton),xpq0(1),x0
c      endif 
      return
      end


      subroutine TMDstart(xstart,Qstart,xpq0)
      Implicit none
      double precision xstart,Qstart,xpq0(-6:11) 
      double precision xglu,xdnv,xupv,xsbar,xubar,xdbar
      external xglu,xdnv,xupv,xsbar,xubar,xdbar

      Double Precision x0
      Double Precision funcSTART
      External funcSTART
      Integer i
      Integer icount,ic
      
      
      x0 = xstart
      Do i=-6,7
        xpq0(i) = 0
      end do
      Qstart = 1.9*1.9
c old QCDnum convention
c      xpq0(0)=func22(0,x0)
cc      xpq0(-1)=func22(5,x0)
c      xpq0(-1)=func22(5,x0) - func22(3,x0)/2.  ! def
cc      xpq0(1)=func22(1,x0)+func22(5,x0)
c      xpq0(1)=func22(1,x0)+func22(5,x0) - func22(3,x0)/2. ! def
c      xpq0(-2)=func22(4,x0)
c      xpq0(2)=func22(2,x0)+func22(4,x0)
c      xpq0(-3)=func22(3,x0)/2.
c      xpq0(3)=func22(3,x0)/2.
c new convention
      xpq0(0)=funcSTART(0,x0)
      xpq0(-1)=funcSTART(-1,x0)
      xpq0(1)=funcSTART(1,x0)
      xpq0(-2)=funcSTART(-2,x0)
      xpq0(2)=funcSTART(2,x0)
      xpq0(-3)=funcSTART(-3,x0)
      xpq0(3)=funcSTART(3,x0)
      
      
      Do i=-6,7
c      write(6,*) ' starting distribution ',i,xpq0(i),pdf(i)
        if(xpq0(i) .ne. xpq0(i)) then
           icount = icount + 1
           ic=1000000
           if (mod(icount,ic).eq.0)  write(6,*) ' problem in starting distribution: NaN-xpq0 set to zero:  ',i,x0,xpq0(i),icount
           xpq0(i) = 0
        endif
        if(xpq0(i) .gt. 1e10) then
           icount = icount + 1
c        if(xpq0(i) .ne. xpq0(i)) then
           ic=100000
           if (mod(icount,ic).eq.0) write(6,*) ' problem in starting distribution: too large-xpq0=1e10: ',i,x0,xpq0(i),icount
           xpq0(i) = 1e10
         endif
      end do
            
      return
      end
      
      
      subroutine write_pbtmd(filenameC, len)
      Implicit None
      character filenameC*132
      character filename*132
c      character*(*) filename
c      Common/updfout/filename
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      character outdir*132
      Common/ myoutdir/ outdir
      Integer n1,n2,n3
      Parameter (n1=51,n2=51,n3=51)
c      Parameter (n1=51,n2=201,n3=51)
      double Precision kt2,xx,px
      DIMENSION kt2(0:n1),xx(0:n2),px(0:n3)
      Double Precision rx,rq2,rp
      Double Precision RKMS6B,RKMS5B,RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0
      Double Precision RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
      Double Precision x,q2,xTMD(-6:11),phot
      Double Precision xkt
      character adum
      character *72 TXT
      LOGICAL ex
      Integer i,j,k,jj
      Integer iflag
      Integer IRR
      Integer ikincut,Ipgg,ns_sel,Iqqbar
      Real Qg0,Qscal,QCDlam
      
      Integer Idebug
      Common/mychecks/Idebug
      
      Integer len
      
      Iflag = 1 ! write only TMDs
      
      filename=filenameC(1:len)
c      filename=trim(filenameC)//char(0)
      
      Write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      Write(6,*) '++  write_pbtmd: Convolute TMD kernel = ',iflag, 'with starting distr  ++'
      Write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

ccc just read one of the files to get the proper binning
      if(iflag.eq.0) then
c         mygridfiles='TMDGrids/tmd-grid-gluon_int.dat'
         mygridfiles_d=trim(trim(outdir)//'/tmd-grid-gluon_int.dat')
        else
c         mygridfiles='TMDGrids/tmd-grid-gluon.dat'
         mygridfiles_d=trim(trim(outdir)//'/tmd-grid-gluon.dat')
      endif
c      write(6,*) ' write_pbtmd: mygridfiles = ',trim(mygridfiles)
      write(6,*) ' write_pbtmd: filename = ',filename
      
      open(30,FILE=trim(mygridfiles_d), FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
200   Read(30,101) TXT
101   Format(A72)
C         WRITE(6,101) '  line ',TXT
      If(TXT(1:4).EQ.'  Qg') then 
        read(txt,1000) adum,Qg0,adum,ikincut
1000    format(A7,f12.8,A10,I6)
c	      read(txt,1000) Qg0,ikincut
c1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         WRITE(6,101) ' 1st line ',TXT
        goto 200
      Endif
      If(TXT(1:4).EQ.' Qg0') then 
c	      read(txt,1000) Qg0,ikincut
c         WRITE(6,101) ' 1st line ',TXT
        read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
c1000        format(A7,f12.8,A10,I6)
        goto 200
      Endif
      If(TXT(1:4).EQ.' Ipg') then 
c	      read(txt,1001) Ipgg,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
        read(txt,1001) adum,Ipgg,adum,ns_sel
1001    format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
        goto 200
      Endif
      If(TXT(1:6).EQ.' Qscal') then 
        read(txt,1002)  adum,Qscal, adum,Iqqbar
1002    format(A9,f7.3,A10,I4)
c	      read(txt,1002) Qscal,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
C         WRITE(6,101) '2nd line ',TXT
        goto 200
      Endif
      If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
        read(txt,1003) adum,QCDLam
1003    format(A9,f12.8)
c            PARU(112)=QCDLam
        goto 200
      Endif
      If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
      Endif
      write(6,*) ' check ',n1,n2,n3
      if(Iflag.eq.0) then 
         do i=1,n2
            do k=1,n3
c for new  structure with W,Z
c                READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                         RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
c for old structure with photon only
                READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                         RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
            xx(i) = rx
            px(k) = rp
            end do
         end do
      elseif(iflag.eq.1) then
         do j=1,n1
          do i=1,n2
            do k=1,n3
c for old structure with photon only
                READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                         RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
            xx(i) = rx
            px(k) = rp
            kt2(j) = rq2
            end do
           end do
         end do
      else
        write(6,*) ' wrong iflag selected: stop ', Iflag
        stop
      endif
      close(30)
      
      inquire(FILE=filename,EXIST=ex)
      if(ex) then
         open(31,FILE=trim(filename), FORM='formatted',STATUS='old',IOSTAT=IRR,ERR=51 )
         write(6,*) '  old file found: overwritten ',filename
         else
         open(31,FILE=trim(filename), FORM='formatted',STATUS='new',IOSTAT=IRR,ERR=51 )
         write(6,*) ' new file created ',filename
      Endif
      write(31,10000) Qg0,ikincut
10000 format(' Qg0 = ',f12.8,' ikincut= ',I6)
      write(31,10001) Ipgg,ns_sel
10001 format(' Ipgg = ',I4,' ns_sel = ',I4)
      write(31,10002) QCDlam
10002 format(' QCDlam = ',f12.8)
      write(31,10100)
10100 format(' ln(xg),  ln(k_t^2), ln(p),  xgx')
      if(iflag.eq.0) then
      do i=1,n2
        do k=1,n3
        x = exp(xx(i))
        q2 = exp(px(k))
        q2=q2*q2
        xkt=-99999.
c        idebug = 1 
        call TMDconv(x,xkt,q2,xTMD)
        do jj=-6,11
          if(xTMD(jj).ne.xTMD(jj)) xTMD(jj) = 0.
        enddo
       
        write(31,*,Err=90 ) xx(i),px(k),(xTMD(jj),jj=-6,11)
        end do
c        write(6,*) ' outer loop:  i = ',i
      end do
      elseif(iflag.eq.1) then
      do j=1,n1
        do i=1,n2
          do k=1,n3
          x = exp(xx(i))
          q2 = exp(px(k))
          q2 = q2*q2
          xkt = exp(kt2(j))
c          idebug = 1 
          call TMDconv(x,xkt,q2,xTMD)
          do jj=-6,11
            xTMD(jj) = xTMD(jj)*xkt
            if(xTMD(jj).ne.xTMD(jj)) xTMD(jj) = 0.
          enddo
          write(31,*,Err=90 ) xx(i),kt2(j),px(k),(xTMD(jj),jj=-6,11) 
c         write(6,*,Err=90 ) ' test TMD ', xx(i),kt2(j),px(k),xTMD(7)
          end do
          end do
c        write(6,*) ' outer loop:  j = ',j
      end do
      
      endif
      close(31)
      
      Return 
   50 write(6,*) ' write_pbtmd: error opening file = ',trim(mygridfiles_d)
      Return
   51 write(6,*) ' write_pbtmd: error opening file = ',trim(filename)
      Return
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS0
      Return
      End

*
* Mathlib gen
*
*
      SUBROUTINE DADAPT(F,A,B,NSEG,RELTOL,ABSTOL,RES,ERR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
C     RES = Estimated Integral of F from A to B,
C     ERR = Estimated absolute error on RES.
C     NSEG  specifies how the adaptation is to be done:
C        =0   means use previous binning,
C        =1   means fully automatic, adapt until tolerance attained.
C        =n>1 means first split interval into n equal segments,
C             then adapt as necessary to attain tolerance.
C     The specified tolerances are:
C            relative: RELTOL ;  absolute: ABSTOL.
C        It stops when one OR the other is satisfied, or number of
C        segments exceeds NDIM.  Either TOLA or TOLR (but not both!)
C        can be set to zero, in which case only the other is used.
 
      EXTERNAL F
 
      PARAMETER (NDIM=100)
      PARAMETER (R1 = 1, HF = R1/2)
 
      DIMENSION XLO(NDIM),XHI(NDIM),TVAL(NDIM),TERS(NDIM)
      SAVE XLO,XHI,TVAL,TERS,NTER
      DATA NTER /0/
 
      IF(NSEG .LE. 0)  THEN
       IF(NTER .EQ. 0) THEN
        NSEGD=1
        GO TO 2
       ENDIF
       TVALS=0
       TERSS=0
       DO 1 I = 1,NTER
       CALL DGS56P(F,XLO(I),XHI(I),TVAL(I),TE)
       TERS(I)=TE**2
       TVALS=TVALS+TVAL(I)
       TERSS=TERSS+TERS(I)
    1  CONTINUE
       ROOT= SQRT(2*TERSS)
       GO TO 9
      ENDIF
      NSEGD=MIN(NSEG,NDIM)
    2 XHIB=A
      BIN=(B-A)/NSEGD
      DO 3 I = 1,NSEGD
      XLO(I)=XHIB
      XLOB=XLO(I)
      XHI(I)=XHIB+BIN
      IF(I .EQ. NSEGD) XHI(I)=B
      XHIB=XHI(I)
      CALL DGS56P(F,XLOB,XHIB,TVAL(I),TE)
      TERS(I)=TE**2
    3 CONTINUE
      NTER=NSEGD
      DO 4 ITER = 1,NDIM
      TVALS=TVAL(1)
      TERSS=TERS(1)
      DO 5 I = 2,NTER
      TVALS=TVALS+TVAL(I)
      TERSS=TERSS+TERS(I)
    5 CONTINUE
      ROOT= SQRT(2*TERSS)
      IF(ROOT .LE. ABSTOL .OR. ROOT .LE. RELTOL*ABS(TVALS)) GO TO 9
      IF(NTER .EQ. NDIM) GO TO 9
      BIGE=TERS(1)
      IBIG=1
      DO 6 I = 2,NTER
      IF(TERS(I) .GT. BIGE) THEN
       BIGE=TERS(I)
       IBIG=I
      ENDIF
    6 CONTINUE
      NTER=NTER+1
      XHI(NTER)=XHI(IBIG)
      XNEW=HF*(XLO(IBIG)+XHI(IBIG))
      XHI(IBIG)=XNEW
      XLO(NTER)=XNEW
      CALL DGS56P(F,XLO(IBIG),XHI(IBIG),TVAL(IBIG),TE)
      TERS(IBIG)=TE**2
      CALL DGS56P(F,XLO(NTER),XHI(NTER),TVAL(NTER),TE)
      TERS(NTER)=TE**2
    4 CONTINUE
    9 RES=TVALS
      ERR=ROOT
      RETURN
      END


*
* Mathlib gen
*
*
      SUBROUTINE DGS56P(F,A,B,RES,ERR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      PARAMETER (R1 = 1, HF = R1/2)
      DIMENSION X5(5),W5(5),X6(6),W6(6)
 
      DATA (X5(I),W5(I),I=1,5)
     1/4.6910077030668004D-02, 1.1846344252809454D-01,
     2 2.3076534494715846D-01, 2.3931433524968324D-01,
     3 5.0000000000000000D-01, 2.8444444444444444D-01,
     4 7.6923465505284154D-01, 2.3931433524968324D-01,
     5 9.5308992296933200D-01, 1.1846344252809454D-01/
 
      DATA (X6(I),W6(I),I=1,6)
     1/3.3765242898423989D-02, 8.5662246189585178D-02,
     2 1.6939530676686775D-01, 1.8038078652406930D-01,
     3 3.8069040695840155D-01, 2.3395696728634552D-01,
     4 6.1930959304159845D-01, 2.3395696728634552D-01,
     5 8.3060469323313225D-01, 1.8038078652406930D-01,
     6 9.6623475710157601D-01, 8.5662246189585178D-02/
 
      RANG=B-A
      E5=0
      E6=0
      DO 1 I = 1,5
      E5=E5+W5(I)*F(A+RANG*X5(I))
      E6=E6+W6(I)*F(A+RANG*X6(I))
    1 CONTINUE
      E6=E6+W6(6)*F(A+RANG*X6(6))
      RES=HF*(E6+E5)*RANG
      ERR=ABS((E6-E5)*RANG)
      RETURN
      END
*
* Mathlib gen
*
*

      SUBROUTINE RADAPT(F,A,B,NSEG,RELTOL,ABSTOL,RES,ERR)
 
C     RES = Estimated Integral of F from A to B,
C     ERR = Estimated absolute error on RES.
C     NSEG  specifies how the adaptation is to be done:
C        =0   means use previous binning,
C        =1   means fully automatic, adapt until tolerance attained.
C        =n>1 means first split interval into n equal segments,
C             then adapt as necessary to attain tolerance.
C     The specified tolerances are:
C            relative: RELTOL ;  absolute: ABSTOL.
C        It stops when one OR the other is satisfied, or number of
C        segments exceeds NDIM.  Either TOLA or TOLR (but not both!)
C        can be set to zero, in which case only the other is used.
 
      DOUBLE PRECISION TVALS,TERSS
      EXTERNAL F
 
      PARAMETER (NDIM=100)
      PARAMETER (R1 = 1, HF = R1/2)
 
      DIMENSION XLO(NDIM),XHI(NDIM),TVAL(NDIM),TERS(NDIM)
      SAVE XLO,XHI,TVAL,TERS,NTER
      DATA NTER /0/
 
      IF(NSEG .LE. 0)  THEN
       IF(NTER .EQ. 0) THEN
        NSEGD=1
        GO TO 2
       ENDIF
       TVALS=0
       TERSS=0
       DO 1 I = 1,NTER
       CALL RGS56P(F,XLO(I),XHI(I),TVAL(I),TE)
       TERS(I)=TE**2
       TVALS=TVALS+TVAL(I)
       TERSS=TERSS+TERS(I)
    1  CONTINUE
       ROOT= SQRT(2*TERSS)
       GO TO 9
      ENDIF
      NSEGD=MIN(NSEG,NDIM)
    2 XHIB=A
      BIN=(B-A)/NSEGD
      DO 3 I = 1,NSEGD
      XLO(I)=XHIB
      XLOB=XLO(I)
      XHI(I)=XHIB+BIN
      IF(I .EQ. NSEGD) XHI(I)=B
      XHIB=XHI(I)
      CALL RGS56P(F,XLOB,XHIB,TVAL(I),TE)
      TERS(I)=TE**2
    3 CONTINUE
      NTER=NSEGD
      DO 4 ITER = 1,NDIM
      TVALS=TVAL(1)
      TERSS=TERS(1)
      DO 5 I = 2,NTER
      TVALS=TVALS+TVAL(I)
      TERSS=TERSS+TERS(I)
    5 CONTINUE
      ROOT= SQRT(2*TERSS)
      IF(ROOT .LE. ABSTOL .OR. ROOT .LE. RELTOL*ABS(TVALS)) GO TO 9
      IF(NTER .EQ. NDIM) GO TO 9
      BIGE=TERS(1)
      IBIG=1
      DO 6 I = 2,NTER
      IF(TERS(I) .GT. BIGE) THEN
       BIGE=TERS(I)
       IBIG=I
      ENDIF
    6 CONTINUE
      NTER=NTER+1
      XHI(NTER)=XHI(IBIG)
      XNEW=HF*(XLO(IBIG)+XHI(IBIG))
      XHI(IBIG)=XNEW
      XLO(NTER)=XNEW
      CALL RGS56P(F,XLO(IBIG),XHI(IBIG),TVAL(IBIG),TE)
      TERS(IBIG)=TE**2
      CALL RGS56P(F,XLO(NTER),XHI(NTER),TVAL(NTER),TE)
      TERS(NTER)=TE**2
    4 CONTINUE
    9 RES=TVALS
      ERR=ROOT
      RETURN
      END

*
* Mathlib gen
*
*
      SUBROUTINE RGS56P(F,A,B,RES,ERR)
      DOUBLE PRECISION E5,E6
 
      PARAMETER (R1 = 1, HF = R1/2)
      DIMENSION X5(5),W5(5),X6(6),W6(6)
 
      DATA (X5(I),W5(I),I=1,5)
     1/4.6910077030668004D-02, 1.1846344252809454D-01,
     2 2.3076534494715846D-01, 2.3931433524968324D-01,
     3 5.0000000000000000D-01, 2.8444444444444444D-01,
     4 7.6923465505284154D-01, 2.3931433524968324D-01,
     5 9.5308992296933200D-01, 1.1846344252809454D-01/
 
      DATA (X6(I),W6(I),I=1,6)
     1/3.3765242898423989D-02, 8.5662246189585178D-02,
     2 1.6939530676686775D-01, 1.8038078652406930D-01,
     3 3.8069040695840155D-01, 2.3395696728634552D-01,
     4 6.1930959304159845D-01, 2.3395696728634552D-01,
     5 8.3060469323313225D-01, 1.8038078652406930D-01,
     6 9.6623475710157601D-01, 8.5662246189585178D-02/
 
      RANG=B-A
      E5=0
      E6=0
      DO 1 I = 1,5
      E5=E5+W5(I)*F(A+RANG*X5(I))
      E6=E6+W6(I)*F(A+RANG*X6(I))
    1 CONTINUE
      E6=E6+W6(6)*F(A+RANG*X6(6))
      RES=HF*(E6+E5)*RANG
      ERR=ABS((E6-E5)*RANG)
      RETURN
      END

         
           
      SUBROUTINE iTMDgridg(X,P,XPQ)
      Implicit None
      Integer n1,n2
      double precision XPQ(-6:11),X,P
      Parameter (n1=51,n2=51)
      double precision xx,px
      DIMENSION XX(0:n1),PX(0:n2)

      Double Precision XA(3),A(N1+N2)
      Double Precision f_grid0(n1,n2),f_grid1(n1,n2),f_grid2(n1,n2),f_grid3(n1,n2),f_grid4(n1,n2),f_grid5(n1,n2),f_grid6(n1,n2)
      Double Precision f_grid7(n1,n2),f_grid8(n1,n2),f_grid9(n1,n2),f_grid10(n1,n2),f_grid11(n1,n2)
      Double Precision f_grid1m(n1,n2),f_grid2m(n1,n2),f_grid3m(n1,n2),f_grid4m(n1,n2),f_grid5m(n1,n2),f_grid6m(n1,n2)
      INTEGER NA(2)
      DATA NA/n1,n2/
      Double Precision DHFINT
      Double Precision DHPOLINT
      
      Double Precision x1a(n1),x2a(n2),xin,yin

      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision Qscal 
      Integer Iqqbar
      Common/GLUDAT2/Qscal,Iqqbar 
      Integer i,j,k,in,irr,igrid,ip

      Double Precision scal,rx,rq2,rp,glu,up,down,str,usea,dsea,ssea,charm,csea,bot,bsea,top,tsea,phot
      Double Precision z0,wplus,wminus,higgs
      Double Precision RKMS6B,RKMS5B,RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6
      Double Precision RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
      
      Integer kx,ky,mx,my,ndimz,knot
      Double Precision DSPPS2
      External DSPPS2      
      Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
      common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      Real QCDLam
      character *72 TXT
      character adum
      LOGICAL FIRST
      DATA FIRST/.TRUE./
           
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      
      character mygridfile_old*132
      
      Logical lwrite/.true./
            
      Integer inter
      Data inter/1/
      
      
c      write(6,*) ' iTMDgridg  ',mygridfiles,P
c      write(6,*) ' iTMDgridg  inter = ',inter
      IF(first) THEN
      if(mygridfiles_d.ne.'mygridfile_old') then
c         write(6,*) ' iTMDg: files ',mygridfiles,mygridfile_old
         mygridfile_old=mygridfiles_d
         n1min = 0
         n1max = 0
         n2min = 0
         n2max = 0
         n3min = 0
         n3max = 0
         ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
         ns_sel = -1
c        if(lwrite)  write(6,*) ' iTMDgridg ',first,iset
    
c         len = len_trim(mygridfiles)
c         my_string=mygridfiles(1:len)//char(0)
c          'TMD'//trim(temp_wname)//'0'//char(0)
c         mygridfiles=trim(mygridfiles)
           if(lwrite) write(6,*) ' iTMDgridg reading: ',trim(mygridfiles_d)
           open(30,FILE=trim(mygridfiles_d), FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
c           open(30,FILE=trim(my_string), FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
c         if(lwrite) write(6,*) ' iTMDgridg: read ',trim(mygridfiles)                
200      Read(30,101) TXT
  101    Format(A72)
c         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
            read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
            read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
            read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
            read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
C         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) 'iTMDgridg ln line:',TXT, TXT(32:36), TXT(37:40),TXT(41:44)
c         write(6,*) ' end ln line read '
          Endif
          
          if(lwrite) write(6,'(''  iTMDgridg: starting scale Q0 = '',F8.3)'),Qg0
c         if(lwrite) write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
c         if(lwrite) write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,' ns_sel = ',ns_sel
c         if(lwrite) write(6,*) ' QCD_lam used in uPDF: ',QCDlam
c         if(lwrite) write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar

            do i=1,n1
               do k=1,n2

c                  READ(30,*,Err=90,END=777 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,
c     &            RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
c                  write(6,*) ' reading EWK set ' 
c                  goto 888
c777               READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,
c     &            RKMS6                                

                  if(ns_sel.ge.1000) then 
c                     write(6,*) ' reading full EW set '
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
                  else
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6
                  endif 
c                  write(6,*) ' reading old set ' 
c888               continue     
     

                  xx(i) = rx
                  px(k) = rp
c                  write(6,*) ' check read ',i,exp(xx(i)),k,exp(px(k))
                  f_grid0(i,k) = max(0.d0,rkms0)*scal
                  f_grid1(i,k) = max(0.d0,rkms1)*scal
                  f_grid2(i,k) = max(0.d0,rkms2)*scal
                  f_grid3(i,k) = max(0.d0,rkms3)*scal
                  f_grid4(i,k) = max(0.d0,rkms4)*scal
                  f_grid5(i,k) = max(0.d0,rkms5)*scal
                  f_grid6(i,k) = max(0.d0,rkms6)*scal
                  f_grid7(i,k) = max(0.d0,rkms7)*scal
                  f_grid8(i,k) = max(0.d0,rkms8)*scal
                  f_grid9(i,k) = max(0.d0,rkms9)*scal
                  f_grid10(i,k) = max(0.d0,rkms10)*scal
                  f_grid11(i,k) = max(0.d0,rkms11)*scal
                  f_grid1m(i,k) = max(0.d0,rkms1b)*scal
                  f_grid2m(i,k) = max(0.d0,rkms2b)*scal
                  f_grid3m(i,k) = max(0.d0,rkms3b)*scal
                  f_grid4m(i,k) = max(0.d0,rkms4b)*scal
                  f_grid5m(i,k) = max(0.d0,rkms5b)*scal
                  f_grid6m(i,k) = max(0.d0,rkms6b)*scal
c                  if(f_grid3m(i,k).gt.0) write(6,*) 'iTMDgrid ',f_grid3m(i,k),i,k
c                   write(6,*) rx,rp
               enddo
            enddo

c         if(lwrite) write(6,*) ' end of file at ',i,k

         IN=0
         DO I=1,N1
            IN=IN+1
            A(IN) = xx(I)
            x1a(i) = xx(i)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = px(I)
            x2a(i) = px(i)
         ENDDO
         FIRST=.FALSE.
c         write(6,*) '  parton densities read from file unit 30 '
         if (inter.eq.1) then
             kx = 1
             ky = 1
             knot = 1
             ndimz = n1 * n2
             mx = n1 + kx + 1
             my = n2 + ky + 1
         endif

         Close(30)
      endif
      ENDIF

      ncall = ncall +1 


      XA(1) = log(X)
      XA(2) = log(P)
      xin = xa(1)
      yin = xa(2)
      if(xa(2).lt.px(1)) then
         n2min = n2min + 1  
         if(n2min.lt.10) then 
            write(6,*) ' iTMDgridg: p out of range ',p,' min p ',exp(px(1))
            write(6,*) ' iTMDgridg: p out of range ',xa(2),px(1)
         elseif(n2min.eq.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
         endif
         xa(2)=px(1)
      endif
      if(xa(2).gt.px(n2)) then
         n2max = n2max + 1
         if(n2max.lt.10) then 
            write(6,*) ' iTMDgridg: p out of range ',p,' max p ',exp(px(n2))
         elseif(n2max.eq.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
         endif
         xa(2)=px(n2)
      endif
c      write(6,*) xa(1) 
      if(xa(1).ge.xx(n1)) xa(1)=xx(n1)-0.0001
      if(xa(1).lt.xx(1)) then
         n1min = n1min + 1
         if(n1min.lt.10) then 
            write(6,*) ' iTMDgridg:  x out of range ',x,' min ',exp(xx(1))
         elseif(n1min.eq.10) then 
            write(6,*) '  iTMDgridg: x out of range ',x,' min ',exp(xx(1))
            write(6,*) '  iTMDgridg: last message printed: min x'
         endif
         xa(1)=xx(1)
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(inter.eq.0) then
         glu  = DHFINT(2,XA,NA,A,f_grid0)
         up   = DHFINT(2,XA,NA,A,f_grid2)
         down = DHFINT(2,XA,NA,A,f_grid1)
         str  = DHFINT(2,XA,NA,A,f_grid3)
         usea = DHFINT(2,XA,NA,A,f_grid2m)
         dsea = DHFINT(2,XA,NA,A,f_grid1m)
         ssea = DHFINT(2,XA,NA,A,f_grid3m)
         charm= DHFINT(2,XA,NA,A,f_grid4)
         csea = DHFINT(2,XA,NA,A,f_grid4m)
         bot  = DHFINT(2,XA,NA,A,f_grid5)
         bsea = DHFINT(2,XA,NA,A,f_grid5m)
         top  = DHFINT(2,XA,NA,A,f_grid6)
         tsea = DHFINT(2,XA,NA,A,f_grid6m)
         phot = DHFINT(2,XA,NA,A,f_grid7)
         z0 =   DHFINT(2,XA,NA,A,f_grid8)
         wplus = DHFINT(2,XA,NA,A,f_grid9)
         wminus = DHFINT(2,XA,NA,A,f_grid10)
         higgs = DHFINT(2,XA,NA,A,f_grid11)

c      write(6,*) ' DHFINT glu = ', glu
       elseif(inter.eq.1) then
         glu  = dhpolint(2,XA,NA,A,f_grid0)
         up   = dhpolint(2,XA,NA,A,f_grid2)
         down = dhpolint(2,XA,NA,A,f_grid1)
         str  = dhpolint(2,XA,NA,A,f_grid3)
         usea = dhpolint(2,XA,NA,A,f_grid1m)
         dsea = dhpolint(2,XA,NA,A,f_grid2m)
         ssea = dhpolint(2,XA,NA,A,f_grid3m)
         charm= dhpolint(2,XA,NA,A,f_grid4)
         csea = dhpolint(2,XA,NA,A,f_grid4m)
         bot  = dhpolint(2,XA,NA,A,f_grid5)
         bsea = dhpolint(2,XA,NA,A,f_grid5m)
         top  = dhpolint(2,XA,NA,A,f_grid6)
         tsea = dhpolint(2,XA,NA,A,f_grid6m)
         phot = dhpolint(2,XA,NA,A,f_grid7)
         z0 =   dhpolint(2,XA,NA,A,f_grid8)
         wplus= dhpolint(2,XA,NA,A,f_grid9)
         wminus=dhpolint(2,XA,NA,A,f_grid10)
         higgs =dhpolint(2,XA,NA,A,f_grid11)
         
cc         write(6,*) ' glu ',glu,glun
      else
        write(6,*) ' iTMDgridg inter = ',inter,' not implemented '
      endif


      DO  IP=-6,11
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(up.gt.0.) xpq(2) = up
      if(down.gt.0.) xpq(1) = down
      if(usea.gt.0) xpq(-2) = usea
      if(dsea.gt.0) xpq(-1) = dsea
      if(str.gt.0.) xpq(3) = str
      if(ssea.gt.0.) xpq(-3) = ssea
      if(charm.gt.0.) xpq(4) = charm
      if(csea.gt.0.) xpq(-4) = csea
      if(bot.gt.0.) xpq(5) = bot
      if(bsea.gt.0.) xpq(-5) = bsea
      if(top.gt.0.) xpq(6) = top
      if(tsea.gt.0.) xpq(-6) = tsea
      if(phot.gt.0.) xpq(7) = phot
      if(z0.gt.0.) xpq(8) = z0
      if(wplus.gt.0.) xpq(9) = wplus
      if(wminus.gt.0.) xpq(10) = wminus
      if(higgs.gt.0.) xpq(11) = higgs
c      if(z0.gt.0) then
c         write(6,*) ' iTMDgridg x,p,xpq',x,p,xpq(8)
ccc         write(6,*) ' i,j,k ',i,j,k
c      endif
c      do i=-6,7
c        xpq(i) = xpq(i)/(1.-min(x,0.999))
c      end do   
c      write(6,*) ' new iTMDgridg x,q2,p,xpq',x,q2,p,glu
      return
   50 write(6,*) ' iTMDgridg: error in opening file ',mygridfiles_d,' test'
      stop
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS0
      stop
      END
      
      
             
      SUBROUTINE iTMDgridq(X,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      double precision XPQ(-6:11),X,P
*! new
      Parameter (n1=51,n2=51)
      Double Precision xx,px
      DIMENSION XX(0:n1),PX(0:n2)

      Double Precision XA(3),A(N1+N2)
      Double Precision f_grid0(n1,n2),f_grid1(n1,n2),f_grid2(n1,n2),f_grid3(n1,n2),f_grid4(n1,n2),f_grid5(n1,n2),f_grid6(n1,n2)
      Double Precision f_grid7(n1,n2),f_grid8(n1,n2),f_grid9(n1,n2),f_grid10(n1,n2),f_grid11(n1,n2)
      Double Precision f_grid1m(n1,n2),f_grid2m(n1,n2),f_grid3m(n1,n2),f_grid4m(n1,n2),f_grid5m(n1,n2),f_grid6m(n1,n2)
      INTEGER NA(2)
      DATA NA/n1,n2/
      Double Precision DHFINT, dhpolint

      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision Qscal 
      Integer Iqqbar
      Common/GLUDAT2/Qscal,Iqqbar 
      Integer i,j,k,in,irr,igrid,ip
      Double Precision scal,rx,rq2,rp,glu,up,down,str,usea,dsea,ssea,charm,csea,bot,bsea,top,tsea,phot
      Double Precision z0,wplus,wminus,higgs
      Double Precision RKMS6B,RKMS5B,RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
      Double Precision RKMS8,RKMS9,RKMS10,RKMS11
      Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
      common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      Real QCDLam
      character *72 TXT
      character adum
      LOGICAL FIRST
      DATA FIRST/.TRUE./
           
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      
      character mygridfile_old*132
      
      Logical lwrite/.true./
      
      Integer inter
      data inter/1/
      
c      write(6,*) ' in quark file '
c      write(6,*) ' iTMDgrid first ',first,first_meoffsh
c         write(6,*) ' iTMDgridq: files ',mygridfiles,mygridfile_old
      IF(first) THEN
      if(mygridfiles_d.ne.'mygridfile_old') then
c         write(6,*) ' iTMDq: files ',mygridfiles,mygridfile_old
         mygridfile_old=mygridfiles_d
         n1min = 0
         n1max = 0
         n2min = 0
         n2max = 0
         n3min = 0
         n3max = 0
         ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
         ns_sel = -1
c        if(lwrite)  write(6,*) ' iTMDgridq ',first,iset
           if(lwrite) write(6,*) ' iTMDgridq_d reading: ',trim(mygridfiles_d)
           open(30,FILE=mygridfiles_d, FORM='formatted',
     +          STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
c         if(lwrite) write(6,*) ' iTMDgridq_d: read  ',trim(mygridfiles)                
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx') 
200      Read(30,101) TXT
  101    Format(A72)
c         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
            read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
            read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
            read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
            read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
C         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif
         
c         if(lwrite) write(6,*) ' iTMDgridq_d: starting scale Q0 = ',Qg0
         if(lwrite) write(6,'(''  iTMDgridq_d: starting scale Q0 = '',F8.3)'),Qg0
c         if(lwrite) write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
c         if(lwrite) write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,' ns_sel = ',ns_sel
c         if(lwrite) write(6,*) ' QCD_lam used in uPDF: ',QCDlam
c         if(lwrite) write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar

            do i=1,n1
               do k=1,n2
cc                  write(6,*) j,i,k,RX,RQ2,RP,RKMS
cc                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS
cc                  write(6,*) ' in iTMDgrid '

                  if(ns_sel.ge.1000) then 
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
                  else
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6
                  endif 
c                  write(6,*) ' grid output format ',RX,RP 

                  xx(i) = rx
                  px(k) = rp
                  f_grid0(i,k) = max(0.d0,rkms0)*scal
                  f_grid1(i,k) = max(0.d0,rkms1)*scal
                  f_grid2(i,k) = max(0.d0,rkms2)*scal
                  f_grid3(i,k) = max(0.d0,rkms3)*scal
                  f_grid4(i,k) = max(0.d0,rkms4)*scal
                  f_grid5(i,k) = max(0.d0,rkms5)*scal
                  f_grid6(i,k) = max(0.d0,rkms6)*scal
                  f_grid7(i,k) = max(0.d0,rkms7)*scal
                  f_grid8(i,k) = max(0.d0,rkms8)*scal
                  f_grid9(i,k) = max(0.d0,rkms9)*scal
                  f_grid10(i,k) = max(0.d0,rkms10)*scal
                  f_grid11(i,k) = max(0.d0,rkms11)*scal
                  f_grid1m(i,k) = max(0.d0,rkms1b)*scal
                  f_grid2m(i,k) = max(0.d0,rkms2b)*scal
                  f_grid3m(i,k) = max(0.d0,rkms3b)*scal
                  f_grid4m(i,k) = max(0.d0,rkms4b)*scal
                  f_grid5m(i,k) = max(0.d0,rkms5b)*scal
                  f_grid6m(i,k) = max(0.d0,rkms6b)*scal
c                  if(f_grid3m(i,k).gt.0) write(6,*) 'iTMDgrid ',f_grid3m(i,k),i,k
c                   write(6,*) rx,rp
               enddo
            enddo

c         if(lwrite) write(6,*) ' end of file at ',i,k
         IN=0
         DO I=1,N1
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         FIRST=.FALSE.
c         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
      endif
      ENDIF

      ncall = ncall +1   


      XA(1) = log(X)
      XA(2) = log(P)
      if(xa(2).lt.px(1)) then
         n2min = n2min + 1

         if(n2min.lt.10) then 
            write(6,*) ' iTMDgridq-d: p out of range ',p,' min p ',exp(px(1))
            write(6,*) ' iTMDgridq_d: p out of range ',xa(2),px(1)
         elseif(n2min.eq.10) then 
c            write(6,*) ' iTMDgridq_d: p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
         endif
         xa(2)=px(1)
      endif
      if(xa(2).gt.px(n2)) then
         n2max = n2max + 1
         if(n2max.lt.10) then 
            write(6,*) ' iTMDgridq_d: p out of range ',p,' max p ',exp(px(n3))
         elseif(n2max.eq.10) then 
c            write(6,*) ' iTMDgridq_d: p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
         endif
         xa(2)=px(n2)
      endif
      if(xa(1).ge.xx(n1)) xa(1)=xx(n1)-0.0001
      if(xa(1).lt.xx(1)) then
         n1min = n1min + 1
         if(n1min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
         elseif(n1min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
         endif
         xa(1)=xx(1)
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(inter.eq.0) then
         glu  = DHFINT(2,XA,NA,A,f_grid0)
         up   = DHFINT(2,XA,NA,A,f_grid2)
         down = DHFINT(2,XA,NA,A,f_grid1)
         str  = DHFINT(2,XA,NA,A,f_grid3)
         usea = DHFINT(2,XA,NA,A,f_grid2m)
         dsea = DHFINT(2,XA,NA,A,f_grid1m)
         ssea = DHFINT(2,XA,NA,A,f_grid3m)
         charm= DHFINT(2,XA,NA,A,f_grid4)
         csea = DHFINT(2,XA,NA,A,f_grid4m)
         bot  = DHFINT(2,XA,NA,A,f_grid5)
         bsea = DHFINT(2,XA,NA,A,f_grid5m)
         top  = DHFINT(2,XA,NA,A,f_grid6)
         tsea = DHFINT(2,XA,NA,A,f_grid6m)
         phot = DHFINT(2,XA,NA,A,f_grid7)
         z0 =   DHFINT(2,XA,NA,A,f_grid8)
         wplus =DHFINT(2,XA,NA,A,f_grid9)
         wminus=DHFINT(2,XA,NA,A,f_grid10)
         higgs =DHFINT(2,XA,NA,A,f_grid11)
      elseif(inter.eq.1) then
         glu  = dhpolint(2,XA,NA,A,f_grid0)
         up   = dhpolint(2,XA,NA,A,f_grid2)
         down = dhpolint(2,XA,NA,A,f_grid1)
         str  = dhpolint(2,XA,NA,A,f_grid3)
         usea = dhpolint(2,XA,NA,A,f_grid2m)
         dsea = dhpolint(2,XA,NA,A,f_grid1m)
         ssea = dhpolint(2,XA,NA,A,f_grid3m)
         charm= dhpolint(2,XA,NA,A,f_grid4)
         csea = dhpolint(2,XA,NA,A,f_grid4m)
         bot  = dhpolint(2,XA,NA,A,f_grid5)
         bsea = dhpolint(2,XA,NA,A,f_grid5m)
         top  = dhpolint(2,XA,NA,A,f_grid6)
         tsea = dhpolint(2,XA,NA,A,f_grid6m)
         phot = dhpolint(2,XA,NA,A,f_grid7)
         z0 =   dhpolint(2,XA,NA,A,f_grid8)
         wplus =dhpolint(2,XA,NA,A,f_grid9)
         wminus=dhpolint(2,XA,NA,A,f_grid10)
         higgs =dhpolint(2,XA,NA,A,f_grid11)
      else
        write(6,*) ' iTMDgridq_d inter = ',inter,' not implemented '
      endif
      DO  IP=-6,11
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(up.gt.0.) xpq(2) = up
      if(down.gt.0.) xpq(1) = down
      if(usea.gt.0) xpq(-2) = usea
      if(dsea.gt.0) xpq(-1) = dsea
      if(str.gt.0.) xpq(3) = str
      if(ssea.gt.0.) xpq(-3) = ssea
      if(charm.gt.0.) xpq(4) = charm
      if(csea.gt.0.) xpq(-4) = csea
      if(bot.gt.0.) xpq(5) = bot
      if(bsea.gt.0.) xpq(-5) = bsea
      if(top.gt.0.) xpq(6) = top
      if(tsea.gt.0.) xpq(-6) = tsea
      if(phot.gt.0.) xpq(7) = phot
      if(z0.gt.0.) xpq(8) = z0
      if(wplus.gt.0.) xpq(9) = wplus
      if(wminus.gt.0.) xpq(10) = wminus
      if(higgs.gt.0.) xpq(11) = higgs
c      if(z0.gt.0) then
c         write(6,*) ' iTMDgridq x,p,xpq',x,p,xpq(8)
ccc         write(6,*) ' i,j,k ',i,j,k
c      endif
c      do i=-6,7
c        xpq(i) = xpq(i)/(1.-min(x,0.999))
c      end do   
         
c         write(6,*) ' new iTMDgrid x,q2,p,xpq',x,q2,p,glu
      return
   50 write(6,*) ' iTMDgridq_d: error in opening file ', mygridfiles_d
      stop
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS0
      stop
      END
      

      SUBROUTINE iTMDgridq_u(X,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      double precision XPQ(-6:11),X,P
*! new
      Parameter (n1=51,n2=51)
      Double Precision xx,px
      DIMENSION XX(0:n1),PX(0:n2)

      Double Precision XA(3),A(N1+N2)
      Double Precision f_grid0(n1,n2),f_grid1(n1,n2),f_grid2(n1,n2),f_grid3(n1,n2),f_grid4(n1,n2),f_grid5(n1,n2),f_grid6(n1,n2)
      Double Precision f_grid7(n1,n2),f_grid8(n1,n2),f_grid9(n1,n2),f_grid10(n1,n2),f_grid11(n1,n2)
      Double Precision f_grid1m(n1,n2),f_grid2m(n1,n2),f_grid3m(n1,n2),f_grid4m(n1,n2),f_grid5m(n1,n2),f_grid6m(n1,n2)
      INTEGER NA(2)
      DATA NA/n1,n2/
      Double Precision DHFINT, dhpolint

      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision Qscal 
      Integer Iqqbar
      Common/GLUDAT2/Qscal,Iqqbar 
      Integer i,j,k,in,irr,igrid,ip
      Double Precision scal,rx,rq2,rp,glu,up,down,str,usea,dsea,ssea,charm,csea,bot,bsea,top,tsea,phot
      Double Precision z0,wplus,wminus,higgs
      Double Precision RKMS6B,RKMS5B,RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
      Double Precision RKMS8,RKMS9,RKMS10,RKMS11
      Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
      common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      Real QCDLam
      character *72 TXT
      character adum
      LOGICAL FIRST
      DATA FIRST/.TRUE./
           
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      
      character mygridfile_old*132
      
      Logical lwrite/.true./
      
      Integer inter
      data inter/1/
      
c      write(6,*) ' in quark file '
c      write(6,*) ' iTMDgrid first ',first,first_meoffsh
c         write(6,*) ' iTMDgridq: files ',mygridfiles,mygridfile_old
      IF(first) THEN
      if(mygridfiles_u.ne.'mygridfile_old') then
c         write(6,*) ' iTMDq: files ',mygridfiles,mygridfile_old
         mygridfile_old=mygridfiles_u
         n1min = 0
         n1max = 0
         n2min = 0
         n2max = 0
         n3min = 0
         n3max = 0
         ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
         ns_sel = -1
c        if(lwrite)  write(6,*) ' iTMDgridq ',first,iset
           if(lwrite) write(6,*) ' iTMDgridq_u reading: ',trim(mygridfiles_u)
           open(30,FILE=mygridfiles_u, FORM='formatted',
     +          STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
c         if(lwrite) write(6,*) ' iTMDgridq: read  ',trim(mygridfiles_u)                
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx') 
200      Read(30,101) TXT
  101    Format(A72)
c         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
            read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
            read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
            read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
            read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
C         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif

c         if(lwrite) write(6,*) ' iTMDgridq: starting scale Q0 = ',Qg0
         if(lwrite) write(6,'(''  iTMDgridq_u: starting scale Q0 = '',F8.3)'),Qg0
c         if(lwrite) write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
c         if(lwrite) write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,' ns_sel = ',ns_sel
c         if(lwrite) write(6,*) ' QCD_lam used in uPDF: ',QCDlam
c         if(lwrite) write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar

            do i=1,n1
               do k=1,n2
cc                  write(6,*) j,i,k,RX,RQ2,RP,RKMS
cc                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS
cc                  write(6,*) ' in iTMDgrid '


c                  READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11

                  if(ns_sel.ge.1000) then 
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
                  else
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6
                  endif 

                  xx(i) = rx
                  px(k) = rp
                  f_grid0(i,k) = max(0.d0,rkms0)*scal
                  f_grid1(i,k) = max(0.d0,rkms1)*scal
                  f_grid2(i,k) = max(0.d0,rkms2)*scal
                  f_grid3(i,k) = max(0.d0,rkms3)*scal
                  f_grid4(i,k) = max(0.d0,rkms4)*scal
                  f_grid5(i,k) = max(0.d0,rkms5)*scal
                  f_grid6(i,k) = max(0.d0,rkms6)*scal
                  f_grid7(i,k) = max(0.d0,rkms7)*scal
                  f_grid8(i,k) = max(0.d0,rkms8)*scal
                  f_grid9(i,k) = max(0.d0,rkms9)*scal
                  f_grid10(i,k) = max(0.d0,rkms10)*scal
                  f_grid11(i,k) = max(0.d0,rkms11)*scal
                  f_grid1m(i,k) = max(0.d0,rkms1b)*scal
                  f_grid2m(i,k) = max(0.d0,rkms2b)*scal
                  f_grid3m(i,k) = max(0.d0,rkms3b)*scal
                  f_grid4m(i,k) = max(0.d0,rkms4b)*scal
                  f_grid5m(i,k) = max(0.d0,rkms5b)*scal
                  f_grid6m(i,k) = max(0.d0,rkms6b)*scal
c                  if(f_grid3m(i,k).gt.0) write(6,*) 'iTMDgrid ',f_grid3m(i,k),i,k
c                   write(6,*) rx,rp
               enddo
            enddo

c         if(lwrite) write(6,*) ' end of file at ',i,k
         IN=0
         DO I=1,N1
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         FIRST=.FALSE.
c         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
      endif
      ENDIF

      ncall = ncall +1   


      XA(1) = log(X)
      XA(2) = log(P)
      if(xa(2).lt.px(1)) then
         n2min = n2min + 1
         
         if(n2min.lt.10) then 
            write(6,*) ' iTMDgridq_u: p out of range ',p,' min p ',exp(px(1))
            write(6,*) ' iTMDgridq_u: p out of range ',xa(2),px(1)
         elseif(n2min.eq.10) then 
c            write(6,*) ' iTMDgridq_u: p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
         endif
         xa(2)=px(1)
      endif
      if(xa(2).gt.px(n2)) then
         n2max = n2max + 1
         if(n2max.lt.10) then 
            write(6,*) ' iTMDgridqu_u: p out of range ',p,' max p ',exp(px(n3))
         elseif(n2max.eq.10) then 
c            write(6,*) ' iTMDgridqu_u: p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
         endif
         xa(2)=px(n2)
      endif
      if(xa(1).ge.xx(n1)) xa(1)=xx(n1)-0.0001
      if(xa(1).lt.xx(1)) then
         n1min = n1min + 1
         if(n1min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
         elseif(n1min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
         endif
         xa(1)=xx(1)
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(inter.eq.0) then
         glu  = DHFINT(2,XA,NA,A,f_grid0)
         up   = DHFINT(2,XA,NA,A,f_grid2)
         down = DHFINT(2,XA,NA,A,f_grid1)
         str  = DHFINT(2,XA,NA,A,f_grid3)
         usea = DHFINT(2,XA,NA,A,f_grid2m)
         dsea = DHFINT(2,XA,NA,A,f_grid1m)
         ssea = DHFINT(2,XA,NA,A,f_grid3m)
         charm= DHFINT(2,XA,NA,A,f_grid4)
         csea = DHFINT(2,XA,NA,A,f_grid4m)
         bot  = DHFINT(2,XA,NA,A,f_grid5)
         bsea = DHFINT(2,XA,NA,A,f_grid5m)
         top  = DHFINT(2,XA,NA,A,f_grid6)
         tsea = DHFINT(2,XA,NA,A,f_grid6m)
         phot = DHFINT(2,XA,NA,A,f_grid7)
         z0 =   DHFINT(2,XA,NA,A,f_grid8)
         wplus =DHFINT(2,XA,NA,A,f_grid9)
         wminus=DHFINT(2,XA,NA,A,f_grid10)
         higgs =DHFINT(2,XA,NA,A,f_grid11)
      elseif(inter.eq.1) then
         glu  = dhpolint(2,XA,NA,A,f_grid0)
         up   = dhpolint(2,XA,NA,A,f_grid2)
         down = dhpolint(2,XA,NA,A,f_grid1)
         str  = dhpolint(2,XA,NA,A,f_grid3)
         usea = dhpolint(2,XA,NA,A,f_grid2m)
         dsea = dhpolint(2,XA,NA,A,f_grid1m)
         ssea = dhpolint(2,XA,NA,A,f_grid3m)
         charm= dhpolint(2,XA,NA,A,f_grid4)
         csea = dhpolint(2,XA,NA,A,f_grid4m)
         bot  = dhpolint(2,XA,NA,A,f_grid5)
         bsea = dhpolint(2,XA,NA,A,f_grid5m)
         top  = dhpolint(2,XA,NA,A,f_grid6)
         tsea = dhpolint(2,XA,NA,A,f_grid6m)
         phot = dhpolint(2,XA,NA,A,f_grid7)
         z0 =   dhpolint(2,XA,NA,A,f_grid8)
         wplus =dhpolint(2,XA,NA,A,f_grid9)
         wminus=dhpolint(2,XA,NA,A,f_grid10)
         higgs =dhpolint(2,XA,NA,A,f_grid11)
      else
        write(6,*) ' iTMDgridq_u inter = ',inter,' not implemented '
      endif
      DO  IP=-6,11
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(up.gt.0.) xpq(2) = up
      if(down.gt.0.) xpq(1) = down
      if(usea.gt.0) xpq(-2) = usea
      if(dsea.gt.0) xpq(-1) = dsea
      if(str.gt.0.) xpq(3) = str
      if(ssea.gt.0.) xpq(-3) = ssea
      if(charm.gt.0.) xpq(4) = charm
      if(csea.gt.0.) xpq(-4) = csea
      if(bot.gt.0.) xpq(5) = bot
      if(bsea.gt.0.) xpq(-5) = bsea
      if(top.gt.0.) xpq(6) = top
      if(tsea.gt.0.) xpq(-6) = tsea
      if(phot.gt.0.) xpq(7) = phot
      if(z0.gt.0.) xpq(8) = z0
      if(wplus.gt.0.) xpq(9) = wplus
      if(wminus.gt.0.) xpq(10) = wminus
      if(higgs.gt.0.) xpq(11) = higgs
c      if(z0.gt.0) then
c         write(6,*) ' iTMDgridq_u x,p,xpq',x,p,xpq(8)
ccc         write(6,*) ' i,j,k ',i,j,k
c      endif
c      do i=-6,7
c        xpq(i) = xpq(i)/(1.-min(x,0.999))
c      end do   
         
c         write(6,*) ' new iTMDgridu x,q2,p,xpq',x,q2,p,glu
      return
   50 write(6,*) ' iTMDgridq_u: error in opening file ', mygridfiles_u
      stop
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS0
      stop
      END

      SUBROUTINE iTMDgridq_d(X,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      double precision XPQ(-6:11),X,P
*! new
      Parameter (n1=51,n2=51)
      Double Precision xx,px
      DIMENSION XX(0:n1),PX(0:n2)

      Double Precision XA(3),A(N1+N2)
      Double Precision f_grid0(n1,n2),f_grid1(n1,n2),f_grid2(n1,n2),f_grid3(n1,n2),f_grid4(n1,n2),f_grid5(n1,n2),f_grid6(n1,n2)
      Double Precision f_grid7(n1,n2),f_grid8(n1,n2),f_grid9(n1,n2),f_grid10(n1,n2),f_grid11(n1,n2)
      Double Precision f_grid1m(n1,n2),f_grid2m(n1,n2),f_grid3m(n1,n2),f_grid4m(n1,n2),f_grid5m(n1,n2),f_grid6m(n1,n2)
      INTEGER NA(2)
      DATA NA/n1,n2/
      Double Precision DHFINT, dhpolint

      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision Qscal 
      Integer Iqqbar
      Common/GLUDAT2/Qscal,Iqqbar 
      Integer i,j,k,in,irr,igrid,ip
      Double Precision scal,rx,rq2,rp,glu,up,down,str,usea,dsea,ssea,charm,csea,bot,bsea,top,tsea,phot
      Double Precision z0,wplus,wminus,higgs
      Double Precision RKMS6B,RKMS5B,RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
      Double Precision RKMS8,RKMS9,RKMS10,RKMS11
      Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
      common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      Real QCDLam
      character *72 TXT
      character adum
      LOGICAL FIRST
      DATA FIRST/.TRUE./
           
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      
      
      character mygridfile_old*132
      
      Logical lwrite/.true./
      
      Integer inter
      data inter/1/
      
c      write(6,*) ' in quark file '
c      write(6,*) ' iTMDgrid first ',first,first_meoffsh
c         write(6,*) ' iTMDgridq: files ',mygridfiles,mygridfile_old
      IF(first) THEN
      if(mygridfiles_d.ne.'mygridfile_old') then
c         write(6,*) ' iTMDq: files ',mygridfiles,mygridfile_old
         mygridfile_old=mygridfiles_d
         n1min = 0
         n1max = 0
         n2min = 0
         n2max = 0
         n3min = 0
         n3max = 0
         ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
         ns_sel = -1
c        if(lwrite)  write(6,*) ' iTMDgridq ',first,iset
           if(lwrite) write(6,*) ' iTMDgridq_d reading: ',trim(mygridfiles_d)
           open(30,FILE=mygridfiles_d, FORM='formatted',
     +          STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
c         if(lwrite) write(6,*) ' iTMDgridq: read  ',trim(mygridfiles_u)                
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx') 
200      Read(30,101) TXT
  101    Format(A72)
c         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
            read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
            read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
            read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
            read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
C         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif

c         if(lwrite) write(6,*) ' iTMDgridq: starting scale Q0 = ',Qg0
         if(lwrite) write(6,'(''  iTMDgridq_d: starting scale Q0 = '',F8.3)'),Qg0
c         if(lwrite) write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
c         if(lwrite) write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,' ns_sel = ',ns_sel
c         if(lwrite) write(6,*) ' QCD_lam used in uPDF: ',QCDlam
c         if(lwrite) write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar

            do i=1,n1
               do k=1,n2
cc                  write(6,*) j,i,k,RX,RQ2,RP,RKMS
cc                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS
cc                  write(6,*) ' in iTMDgrid '


c                  READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11

                  if(ns_sel.ge.1000) then 
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
                  else
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6
                  endif 

                  xx(i) = rx
                  px(k) = rp
                  f_grid0(i,k) = max(0.d0,rkms0)*scal
                  f_grid1(i,k) = max(0.d0,rkms1)*scal
                  f_grid2(i,k) = max(0.d0,rkms2)*scal
                  f_grid3(i,k) = max(0.d0,rkms3)*scal
                  f_grid4(i,k) = max(0.d0,rkms4)*scal
                  f_grid5(i,k) = max(0.d0,rkms5)*scal
                  f_grid6(i,k) = max(0.d0,rkms6)*scal
                  f_grid7(i,k) = max(0.d0,rkms7)*scal
                  f_grid8(i,k) = max(0.d0,rkms8)*scal
                  f_grid9(i,k) = max(0.d0,rkms9)*scal
                  f_grid10(i,k) = max(0.d0,rkms10)*scal
                  f_grid11(i,k) = max(0.d0,rkms11)*scal
                  f_grid1m(i,k) = max(0.d0,rkms1b)*scal
                  f_grid2m(i,k) = max(0.d0,rkms2b)*scal
                  f_grid3m(i,k) = max(0.d0,rkms3b)*scal
                  f_grid4m(i,k) = max(0.d0,rkms4b)*scal
                  f_grid5m(i,k) = max(0.d0,rkms5b)*scal
                  f_grid6m(i,k) = max(0.d0,rkms6b)*scal
c                  if(f_grid3m(i,k).gt.0) write(6,*) 'iTMDgrid ',f_grid3m(i,k),i,k
c                   write(6,*) rx,rp
               enddo
            enddo

c         if(lwrite) write(6,*) ' end of file at ',i,k
         IN=0
         DO I=1,N1
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         FIRST=.FALSE.
c         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
      endif
      ENDIF

      ncall = ncall +1   


      XA(1) = log(X)
      XA(2) = log(P)
      if(xa(2).lt.px(1)) then
         n2min = n2min + 1
         
         if(n2min.lt.10) then 
            write(6,*) ' iTMDgridq_d: p out of range ',p,' min p ',exp(px(1))
            write(6,*) ' iTMDgridq_d: p out of range ',xa(2),px(1)
         elseif(n2min.eq.10) then 
c            write(6,*) ' iTMDgridq_d: p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
         endif
         xa(2)=px(1)
      endif
      if(xa(2).gt.px(n2)) then
         n2max = n2max + 1
         if(n2max.lt.10) then 
            write(6,*) ' iTMDgridqu_u: p out of range ',p,' max p ',exp(px(n3))
         elseif(n2max.eq.10) then 
c            write(6,*) ' iTMDgridqu_u: p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
         endif
         xa(2)=px(n2)
      endif
      if(xa(1).ge.xx(n1)) xa(1)=xx(n1)-0.0001
      if(xa(1).lt.xx(1)) then
         n1min = n1min + 1
         if(n1min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
         elseif(n1min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
         endif
         xa(1)=xx(1)
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(inter.eq.0) then
         glu  = DHFINT(2,XA,NA,A,f_grid0)
         up   = DHFINT(2,XA,NA,A,f_grid2)
         down = DHFINT(2,XA,NA,A,f_grid1)
         str  = DHFINT(2,XA,NA,A,f_grid3)
         usea = DHFINT(2,XA,NA,A,f_grid2m)
         dsea = DHFINT(2,XA,NA,A,f_grid1m)
         ssea = DHFINT(2,XA,NA,A,f_grid3m)
         charm= DHFINT(2,XA,NA,A,f_grid4)
         csea = DHFINT(2,XA,NA,A,f_grid4m)
         bot  = DHFINT(2,XA,NA,A,f_grid5)
         bsea = DHFINT(2,XA,NA,A,f_grid5m)
         top  = DHFINT(2,XA,NA,A,f_grid6)
         tsea = DHFINT(2,XA,NA,A,f_grid6m)
         phot = DHFINT(2,XA,NA,A,f_grid7)
         z0 =   DHFINT(2,XA,NA,A,f_grid8)
         wplus =DHFINT(2,XA,NA,A,f_grid9)
         wminus=DHFINT(2,XA,NA,A,f_grid10)
         higgs =DHFINT(2,XA,NA,A,f_grid11)
      elseif(inter.eq.1) then
         glu  = dhpolint(2,XA,NA,A,f_grid0)
         up   = dhpolint(2,XA,NA,A,f_grid2)
         down = dhpolint(2,XA,NA,A,f_grid1)
         str  = dhpolint(2,XA,NA,A,f_grid3)
         usea = dhpolint(2,XA,NA,A,f_grid2m)
         dsea = dhpolint(2,XA,NA,A,f_grid1m)
         ssea = dhpolint(2,XA,NA,A,f_grid3m)
         charm= dhpolint(2,XA,NA,A,f_grid4)
         csea = dhpolint(2,XA,NA,A,f_grid4m)
         bot  = dhpolint(2,XA,NA,A,f_grid5)
         bsea = dhpolint(2,XA,NA,A,f_grid5m)
         top  = dhpolint(2,XA,NA,A,f_grid6)
         tsea = dhpolint(2,XA,NA,A,f_grid6m)
         phot = dhpolint(2,XA,NA,A,f_grid7)
         z0 =   dhpolint(2,XA,NA,A,f_grid8)
         wplus =dhpolint(2,XA,NA,A,f_grid9)
         wminus=dhpolint(2,XA,NA,A,f_grid10)
         higgs =dhpolint(2,XA,NA,A,f_grid11)
      else
        write(6,*) ' iTMDgridq_d inter = ',inter,' not implemented '
      endif
      DO  IP=-6,11
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(up.gt.0.) xpq(2) = up
      if(down.gt.0.) xpq(1) = down
      if(usea.gt.0) xpq(-2) = usea
      if(dsea.gt.0) xpq(-1) = dsea
      if(str.gt.0.) xpq(3) = str
      if(ssea.gt.0.) xpq(-3) = ssea
      if(charm.gt.0.) xpq(4) = charm
      if(csea.gt.0.) xpq(-4) = csea
      if(bot.gt.0.) xpq(5) = bot
      if(bsea.gt.0.) xpq(-5) = bsea
      if(top.gt.0.) xpq(6) = top
      if(tsea.gt.0.) xpq(-6) = tsea
      if(phot.gt.0.) xpq(7) = phot
      if(z0.gt.0.) xpq(8) = z0
      if(wplus.gt.0.) xpq(9) = wplus
      if(wminus.gt.0.) xpq(10) = wminus
      if(higgs.gt.0.) xpq(11) = higgs
c      if(z0.gt.0) then
c         write(6,*) ' iTMDgridq_d x,p,xpq',x,p,xpq(8)
ccc         write(6,*) ' i,j,k ',i,j,k
c      endif
c      do i=-6,7
c        xpq(i) = xpq(i)/(1.-min(x,0.999))
c      end do   
         
c         write(6,*) ' new iTMDgridu x,q2,p,xpq',x,q2,p,glu
      return
   50 write(6,*) ' iTMDgridq_d: error in opening file ', mygridfiles_d
      stop
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS0
      stop
      END

      SUBROUTINE iTMDgridq_dbar(X,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      double precision XPQ(-6:11),X,P
*! new
      Parameter (n1=51,n2=51)
      Double Precision xx,px
      DIMENSION XX(0:n1),PX(0:n2)

      Double Precision XA(3),A(N1+N2)
      Double Precision f_grid0(n1,n2),f_grid1(n1,n2),f_grid2(n1,n2),f_grid3(n1,n2),f_grid4(n1,n2),f_grid5(n1,n2),f_grid6(n1,n2)
      Double Precision f_grid7(n1,n2),f_grid8(n1,n2),f_grid9(n1,n2),f_grid10(n1,n2),f_grid11(n1,n2)
      Double Precision f_grid1m(n1,n2),f_grid2m(n1,n2),f_grid3m(n1,n2),f_grid4m(n1,n2),f_grid5m(n1,n2),f_grid6m(n1,n2)
      INTEGER NA(2)
      DATA NA/n1,n2/
      Double Precision DHFINT, dhpolint

      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision Qscal 
      Integer Iqqbar
      Common/GLUDAT2/Qscal,Iqqbar 
      Integer i,j,k,in,irr,igrid,ip
      Double Precision scal,rx,rq2,rp,glu,up,down,str,usea,dsea,ssea,charm,csea,bot,bsea,top,tsea,phot
      Double Precision z0,wplus,wminus,higgs
      Double Precision RKMS6B,RKMS5B,RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
      Double Precision RKMS8,RKMS9,RKMS10,RKMS11
      Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
      common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      Real QCDLam
      character *72 TXT
      character adum
      LOGICAL FIRST
      DATA FIRST/.TRUE./
           
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      
      
      character mygridfile_old*132
      
      Logical lwrite/.true./
      
      Integer inter
      data inter/1/
      
c      write(6,*) ' in quark file '
c      write(6,*) ' iTMDgrid first ',first,first_meoffsh
c         write(6,*) ' iTMDgridq: files ',mygridfiles,mygridfile_old
      IF(first) THEN
      if(mygridfiles_dbar.ne.'mygridfile_old') then
c         write(6,*) ' iTMDq: files ',mygridfiles,mygridfile_old
         mygridfile_old=mygridfiles_dbar
         n1min = 0
         n1max = 0
         n2min = 0
         n2max = 0
         n3min = 0
         n3max = 0
         ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
         ns_sel = -1
c        if(lwrite)  write(6,*) ' iTMDgridq ',first,iset
           if(lwrite) write(6,*) ' iTMDgridq_dbar reading: ',trim(mygridfiles_dbar)
           open(30,FILE=mygridfiles_dbar, FORM='formatted',
     +          STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
c         if(lwrite) write(6,*) ' iTMDgridq: read  ',trim(mygridfiles_u)                
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx') 
200      Read(30,101) TXT
  101    Format(A72)
c         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
            read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
            read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
            read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
            read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
C         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif

c         if(lwrite) write(6,*) ' iTMDgridq: starting scale Q0 = ',Qg0
         if(lwrite) write(6,'(''  iTMDgridq_dbar: starting scale Q0 = '',F8.3)'),Qg0
c         if(lwrite) write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
c         if(lwrite) write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,' ns_sel = ',ns_sel
c         if(lwrite) write(6,*) ' QCD_lam used in uPDF: ',QCDlam
c         if(lwrite) write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar

            do i=1,n1
               do k=1,n2
cc                  write(6,*) j,i,k,RX,RQ2,RP,RKMS
cc                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS
cc                  write(6,*) ' in iTMDgrid '


c                  READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11

                  if(ns_sel.ge.1000) then 
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
                  else
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6
                  endif 

                  xx(i) = rx
                  px(k) = rp
                  f_grid0(i,k) = max(0.d0,rkms0)*scal
                  f_grid1(i,k) = max(0.d0,rkms1)*scal
                  f_grid2(i,k) = max(0.d0,rkms2)*scal
                  f_grid3(i,k) = max(0.d0,rkms3)*scal
                  f_grid4(i,k) = max(0.d0,rkms4)*scal
                  f_grid5(i,k) = max(0.d0,rkms5)*scal
                  f_grid6(i,k) = max(0.d0,rkms6)*scal
                  f_grid7(i,k) = max(0.d0,rkms7)*scal
                  f_grid8(i,k) = max(0.d0,rkms8)*scal
                  f_grid9(i,k) = max(0.d0,rkms9)*scal
                  f_grid10(i,k) = max(0.d0,rkms10)*scal
                  f_grid11(i,k) = max(0.d0,rkms11)*scal
                  f_grid1m(i,k) = max(0.d0,rkms1b)*scal
                  f_grid2m(i,k) = max(0.d0,rkms2b)*scal
                  f_grid3m(i,k) = max(0.d0,rkms3b)*scal
                  f_grid4m(i,k) = max(0.d0,rkms4b)*scal
                  f_grid5m(i,k) = max(0.d0,rkms5b)*scal
                  f_grid6m(i,k) = max(0.d0,rkms6b)*scal
c                  if(f_grid3m(i,k).gt.0) write(6,*) 'iTMDgrid ',f_grid3m(i,k),i,k
c                   write(6,*) rx,rp
               enddo
            enddo

c         if(lwrite) write(6,*) ' end of file at ',i,k
         IN=0
         DO I=1,N1
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         FIRST=.FALSE.
c         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
      endif
      ENDIF

      ncall = ncall +1   


      XA(1) = log(X)
      XA(2) = log(P)
      if(xa(2).lt.px(1)) then
         n2min = n2min + 1
         
         if(n2min.lt.10) then 
            write(6,*) ' iTMDgridq_dbar: p out of range ',p,' min p ',exp(px(1))
            write(6,*) ' iTMDgridq_dbar: p out of range ',xa(2),px(1)
         elseif(n2min.eq.10) then 
c            write(6,*) ' iTMDgridq_dbar: p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
         endif
         xa(2)=px(1)
      endif
      if(xa(2).gt.px(n2)) then
         n2max = n2max + 1
         if(n2max.lt.10) then 
            write(6,*) ' iTMDgridqu_u: p out of range ',p,' max p ',exp(px(n3))
         elseif(n2max.eq.10) then 
c            write(6,*) ' iTMDgridqu_u: p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
         endif
         xa(2)=px(n2)
      endif
      if(xa(1).ge.xx(n1)) xa(1)=xx(n1)-0.0001
      if(xa(1).lt.xx(1)) then
         n1min = n1min + 1
         if(n1min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
         elseif(n1min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
         endif
         xa(1)=xx(1)
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(inter.eq.0) then
         glu  = DHFINT(2,XA,NA,A,f_grid0)
         up   = DHFINT(2,XA,NA,A,f_grid2)
         down = DHFINT(2,XA,NA,A,f_grid1)
         str  = DHFINT(2,XA,NA,A,f_grid3)
         usea = DHFINT(2,XA,NA,A,f_grid2m)
         dsea = DHFINT(2,XA,NA,A,f_grid1m)
         ssea = DHFINT(2,XA,NA,A,f_grid3m)
         charm= DHFINT(2,XA,NA,A,f_grid4)
         csea = DHFINT(2,XA,NA,A,f_grid4m)
         bot  = DHFINT(2,XA,NA,A,f_grid5)
         bsea = DHFINT(2,XA,NA,A,f_grid5m)
         top  = DHFINT(2,XA,NA,A,f_grid6)
         tsea = DHFINT(2,XA,NA,A,f_grid6m)
         phot = DHFINT(2,XA,NA,A,f_grid7)
         z0 =   DHFINT(2,XA,NA,A,f_grid8)
         wplus =DHFINT(2,XA,NA,A,f_grid9)
         wminus=DHFINT(2,XA,NA,A,f_grid10)
         higgs =DHFINT(2,XA,NA,A,f_grid11)
      elseif(inter.eq.1) then
         glu  = dhpolint(2,XA,NA,A,f_grid0)
         up   = dhpolint(2,XA,NA,A,f_grid2)
         down = dhpolint(2,XA,NA,A,f_grid1)
         str  = dhpolint(2,XA,NA,A,f_grid3)
         usea = dhpolint(2,XA,NA,A,f_grid2m)
         dsea = dhpolint(2,XA,NA,A,f_grid1m)
         ssea = dhpolint(2,XA,NA,A,f_grid3m)
         charm= dhpolint(2,XA,NA,A,f_grid4)
         csea = dhpolint(2,XA,NA,A,f_grid4m)
         bot  = dhpolint(2,XA,NA,A,f_grid5)
         bsea = dhpolint(2,XA,NA,A,f_grid5m)
         top  = dhpolint(2,XA,NA,A,f_grid6)
         tsea = dhpolint(2,XA,NA,A,f_grid6m)
         phot = dhpolint(2,XA,NA,A,f_grid7)
         z0 =   dhpolint(2,XA,NA,A,f_grid8)
         wplus =dhpolint(2,XA,NA,A,f_grid9)
         wminus=dhpolint(2,XA,NA,A,f_grid10)
         higgs =dhpolint(2,XA,NA,A,f_grid11)
      else
        write(6,*) ' iTMDgridq_dbar inter = ',inter,' not implemented '
      endif
      DO  IP=-6,11
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(up.gt.0.) xpq(2) = up
      if(down.gt.0.) xpq(1) = down
      if(usea.gt.0) xpq(-2) = usea
      if(dsea.gt.0) xpq(-1) = dsea
      if(str.gt.0.) xpq(3) = str
      if(ssea.gt.0.) xpq(-3) = ssea
      if(charm.gt.0.) xpq(4) = charm
      if(csea.gt.0.) xpq(-4) = csea
      if(bot.gt.0.) xpq(5) = bot
      if(bsea.gt.0.) xpq(-5) = bsea
      if(top.gt.0.) xpq(6) = top
      if(tsea.gt.0.) xpq(-6) = tsea
      if(phot.gt.0.) xpq(7) = phot
      if(z0.gt.0.) xpq(8) = z0
      if(wplus.gt.0.) xpq(9) = wplus
      if(wminus.gt.0.) xpq(10) = wminus
      if(higgs.gt.0.) xpq(11) = higgs
c      if(z0.gt.0) then
c         write(6,*) ' iTMDgridq_dbar x,p,xpq',x,p,xpq(8)
ccc         write(6,*) ' i,j,k ',i,j,k
c      endif
c      do i=-6,7
c        xpq(i) = xpq(i)/(1.-min(x,0.999))
c      end do   
         
c         write(6,*) ' new iTMDgridu x,q2,p,xpq',x,q2,p,glu
      return
   50 write(6,*) ' iTMDgridq_dbar: error in opening file ', mygridfiles_dbar
      stop
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS0
      stop
      END


      SUBROUTINE iTMDgridq_ubar(X,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      double precision XPQ(-6:11),X,P
*! new
      Parameter (n1=51,n2=51)
      Double Precision xx,px
      DIMENSION XX(0:n1),PX(0:n2)

      Double Precision XA(3),A(N1+N2)
      Double Precision f_grid0(n1,n2),f_grid1(n1,n2),f_grid2(n1,n2),f_grid3(n1,n2),f_grid4(n1,n2),f_grid5(n1,n2),f_grid6(n1,n2)
      Double Precision f_grid7(n1,n2),f_grid8(n1,n2),f_grid9(n1,n2),f_grid10(n1,n2),f_grid11(n1,n2)
      Double Precision f_grid1m(n1,n2),f_grid2m(n1,n2),f_grid3m(n1,n2),f_grid4m(n1,n2),f_grid5m(n1,n2),f_grid6m(n1,n2)
      INTEGER NA(2)
      DATA NA/n1,n2/
      Double Precision DHFINT, dhpolint

      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision Qscal 
      Integer Iqqbar
      Common/GLUDAT2/Qscal,Iqqbar 
      Integer i,j,k,in,irr,igrid,ip
      Double Precision scal,rx,rq2,rp,glu,up,down,str,usea,dsea,ssea,charm,csea,bot,bsea,top,tsea,phot
      Double Precision z0,wplus,wminus,higgs
      Double Precision RKMS6B,RKMS5B,RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
      Double Precision RKMS8,RKMS9,RKMS10,RKMS11
      Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
      common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      Real QCDLam
      character *72 TXT
      character adum
      LOGICAL FIRST
      DATA FIRST/.TRUE./
           
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      
      character mygridfile_old*132
      
      Logical lwrite/.true./
      
      Integer inter
      data inter/1/
      
c      write(6,*) ' in quark file '
c      write(6,*) ' iTMDgrid first ',first,first_meoffsh
c         write(6,*) ' iTMDgridq: files ',mygridfiles,mygridfile_old
      IF(first) THEN
      if(mygridfiles_ubar.ne.'mygridfile_old') then
c         write(6,*) ' iTMDq: files ',mygridfiles,mygridfile_old
         mygridfile_old=mygridfiles_ubar
         n1min = 0
         n1max = 0
         n2min = 0
         n2max = 0
         n3min = 0
         n3max = 0
         ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
         ns_sel = -1
c        if(lwrite)  write(6,*) ' iTMDgridq ',first,iset
           if(lwrite) write(6,*) ' iTMDgridq_ubar reading: ',trim(mygridfiles_ubar)
           open(30,FILE=mygridfiles_ubar, FORM='formatted',
     +          STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
c         if(lwrite) write(6,*) ' iTMDgridq: read  ',trim(mygridfiles_u)                
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx') 
200      Read(30,101) TXT
  101    Format(A72)
c         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
            read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
            read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
            read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
            read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
C         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif

c         if(lwrite) write(6,*) ' iTMDgridq: starting scale Q0 = ',Qg0
         if(lwrite) write(6,'(''  iTMDgridq_ubar: starting scale Q0 = '',F8.3)'),Qg0
c         if(lwrite) write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
c         if(lwrite) write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,' ns_sel = ',ns_sel
c         if(lwrite) write(6,*) ' QCD_lam used in uPDF: ',QCDlam
c         if(lwrite) write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar

            do i=1,n1
               do k=1,n2
cc                  write(6,*) j,i,k,RX,RQ2,RP,RKMS
cc                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS
cc                  write(6,*) ' in iTMDgrid '


c                  READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11

                  if(ns_sel.ge.1000) then 
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
                  else
                     READ(30,*,Err=90 ) RX,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6
                  endif 

                  xx(i) = rx
                  px(k) = rp
                  f_grid0(i,k) = max(0.d0,rkms0)*scal
                  f_grid1(i,k) = max(0.d0,rkms1)*scal
                  f_grid2(i,k) = max(0.d0,rkms2)*scal
                  f_grid3(i,k) = max(0.d0,rkms3)*scal
                  f_grid4(i,k) = max(0.d0,rkms4)*scal
                  f_grid5(i,k) = max(0.d0,rkms5)*scal
                  f_grid6(i,k) = max(0.d0,rkms6)*scal
                  f_grid7(i,k) = max(0.d0,rkms7)*scal
                  f_grid8(i,k) = max(0.d0,rkms8)*scal
                  f_grid9(i,k) = max(0.d0,rkms9)*scal
                  f_grid10(i,k) = max(0.d0,rkms10)*scal
                  f_grid11(i,k) = max(0.d0,rkms11)*scal
                  f_grid1m(i,k) = max(0.d0,rkms1b)*scal
                  f_grid2m(i,k) = max(0.d0,rkms2b)*scal
                  f_grid3m(i,k) = max(0.d0,rkms3b)*scal
                  f_grid4m(i,k) = max(0.d0,rkms4b)*scal
                  f_grid5m(i,k) = max(0.d0,rkms5b)*scal
                  f_grid6m(i,k) = max(0.d0,rkms6b)*scal
c                  if(f_grid3m(i,k).gt.0) write(6,*) 'iTMDgrid ',f_grid3m(i,k),i,k
c                   write(6,*) rx,rp
               enddo
            enddo

c         if(lwrite) write(6,*) ' end of file at ',i,k
         IN=0
         DO I=1,N1
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         FIRST=.FALSE.
c         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
      endif
      ENDIF

      ncall = ncall +1   


      XA(1) = log(X)
      XA(2) = log(P)
      if(xa(2).lt.px(1)) then
         n2min = n2min + 1
         
         if(n2min.lt.10) then 
            write(6,*) ' iTMDgridq_ubar: p out of range ',p,' min p ',exp(px(1))
            write(6,*) ' iTMDgridq_ubar: p out of range ',xa(2),px(1)
         elseif(n2min.eq.10) then 
c            write(6,*) ' iTMDgridq_ubar: p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
         endif
         xa(2)=px(1)
      endif
      if(xa(2).gt.px(n2)) then
         n2max = n2max + 1
         if(n2max.lt.10) then 
            write(6,*) ' iTMDgridqu_u: p out of range ',p,' max p ',exp(px(n3))
         elseif(n2max.eq.10) then 
c            write(6,*) ' iTMDgridqu_u: p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
         endif
         xa(2)=px(n2)
      endif
      if(xa(1).ge.xx(n1)) xa(1)=xx(n1)-0.0001
      if(xa(1).lt.xx(1)) then
         n1min = n1min + 1
         if(n1min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
         elseif(n1min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
         endif
         xa(1)=xx(1)
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(inter.eq.0) then
         glu  = DHFINT(2,XA,NA,A,f_grid0)
         up   = DHFINT(2,XA,NA,A,f_grid2)
         down = DHFINT(2,XA,NA,A,f_grid1)
         str  = DHFINT(2,XA,NA,A,f_grid3)
         usea = DHFINT(2,XA,NA,A,f_grid2m)
         dsea = DHFINT(2,XA,NA,A,f_grid1m)
         ssea = DHFINT(2,XA,NA,A,f_grid3m)
         charm= DHFINT(2,XA,NA,A,f_grid4)
         csea = DHFINT(2,XA,NA,A,f_grid4m)
         bot  = DHFINT(2,XA,NA,A,f_grid5)
         bsea = DHFINT(2,XA,NA,A,f_grid5m)
         top  = DHFINT(2,XA,NA,A,f_grid6)
         tsea = DHFINT(2,XA,NA,A,f_grid6m)
         phot = DHFINT(2,XA,NA,A,f_grid7)
         z0 =   DHFINT(2,XA,NA,A,f_grid8)
         wplus =DHFINT(2,XA,NA,A,f_grid9)
         wminus=DHFINT(2,XA,NA,A,f_grid10)
         higgs =DHFINT(2,XA,NA,A,f_grid11)
      elseif(inter.eq.1) then
         glu  = dhpolint(2,XA,NA,A,f_grid0)
         up   = dhpolint(2,XA,NA,A,f_grid2)
         down = dhpolint(2,XA,NA,A,f_grid1)
         str  = dhpolint(2,XA,NA,A,f_grid3)
         usea = dhpolint(2,XA,NA,A,f_grid2m)
         dsea = dhpolint(2,XA,NA,A,f_grid1m)
         ssea = dhpolint(2,XA,NA,A,f_grid3m)
         charm= dhpolint(2,XA,NA,A,f_grid4)
         csea = dhpolint(2,XA,NA,A,f_grid4m)
         bot  = dhpolint(2,XA,NA,A,f_grid5)
         bsea = dhpolint(2,XA,NA,A,f_grid5m)
         top  = dhpolint(2,XA,NA,A,f_grid6)
         tsea = dhpolint(2,XA,NA,A,f_grid6m)
         phot = dhpolint(2,XA,NA,A,f_grid7)
         z0 =   dhpolint(2,XA,NA,A,f_grid8)
         wplus =dhpolint(2,XA,NA,A,f_grid9)
         wminus=dhpolint(2,XA,NA,A,f_grid10)
         higgs =dhpolint(2,XA,NA,A,f_grid11)
      else
        write(6,*) ' iTMDgridq_ubar inter = ',inter,' not implemented '
      endif
      DO  IP=-6,11
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(up.gt.0.) xpq(2) = up
      if(down.gt.0.) xpq(1) = down
      if(usea.gt.0) xpq(-2) = usea
      if(dsea.gt.0) xpq(-1) = dsea
      if(str.gt.0.) xpq(3) = str
      if(ssea.gt.0.) xpq(-3) = ssea
      if(charm.gt.0.) xpq(4) = charm
      if(csea.gt.0.) xpq(-4) = csea
      if(bot.gt.0.) xpq(5) = bot
      if(bsea.gt.0.) xpq(-5) = bsea
      if(top.gt.0.) xpq(6) = top
      if(tsea.gt.0.) xpq(-6) = tsea
      if(phot.gt.0.) xpq(7) = phot
      if(z0.gt.0.) xpq(8) = z0
      if(wplus.gt.0.) xpq(9) = wplus
      if(wminus.gt.0.) xpq(10) = wminus
      if(higgs.gt.0.) xpq(11) = higgs
c      if(z0.gt.0) then
c         write(6,*) ' iTMDgridq_ubar x,p,xpq',x,p,xpq(8)
ccc         write(6,*) ' i,j,k ',i,j,k
c      endif
c      do i=-6,7
c        xpq(i) = xpq(i)/(1.-min(x,0.999))
c      end do   
         
c         write(6,*) ' new iTMDgridu x,q2,p,xpq',x,q2,p,glu
      return
   50 write(6,*) ' iTMDgridq_ubar: error in opening file ', mygridfiles_ubar
      stop
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS0
      stop
      END






*
* Kernlib
*
*
          FUNCTION DHFINT(NARG,ARG,NENT,ENT,TABLE)
          Implicit None
C
C   INTERPOLATION ROUTINE. AUTHOR C. LETERTRE.
C   MODIFIED BY B. SCHORR, 1.07.1982.
C

          INTEGER   NENT(500)
          Double Precision     ARG(500),   ENT(500),TABLE(1000000)
          INTEGER   INDEX(32)
          Double Precision       WEIGHT(32)
          Double Precision DHFINT
          Integer NARG,LMAX,ISTEP,KNOTS,N,NDIM,ISHIFT,I,K
          Integer LMIN,LOCA,LOCB,LOCC,LGFILE
          Double Precision X,H,ETA
          LOGICAL   MFLAG,    RFLAG
          DHFINT  =  0.
          IF(NARG .LT. 1  .OR.  NARG .GT. 5)  GOTO 300
          LMAX      =  0
          ISTEP     =  1
          KNOTS     =  1
          INDEX(1)  =  1
          WEIGHT(1) =  1.
          DO 100    N  =  1, NARG
             X     =  ARG(N)
             NDIM  =  NENT(N)
             LOCA  =  LMAX
             LMIN  =  LMAX + 1
             LMAX  =  LMAX + NDIM
             IF(NDIM .GT. 2)  GOTO 10
             IF(NDIM .EQ. 1)  GOTO 100
             H  =  X - ENT(LMIN)
             IF(H .EQ. 0.)  GOTO 90
             ISHIFT  =  ISTEP
             IF(X-ENT(LMIN+1) .EQ. 0.)  GOTO 21
             ISHIFT  =  0
             ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
             GOTO 30
  10         LOCB  =  LMAX + 1
  11         LOCC  =  (LOCA+LOCB) / 2
             IF(X-ENT(LOCC))  12, 20, 13
  12         LOCB  =  LOCC
             GOTO 14
  13         LOCA  =  LOCC
  14         IF(LOCB-LOCA .GT. 1)  GOTO 11
             LOCA    =  MIN0( MAX0(LOCA,LMIN), LMAX-1 )
             ISHIFT  =  (LOCA - LMIN) * ISTEP
             ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
             GOTO 30
  20         ISHIFT  =  (LOCC - LMIN) * ISTEP
  21         DO 22  K  =  1, KNOTS
                INDEX(K)  =  INDEX(K) + ISHIFT
  22            CONTINUE
             GOTO 90
  30         DO 31  K  =  1, KNOTS
                INDEX(K)         =  INDEX(K) + ISHIFT
                INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
                WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
                WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
  31            CONTINUE
             KNOTS  =  2*KNOTS
  90         ISTEP  =  ISTEP * NDIM
 100         CONTINUE
          DO 200    K  =  1, KNOTS
             I  =  INDEX(K)
             DHFINT  =  DHFINT + WEIGHT(K) * TABLE(I)
 200         CONTINUE
          RETURN
 300      Continue
          IF(MFLAG) THEN
             IF(LGFILE .EQ. 0) THEN
                WRITE(*,1000) NARG
             ELSE
                WRITE(LGFILE,1000) NARG
             ENDIF
          ENDIF
          IF(.NOT. RFLAG) CALL ABEND
          RETURN
1000      FORMAT('  HFUNCTION DHFINT ... NARG =',I6,
     +              '  NOT WITHIN RANGE')
          END
      
      
      SUBROUTINE TMDgridg(X,Q2,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      double precision XPQ(-6:11),X,Q2,P
*! new
      Parameter (n1=51,n2=51,n3=51)
      double precision Q2x,xx,px
      DIMENSION Q2X(0:n1),XX(0:n2),PX(0:n3)
      Double Precision XA(3),A(N1+N2+N3)
      Double Precision f_grid0(n1,n2,n3),f_grid1(n1,n2,n3),f_grid2(n1,n2,n3),
     & f_grid3(n1,n2,n3),f_grid4(n1,n2,n3),f_grid5(n1,n2,n3),f_grid6(n1,n2,n3),f_grid7(n1,n2,n3)
      Double Precision f_grid1m(n1,n2,n3),f_grid2m(n1,n2,n3),
     & f_grid3m(n1,n2,n3),f_grid4m(n1,n2,n3),f_grid5m(n1,n2,n3),f_grid6m(n1,n2,n3)
      double precision f_grid8(n1,n2,n3),f_grid9(n1,n2,n3),f_grid10(n1,n2,n3),f_grid11(n1,n2,n3)
      INTEGER NA(3)
      DATA NA/n1,n2,n3/
      Double Precision DHFINT
      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision Qscal
      Integer Iqqbar
      Common/GLUDAT2/Qscal,Iqqbar
      Integer i,j,k,in,irr,igrid,ip
      Double Precision scal,rx,rq2,rp,glu,up,usea,down,dsea,str,ssea,charm,csea,bot,bsea,top,tsea,phot
      double precision z0,wplus,wminus,higgs
      Double Precision RKMS6B,RKMS5B,RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
      double precision RKMS8,RKMS9,RKMS10,RKMS11
      Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
      common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      character *72 TXT
      character adum
      Real QCDlam
      LOGICAL FIRST
      DATA FIRST/.TRUE./

            
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      
      character mygridfile_old*132
      
      Logical lwrite/.true./
     

c      write(6,*) ' TMDg: files ',mygridfiles,p
      IF(FIRST) THEN
      if(mygridfiles_d.ne.'mygridfile_old') then
         write(6,*) ' TMDg: files ',mygridfiles_d,mygridfile_old
         mygridfile_old=mygridfiles_d
         n1min = 0
         n1max = 0
         n2min = 0
         n2max = 0
         n3min = 0
         n3max = 0
         ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
         ns_sel = -1
        if(lwrite)  write(6,*) ' TMDgridg ',first 
           if(lwrite) write(6,*) ' TMDgridg reading: ',mygridfiles_d
           open(30,FILE=mygridfiles_d, FORM='formatted',
     +          STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
         if(lwrite) write(6,*) ' TMDgridg: read ',mygridfiles_d                 
200      Read(30,101) TXT
  101    Format(A72)
C         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
c	      read(txt,1000) Qg0,ikincut
c1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
            read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
c	      read(txt,1000) Qg0,ikincut
            read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
c	      read(txt,1001) Ipgg,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
cc         WRITE(6,101) '2nd line ',TXT
           read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
c            PARU(112)=QCDLam
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
c	      read(txt,1002) Qscal,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
C         WRITE(6,101) '2nd line ',TXT
            read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif

c         write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
c         write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,
c     &	   ' ns_sel = ',ns_sel
c         write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar

         do j=1,n1
            do i=1,n2
               do k=1,n3
c                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
c                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6 

                  if(ns_sel.ge.1000) then 
                     READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
                  else
                     READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6 
                  endif 
                  xx(i) = rx
                  q2x(j) = rq2
                  px(k) = rp
                  f_grid0(j,i,k) = max(0.d0,rkms0)*scal
                  f_grid1(j,i,k) = max(0.d0,rkms1)*scal
                  f_grid2(j,i,k) = max(0.d0,rkms2)*scal
                  f_grid3(j,i,k) = max(0.d0,rkms3)*scal
                  f_grid4(j,i,k) = max(0.d0,rkms4)*scal
                  f_grid5(j,i,k) = max(0.d0,rkms5)*scal
                  f_grid6(j,i,k) = max(0.d0,rkms6)*scal
                  f_grid7(j,i,k) = max(0.d0,rkms7)*scal
                  f_grid8(j,i,k) = max(0.d0,rkms8)*scal
                  f_grid9(j,i,k) = max(0.d0,rkms9)*scal
                  f_grid10(j,i,k) = max(0.d0,rkms10)*scal
                  f_grid11(j,i,k) = max(0.d0,rkms11)*scal
                  f_grid1m(j,i,k) = max(0.d0,rkms1b)*scal
                  f_grid2m(j,i,k) = max(0.d0,rkms2b)*scal
                  f_grid3m(j,i,k) = max(0.d0,rkms3b)*scal
                  f_grid4m(j,i,k) = max(0.d0,rkms4b)*scal
                  f_grid5m(j,i,k) = max(0.d0,rkms5b)*scal
                  f_grid6m(j,i,k) = max(0.d0,rkms6b)*scal
c                  write(6,*)f_grid(j,i,k) 
               enddo
            enddo
         enddo
         write(6,*) ' TMDgridg: after reading ',mygridfiles_d,j,i,k
c         write(6,*) ' after reading ccfm-kernel.dat ',j,i,k      
         IN=0
         DO I=1,n1
            IN=IN+1
            A(IN) = q2x(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N3
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         FIRST=.FALSE.
c         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
      endif
      ENDIF

      ncall = ncall +1 


      XA(1) = log(Q2)
      XA(2) = log(X)
      XA(3) = log(P)
      if(xa(3).lt.px(1)) then
         n3min = n3min + 1

         if(n3min.lt.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' p out of range ',xa(3),px(1)
         elseif(n3min.eq.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
         endif
         xa(3)=px(1)
      endif
      if(xa(3).gt.px(n3)) then
         n3max = n3max + 1
         if(n3max.lt.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
         elseif(n3max.eq.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
         endif
         xa(3)=px(n3)
      endif
      if(xa(2).ge.xx(n2)) xa(2)=xx(n2)-0.0001
      if(xa(2).lt.xx(1)) then
         n2min = n2min + 1
         if(n2min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
         elseif(n2min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
         endif
         xa(2)=xx(1)
      endif 
      if(xa(1).lt.q2x(1)) then
         n1min = n1min + 1
         if(n1min.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
         elseif(n1min.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
c            write(6,*) ' last message printed: min k2'
         endif
         xa(1)=q2x(1)
      endif
      if(xa(1).ge.q2x(n1)) then
         n1max = n1max + 1
         if(n1max.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
         elseif(n1max.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
c            write(6,*) ' last message printed: max k2'
         endif
         xa(1)=q2x(n1)-0.1
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(igrid.eq.1) then
         if(xa(1).lt.q2x(1)) then
            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
            xa(1)=q2x(1)
         endif
         if(xa(1).ge.q2x(n1)) then
            write(6,*) '  k2 out of range: x = ', x,n1,q2x(n1)
            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
            xa(1)=q2x(n1)-0.1
         endif
         i=0
   20    i=i+1
         if(xa(1).gt.a(na(1))) then
c       write(6,*) ' q2  ',xa(1),a(na(1))
            write(6,*) ' q2 not found ',q2,a(na(1)),q2x(n1),xa(1)
            i=na(1)
         else
            if(xa(1).ge.A(i).and.xa(1).lt.a(i+1)) Then
            else
               if(i.le.na(1)) then
                  goto 20
               else
                  write(6,*) ' at 20: q2 not found ',i,q2
               endif
            endif
         endif
         j=0
   30    j=j+1
         if(xa(2).ge.A(na(1)+j).and.xa(2).lt.a(na(1)+j+1)) Then
         else
            if(j.le.na(2)) then
               goto 30
            else
               write(6,*) ' at 30: x not found ',x,xa(2),j
            endif
         endif
         k=0
   40    k=k+1
         if(xa(3).ge.a(na(1)+na(2)+na(3))) then
            k=na(3)
c       write(6,*) ' p  ',xa(3),a(na(1)+na(2)+na(3))
         else
            if(xa(3).ge.A(na(1)+na(2)+k).and. xa(3).lt.a(na(1)+na(2)+k+1)) Then
            else
               if(k.le.na(3)) then
                  goto 40
               else
                  write(6,*) ' at 40: p not found ',k,p
               endif
            endif
         endif

         glu  = DHFINT(3,XA,NA,A,f_grid0)/q2
         up = DHFINT(3,XA,NA,A,f_grid2)/q2
         down = DHFINT(3,XA,NA,A,f_grid1)/q2
         str  = DHFINT(3,XA,NA,A,f_grid3)/q2
         usea = DHFINT(3,XA,NA,A,f_grid2m)/q2
         dsea = DHFINT(3,XA,NA,A,f_grid1m)/q2
         ssea = DHFINT(3,XA,NA,A,f_grid3m)/q2
         charm= DHFINT(3,XA,NA,A,f_grid4)/q2
         csea = DHFINT(3,XA,NA,A,f_grid4m)/q2
         bot= DHFINT(3,XA,NA,A,f_grid5)/q2
         bsea= DHFINT(3,XA,NA,A,f_grid5m)/q2
         top= DHFINT(3,XA,NA,A,f_grid6)/q2
         tsea= DHFINT(3,XA,NA,A,f_grid6m)/q2
         phot= DHFINT(3,XA,NA,A,f_grid7)/q2
         z0= DHFINT(3,XA,NA,A,f_grid8)/q2
         wplus= DHFINT(3,XA,NA,A,f_grid9)/q2
         wminus= DHFINT(3,XA,NA,A,f_grid10)/q2
         higgs= DHFINT(3,XA,NA,A,f_grid11)/q2
         
      else
         glu  = DHFINT(3,XA,NA,A,f_grid0)/q2
         up = DHFINT(3,XA,NA,A,f_grid2)/q2
         down = DHFINT(3,XA,NA,A,f_grid1)/q2
         str  = DHFINT(3,XA,NA,A,f_grid3)/q2
         usea = DHFINT(3,XA,NA,A,f_grid2m)/q2
         dsea = DHFINT(3,XA,NA,A,f_grid1m)/q2
         ssea = DHFINT(3,XA,NA,A,f_grid3m)/q2
         charm= DHFINT(3,XA,NA,A,f_grid4)/q2
         csea = DHFINT(3,XA,NA,A,f_grid4m)/q2
         bot= DHFINT(3,XA,NA,A,f_grid5)/q2
         bsea= DHFINT(3,XA,NA,A,f_grid5m)/q2
         top= DHFINT(3,XA,NA,A,f_grid6)/q2
         tsea= DHFINT(3,XA,NA,A,f_grid6m)/q2
         phot= DHFINT(3,XA,NA,A,f_grid7)/q2
         z0= DHFINT(3,XA,NA,A,f_grid8)/q2
         wplus= DHFINT(3,XA,NA,A,f_grid9)/q2
         wminus= DHFINT(3,XA,NA,A,f_grid10)/q2
         higgs= DHFINT(3,XA,NA,A,f_grid11)/q2
      endif

c           write(6,*) ' new TMDgridg x,q2,p,xpq',x,q2,p,glu
    
      DO  IP=-6,11
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(up.gt.0.) xpq(2) = up
      if(usea.gt.0.) xpq(-2) = usea
      if(down.gt.0.) xpq(1) = down
      if(dsea.gt.0.) xpq(-1) = dsea
      if(str.gt.0.) xpq(3) = str
      if(ssea.gt.0.) xpq(-3) = ssea
      if(charm.gt.0.) xpq(4) = charm
      if(csea.gt.0.) xpq(-4) = csea
      if(bot.gt.0.) xpq(5) = bot
      if(bsea.gt.0.) xpq(-5) = bsea
      if(top.gt.0.) xpq(6) = top
      if(tsea.gt.0.) xpq(-6) = tsea
      if(phot.gt.0.) xpq(7) = phot
      if(z0.gt.0.) xpq(8) = z0
      if(wplus.gt.0.) xpq(9) = wplus
      if(wminus.gt.0.) xpq(10) = wminus
      if(higgs.gt.0.) xpq(11) = higgs
      if(glu.gt.5E6) then
         write(6,*) ' new TMDgridg gluon x,q2,p,xpq',x,q2,p,glu
c         write(6,*) ' i,j,k ',i,j,k
      endif
      if(ssea.gt.5E6) then
         write(6,*) ' new TMDgridg sea x,q2,p,xpq',x,q2,p,ssea
c         write(6,*) ' i,j,k ',i,j,k
      endif
      return


   50 write(6,*) ' TMDgridg: error in opening file '
      stop
   90 continue
      write(6,*) ' TMDgridg: end of file at ',i,j,k,RX,RQ2,RP,RKMS1,RKMS2,RKMS3
      stop
      END

      SUBROUTINE TMDgridq(X,Q2,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      double precision XPQ(-6:11),X,Q2,P
*! new
      Parameter (n1=51,n2=51,n3=51)

      double precision Q2x,xx,px
      DIMENSION Q2X(0:n1),XX(0:n2),PX(0:n3)
      Double Precision XA(3),A(N1+N2+N3)
      Double Precision f_grid0(n1,n2,n3),f_grid1(n1,n2,n3),f_grid2(n1,n2,n3),
     & f_grid3(n1,n2,n3),f_grid4(n1,n2,n3),f_grid5(n1,n2,n3),f_grid6(n1,n2,n3),f_grid7(n1,n2,n3)
      Double Precision f_grid1m(n1,n2,n3),f_grid2m(n1,n2,n3),
     & f_grid3m(n1,n2,n3),f_grid4m(n1,n2,n3),f_grid5m(n1,n2,n3),f_grid6m(n1,n2,n3)
      Double Precision f_grid8(n1,n2,n3),f_grid9(n1,n2,n3),f_grid10(n1,n2,n3),f_grid11(n1,n2,n3)
      INTEGER NA(3)
      DATA NA/n1,n2,n3/
      Double Precision DHFINT
      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision Qscal
      Integer Iqqbar
      Common/GLUDAT2/Qscal,Iqqbar
      Integer i,j,k,in,irr,igrid,ip

      Double Precision scal,rx,rq2,rp,glu,up,usea,down,dsea,str,ssea,charm,csea,bot,bsea,top,tsea,phot
      Double Precision z0,wplus,wminus,higgs
      Double Precision RKMS6B,RKMS5B,RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
      Double Precision RKMS8,RKMS9,RKMS10,RKMS11
      Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
      common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      character *72 TXT
      character adum
      Real QCDlam
      LOGICAL FIRST
      DATA FIRST/.TRUE./

            
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      
      character mygridfile_old*132
      
      Logical lwrite/.true./

c         write(6,*) ' TMDgridq: files ',mygridfiles,mygridfile_old

      IF(FIRST) THEN
      if(mygridfiles_d.ne.'mygridfile_old') then
         write(6,*) ' TMDgridqd: files ',mygridfiles_d,mygridfile_old
         mygridfile_old=mygridfiles_d
         n1min = 0
         n1max = 0
         n2min = 0
         n2max = 0
         n3min = 0
         n3max = 0
         ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
         ns_sel = -1
         if(lwrite)  write(6,*) ' TMDgridq_d ',first 
           if(lwrite) write(6,*) ' TMDgridq_d reading: ',mygridfiles_d
           open(30,FILE=mygridfiles_d, FORM='formatted',
     +          STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
         if(lwrite) write(6,*) ' TMDgridq_d: read ',mygridfiles_d                 
200      Read(30,101) TXT
  101    Format(A72)
C         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
c	      read(txt,1000) Qg0,ikincut
c1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
            read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
c	      read(txt,1000) Qg0,ikincut
            read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
c	      read(txt,1001) Ipgg,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
cc         WRITE(6,101) '2nd line ',TXT
            read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
c            PARU(112)=QCDLam
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
c	      read(txt,1002) Qscal,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
C         WRITE(6,101) '2nd line ',TXT
            read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif


         do j=1,n1
            do i=1,n2
               do k=1,n3
c                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
c                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6
                  if(ns_sel.ge.1000) then 
                     READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
                  else
                     READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6 
                  endif 
                  xx(i) = rx
                  q2x(j) = rq2
                  px(k) = rp
                  f_grid0(j,i,k) = max(0.d0,rkms0)*scal
                  f_grid1(j,i,k) = max(0.d0,rkms1)*scal
                  f_grid2(j,i,k) = max(0.d0,rkms2)*scal
                  f_grid3(j,i,k) = max(0.d0,rkms3)*scal
                  f_grid4(j,i,k) = max(0.d0,rkms4)*scal
                  f_grid5(j,i,k) = max(0.d0,rkms5)*scal
                  f_grid6(j,i,k) = max(0.d0,rkms6)*scal
                  f_grid7(j,i,k) = max(0.d0,rkms7)*scal
                  f_grid8(j,i,k) = max(0.d0,rkms8)*scal
                  f_grid9(j,i,k) = max(0.d0,rkms9)*scal
                  f_grid10(j,i,k) = max(0.d0,rkms10)*scal
                  f_grid11(j,i,k) = max(0.d0,rkms11)*scal
                  f_grid1m(j,i,k) = max(0.d0,rkms1b)*scal
                  f_grid2m(j,i,k) = max(0.d0,rkms2b)*scal
                  f_grid3m(j,i,k) = max(0.d0,rkms3b)*scal
                  f_grid4m(j,i,k) = max(0.d0,rkms4b)*scal
                  f_grid5m(j,i,k) = max(0.d0,rkms5b)*scal
                  f_grid6m(j,i,k) = max(0.d0,rkms6b)*scal
               enddo
            enddo
         enddo
         write(6,*) ' TMDgridq_d: after reading ',mygridfiles_d,j,i,k
         IN=0
         DO I=1,n1
            IN=IN+1
            A(IN) = q2x(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N3
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         FIRST=.FALSE.
c         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
      endif
      ENDIF

      ncall = ncall +1 


      XA(1) = log(Q2)
      XA(2) = log(X)
      XA(3) = log(P)
      if(xa(3).lt.px(1)) then
         n3min = n3min + 1

         if(n3min.lt.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' p out of range ',xa(3),px(1)
         elseif(n3min.eq.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
         endif
         xa(3)=px(1)
      endif
      if(xa(3).gt.px(n3)) then
         n3max = n3max + 1
         if(n3max.lt.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
         elseif(n3max.eq.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
         endif
         xa(3)=px(n3)
      endif
      if(xa(2).ge.xx(n2)) xa(2)=xx(n2)-0.0001
      if(xa(2).lt.xx(1)) then
         n2min = n2min + 1
         if(n2min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
         elseif(n2min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
         endif
         xa(2)=xx(1)
      endif 
      if(xa(1).lt.q2x(1)) then
         n1min = n1min + 1
         if(n1min.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
         elseif(n1min.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
c            write(6,*) ' last message printed: min k2'
         endif
         xa(1)=q2x(1)
      endif
      if(xa(1).ge.q2x(n1)) then
         n1max = n1max + 1
         if(n1max.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
         elseif(n1max.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
c            write(6,*) ' last message printed: max k2'
         endif
         xa(1)=q2x(n1)-0.1
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(igrid.eq.1) then
         if(xa(1).lt.q2x(1)) then
            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
            xa(1)=q2x(1)
         endif
         if(xa(1).ge.q2x(n1)) then
            write(6,*) '  k2 out of range: x = ', x,n1,q2x(n1)
            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
            xa(1)=q2x(n1)-0.1
         endif
         i=0
   20    i=i+1
         if(xa(1).gt.a(na(1))) then
c       write(6,*) ' q2  ',xa(1),a(na(1))
            write(6,*) ' q2 not found ',q2,a(na(1)),q2x(n1),xa(1)
            i=na(1)
         else
            if(xa(1).ge.A(i).and.xa(1).lt.a(i+1)) Then
            else
               if(i.le.na(1)) then
                  goto 20
               else
                  write(6,*) ' at 20: q2 not found ',i,q2
               endif
            endif
         endif
         j=0
   30    j=j+1
         if(xa(2).ge.A(na(1)+j).and.xa(2).lt.a(na(1)+j+1)) Then
         else
            if(j.le.na(2)) then
               goto 30
            else
               write(6,*) ' at 30: x not found ',x,xa(2),j
            endif
         endif
         k=0
   40    k=k+1
         if(xa(3).ge.a(na(1)+na(2)+na(3))) then
            k=na(3)
c       write(6,*) ' p  ',xa(3),a(na(1)+na(2)+na(3))
         else
            if(xa(3).ge.A(na(1)+na(2)+k).and. xa(3).lt.a(na(1)+na(2)+k+1)) Then
            else
               if(k.le.na(3)) then
                  goto 40
               else
                  write(6,*) ' at 40: p not found ',k,p
               endif
            endif
         endif

         glu  = DHFINT(3,XA,NA,A,f_grid0)/q2
         up = DHFINT(3,XA,NA,A,f_grid2)/q2
         down = DHFINT(3,XA,NA,A,f_grid1)/q2
         str  = DHFINT(3,XA,NA,A,f_grid3)/q2
         usea = DHFINT(3,XA,NA,A,f_grid2m)/q2
         dsea = DHFINT(3,XA,NA,A,f_grid1m)/q2
         ssea = DHFINT(3,XA,NA,A,f_grid3m)/q2
         charm= DHFINT(3,XA,NA,A,f_grid4)/q2
         csea = DHFINT(3,XA,NA,A,f_grid4m)/q2
         bot= DHFINT(3,XA,NA,A,f_grid5)/q2
         bsea= DHFINT(3,XA,NA,A,f_grid5m)/q2
         top= DHFINT(3,XA,NA,A,f_grid6)/q2
         tsea= DHFINT(3,XA,NA,A,f_grid6m)/q2
         phot= DHFINT(3,XA,NA,A,f_grid7)/q2
         z0= DHFINT(3,XA,NA,A,f_grid8)/q2
         wplus= DHFINT(3,XA,NA,A,f_grid9)/q2
         wminus= DHFINT(3,XA,NA,A,f_grid10)/q2
         higgs= DHFINT(3,XA,NA,A,f_grid11)/q2
         
      else
         glu  = DHFINT(3,XA,NA,A,f_grid0)/q2
         up = DHFINT(3,XA,NA,A,f_grid2)/q2
         down = DHFINT(3,XA,NA,A,f_grid1)/q2
         str  = DHFINT(3,XA,NA,A,f_grid3)/q2
         usea = DHFINT(3,XA,NA,A,f_grid2m)/q2
         dsea = DHFINT(3,XA,NA,A,f_grid1m)/q2
         ssea = DHFINT(3,XA,NA,A,f_grid3m)/q2
         charm= DHFINT(3,XA,NA,A,f_grid4)/q2
         csea = DHFINT(3,XA,NA,A,f_grid4m)/q2
         bot= DHFINT(3,XA,NA,A,f_grid5)/q2
         bsea= DHFINT(3,XA,NA,A,f_grid5m)/q2
         top= DHFINT(3,XA,NA,A,f_grid6)/q2
         tsea= DHFINT(3,XA,NA,A,f_grid6m)/q2
         phot= DHFINT(3,XA,NA,A,f_grid7)/q2
         z0= DHFINT(3,XA,NA,A,f_grid8)/q2
         wplus= DHFINT(3,XA,NA,A,f_grid9)/q2
         wminus= DHFINT(3,XA,NA,A,f_grid10)/q2
         higgs= DHFINT(3,XA,NA,A,f_grid11)/q2
c         write(6,*) ' new TMDgridqd x,q2,p ',x,q2,p,glu,uval,dval,usea,dsea,ssea
      endif

    
      DO  IP=-6,11
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(up.gt.0.) xpq(2) = up
      if(usea.gt.0.) xpq(-2) = usea
      if(down.gt.0.) xpq(1) = down
      if(dsea.gt.0.) xpq(-1) = dsea
      if(str.gt.0.) xpq(3) = str
      if(ssea.gt.0.) xpq(-3) = ssea
      if(charm.gt.0.) xpq(4) = charm
      if(csea.gt.0.) xpq(-4) = csea
      if(bot.gt.0.) xpq(5) = bot
      if(bsea.gt.0.) xpq(-5) = bsea
      if(top.gt.0.) xpq(6) = top
      if(tsea.gt.0.) xpq(-6) = tsea
      if(phot.gt.0.) xpq(7) = phot
      if(z0.gt.0.) xpq(8) = z0
      if(wplus.gt.0.) xpq(9) = wplus
      if(wminus.gt.0.) xpq(10) = wminus
      if(higgs.gt.0.) xpq(11) = higgs
      if(glu.gt.5E6) then
         write(6,*) ' new TMDgridq_d gluon x,q2,p,xpq',x,q2,p,glu
c         write(6,*) ' i,j,k ',i,j,k
      endif
      if(ssea.gt.5E6) then
         write(6,*) ' new TMDgridq_d sea x,q2,p,xpq',x,q2,p,ssea
c         write(6,*) ' i,j,k ',i,j,k
      endif
      return


   50 write(6,*) ' TMDgridq_d: error in opening file '
      stop
   90 continue
      write(6,*) ' TMDgridq_d: end of file at ',i,j,k,RX,RQ2,RP,RKMS1,RKMS2,RKMS3
      stop
      END

      SUBROUTINE TMDgridq_u(X,Q2,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      double precision XPQ(-6:11),X,Q2,P
*! new
      Parameter (n1=51,n2=51,n3=51)

      double precision Q2x,xx,px
      DIMENSION Q2X(0:n1),XX(0:n2),PX(0:n3)
      Double Precision XA(3),A(N1+N2+N3)
      Double Precision f_grid0(n1,n2,n3),f_grid1(n1,n2,n3),f_grid2(n1,n2,n3),
     & f_grid3(n1,n2,n3),f_grid4(n1,n2,n3),f_grid5(n1,n2,n3),f_grid6(n1,n2,n3),f_grid7(n1,n2,n3)
      Double Precision f_grid1m(n1,n2,n3),f_grid2m(n1,n2,n3),
     & f_grid3m(n1,n2,n3),f_grid4m(n1,n2,n3),f_grid5m(n1,n2,n3),f_grid6m(n1,n2,n3)
      Double Precision f_grid8(n1,n2,n3),f_grid9(n1,n2,n3),f_grid10(n1,n2,n3),f_grid11(n1,n2,n3)
      INTEGER NA(3)
      DATA NA/n1,n2,n3/
      Double Precision DHFINT
      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision Qscal
      Integer Iqqbar
      Common/GLUDAT2/Qscal,Iqqbar
      Integer i,j,k,in,irr,igrid,ip

      Double Precision scal,rx,rq2,rp,glu,up,usea,down,dsea,str,ssea,charm,csea,bot,bsea,top,tsea,phot
      Double Precision z0,wplus,wminus,higgs
      Double Precision RKMS6B,RKMS5B,RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7
      Double Precision RKMS8,RKMS9,RKMS10,RKMS11
      Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
      common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      character *72 TXT
      character adum
      Real QCDlam
      LOGICAL FIRST
      DATA FIRST/.TRUE./

            
      character mygridfiles_d*132 , mygridfiles_u*132, mygridfiles_ubar*132, mygridfiles_dbar*132 
      Common/ myfiles/ mygridfiles_d, mygridfiles_u, mygridfiles_ubar, mygridfiles_dbar
      
      character mygridfile_old*132
      
      Logical lwrite/.true./

c         write(6,*) ' TMDgridqu: files ',mygridfiles,mygridfile_old
c      goto 555
      IF(FIRST) THEN
      if(mygridfiles_u.ne.'mygridfile_old') then
         write(6,*) ' TMDgridq_u: files ',mygridfiles_u,mygridfile_old
         mygridfile_old=mygridfiles_u
         n1min = 0
         n1max = 0
         n2min = 0
         n2max = 0
         n3min = 0
         n3max = 0
         ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
         ns_sel = -1
         if(lwrite) then 
           write(6,*) ' TMDgridq_u ',first 
           write(6,*) ' TMDgridq_u reading: ',mygridfiles_u
         endif
         open(30,FILE=mygridfiles_u, FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
         if(lwrite) write(6,*) ' TMDgridq_u: read ',mygridfiles_u      
                             
200      Read(30,101) TXT
  101    Format(A72)
C         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
c	      read(txt,1000) Qg0,ikincut
c1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
            read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
c	      read(txt,1000) Qg0,ikincut
            read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
c	      read(txt,1001) Ipgg,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
cc         WRITE(6,101) '2nd line ',TXT
            read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
c            PARU(112)=QCDLam
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
c	      read(txt,1002) Qscal,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
C         WRITE(6,101) '2nd line ',TXT
	      read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif
        
         do j=1,n1
            do i=1,n2
               do k=1,n3
c                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
c                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
c     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6
                  if(ns_sel.ge.1000) then 
                     READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6,RKMS7,RKMS8,RKMS9,RKMS10,RKMS11
                  else
                     READ(30,*,Err=90 ) RX,RQ2,RP,RKMS6b,RKMS5b,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,
     &                               RKMS1,RKMS2,RKMS3,RKMS4,RKMS5,RKMS6 
                  endif 
                  xx(i) = rx
                  q2x(j) = rq2
                  px(k) = rp
                  f_grid0(j,i,k) = max(0.d0,rkms0)*scal
                  f_grid1(j,i,k) = max(0.d0,rkms1)*scal
                  f_grid2(j,i,k) = max(0.d0,rkms2)*scal
                  f_grid3(j,i,k) = max(0.d0,rkms3)*scal
                  f_grid4(j,i,k) = max(0.d0,rkms4)*scal
                  f_grid5(j,i,k) = max(0.d0,rkms5)*scal
                  f_grid6(j,i,k) = max(0.d0,rkms6)*scal
                  f_grid7(j,i,k) = max(0.d0,rkms7)*scal
                  f_grid8(j,i,k) = max(0.d0,rkms8)*scal
                  f_grid9(j,i,k) = max(0.d0,rkms9)*scal
                  f_grid10(j,i,k) = max(0.d0,rkms10)*scal
                  f_grid11(j,i,k) = max(0.d0,rkms11)*scal
                  f_grid1m(j,i,k) = max(0.d0,rkms1b)*scal
                  f_grid2m(j,i,k) = max(0.d0,rkms2b)*scal
                  f_grid3m(j,i,k) = max(0.d0,rkms3b)*scal
                  f_grid4m(j,i,k) = max(0.d0,rkms4b)*scal
                  f_grid5m(j,i,k) = max(0.d0,rkms5b)*scal
                  f_grid6m(j,i,k) = max(0.d0,rkms6b)*scal
               enddo
            enddo
         enddo
         write(6,*) ' TMDgridq_u: after reading ',mygridfiles_u,j,i,k
         IN=0
         DO I=1,n1
            IN=IN+1
            A(IN) = q2x(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N3
            IN=IN+1
            A(IN) = px(I)
         ENDDO
345      continue         
         FIRST=.FALSE.
c         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
         
      endif
      ENDIF

      ncall = ncall +1 
 
      XA(1) = log(Q2)
      XA(2) = log(X)
      XA(3) = log(P)
      if(xa(3).lt.px(1)) then
         n3min = n3min + 1
        
         if(n3min.lt.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' p out of range ',xa(3),px(1)
         elseif(n3min.eq.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
         endif
         xa(3)=px(1)
      endif
      if(xa(3).gt.px(n3)) then
         n3max = n3max + 1
         if(n3max.lt.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
         elseif(n3max.eq.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
         endif
         xa(3)=px(n3)
      endif
      if(xa(2).ge.xx(n2)) xa(2)=xx(n2)-0.0001
      if(xa(2).lt.xx(1)) then
         n2min = n2min + 1
         if(n2min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
         elseif(n2min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
         endif
         xa(2)=xx(1)
         endif 
      if(xa(1).lt.q2x(1)) then
         n1min = n1min + 1
         if(n1min.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
         elseif(n1min.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
c            write(6,*) ' last message printed: min k2'
         endif
         xa(1)=q2x(1)
      endif
      if(xa(1).ge.q2x(n1)) then
         n1max = n1max + 1
         if(n1max.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
         elseif(n1max.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
c            write(6,*) ' last message printed: max k2'
         endif
         xa(1)=q2x(n1)-0.1
      endif
c check if interpolation or grid wanted
      igrid = 0
      
      
      if(igrid.eq.1) then
         if(xa(1).lt.q2x(1)) then
            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
            xa(1)=q2x(1)
         endif
         if(xa(1).ge.q2x(n1)) then
            write(6,*) '  k2 out of range: x = ', x,n1,q2x(n1)
            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
            xa(1)=q2x(n1)-0.1
         endif
         i=0
   20    i=i+1
         if(xa(1).gt.a(na(1))) then
c       write(6,*) ' q2  ',xa(1),a(na(1))
            write(6,*) ' q2 not found ',q2,a(na(1)),q2x(n1),xa(1)
            i=na(1)
         else
            if(xa(1).ge.A(i).and.xa(1).lt.a(i+1)) Then
            else
               if(i.le.na(1)) then
                  goto 20
               else
                  write(6,*) ' at 20: q2 not found ',i,q2
               endif
            endif
         endif
         j=0
   30    j=j+1
         if(xa(2).ge.A(na(1)+j).and.xa(2).lt.a(na(1)+j+1)) Then
         else
            if(j.le.na(2)) then
               goto 30
            else
               write(6,*) ' at 30: x not found ',x,xa(2),j
            endif
         endif
         k=0
   40    k=k+1
         if(xa(3).ge.a(na(1)+na(2)+na(3))) then
            k=na(3)
c       write(6,*) ' p  ',xa(3),a(na(1)+na(2)+na(3))
         else
            if(xa(3).ge.A(na(1)+na(2)+k).and. xa(3).lt.a(na(1)+na(2)+k+1)) Then
            else
               if(k.le.na(3)) then
                  goto 40
               else
                  write(6,*) ' at 40: p not found ',k,p
               endif
            endif
         endif

         glu  = DHFINT(3,XA,NA,A,f_grid0)/q2
         up = DHFINT(3,XA,NA,A,f_grid2)/q2
         down = DHFINT(3,XA,NA,A,f_grid1)/q2
         str  = DHFINT(3,XA,NA,A,f_grid3)/q2
         usea = DHFINT(3,XA,NA,A,f_grid2m)/q2
         dsea = DHFINT(3,XA,NA,A,f_grid1m)/q2
         ssea = DHFINT(3,XA,NA,A,f_grid3m)/q2
         charm= DHFINT(3,XA,NA,A,f_grid4)/q2
         csea = DHFINT(3,XA,NA,A,f_grid4m)/q2
         bot= DHFINT(3,XA,NA,A,f_grid5)/q2
         bsea= DHFINT(3,XA,NA,A,f_grid5m)/q2
         top= DHFINT(3,XA,NA,A,f_grid6)/q2
         tsea= DHFINT(3,XA,NA,A,f_grid6m)/q2
         phot= DHFINT(3,XA,NA,A,f_grid7)/q2
         z0= DHFINT(3,XA,NA,A,f_grid8)/q2
         wplus= DHFINT(3,XA,NA,A,f_grid9)/q2
         wminus= DHFINT(3,XA,NA,A,f_grid10)/q2
         higgs= DHFINT(3,XA,NA,A,f_grid11)/q2
         
      else
         glu  = DHFINT(3,XA,NA,A,f_grid0)/q2
         up = DHFINT(3,XA,NA,A,f_grid2)/q2
         down = DHFINT(3,XA,NA,A,f_grid1)/q2
         str  = DHFINT(3,XA,NA,A,f_grid3)/q2
         usea = DHFINT(3,XA,NA,A,f_grid2m)/q2
         dsea = DHFINT(3,XA,NA,A,f_grid1m)/q2
         ssea = DHFINT(3,XA,NA,A,f_grid3m)/q2
         charm= DHFINT(3,XA,NA,A,f_grid4)/q2
         csea = DHFINT(3,XA,NA,A,f_grid4m)/q2
         bot= DHFINT(3,XA,NA,A,f_grid5)/q2
         bsea= DHFINT(3,XA,NA,A,f_grid5m)/q2
         top= DHFINT(3,XA,NA,A,f_grid6)/q2
         tsea= DHFINT(3,XA,NA,A,f_grid6m)/q2
         phot= DHFINT(3,XA,NA,A,f_grid7)/q2
         z0= DHFINT(3,XA,NA,A,f_grid8)/q2
         wplus= DHFINT(3,XA,NA,A,f_grid9)/q2
         wminus= DHFINT(3,XA,NA,A,f_grid10)/q2
         higgs= DHFINT(3,XA,NA,A,f_grid11)/q2
c         write(6,*) ' new TMDgridqu x,q2,p ',x,q2,p,glu,uval,dval,usea,dsea,ssea
      endif

   
      DO  IP=-6,11
         XPQ(IP)=0.0
      ENDDO
            
      if(glu.gt.0.) xpq(0) = glu
      if(up.gt.0.) xpq(2) = up
      if(usea.gt.0.) xpq(-2) = usea
      if(down.gt.0.) xpq(1) = down
      if(dsea.gt.0.) xpq(-1) = dsea
      if(str.gt.0.) xpq(3) = str
      if(ssea.gt.0.) xpq(-3) = ssea
      if(charm.gt.0.) xpq(4) = charm
      if(csea.gt.0.) xpq(-4) = csea
      if(bot.gt.0.) xpq(5) = bot
      if(bsea.gt.0.) xpq(-5) = bsea
      if(top.gt.0.) xpq(6) = top
      if(tsea.gt.0.) xpq(-6) = tsea
      if(phot.gt.0.) xpq(7) = phot
      if(z0.gt.0.) xpq(8) = z0
      if(wplus.gt.0.) xpq(9) = wplus
      if(wminus.gt.0.) xpq(10) = wminus
      if(higgs.gt.0.) xpq(11) = higgs
      if(glu.gt.5E6) then
         write(6,*) ' new TMDgridq_u gluon x,q2,p,xpq',x,q2,p,glu
c         write(6,*) ' i,j,k ',i,j,k
      endif
      if(ssea.gt.5E6) then
         write(6,*) ' new TMDgridq_u sea x,q2,p,xpq',x,q2,p,ssea
c         write(6,*) ' i,j,k ',i,j,k
      endif
555   return


   50 write(6,*) ' TMDgridq_u: error in opening file '
      stop
   90 continue
      write(6,*) ' TMDgridqu: end of file at ',i,j,k,RX,RQ2,RP,RKMS1,RKMS2,RKMS3
      stop
      END


      Subroutine iTMDgrids(x,q2,xpq)
      Implicit none
      Real x,q2,xpq(-6:6)
      Integer i
      if(x.gt.0.9) then
        do i=-6,6
        xpq(i) = x/0.1
        end do
        else
         do i=-6,6
        xpq(i) = 0
        end do
      endif
      Return
      end

 
      
      
      
      function dhpolint(NARG,ARG,NENT,ENT,TABLE)
      Implicit None
C
C
      Integer NARG
      INTEGER   NENT(500)
      Double Precision     ARG(500),   ENT(100000)
c      Double Precision TABLE(NENT(1),NENT(2),NENT(3))
      Double Precision TABLE(NENT(1),NENT(2))
c      Double Precision TABLE(1000000)
      Double Precision DHPOLINT
      Integer N,NDIM,LMAX,LMIN
      Double Precision X
      Integer N1,N2,N3
      Integer jj(narg)
      Parameter (N1=3, N2=3, N3=3)
      Double precision x1a(N1),x2a(N2),x3a(N3),ya(N1,N2)
      Integer i,j,i1,j1,nn,ii,ji
      Double Precision y,dy
           
     
      dhpolint = 0.
      LMAX = 0
      
      DO N  =  1, NARG
         X     =  ARG(N)
         NDIM  = NENT(N)
c         write(6,*) ' x, ndim ',x,ndim,n
         LMIN  = LMAX + 1
         LMAX  = LMAX + NDIM
         call locate_sm(ENT(LMIN),NDIM,X,jj(N))
c         write(6,*) ' after locate_sm ', x,jj(n)
c         stop
         if(n.eq.1) then 
           nn = n1
           elseif(n.eq.2) then
           nn = n2
           elseif(n.eq.3) then
           nn = n3
           endif 
c         write(6,*) ' check ',nn,jj(n)
         i1=0
         do i=jj(n),jj(n)+nn-1
           j=lmin-1 + i
           i1=i1+1
c           write(6,*) ' index i,i1,j = ',i,i1,j,lmin
           if(n.eq.1) then 
           x1a(i1) = ent(j)
           elseif(n.eq.2) then
           x2a(i1) = ent(j)
           elseif(n.eq.3) then
           x3a(i1) = ent(j)
           endif 
         end do 
      End Do
      i1=0
c      write(6,*) ' jj ', jj(1),jj(1)+n1-1,jj(2),jj(2)+n2-1
      do i=jj(1),jj(1)+n1-1
        i1=i1+1
        j1=0
        ii=i
c        if(i.ge.51) ii=51
        if(i.ge.NENT(1)) ii=NENT(1)
        do j=jj(2),jj(2)+n2-1
          j1=j1+1
          ji=j
c          if(j.ge.51) ji=51
          if(j.ge.NENT(2)) ji=NENT(2)
c          write(6,*) ' table points ',ii,ji,i1,j1,n1,n2
          ya(i1,j1) = table(ii,ji)
        end do
      end do
c      write(6,*) ' x1a ',(x1a(i),i=1,n1)
c      write(6,*) ' x2a ',(x2a(i),i=1,n2)
c      write(6,*) ' ya1 ',(ya(i,1),i=1,3)
c      write(6,*) ' ya2 ',(ya(i,2),i=1,3)
c      write(6,*) ' ya3 ',(ya(i,3),i=1,3)
c      write(6,*) ' arg ',arg(1),arg(2)
      if(narg.eq.2) then
         call polin2(x1a,x2a,ya,n1,n2,arg(1),arg(2),y,dy)
      elseif(narg.eq.3) then
c         call polin3(x1a,x2a,x3a,yb,n1,n2,n3,arg(1),arg(2),arg(3),y,dy)
      else
         write(6,*) ' dhpolint narg = ',narg,' not implemented '
         stop
      endif
c      write(6,*) ' result: y = ',y,' dy = ',dy
      dhpolint = y
c      dhpolint = ya(1,1)
      RETURN
1000  FORMAT('  HFUNCTION DPOLINT ... NARG =',I6,
     +       '  NOT WITHIN RANGE')
      END
      SUBROUTINE locate_sm(xx,n,x,j)

c     From Numerical recipes:
c     -----------------------
c     Given an array XX of length N, and a given value of X, returns a
c     value of J such that X is between XX(J) and XX(J+1).  XX must be
c     monotonic, either increasing or decreasing. J=0 or J=N is
c     returned to indicate that X is out of range. 
      Implicit None     
      Integer j,n
      Double Precision x, xx(n)
      Integer jl,jm,ju
      jl=0
      ju=n+1
 10   IF(ju-jl.gt.1)THEN
          jm=(ju+jl)/2
          IF((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))THEN
              jl=jm
          ELSE
              ju=jm
          ENDIF
          go to 10
      ENDIF
      if(x.eq.xx(1)) then
        j=1
      else if(x.eq.xx(n)) then
        j=n-1
      else
        j=jl
      endif
      RETURN
      END
      SUBROUTINE polin3(x1a,x2a,x3a,ya,n1,n2,n3,x1,x2,x3,y,dy)
      Double Precision dy,x1,x2,x3,y,x1a(n1),x2a(n2),x3a(n3)
      Double Precision ya(n1,n2,n3)
      PARAMETER (N1MAX=20,N2MAX=20,N3MAX=20)
CU    USES polint
      Double Precision yn1tmp(N1MAX),yn2tmp(N2MAX),yn3tmp(N3MAX)
      do 12 j1=1,n1
        do 11 j2=1,n2
          do 10 j3=1,n3
             yn3tmp(j3)=ya(j1,j2,j3)
10        continue
          call polint(x3a,yn3tmp,n3,x3,yn2tmp(j2),dy)
11      continue
        call polint(x2a,yn2tmp,n2,x2,yn1tmp(j1),dy)
12    continue
      call polint(x1a,yn1tmp,n1,x1,y,dy)
      return
      END

      SUBROUTINE polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
      INTEGER m,n,NMAX,MMAX
      Double Precision dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      PARAMETER (NMAX=20,MMAX=20)
CU    USES polint
      INTEGER j,k
      Double Precision ymtmp(MMAX),yntmp(NMAX)
      do 12 j=1,m
        do 11 k=1,n
          yntmp(k)=ya(j,k)
11      continue
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
c        call ratint(x2a,yntmp,n,x2,ymtmp(j),dy)
12    continue
      call polint(x1a,ymtmp,m,x1,y,dy)
c      call ratint(x1a,ymtmp,m,x1,y,dy)
      return
      END
      
      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
C     Adapted from "Numerical Recipes" 
      IMPLICIT NONE
      INTEGER NS,I,M,N
      DOUBLE PRECISION XA(N),YA(N),X,Y,DY
      Integer NMAX
      Parameter (NMAX = 10)
      DOUBLE PRECISION C(NMAX),DIF,DIFT,HO,HP,W,D(NMAX),DEN

      NS=1
      DIF=DABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) then
             write(6,*) ' failure in polint '
          endif 
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
      
      SUBROUTINE ratint(xa,ya,n,x,y,dy)
C     Adapted from "Numerical Recipes" 
      IMPLICIT NONE
      INTEGER n,NMAX
      double precision dy,x,y,xa(n),ya(n),TINY
      PARAMETER (NMAX=10,TINY=1.e-25) 
      INTEGER i,m,ns
      double precision dd,h,hh,t,w,c(NMAX),d(NMAX) 
      ns=1
      hh=abs(x-xa(1)) 

      do i=1,n 
      	h=abs(x-xa(i)) 
            if (h.eq.0.)then 
                  y=ya(i)
                  dy=0.0
                  return
            else if (h.lt.hh) then 
                  ns=i 
                  hh=h 
            endif 
            c(i)=ya(i) 
            d(i)=ya(i)+TINY 
      end do 

      y=ya(ns) 
      ns=ns-1
      do m=1,n-1 
           do i=1,n-m 
               w=c(i+1)-d(i)
               h=xa(i+m)-x 
               t=(xa(i)-x)*d(i)/h
               dd=t-c(i+1)
               if(dd.eq.0.) then
                     write(6,*) 'failure in ratint',x,xa(i),i
                     dd=tiny
                     endif
               dd=w/dd 
               d(i)=c(i+1)*dd 
               c(i)=t*dd 
            end do 
            if (2*ns.lt.n-m)then 
               dy=c(ns+1) 
            else 
               dy=d(ns) 
               ns=ns-1
            endif
            y=y+dy 
      end do 
      return 
      END 

      SUBROUTINE ABEND
C
C CERN PROGLIB# Z035    ABEND           .VERSION KERNFOR  4.31  911111
C ORIG.  8/02/88  JZ
C

      STOP  7
      END
