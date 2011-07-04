
      subroutine Error_Bands_Pumplin

      implicit none

      include 'steering.inc'
      include 'pdfparam.inc'
      include 'endmini.inc'
      include 'alphas.inc'
      
      double precision hquad, umat, vmat,
     &     tvec, errmat
      common /umatco/ hquad(MNI,MNI), umat(MNI,MNI), vmat(MNI,MNI),
     &     tvec(MNI), errmat(MNI,MNI)
      integer sign
      double precision a
      dimension a(MNI)

      integer i,j,npar,idx,idx2,kflag
      character*25 name,name2
      character*20 base,base2
      character tag(14)*3
      data (tag(i),i=1,14) /'s01','s02','s03','s04','s05',
     +     's06','s07','s08','s09','s10',
     +     's11','s12','s13','s14'/


      include 'couplings.inc'

C SG: x-dependent fs:
      double precision fs0
      double precision fshermes

      
      character*20 parname
      integer  ind,ind2,jext,iint

      integer nparmax
      parameter (nparmax=1000)

      integer  iunint(nparmax)  ! internal param. number
      integer  iexint(nparmax)  ! external param. number
      double precision 
     $     parval, parerr,parlolim,parhilim

C---------------------------------------------------------------
      
C
C 22 Dec 2010: SG: Fix relation between internal and external params.
C
      do ind=1,mpar0
         call mnpout(ind,parname,parval,parerr,parlolim,
     $        parhilim,iunint(ind))
         write (6,*) 'Parameter',ind,' name=',parname
         write (6,*) 'Internal index=',iunint(ind)
         write (6,*) ' '
      enddo

      do ind=1,mpar
         do ind2=1,mpar0
            if (iunint(ind2).eq.ind) then
               iexint(ind) = ind2
            endif
         enddo
      enddo
      print *,'Relation internal->extrnal index:'
      do ind=1,mpar
         write (6,*) 'internal=',ind,' external=',iexint(ind)
      enddo
      

C
C 31 Oct 2009: SG: add x-dependent strange 
C
      if (ifsttype.eq.0) then
         fs0 = fstrange
      else
         fs0 = fshermes(0.D0)
      endif


      write(6,*) 'mpar = ',mpar,mpar0

      npar = mpar0

      do i=1,npar
         write(6,*) 'pkeep(i) = ',pkeep(i)
      enddo


      do j=1,mpar
  
         jext = iexint(j)

         base = 'pdfs_q2val_'//tag(j)
         idx = index(base,' ')-1
         
         base2 = 'pdfs_'//tag(j)
         idx2  = index(base2,' ')-1

         do sign=-1,1,2
            if (sign.eq.-1) then
               if (idx.gt.0) then
                  name  = base(1:idx)//'m_'
                  name2 = base2(1:idx2)//'m.lhgrid'
               else
                  name = base//'m_'
                  name2 = base2//'m.lhgrid'
               endif
            else
               if (idx.gt.0) then
                  name = base(1:idx)//'p_'
                  name2 = base2(1:idx2)//'p.lhgrid'
               else
                  name = base//'p_'
                  name2 = base2//'m.lhgrid'
               endif
            endif
            write(6,*) 'name = ',name


            do i=1,npar

               a(i) = pkeep(i) 

               iint = iunint(i)
               if (iint.gt.0) then
                  a(i) = a(i) + sign * umat(iint,j)
               endif
               write(6,*) 'i a(i) = ',i,a(i)

            enddo

            if (iparam.eq.1) then

               Bg  = a(1)
               Cg  = a(2)
               Dg  = a(3)

               Bu  = a(4)
               Cu  = a(5)
               Fu  = a(6)

               Ad  = a(7)
               Cd  = a(8)

               Cubar = a(9)
               Cdbar = a(10)
              
               Bd = Bu
               Bubar = Bu
               Bdbar = Bu
               Adbar = Ad
               aU = aD * (1.-fs0)/(1.-fcharm)
               aUbar = aU

* fixed for H1param
               Alphas = a(11)
               Eu = a(12)

            elseif (iparam.eq.11) then

               Bg  = a(1)
               Cg  = a(2)

               Bu  = a(3)
               Cu  = a(4)
               Eu = a(5)
               Fu  = a(6)

               Ad  = a(7)
               Cd  = a(8)

               Cubar = a(9)
               Cdbar = a(10)

               Bd = Bu
               Bubar = Bu
               Bdbar = Bu
               Adbar = Ad
               aU = aD * (1.-fs0)/(1.-fcharm)
               aUbar = aU

* fixed for optimized H1param
               Dg = a(11)
               Alphas = a(12)

            elseif (iparam.eq.2) then

               Bg = a(1)
               Cg = a(2)
               Dg = a(3)

               Buv = a(4)
               Bdv = Buv
               Cuv = a(5)
               Duv = a(6)

               Cdv = a(7)
               Ddv = a(8)

               Adbar = a(9)
               Aubar = Adbar * (1.-fs0)/(1.-fcharm)

               Bdbar = a(10)
               Bubar = Bdbar

               Cdbar = a(11)
               Cubar = a(12)

* fixed for inbetween
               Alphas = a(13)
               Euv = a(14)

            elseif (iparam.eq.22.or.iparam.eq.222.or.iparam.eq.221) then

               Bg = a(1)
               Cg = a(2)

               Buv = a(3)
               if (iparam.eq.221) then
                  Bdv=a(17)
		else
               Bdv = Buv
               endif
               Cuv = a(4)
               Duv = a(5)

               Cdv = a(6)

               Adbar = a(7)
               Aubar = Adbar * (1.-fs0)/(1.-fcharm)

               Bdbar = a(8)
               Bubar = Bdbar

               Cdbar = a(9)
               Cubar = a(10)

               Euv = a(11)

* fixed for optimized in between
               Dg = a(12)
               Ddv = a(13)
               Ddbar=a(15)
               Dubar=a(16)	
               Apg=a(18)
               Bpg=a(19)
               Rfudge=a(20)
               afudge=a(21)
               f2ht1=a(22)
               f2ht2=a(23)
               if (iparam.eq.222) then
                  Cpg=25.
               endif
               Alphas = a(14)
            elseif (iparam.eq.225) then

               Bg = a(1)
               Cg = a(2)

               Buv = a(3)
               Bdv = Buv
               Cuv = a(4)
               Duv = a(5)

               Cdv = a(6)

               Adbar = a(7)
               Aubar = a(15)

               Bdbar = a(8)
               Bubar = a(16)

               Cdbar = a(9)
               Cubar = a(10)

               Euv = a(11)

* fixed for optimized in between
               Dg = a(12)
               Ddv = a(13)
               Alphas = a(14)

            elseif (iparam.eq.3) then ! g,uval,dval,sea as in ZEUS-S 2002 fit
               
               Bg = a(1)
               Cg = a(2)
               Dg = 0.          ! Different from H1PDF2k
               
               Cuv = a(3)
               Buv = 0.5
               Duv = a(4)
               Fuv=0.
               
               Cdv = a(5)
               Bdv = 0.5
               Ddv = a(6)
               Fdv=0.
               
               Asea = a(7)
               Bsea = a(8)
               Csea = a(9)
               Dsea = a(10)
               
               Adel = a(11)
               Bdel = 0.5
               Cdel = Csea +2.

            elseif (iparam.eq.229) then
         
               Bg = a(1)
               Cg = a(2)
         
               Buv = a(3)
               
               
               Bdv=a(17)
               
               Cuv = a(4)
               Duv = a(5)
         
               Cdv = a(6)
               
               Adbar = a(7)
               Aubar = Adbar * (1.-fstrange)/(1.-fcharm)
               
               Bdbar = a(8)
               Bubar = Bdbar
               
               Cdbar = a(9)
               Cubar = a(10)
               
               Euv = a(11)
               
* fixed for optimized in between
               Dg = a(12)
               Ddv = a(13)
               Alphas = a(14)
               
               Ddbar=a(15)
               Dubar=a(16)	
               
               Apg=a(18)
               Bpg=a(19)
               
               Rfudge=a(20)
               afudge=a(21)
               Cpg=25.
               


            elseif (iparam.eq.4) then ! g,uval,dval,sea as in ZEUS-JET fit

               Bg = a(1)
               Cg = a(2)
               Dg = a(3)         
               
               Buv = a(4)
               Cuv = a(5)
               Duv = a(6)
               
               Bdv = Buv
               Cdv = a(7)
               Ddv = a(8)
               
               Asea = a(9)
               Bsea = a(10)
               Csea = a(11)

*  dbar-ubar (not Ubar - Dbar), Adel fixed to output of ZEUS-S fit   
 
               Adel = 0.27          
               Bdel = 0.5
               Cdel = Csea +2.

* fixed for ZEUS-JETS
               Alphas = a(12)
               Euv = a(13)

            elseif (iparam.eq.24) then ! g,uval,dval,sea as in ZEUS-JET fit

               Bg = a(1)
               Cg = a(2)
               Dg = a(3)         
               
               Buv = a(4)
               Cuv = a(5)
               Duv = a(6)
               
               Bdv = Buv
               Cdv = a(7)

               Euv = a(8)
               
               Asea = a(9)
               Bsea = a(10)
               Csea = a(11)

*  dbar-ubar (not Ubar - Dbar), Adel fixed to output of ZEUS-S fit   
 
               Adel = 0.27          
               Bdel = 0.5
               Cdel = Csea +2.

* fixed for ZEUS-JETS optimized
              
               Ddv = a(12)
               Alphas = a(13)

            endif 

            kflag = 0
            call SumRules(kflag)
            call Evolution

            print *,name2
            open (76,file=name2,status='unknown')
            call store_pdfs(name)
            close (76)

         enddo

      enddo

      return
      end


