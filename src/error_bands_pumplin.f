
      subroutine Error_Bands_Pumplin

      implicit none

      include 'steering.inc'
      include 'pdfparam.inc'
      include 'endmini.inc'
      include 'alphas.inc'
      include 'thresholds.inc'
      
      double precision hquad, umat, vmat,
     &     tvec, errmat
      common /umatco/ hquad(MNI,MNI), umat(MNI,MNI), vmat(MNI,MNI),
     &     tvec(MNI), errmat(MNI,MNI)
      integer sign
      double precision a
      dimension a(MNI)

      integer i,j,npar,idx,idx2,kflag,ii
      character*48 name,name2
      character*48 base,base2
      character tag(40)*3
      data (tag(i),i=1,40) /'s01','s02','s03','s04','s05',
     +     's06','s07','s08','s09','s10',
     +     's11','s12','s13','s14','s15',
     +     's16','s17','s18','s19','s20',
     +     's21','s22','s23','s24','s25',
     +     's26','s27','s28','s29','s30',
     +     's31','s32','s33','s34','s35',
     +     's36','s37','s38','s39','s40'/


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
C  Fix relation between internal and external params.
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
      print *,'Relation internal->external index:'
      do ind=1,mpar
         write (6,*) 'internal=',ind,' external=',iexint(ind)
      enddo

      write(6,*) 'mpar = ',mpar,mpar0

      npar = mpar0

      do i=1,npar
         write(6,*) 'pkeep(i) = ',pkeep(i)
      enddo


      do j=1,mpar
  
         jext = iexint(j)

         base = 'output/pdfs_q2val_'//tag(j)
         idx = index(base,' ')-1
         
         base2 = 'output/pdfs_'//tag(j)
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

C
C Initialise parameters for standard Parametrisation
C
            if (iparam.ne.171717) then
               Ag=a(1)
               Bg=a(2)
               Cg=a(3)
               Dg=a(4)
               Eg=a(5)
               Fg=a(6)
               Apg=a(7)
               Bpg=a(8)
               Cpg=a(9)

               Auv=a(11)
               Buv=a(12)
               Cuv=a(13)
               Duv=a(14)
               Euv=a(15)
               Fuv=a(16)
            
               Adv=a(21)
               Bdv=a(22)
               Cdv=a(23)
               Ddv=a(24)
               Edv=a(25)
               Fdv=a(26)
               
               Aubar=a(31)
               Bubar=a(32)
               Cubar=a(33)
               Dubar=a(34)
               
               AU=a(51)
               BU=a(52)
               CU=a(53)
               DU=a(54)
               
               Adbar=a(41)
               Bdbar=a(42)
               Cdbar=a(43)
               Ddbar=a(44)
               
               AD=a(61)
               BD=a(62)
               CD=a(63)
               DD=a(64)
            
               
               Asea=a(71)
               Bsea=a(72)
               Csea=a(73)
               Dsea=a(74)
               
               Adel=a(81)
               Bdel=a(82)
               Cdel=a(83)
               Ddel=a(84)
      


               alphas=a(95)
               fstrange=a(96)

               if (q0.ge.qc) then
                  fcharm=a(97)
               else
                  fcharm=0.
               endif
              
               
C
C  add x-dependent strange 
C     
               if (ifsttype.eq.0) then
                  fs0 = fstrange
               else
                  fs0 = fshermes(0.D0)
               endif
            
            else
! get the ctpara            

               do ii=1,6
                  ctglue(ii) = a(ii)
                  ctuval(ii) = a(10+ii)
                  ctdval(ii) = a(20+ii)
                  ctubar(ii) = a(30+ii)
                  ctdbar(ii) = a(40+ii)
                  ctother(ii)= a(94+ii)

               enddo
            endif

            if ((iparam.eq.1).or.(iparam.eq.11)) then
               Bd = Bu
               Bubar = Bu
               Bdbar = Bu
               Adbar = Ad
               aU = aD * (1.-fs0)/(1.-fcharm)
               aUbar = aU

* fixed for H1param
            elseif (iparam.eq.11) then

               Bd = Bu
               Bubar = Bu
               Bdbar = Bu
               Adbar = Ad
               aU = aD * (1.-fs0)/(1.-fcharm)
               aUbar = aU

* fixed for optimized H1param

            elseif (iparam.eq.2) then
               Aubar = Adbar * (1.-fs0)/(1.-fcharm)
               Bubar = Bdbar
            elseif (iparam.eq.22.or.iparam.eq.222.or.iparam.eq.221) then


               if (iparam.ne.221) then
                  Bdv = Buv
               endif

               Aubar = Adbar * (1.-fs0)/(1.-fcharm)
               Bubar = Bdbar

            elseif (iparam.eq.3) then ! g,uval,dval,sea as in ZEUS-S 2002 fit
               

               Buv = 0.5
               Bdv = 0.5
               Bdel = 0.5
               Cdel = Csea +2.

            elseif (iparam.eq.229) then
         

               Aubar = Adbar * (1.-fstrange)/(1.-fcharm)
               Bubar = Bdbar

            elseif ((iparam.eq.4).or.(iparam.eq.24)) then ! g,uval,dval,sea as in ZEUS-JET fit
               
               Bdv = Buv
*  dbar-ubar (not Ubar - Dbar), Adel fixed to output of ZEUS-S fit   
 
               Adel = 0.27          
               Bdel = 0.5
               Cdel = Csea +2.

            endif 

            kflag = 0
            call SumRules(kflag)
            call Evolution

C            print *,name2
            open (76,file=name2,status='unknown')
            call store_pdfs(name)
            close (76)

         enddo

      enddo

      return
      end


