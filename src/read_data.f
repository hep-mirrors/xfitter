      subroutine read_data
      
c     
c     Read data
c     


      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'

      
      include 'ntot.inc'
      include 'systematics.inc'
      include 'hadcor.inc'
      include 'indata.inc'
      include 'couplings.inc'


      character*10  cdummy
      character   adum 
      character*1   a1dum 

      double precision weak,pi,twopi,alphaem,Gf,convfac
C     -jf
      double precision ECM
c     
      dimension WEAK(NMAX)

      double precision fac

      logical FailCuts
      logical FIRST             !  true : cov matrix recalculated
      logical GNORM             !  correlated part for the luminosity errors

      real*4 q2,x,y,s,stat,tot,dum,cor,sta,sys
      real*4 e1,e2,e3,e4,e5,e6

      real*4 new1,new2,new3,new4,new5,new6
      real*4 new7,new8,new9,new10,new11,new12
      real*4 new13,new14
      integer insys
      real*4 phi
      real*4 dall,df2,df3,dfl
      real*4 stheo
      real*4 hcor,z0cor,uncp,uncm,e4p,e4m,totp,totm,hcorerr
      real*4 ehcorp,ehcorm
      real*4 sysp,sysm

      real*4 proc

      double precision Q2mean(6)
      double precision ETmean(24)
      double precision Xmean(24)

C      // Gouzevitch 2009/03/31

      double precision Width_q2(7),Width_e(5)
      double precision q2bin(7),ebin(5)
      

cv ///// study of PDF uncertainties
cv
      double precision alnorm
      external alnorm
      real logshift
      external logshift
      real lsig, lmu, lrunif
      real dummy, dummy_st
      integer vi,icount
      double precision  voica, voica_fl, voica_un
      real ranmflat
      double precision rand_shift(NSYS)
      double precision r_sh_fl(NSYS)
      integer numsys
      real rndsh, ranflat
      integer num,iseedrand, idate,is,ntime, ndate
      integer*4  timeArray(3)
      real voica_def, voica_deh, voica_lo
      real f_un
      COMMON/SLATE/IS(40)
      real alumlognorm(300)
      real alumierr(300)
      data alumierr/300*0.0/
cv////////
 

      integer i,j,k,iset,n0,isys,iq2bin,iebin,jsys,idx
      integer ixbin, ietabin, ietbin, used
      double precision xlum0,corlum,xlum,dp
      double precision newunc,unc,ee,ep,sep,epolar,syst,phicc
      double precision eres,eeff,edmul,ebmul,efrag,emod,euds,ephi
      double precision f2_dglap,stat_dglap,sysu_dglap,sysd_dglap
      double precision f2_ccfm,stat_ccfm,sysu_ccfm,sysd_ccfm
      double precision sysu,sysd,kfac
      
      double precision ed1,ed2,ed3,ed4,ed5
      double precision ed6,ed7,ed8,ed9,ed10,ed11
      double precision ed12,ed13,ed14,ed15,ed16,ed17
      double precision ed18,ed19,ed20, ed21,ed22,ed23,ed24,ed25,ed26
      double precision fl, corr

      double precision dum1, dum2, dum3, dum4
      double precision dum5, dum6, dum7, dum8
cv fltoy
      double precision col1, col2, col3
cv ZEUS, H1:
      double precision cor1, cor2


******
      integer iss,n1,n2,n3,n4,jj,icombo1
      double precision dumsys(nsys),factorCC
      double precision sys35, sysel, sysmu,coel,comu

*new2009 jf
* nocorrel = 0, all  (=nsyshz) correl kept in data sets 28 to 31
* nocorrel = 1, all correl removed in data sets 28 to 31
* nocorrel = 2, all correl except the last four or two (i.e., (nsyshz - nsyshzm4))
* nocorrel = 3 and nfirst = number of correlated errors, ranked by decreasing size:
* which are kept  while (nsyshz-nfirst) errors become uncorrelated 
*     if nfirst =nsyshz => nocorrel 3=0
*     if nfirst =0 => nocorrel 3=1
*
*new4 jf
* nlowq2 = 1, new average of H1-ZEUS data sets includung new H1-MB
* number of systematics nsyshz increases from 47 to 96 -> 102
* number of data points in set 29 increases 
* 
      integer isup,nsyshz, nsyshzp150
      integer nsyshzp100,nsyshzm3,nsyshzm4,newdat62
      double precision totot,totot2,tototm
      real aaaa
* 22 Jan jf : unc2 double precision and nfirst integer introduced for the 
* nocorrel = 3 option
      double precision unc2
      
cv+sh **********
* Add sys info for FL data (low energy)
* 25 experimental errors, 3 procedural
* nsysfl=28, nsysflm3=28-3, nsysflp300 (max slot!)
      integer nsysfl, nsysflm3, nsysflp300
cv+sh **********
cv new herai+ii data
      integer ntothz,ntothzm3, ntothzp400
cv=highq2
      integer nsyshiq, nsyshiq250
cv=lowq2
      integer nsysloq, nsysloq200
cv=f2cc h1+zeus
      integer nsysfc, nsysfcm1, nsysfcp340
*      integer nfirst

cv      double precision f2, norm, whad2,statqcd, sysqcd, f2qcd, shift


      double precision f2, norm, whad2,statqcd, sysqcd, f2qcd,lambda, shift
cv atlas
      double precision ymin, ymax, mcsta,kfactor, kfactorm, kfactorp

cv ... THIS FOR MC STUDY!
cv saving the lumi uncertainties in an array to be used for lognormal random shifts
      alumierr(1) = 1.7
      alumierr(7) = 1.7
      alumierr(8) = 1.5
      alumierr(14) =1.8
      alumierr(18) =1.5
      alumierr(19) =2.2
      alumierr(32) = 0.005


      do i=1,31
         if (alumierr(i).gt.0) then
            alumierr(i) = sqrt((0.01*alumierr(i))**2-alumierr(32)**2)
         endif
      enddo

cv end of saving the lumi ...
! new data set for nlowq1=1, otherwis =0
cv NOW IN THE STEERING!
c      nlowq2 = 1
c: default values correspond to the PRELIM DIS data



cv+sh  comment: p100(300) means allocating the slot for
cv+sh  maximum number of sys in the total sys array
      nsysfl = 28 ! total number of sys
      nsysflm3=nsysfl-3 ! remove procedural errors
      nsysflp300 = nsysfl + 300


cv+herai+ii data
      ntothz = 134
      ntothzm3 =ntothz-3
      ntothzp400=ntothz+400

cv=highq2

      nsyshiq=32
      nsyshiq250=nsyshiq+250

cv=lowq2

      nsysloq=46
      nsysloq200=nsysloq+200



      if (nlowq2.eq.0) then
         nsyshz = 47
         newdat62 = 383


         nsyshzp100 = nsyshz + 100
         nsyshzm3 = nsyshz - 3      
         nsyshzm4 = nsyshzm3 - 1
        

      elseif (nlowq2.eq.1) then
cv new data set has 100+2 systematic uncertainties, and the new data has 527 points..
         nsyshz = 102
         newdat62 = 527
cv  last two are the procedural errors, this is used with nocorrel=2 case 
cv (errors added in quadrature, up to last two procedural errors)
         nsyshzm4=nsyshz-2
      elseif (nlowq2.eq.2) then
cv new data set has 110+3 systematic uncertainties, and the new data has 528 points..
         nsyshz = 113
         newdat62 = 528
cv  last three are the procedural errors, this is used with nocorrel=2 case 
cv (errors added in quadrature, up to last two procedural errors)
         nsyshzm4=nsyshz-3
      endif


cv=f2cc h1+zeus
      nsysfc = 2 ! total number of sys
      nsysfcm1=nsysfc-1 ! remove procedural errors
      nsysfcp340 = nsysfc + 340


*18 Jan 2009 jf
      nsyshzp100 = nsyshz + 100
*     end special initialisation

      DEBUG = lDEBUG
      ONLINE = lONLINE
      FIRST = lFIRST

      GNORM = .true.

      ISCHARMANY  = .false.
      ISBOTTOMANY = .false.


      pi = 4.d0 * datan(1.d0)
      twopi = 2.d0 * pi
      alphaem = 1.d0 / 137.035999d0


C-2- 05/08/2010: end of the addition --------------


      Gf = 1.16637d-5
      convfac = 3.893796D8
      Mw = 80.41d0
      factorCC=0.d0
      do iset=1,nset
         POLAR(iset) = 0.d0
      enddo

      do iset=1,nset
         ISJET(iset) = .false.
         ISCHARM(iset) = .false.
         ISBOTTOM(iset) = .false.
         ISDYZ(iset) = .false.
         ISDYW(iset) = .false.
      enddo

      do i=1,nsys
         do j=1,ntot
            BETA(i,j) = 0.d0
         enddo
      enddo

      do i=1,nset
         do j=1,nmax
            VAL_Q2(i,j) = 0.d0
            VAL_X(I,J)  = 0.D0
            VAL_Y(I,J)  = 0.D0
            VAL_S(I,J)  = 0.D0
            VAL_E(I,J)  = 0.D0
            VAL_STA(I,J) = 0.D0
            VAL_UNC(I,J) = 0.D0
         enddo
      enddo


      do i=1,nset
         NDATAPOINTS(i) = 0
      enddo
      npoints = 0
      ndis=0
      n0 = 0
      n1 = 0
      n2 = 0
      n3 = 0
      n4 = 0

cv============
cv RANDOM SHIFTS


      if (lRAND) then
         f_un = 2.
c         sys=3
c         sta=1

cv         call system_clock(icount)
cv         call itime(now)
cv         icount=now(1)+now(2)+now(3)
         icount= time()
         print*,' clock = ', icount
         
         call datime(ndate,ntime)
         ntime = ntime*1000000+is(6)
         icount=ntime

cv initialise the random shifts
         do numsys=1,nsys
            rand_shift(numsys)=0.
            r_sh_fl(numsys)=0.
         enddo

cv initialise the random seed gener
         if (iseedmc.ne.0) then
C SG: Overwrite seed to one selected in the steering:
            icount = iseedmc
         endif

         call rmarin(icount,0,0)
         call rluxgo(3,icount,0,0)
         print*,'initialize smeering with a seed isdrn = ',icount
         
c         stop
         do numsys=1,nsys
            call rnorml(rndsh,1)    ! gauss random number
            call ranlux(ranflat,1)   ! uniform random number

            rand_shift(numsys) = rndsh
            r_sh_fl(numsys) = ranflat

            print*,'random numbers: sys, gauss, flat ',
     $           numsys,rand_shift(numsys),
     $           r_sh_fl(numsys)

cv save shifts for lumi uncertainties only (but we are not using it)            
            if ((numsys.eq.1).or.(numsys.eq.7).or. (numsys.eq.14)
     $           .or.(numsys.eq.18).or.(numsys.eq.19)
     $           .or.(numsys.eq.25).or.(numsys.eq.26).or.(numsys.eq.27)
     $           .or.(numsys.eq.32).or.(numsys.eq.8)) then
               alumlognorm(numsys) = alnorm(1.,alumierr(numsys))
            else
               
            endif
         enddo

      endif

cv ============done with initializing=========

      xlum0 = 1.8d0
                                !
      corlum = 0.5d0 
                                ! common for all data sets (theoretical uncertainties on BH xsec)



C
C Read data from namelists:
C
      do i=1,NInputFiles
         call ReadDataFile(InputFileNames(i))
      enddo

C-----------------------------------------
      print*,'number of points', npoints


*     ------------------------------------------------------------------
*     -- with this the normalisations are minuit parameters instead...

      if (lNORMA) then

         do jsys=1,nsys
            do k=1,npoints
               if (jsys.eq.1.or.jsys.eq.7.or.jsys.eq.8.or.
     +              jsys.eq.14.or.jsys.eq.18.or.jsys.eq.19.or.
     +              jsys.eq.25.or.jsys.eq.26.or.jsys.eq.27) then
                  beta(jsys,k) = 0.d0
               endif
            enddo
         enddo

      endif
*     ------------------------------------------------------------------




      do i=1,nsys
         do k=1,npoints
            beta(i,k) = beta(i,k) / 100.
            if (.not.lCORR) beta(i,k)=0.d0
         enddo
      enddo


      do i=1,npoints
         ALPHA(i) = ALPHA(i) / 100.
         if (alpha(i).le.0) write(6,*) 'alpha(i) = 0 for point ',i
      enddo





*     ------------------------------------------------------------------
*     -- Calculate or read the full covariance matrix

      if (ICHI2.eq.3) then

         if (FIRST) then
            call CoVarMatrix
         else
            call Read_CoVarMatrix
         endif

      endif


*     ----------------------------------------------------------------------
*     -- CTEQ-like chi2 :
*     -- compute the matrix (A) as given in Eq. (B.4) of JHEP07 (2002) 012
*     -- This is the matrix "sysa" in common/systema/

      if (ICHI2.eq.2) then
         call Systematics
      endif



      return
      end



      subroutine ReadDataFile(CFile)
C------------------------------------------------------------------------
C
C  Created 20 May 2011 by SG.
C   Read data set using namelist format
C
C------------------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'systematics.inc'

      character *(*) CFile
C Namelist  variables:    
      integer ndataMax,ninfomax,nsystMax
      parameter (ndataMax=1000)
      parameter (ninfoMax=100)
      parameter (nsystMax=500)

      character *80 Name
      integer  NData
      integer  NSyst
      integer  NInfo
      integer  NBinDimension

      
      character *80 BinName(NBinDimensionMax)
      double precision datainfo(ninfoMax)
      character *80 CInfo(ninfoMax)
      character *80 Reaction
      logical Percent(0:nsystMax)
      integer SystematicType(nsystMax)
      integer IndexDataset
      double precision SystScales(nsystMax)
C Extra info about k-factors, applegrid file(s):
      character*80 TheoryInfoFile,TheoryType
      character*80 KFactorNames(NKFactMax)
      integer      NKFactor
C Namelist definition:
      namelist/Data/Name,NData,NSyst,NBinDimension
     $     ,BinName,NInfo,datainfo,CInfo,Reaction,Percent
     $     ,SystematicType, SystScales, IndexDataset
     $     ,TheoryInfoFile,TheoryType,KFactorNames,NKFactor

      double precision XSections(ndataMax)
      double precision StatErrors(ndataMax)
      double precision AllBins(10,ndataMax)
      double precision Syst(nsystmax)

      double precision Akfact(NKFactMax)

      double precision StatError

      double precision UncorError  ! uncorrelated systematics
      double precision TotalError  ! total uncertainty

      integer i,j
      logical LReadKFactor

C Function to check cuts
      logical FailSelectionCuts
      
C-------------------------------------------------------      

C Reset to default:
      NSYST = 0
      NData = 0
      NInfo = 0
      NBinDimension = 0
      Reaction = ' '
      Name     = ' '
      IndexDataSet = 0
      NKFactor = 0
      TheoryInfoFile = ' '
      LReadKFactor = .false.

      open(51,file=CFile,status='old',err=99)

      print *,'Reading data file ...'
      print *,CFile
      read(51,NML=Data,err=98)
C
C Store 
C
      NDATASETS = NDATASETS + 1
      DATASETNUMBER(NDATASETS)   = 10000+NDATASETS
      DATASETLABEL(NDATASETS)    = Name
      DATASETNUMBER(NDATASETS)   = IndexDataset   !!!  
C Reaction info:
      DATASETREACTION(NDATASETS) = Reaction
C Binning info:
      DATASETBinningDimension(NDATASETS) = NBinDimension
      do i=1,NBinDimension
         DATASETBinNames(i,NDATASETS) = BinName(i)
      enddo
C Extra info:
      DATASETInfoDimension(NDATASETS) = NInfo
      do i=1,NInfo
         DATASETInfoNames(i,NDATASETS) = CInfo(i)
         DATASETInfo(i,NDATASETS) =      DataInfo(i)
      enddo
      NQ2BINS(NDATASETS) = 0            ! uff

C Check if we need to read kfactor file:
      if (TheoryInfoFile.ne.' ') then
         if (TheoryType.eq.'kfactor') then
            open (52,file=TheoryInfoFile,status='old',err=100)
            lreadkfactor = .true.
         endif
      endif
      DATASETNKfactors(NDATASETS) = NKFactor
      do i=1,NKFactor
         DATASETKFactorNames(i,NDATASETS) = KFactorNames(i)
      enddo


C Read data info:
      do j=1,NData
         read(51,*)(allbins(i,j),i=1,NBinDimension),XSections(j)
     $        ,StatErrors(j),(syst(i),i=1,NSyst)

         if (lreadkfactor) then
            read (52,*) (akfact(i),i=1,NKFactor)
         endif

C Apply cuts:
         if (FailSelectionCuts(Reaction,NBinDimension,allbins(1,j),BinName)) then
            goto 1717
         endif

C Add a point:
         npoints = npoints+1
         
         if (npoints.ge.NTOT) then
            print '('' ReadDataFile Error, increase NTOT value inside ntot.inc'')'
            stop
         endif

         NDATAPOINTS(NDATASETS) = NDATAPOINTS(NDATASETS) + 1
         DATASETIDX(NDATASETS,NDATAPOINTS(NDATASETS)) = npoints

C Translate errors in %:
         TotalError = 0.
         UncorError = 0.

         if (.not.Percent(0)) then
            StatErrors(j) = StatErrors(j)/XSections(j)*100.
         endif
         TotalError = TotalError + StatErrors(j)**2

         do i=1,NSyst
            if (.not.Percent(i)) then
               syst(i) = syst(i)/XSections(j)*100.
            endif
   
            TotalError = TotalError + Syst(i)**2
            if (SystematicType(i).eq.0) then
C Uncorrelated error:
               UncorError = UncorError +  Syst(i)**2
            endif
         enddo

         UncorError = sqrt(UncorError)
         TotalError = sqrt(TotalError)

         DATEN(npoints) = XSections(j)
         E_UNC(npoints)  = UncorError
         E_TOT(npoints)  = TotalError
         E_STA(npoints)  = StatErrors(j)

         do i=1,NBinDimension
            AbstractBins(i,npoints) = allbins(i,j)
         enddo

         ALPHA(npoints) = sqrt(UncorError**2+StatErrors(j)**2)*DATEN(npoints)
         do i=1,NSyst
            if (SystematicType(i).gt.0) then
               if (SystematicType(i).gt. NSys) then
                  print '(''ReadDataFile Error: requested error source'',i6,'' larger than NSYST='',i6)'
     $                 ,SystematicType(i),NSys
                  print '(''Check SystematicType or increase NSys in systematics.inc'')'
                  stop
               endif

               BETA(SystematicType(i),npoints) = syst(i)
            endif
         enddo

C         print *,'hhhh',alpha(npoints),e_sta(npoints)

         JSET(npoints) = IndexDataset  ! XXXXXXXXXXXX

C Store k-factors:
         if (lreadkfactor) then
            do i=1,nkfactor
               kfactors(i,npoints) = akfact(i)
            enddo
         endif
         
 1717 enddo

      close (51)
      if (lreadkfactor) then
         close (52)
      endif

      print '(''Read'',i8,'' data points for '',A80)',NData,Name
      print '(''Printing first'',i5,'' data points'')',min(Ndata,5)
      print '(20A14)',(BinName(i),i=1,NBinDimension),' sigma'
     $     ,' stat. err'
      do j=1,min(NData,5)
         print '(20E14.4)',(Allbins(i,j),i=1,NBinDimension),XSections(j)
     $        ,StatErrors(j)
      enddo
      return

 98   continue
      print '(''Error reading namelist Data'')'
      print *,CFile
      stop

 99   continue
      print '(''Can not open file '')'
      print *,CFile
      stop
100   continue
      print '(''Can not open file '')'
      print *,TheoryInfoFile
      stop
      end


