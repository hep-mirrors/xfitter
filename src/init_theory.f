      subroutine init_theory_modules
*     ------------------------------------------------

      implicit none 
#include "ntot.inc"
#include "steering.inc"

*     ------------------------------------------------
*     Initialise EW parameters
*     ------------------------------------------------

      call Init_EW_parameters

*     ------------------------------------------------
*     Initialise qcdnum and APFEL
*     ------------------------------------------------
      if(itheory.eq.0.or.itheory.eq.10.or.itheory.eq.11
     $.or.itheory.eq.35) then
C Init evolution code:
         call qcdnum_ini
C Init APFEL if needed
         if(itheory.eq.10.or.itheory.eq.35) call apfel_ini
C Init QEDEVOL if needed
         if(itheory.eq.11) call qedevol_ini

         call Init_heavy_flavours

         if (ewfit.gt.0) call eprc_init(.true.)
      elseif(itheory.ge.100) then       
cc           write(6,*) ' in ini_theory for itheory =',itheory
      endif

*     ------------------------------------------------
*     Initialise calculations for each dataset:
*     ------------------------------------------------
      if(Itheory.ge.100) then
ccc        write(6,*) ' ini_theory: no data sets initialised for theory ',itheory
      else
c      call Init_theory_datasets
      endif

      return
      end



      Subroutine Init_heavy_flavours()
*-----------------------------------------------------
*
*   Init heavy flavour modules
*
*----------------------- ------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
C----------------------------------------


c Extra init for QCDNUM schemes:
      Call ZMVFNS_init()
      end


      subroutine qcdnum_ini
*     ------------------------------------------------
*     QCDNUM initialisation
*     ------------------------------------------------

      implicit none

#include "ntot.inc"
#include "steering.inc"
#include "thresholds.inc"
#include "couplings.inc"
#include "datasets.inc"

         !> QCDNUM grid definitions
c     xmin_grid defined in steering.inc
      integer  iwt_xgrid(NMXGRID)         !> X-grid population

      integer NQGRID                     !> Min. number of Q2 sub-ranges
      PARAMETER (NQGRID=2)

      integer NQGridHF  !> Increase number of Q2 grid points to ensure that HF thresholds and the starting
                        !> Q2 scale are in.
      parameter (NQGridHF=NQGrid+5)  


c set-up of the constants
      integer iosp,nqout
 
      integer NQall     !> Actual number of grid points
      double precision Q2Grid(NQGridHF)
      double precision WQGrid(NQGridHF)

      double precision  QARR(NQGRID),WGT_Q2(NQGRID)
      data iosp/2/                   !x grid, lin/qua/spli

      double precision  tmp

      double precision as0,r20
      data as0/0.364/, r20/2.D0/!, nfin/0/ !alphas, NNLO, VFNS
      integer I,ndum,ierr,j

      double precision a,b, qmz
      
      integer id1,id2
      integer nw,nwords,nx
      integer iqb,iqc,iqt
      double precision tiny
      parameter(tiny=1d-3)


      integer NQ2bins, NXbins !> requested number of x,q2 bins
      logical ReportXGrid, ReadXGrid

      namelist/qcdnum/xmin_grid, iwt_xgrid, iosp, wgt_q2, QARR,
     $     NQ2bins, NXbins, Read_QCDNUM_Tables, ICheck_QCDNUM,
     $     ReportXGrid, ReadXGrid, kmuc, kmub, kmut

      double precision xgrid(NXGridMax)
      integer NXGridAct

C Functions:
      integer iqfrmq

      common/ext_xgrid/xgrid,NXGridAct,ReadXGrid

C Memory parameters of QCDNUM
      integer ver,mxg,mxx,mqq,msr,mce,mbf,mky,nwf

C---------------------------------------------------------------------------------------

C-----  DEFAULTS -----------

C More detailed at high x:
      xmin_grid(1) = 9.9D-7
      xmin_grid(2) = 0.01D0
      xmin_grid(3) = 0.10D0
      xmin_grid(4) = 0.4D0
      xmin_grid(5) = 0.7D0

C Increasingly more dense grid:
      iwt_xgrid(1) = 1
      iwt_xgrid(2) = 2
      iwt_xgrid(3) = 4
      iwt_xgrid(4) = 8
      iwt_xgrid(5) = 16
C Q2 grid weights 
      WGT_q2(1) = 1.d0
      WGT_q2(2) = 1.d0
C Basic Q2 grid:
      QARR(1) = 1.
      QARR(2) = 2.05D8 ! needed for lhapdf grid
C Reduce the Q2 interval if small-x resummation through APFEL is included.
      if(HFSCHEME.eq.3005.or.
     1   HFSCHEME.eq.3055.or.
     2   HFSCHEME.eq.3555)then
         QARR(1) = starting_scale
         QARR(2) = 2.025D7      ! needed for lhapdf grid  
      endif
c      QARR(2) =  64000000.      ! enough for 8 TeV LHC.

C Default sizes
      NQ2bins = 120
      NXbins  = 200

C Default QCDNUM check
      ICheck_QCDNUM = 0

C Default extra info:
      ReportXGrid = .false.
      ReadXGrid= .false.

C----- Rescaling factors for the heavy-quark threholds (used by APFEL)
      kmuc = 1d0
      kmub = 1d0
      kmut = 1d0


C----- Read grid definitions from the steering.txt
      Read_QCDNUM_Tables = .false.

      open (51,file='steering.txt',status='old')
      read (51,NML=QCDNUM,ERR=7117,END=1771)

 1771 continue

      close (51)

C     If any of the rescaling factors kmuc, kmub, or kmut is different from one,
C     check that APFEL is used for the evolution otherwise stop the code.
      if(kmuc.ne.1d0.or.kmub.ne.1d0.or.kmut.ne.1d0)then
         if(itheory.ne.10.and.itheory.ne.35)then
            call HF_errlog(28042016, 'F: '//
     1           'When using displaced heavy quark thresholds, '//
     2           'APFEL must be used for the evolution. '//
     3           'Please set TheoryType = "DGLAP_APFEL" or '//
     4           'TheoryType = "DGLAP_APFEL_QED" in the '//
     5           'steering.txt card.')
         endif
      endif

C--------------------- kmuc/kmub/kmut can be re-defined  in the namelist file:

      q0 = starting_scale
      qc = ( kmuc * HF_MASS(1) )**2
      qb = ( kmub * HF_MASS(2) )**2
      qt = ( kmut * HF_MASS(3) )**2
      qmz= Mz**2



C Check that scales are in proper order:
      if (q0.gt.qc) then
c In fixed-flavour scheme starting scale can be above charm threshold
         if(HFSCHEME.ne.3.and.HFSCHEME.ne.4.and.HFSCHEME.ne.444) then
            print *,'Starting scale must be below charm threshold, stop'
            call HF_stop
         endif
      endif
      
      if (qc.gt.qb) then
         print *,'Charm mass must be below bottom mass, stop'
         call HF_stop
      endif
      
      if (qb.gt.qt) then
         print *,'Bottom mass must be below top mass, stop'
         call HF_stop
      endif

      do i=1,NQGrid 
         Q2Grid(i) = QARR(i)
         WQGrid(i) = WGT_Q2(i)
      enddo
C Add extra points:
      Q2Grid(NQGrid+1) = q0
      Q2Grid(NQGrid+2) = qc
      Q2Grid(NQGrid+3) = qb
      Q2Grid(NQGrid+4) = qmz
      WQGrid(NQGrid+1) = 4.D0
      WQGrid(NQGrid+2) = 2.D0
      WQGrid(NQGrid+3) = 1.3D0
      WQGrid(NQGrid+4) = 1.1D0

      if (qt.lt.qarr(2)) then
         WQGrid(NQGrid+5) = 1.D0
         Q2Grid(NQGrid+5) = qt
      else
         ! no top
         WQGrid(NQGrid+5) = 0.
         Q2Grid(NQGrid+5) = qb
      endif

C Sort the Q2Grid:
      do i=1,NQGridHF
         do j=i+1,NQGridHF
            if ( Q2Grid(j) .lt. Q2GRID(i) ) then
               tmp =  Q2GRID(j)
               Q2Grid(j) = Q2GRID(i)
               Q2GRID(i) = tmp
               tmp =  WQGRID(j)
               WQGrid(j) = WQGRID(i)
               WQGRID(i) = tmp               
            endif
         enddo
      enddo

      NQAll = NQGridHF


C Remove duplicates:
      i = 1
      do while (i.lt.NQAll)
         if (  abs(Q2Grid(i)-Q2Grid(i+1)).lt.1.D-5 .or. 
     $        ( Q2Grid(i).eq.0 ) ) then
            do j=i+1,NQAll-1               
               Q2Grid(j) = Q2Grid(j+1)
               WQGrid(j) = WQGrid(j+1)
            enddo
            NQAll = NQAll-1
         else
            i = i + 1
         endif
      enddo
 
      print *,' '
      print *,'Info FROM QCDNUM_INI'
      print '('' Init Q2 grid with number of nodes='',i5)',NQALL      
      print '('' Q2 values at:'',20F14.1)',(Q2grid(i),i=1,NQALL)
      print '('' Weights are :'',20F14.1)',(WQgrid(i),i=1,NQALL)
      print *,' '
      print '('' Init X  grid with number of nodes='',i5)',NMXGRID
      print '('' X  values at:'',20E11.2)',(xmin_grid(i),i=1,NMXGRID),1.0
      print '('' Weights are :'',20I11)',(iwt_xgrid(i),i=1,NMXGRID)
     $     ,iwt_xgrid(NMXGRID)
      print *,' '

      call qcinit(6,' ')        !initialize
*
*     After the initialization of QCDNUM, the code checks that the memory parameters
*     used to compile QCDNUM are appropriate.
*
      call getint('vers',ver)
      call getint('mxg0',mxg)
      call getint('mxx0',mxx)
      call getint('mqq0',mqq)
      call getint('mst0',msr)
      call getint('mce0',mce)
      call getint('mbf0',mbf)
      call getint('mky0',mky)
      call getint('nwf0',nwf)
*
      if(ver.lt.170113 )then
         call HF_errlog(231020151, 'F: '//
     1   'Obsolete version of QCDNUM. '//
     2   'Install version 17.01.13 or later.')
      endif
      if(mxg.lt.5      )then
         call HF_errlog(231020152, 'F: '//
     1   'QCDNUM memory allocation insufficient. '//
     2   'Recompile QCDNUM with at least mxg0 = 5 in qcdnum.inc.')
      endif
      if(mxx.lt.300    )then
         call HF_errlog(231020153, 'F: '//
     1   'QCDNUM memory allocation insufficient. '//
     2   'Recompile QCDNUM with at least mxx0 = 300 in qcdnum.inc.')
      endif
      if(mqq.lt.150    )then
         call HF_errlog(231020154, 'F: '//
     1   'QCDNUM memory allocation insufficient. '//
     2   'Recompile QCDNUM with at least mqq0 = 150 in qcdnum.inc.')
      endif
      if(msr.lt.30     )then
         call HF_errlog(231020155, 'F: '//
     1   'QCDNUM memory allocation insufficient. '//
     2   'Recompile QCDNUM with at least mst0 = 30 in qcdnum.inc.')
      endif
      if(mce.lt.20     )then
         call HF_errlog(231020156, 'F: '//
     1   'QCDNUM memory allocation insufficient. '//
     2   'Recompile QCDNUM with at least mce0 = 20 in qcdnum.inc.')
      endif
      if(mbf.lt.10     )then
         call HF_errlog(231020157, 'F: '//
     1   'QCDNUM memory allocation insufficient. '//
     2   'Recompile QCDNUM with at least mbf0 = 10 in qcdnum.inc.')
      endif
      if(mky.lt.50     )then
         call HF_errlog(231020158, 'F: '//
     1   'QCDNUM memory allocation insufficient. '//
     2   'Recompile QCDNUM with at least mky0 = 50 in qcdnum.inc.')
      endif
      if(nwf.lt.1200000)then
         call HF_errlog(231020159, 'F: '//
     1   'QCDNUM memory allocation insufficient. '//
     2   'Recompile QCDNUM with at least nwf0 = 1200000 in qcdnum.inc.')
      endif
*
      call setord(I_FIT_ORDER)         !LO, NLO, NNLO

      call gxmake(xmin_grid,iwt_xgrid,nmxgrid,NXBINS,nx,iosp) !x-grid
      call gqmake(Q2Grid,WQGrid,NQAll,NQ2bins,nqout)             !mu2-grid

      print '(''Requested, actual number of x bins are: '',2i5)', 
     $     nXbins,nx
      print '(''Requested, actual number of Q2 bins are: '',2i5)', 
     $     nQ2bins,nqout

      if (ReportXGrid) then
         open (51, file='xgrid.nml',status='unknown')
         call GXCOPY(xgrid, NXGridMax, NXGridAct)
         write (51,'(''&XGrid'')')
         write (51,'(''  NXgrid = '',I5)') NXGridAct
         write (51,'(''  Q20 = '',F10.2)') starting_scale
         write (51,'(''&End'')')
         write (51,'(10E26.18)') (xgrid(j),j=1,NXGridAct)
         close (51)
C         stop
      endif

      if(ReadXGrid) then 
        call ReadXGridNML()
      endif

      iqc =iqfrmq(qc+tiny)  !> Charm                                                                                                                                                                    
      iqb =iqfrmq(qb+tiny)  !> Bottom  
      iqt =iqfrmq(qt+tiny)  !> Top            

C 7/10/16 Reset top threshold if beyond kinematic limit: 
      if (qt.gt.QARR(2)) then
         iqt = 0.
         call hf_errlog(2016100701,
     $  'I: Top threshold beyond kinematic limit: turn off top PDF')
      endif
      

      if ((mod(HFSCHEME,10).eq.3).or.HFSCHEME.eq.4.or.HFSCHEME.eq.444) then
         call setcbt(3,iqc,iqb,iqt) !thesholds in the ffns
         print *,'Fixed Flavour Number Scheme set with nf=3'
      else
        call setcbt(0,iqc,iqb,iqt) !thesholds in the ffns
      endif

      ierr = 1
      if (Read_QCDNUM_Tables) then
         call readwt(22,'unpolarised.wgt',id1,id2,nw,ierr)
      endif
      if(ierr.ne.0) then
         call fillwt(0,id1,id2,nw) !calculate weights
         call dmpwgt(1,22,'unpolarised.wgt')
      else 
         print*,' Read unpolarised weight file'
      endif
      write(6,'(/'' weight: words used ='',I10)') nw     


! setting of the evolution parameters


c      call SETABR(1.D0,0.D0)  ! mur scale variation
c      call ZMDEFQ2(1.D0,0.D0) ! muf scale variation

      print*,'exit qcdnum_Ini'
      return
 7117 continue
      print '(''Error reading QCDNUM namelist. Stop'')'
      call hf_stop
      end


      subroutine ZMVFNS_init()
C-------------------------------------------------------------
C
C Created 31/12/12. Split FF code initialisation from QCDNUM
C
C-------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "couplings.inc"
      integer nwords, ierr
C-------------------------------------------------------------
      ierr = 1
      if (Read_QCDNUM_Tables) then
         call zmreadw(22,'zmstf.wgt',nwords,ierr)
      endif
      if(MASSH.eq.1) then
      hqscale2inmass=bq2*mch*mch
      elseif(MASSH.eq.2) then
      hqscale2inmass=bq2*mbt*mbt
      endif
      if(mod(HFSCHEME,10).eq.0.and.MASSH.eq.1) then
      print*,'factorisation scale for heavy quarks is set to  sqrt(', hqscale1in,'*Q^2 + ',hqscale2in , '* 4m_c^2 )'   
      elseif(mod(HFSCHEME,10).eq.0.and.MASSH.eq.2) then
      print*,'factorisation scale for heavy quarks is set to  sqrt(', hqscale1in,'*Q^2 + ',hqscale2in , '* 4m_b^2 )'   
      endif
c    
      if(ierr.ne.0) then
         call zmfillw(nwords)
      if ( mod(HFSCHEME,10).eq.0) then
         call ZMDEFQ2(aq2,hqscale2inmass) ! muf scale variation     
      endif 
        call zmdumpw(22,'zmstf.wgt')
      else 
         print*,'Read zmstf weight file'
      endif      
      write(6,'(/'' ZMSTF: words used ='',I10)') nwords      
      call zswitch(IPDFSET)

C-------------------------------------------------------------
      end



      Subroutine Init_EW_parameters
C-----------------------------------------------------
C
C  Initialise electroweak parameters. Created 5 June 2011
C
C  Comments:
C    12/08/2011: Move reading of the namelist to read_steer
C
C-----------------------------------------------------
      implicit none
#include "couplings.inc"
#include "steering.inc"
C     
C set derived values
C
      cos2thw = 1.d0 - sin2thw

       ! move initialzation of the thresholds from read_steer:
      HF_MASS(1) = mch 
      HF_MASS(2) = mbt
      HF_MASS(3) = mtp

C-----------------------------------------------------
      end

      Subroutine ReadXGridNML
#include "ntot.inc"
        double precision grid(NXGridMax)
        double precision Q20
        integer NXgrid
        logical ReadXGrid

        common/ext_xgrid/grid,nxgrid,ReadXGrid
        namelist/XGrid/NXgrid,Q20

        open (51,file='xgrid.nml',status='old')
        read (51,NML=XGrid,ERR=8118,END=1881)

 1881   continue
        
        read (51,'(10E26.18)') (grid(idx),idx=1,NXgrid);
        close (51)
        return 
 8118   continue
        print '(''Error reading XGrid namelist. Stop'')'
        call hf_stop
 8128   continue
        print '(''Error reading xgrid. Stop'')'
        call hf_stop
      end
