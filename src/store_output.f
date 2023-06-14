C--------------------------------------------------------------
client.focused "#4c7899" "#285577" "#ffffff" "#2e9ef4"
C> Function to store outputs
C> \param base name for an output directory
C--------------------------------------------------------------
* --------------------------------------------------------------
      subroutine store_pdfs(base)
* --------------------------------------------------------------

      implicit none

#include "steering.inc"

      integer i,ix,idx,iq2,iflag
      double precision q2,x,gval,sing,umin,dmin
*new jf
      double precision QPDFXQ,cplus,splus,bplus,uplus,dplus,U,D,sea,DbmUb
      double precision d_Ubar,d_Dbar,u_sea,d_sea,str,strbar,chm,bot,photon
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

      character*48 name
      character*48 h1name
      character*(*) base
      character*25 fsfc
      character*25 namefsfc

* base = pdfs_q2val_
      character tag(NBANDS)*2
      data (tag(i),i=1,NBANDS) /'01','02','03','04','05',
     $     '06','07','08','09','10',
     $     '11','12','13','14','15','16','17','18','19','20',
     $     '21','22','23','24','25','26','27','28','29','30',
     $     '31','32','33','34','35','36','37','38','39','40' /

#include "thresholds.inc"

      real*4 glun00h100,ubar00h100,uval00h100,dbar00h100
      real*4 dval00h100,cbar00h100,sbar00h100,bbar00h100

      double precision nf

      double precision pdf
      dimension pdf(-N_CHARGE_PDF:N_CHARGE_PDF+N_NEUTRAL_PDF )


      double precision pdfl,hf_get_alphas,alphas(0:160),q2valpdf(0:160)
     $     ,xvalpdf(0:160)
      dimension pdfl(-N_CHARGE_PDF:N_CHARGE_PDF+N_NEUTRAL_PDF )
      integer iq,jx,j

  ! Store how many PDFs are written out:
      integer NPdfs
      parameter (NPdfs = 15)

C---------------------------------------------------------------

      write(6,*) '--------- in store-pdfs -------'

      idx = index(base,' ')-1
      do 999 iq2=1,NBANDS
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
         ! Write basic info on the table:
         write (81,*) q2val(iq2),outnx, NPdfs, outxrange(1), outxrange(2)

         ! Write the names of PDFs
         write (81,'(16(2x,A12))')
     $        ' x ',' g    ',' U    ',' D    ',' Ubar    ', ' Dbar    ',
     $        ' u_val    ', ' d_val    ', ' sea    ' ,' u_sea    ',
     $        ' d_sea    ', ' str    ',' chm    ',' bot    ', '  ph ',
     $        'strbar'

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
     $           * (ix-1.d0)/(outnx-1.d0)
            x = dexp(x)
*     for integral calculation
            if(ix.gt.1) then
               delx = x - xbelow
            endif
            call  hf_get_pdfs(x,q2,pdf)


            gval=pdf(0)

            if (q2.gt.qc) then
               U=pdf(2)+pdf(-4)
            else
               U=pdf(2)
            endif

            if (q2.gt.qb) then
               D=pdf(1)+pdf(-3)+pdf(-5)
            else
               D=pdf(1)+pdf(-3)
            endif

            umin=pdf(2)-pdf(-2)
            dmin=pdf(1)-pdf(-1)

            if (q2.gt.qc) then
               d_Ubar=pdf(-2)+pdf(-4)
            else
               d_Ubar=pdf(-2)
            endif

            if (q2.gt.qb) then
               d_Dbar=pdf(-1)+pdf(-3)+pdf(-5)
            else
               d_Dbar=pdf(-1)+pdf(-3)
            endif

            sea=d_Ubar+d_Dbar

*      DbmUb=d_Dbar-d_Ubar
            u_sea=pdf(-2)
            d_sea=pdf(-1)
            str = (pdf(-3)+pdf(3))/2.d0
            strbar = pdf(-3)


            chm = 0.0d0
            if (q2.gt.qc.and.abs(pdf(-4)).gt.1D-50) then
               chm=pdf(-4)
            endif

            bot = 0.d0
            if (q2.gt.qb.and.abs(pdf(-5)).gt.1D-50) then
               bot=pdf(-5)
            endif

            photon = pdf(7)

* integral calculation to estimate fs and fc
            totstr = totstr + str*delx
            totDbar = totDbar + d_Dbar*delx
            totcha = totcha + chm*delx
            totUbar = totUbar +d_Ubar*delx
            totusea = totusea + u_sea*delx
            totdsea = totdsea + d_sea*delx

            write(81,810)
     +           x,gval,U,D,d_Ubar,d_Dbar,umin,dmin,sea,u_sea,d_sea,str,
     $           chm,bot,photon,strbar
 810        format(16(2x,G12.6))
 811        format(I3,2x,23(2x,G12.6))


         enddo

         close(81)

 999  continue


      if ( WriteLHAPDF5 ) then

cv store for LHAPDF
cv HERAPDF in LHAPDF5 format
cv  PDFs are Glue Uval Dval Ubar Dbar Str Chrm  Bot
         DO Iq=0,160
            Q2=10**(8.30103/160D0*Iq )
            q2valpdf(iq) = q2
            alphas(iq)=hf_get_alphas(q2)
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
c     v         grid(162+iq)=q2

            DO jx=0,160
               IF(Jx.LE.80)THEN
                  X=10**(6D0/120D0*Jx-6D0)
               ELSE
                  X=10**(2D0/80D0*(Jx-80)-201D-2)
               ENDIF
c     v       grid(1+jx)=x
               call  hf_get_pdfs(x,q2,pdfl)
               write(76,666) PDFl(0), PDFl(2)-PDFl(-2), PDFl(1)-PDFl(-1),
     $              PDFl(-2), PDFl(-1),
     $              PDFl(3), PDFl(4), PDFl(5)


            enddo
         ENDDO

         close (76)

      endif

 666   format(8(2x,G12.6))




      return
      end

C--------------------------------------------------
C> \brief Write results of fit
C> \details Write to a text file binning,
C> data points with uncorrelated and total uncertainties,
C> fitted points and pulls.
C--------------------------------------------------

      subroutine WriteFittedPoints

      implicit none

#include "steering.inc"
#include "ntot.inc"
#include "datasets.inc"
#include "indata.inc"
#include "systematics.inc"
#include "theo.inc"

      integer i,j,index,PlotVarColIdx,PreviousPlots
      double precision PlotVar,PullVar

      integer GetBinIndex

      open(90,file=TRIM(OutDirName)//'/fittedresults.txt')
      write(90,*)ndatasets

      PreviousPlots = 0
C Update theory errors, sum up what is already in and theory sources
      do i=1,NPoints
         THEO_TOT_UP(i) = 0
         THEO_TOT_DOWN(i) = 0
      enddo
      do i=1,NPoints
         do j=1,NSys
            if (ISystType(j).eq.iTheorySyst) then
                if (LAsymSyst(j)) then
                  THEO_TOT_DOWN(i) = sqrt(THEO_TOT_DOWN(i)**2 +
     +                    (THEO(i)*MAX(MAX(BetaAsym(j,1,i),
     +                    BetaAsym(j,2,i)),
     +                    0d0)) ** 2)
                  THEO_TOT_UP(i) = sqrt(THEO_TOT_UP(i)**2 +
     +                    (THEO(i)*MAX(MAX(-BetaAsym(j,1,i),
     +                    -BetaAsym(j,2,i)),
     +                    0d0)) ** 2)
               else          !Symmetric errors
                  THEO_TOT_UP(i) = sqrt(THEO_TOT_UP(i)**2
     +                    + (THEO(i)*Beta(j, i)) ** 2)
                  THEO_TOT_DOWN(i) = sqrt(THEO_TOT_DOWN(i)**2
     +                    + (THEO(i)*Beta(j, i)) ** 2)
               endif
            endif
         enddo
      enddo
      do i=1,ndatasets
         write(90,*)DATASETNUMBER(i)
         write(90,*) DATASETLABEL(i)

         do j=1,GNPlots(i)
            PreviousPlots = PreviousPlots + 1
            write(90,16) 'Plot',j,'@',TRIM(GPlotOptions(PreviousPlots))
         enddo
 16      format(A4,i0,A1,A)

c         write(90,*) '     q2          x        y    data     +- uncorr.err'//
c     &        '   +-toterr      theory      pull     dataset  '

         write (90,17) (DATASETBinNames(j,i),j=1,3),'data    '
     $        ,' +- uncor  ',' +- tot   ',' th orig   ','th mod'
     $        , ' therr+   ', 'therr-'
     $        , ' pull   ', 'iset', 'iplot'
 17      format(1X,11(A11,1X),A4,A12)

         do j=1,NDATAPOINTS(i)
            index = DATASETIDX(i,j)

            PlotVarColIdx = GetBinIndex(i,TRIM(Gplotvarcol(i)))

            if(PlotVarColIdx.eq.0.and.GNPlots(i).eq.0) then
               if(Gplotvarcol(i).eq.'undefined') then
                  call HF_Errlog(13021000,
     $                  'W: Plotting options not set for data set: '
     $                   //DATASETLABEL(i))
               else
                  call HF_Errlog(13012901,
     $                 'W: Plotting: Can not find one of the columns')
               endif
               PlotVar = 0.
            else
               if ( PlotVarColIdx.eq.0) then
                  PlotVar = 0
               else
                  PlotVar = AbstractBins(PlotVarColIdx,index)
               endif
            endif

c set pull to zero if no unc error
            if(ALPHA_MOD(index).gt.0d0) then
               PullVar = (DATEN(index)-THEO_MOD(index))/ALPHA_MOD(index)
            else
               PullVar = 0d0
            endif

            write(90,'(1X,11(e11.5,1X),i12,i4,A1,E11.5)')
     $              AbstractBins(1,index),
     $              AbstractBins(2,index),AbstractBins(3,index),
     &           DATEN(index),ALPHA_MOD(index),
     &           E_TOT(index)/100.*DATEN(index),THEO(index), THEO_MOD(index),
     &           THEO_TOT_UP(index),THEO_TOT_DOWN(index),
     &           PullVar,DATASETNUMBER(i), JPLOT(index), '/',PlotVar

cv
c            write(44,111) VQ2(index),VX(index), f2sh(index),flsh(index),
c     &           xf3sh(index)
         enddo
cv         write(34,*), index,i,DATASETNUMBER(i)
      enddo
  111  format(1X, F10.3, 2X, F12.6, 2X, 3(F12.6,2X))
      close(90)


      RETURN
      END

      subroutine my_system(Command)
C
C Hack to resolve "system" name conflict
C
      character *(*) command
      call system(Command)
      end

C-------------------------------------------------------------
C> Write fit results for free parameters
C-------------------------------------------------------------
      subroutine write_pars(ifcn3)
C-------------------------------------------------------------
C Extra output of PDF parameters
C-------------------------------------------------------------
      implicit none
#include "fcn.inc"
#include "endmini.inc"
#include "steering.inc"
      integer ifcn3

      integer i
      double precision val,err,xlo,xhi
      integer ipar
      character*32 parname
      character*300 fname

      double precision, allocatable :: errIterate(:,:)


C-------------------------------------------------------------
      if (ifcn3.lt.10) then
c RP         write (fname,'(''output/parsout_'',i1)') ifcn3
         write (fname,'( a,''/parsout_'',i1)') TRIM(OutDirName),ifcn3
      elseif (ifcn3.lt.100) then
         write (fname,'( a,''/parsout_'',i2)') TRIM(OutDirName),ifcn3
      elseif (ifcn3.lt.1000) then
         write (fname,'( a,''/parsout_'',i3)') TRIM(OutDirName),ifcn3
      elseif (ifcn3.lt.10000) then
         write (fname,'( a,''/parsout_'',i4)') TRIM(OutDirName),ifcn3
      endif

      open (71,file=fname,status='unknown')
! Broken since 2.2.0
!     if (DoBands .and. ifcn3.eq.0) then
!        Allocate(errIterate(MNE,MNE))
!        call GetErrMatScaled(errIterate)
!     endif

      do i=1,mne
         call mnpout(i,parname,val,err,xlo,xhi,ipar)

C
C For bands, replace by "iterate" estimate, if present
C     Broken since 2.2.0
!        if ( Dobands .and. ipar.gt.0 .and. ifcn3.eq.0 ) then
!           if ( errIterate(ipar,ipar).gt.0 ) then
!              err = sqrt(errIterate(ipar,ipar))
!              val = pkeep(i)
!              call hf_errlog(1060402016,
!    $ 'I: Write uncertainties to parsout_0 using Iterate method')
!           endif
!        endif

         if (Trim(parname).ne.'undefined') then
            if (xlo.eq.0.and.xhi.eq.0) then
               write (71,72) i, Trim(parname), val,err
            else
               write (71,72) i, Trim(parname), val,err,xlo,xhi
            endif
         endif
      enddo
! Broken since 2.2.0
!     if (DoBands .and. ifcn3.eq.0) then
!        deallocate(errIterate)
!     endif
 72   format (I5,'   ','''',A,'''',4F12.6)
      close(71)

C-------------------------------------------------------------
      end

C-------------------------------------------------------------------
C
C> Find minuit train which has the best chi2 for the control sample.
C
C--------------------------------------------------------------------
      Subroutine FindBestFCN3

      implicit none
#include "endmini.inc"
#include "steering.inc"
      integer i,iminCont
      double precision aminCont


      double precision val,err,xlo,xhi
      integer ipar
      character*32 parname

C-------------------------------------------------------------------
      aminCont = 1.D30
      do i=1,nfcn3
         if ( chi2cont3(i).lt. aminCont) then
            aminCont = chi2cont3(i)
            iminCont = i
         endif
      enddo
C-------------------------------------------------------------------
      print *,' '
      print *,' '
      print *,' '
      print *,'======================================================'
      print '(''  Use NNPDF overfitting method.
     $   Prepare output PDF files '')'
      print '(''  Best FCN3 call='',i4,'' out of '',i4,'' calls'')',
     $     iminCont,nfcn3
      print '(''  Chi2 control best='',F10.4)',aminCont
      print *,'======================================================'
      print *,' '
      print *,' '
      print *,' '

      ! Dump PDFs for this:
!Broken since 2.2.0
!     call PDF_param_iteration(pkeep3(1,iminCont),2) !Decode params.
C
C Fix some pars by sum-rules:
C
C     call Evolution !I'm commenting this out because it conflicts with
C     the new evolution interface, but I do not know why it was here
C     --Ivan

C ! Ready to store:
cv      open (76,file='output/lhapdf.block.txt',status='unknown')
cv      call store_pdfs('output/pdfs_q2val_')
C store the optimal values

      open (76,file=TRIM(OutDirName)//'/opt_lhapdf.block.txt',
     &   status='unknown')
      call store_pdfs(TRIM(OutDirName)//'/opt_pdfs_q2val_')
      call print_lhapdf6_opt()


      open (71,file=TRIM(OutDirName)//'/parseout_opt',status='unknown')
      do i=1,mne
         call mnpout(i,parname,val,err,xlo,xhi,ipar)
         if (Trim(parname).ne.'undefined') then
            if (xlo.eq.0.and.xhi.eq.0) then
               write (71,72) i, Trim(parname), pkeep3(i,iminCont),err
            else
               write (71,72) i, Trim(parname), pkeep3(i,iminCont),err,xlo,xhi
            endif
         endif
      enddo
 72   format (I5,'   ','''',A,'''',4F12.6)
      close(71)


      end


C----------------------------------------------------
C> The subroutine writes the theoretical values of
C> the cross sections according fo the input bins.
C> The subroutine is used for integration with HERAverager
C----------------------------------------------------
C The functionality is constrained to H1 ZEUS data
C----------------------------------------------------
C> \author P.Belov
c> \date   19/11/2012
C----------------------------------------------------
      subroutine WriteCSforAverager

      implicit none

#include "steering.inc"
#include "ntot.inc"
#include "datasets.inc"
#include "indata.inc"
#include "systematics.inc"
#include "theo.inc"

      integer i,j,index,k,reacindx

      double precision currEcharge

      open(90,file='./'//TRIM(OutDirName)//'/heraverager.dat')
C      write(90,*)ndatasets




        write(90,*) '!* '
        write(90,*)
     $  '!* Swimming set from XFITTER for the HERAverager'
        write(90,*) '&Data'
        write(90,*) '  Name = ''Swimming'' '
        write(90,*) '  NData = ',NPOINTS
        write(90,*) '  NColumn = 5'
        write(90,*) '  ColumnType = 4*''Bin'',''Sigma'' '
        write(90,*) '  ColumnName = ''reaction index'', ''x'', ''Q2'', ''y '',' //
     & ' ''reduced x-section '' '
        write(90,*) '                                                    '
        write(90,*) '  IndexDataset = 666'
        write(90,*) '  Reaction = ''NC e+-p'' '
        write(90,*) '                                                    '
        write(90,*) '&END'

      do i=1,ndatasets
          reacindx = 0
          currEcharge = 0
          do k=1,DATASETInfoDimension(i)
            if(DATASETInfoNames(k,i).eq.'e charge') then
              currEcharge = DATASETInfo(k,i)
            endif
          enddo
          if ((DATASETREACTION(i).eq.'NC e+-p').or.
     $      (DATASETREACTION(i).eq.'NC e+-p Dummy')) then
            if (currEcharge.gt.0) then
              reacindx = 615
            else if (currEcharge.lt.0) then
              reacindx = 515
            endif
          else if ((DATASETREACTION(i).eq.'CC e+-p').or.
     $      (DATASETREACTION(i).eq.'CC e+-p Dummy')) then
            if (currEcharge.gt.0) then
              reacindx = 3615
            else if (currEcharge.lt.0) then
              reacindx = 3515
            endif
          endif
          do j=1,NDATAPOINTS(i)
             index = DATASETIDX(i,j)

             write(90,'(1X,i5,1X,4(e11.5,1X),i4)')
     $              reacindx,
     $              AbstractBins(1,index),
     $              AbstractBins(2,index),
     $              AbstractBins(3,index),
     $              THEO(index)
          enddo

      enddo
  111  format(1X, F10.3, 2X, F12.6, 2X, 3(F12.6,2X))
      close(90)

C      RETURN
      end subroutine


C----------------------------------------------------------------
C> \brief Write theory prediction in format of input data tables.
C> \param NNuisance number of error sets
C> \param Theo_cent central value of theory predction
C> \param SymmetricPDFErr use symmetric or assymmetric errros (beta vs betaasym)
C----------------------------------------------------------------
      Subroutine WriteTheoryFiles(NNuisance,Theo_cent,SymmetricPDFErr)
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "systematics.inc"
#include "datasets.inc"
#include "indata.inc"

      integer NNuisance
      double precision Theo_cent(Ntot)
      logical SymmetricPDFErr
      integer iset, ipoint, j, i

      character*2 c

C---------------------------------------------------------

      ! Loop over data sets

      do iset=1, NDataSets
         if (iset.lt.10) then
            write (c,'(''0'',I1)') iset
         else
            write (c,'(I2)') iset
         endif
         open (51
     $        ,file=Trim(OutDirName)//'/theo_'//c//'.dat'
     $        ,status='unknown')

         ! Write  a header
         write (51,'(''* Theory file for '',A)')
     $        Trim(DATASETLABEL(iset))

         write (51,'(''&Data '')')
         write (51,'(''   Name = "Theory for '',A,''"'')')
     $        Trim(DATASETLABEL(iset))
         write (51,'(''   NData = '',I5)') NDATAPOINTS(iset)

         if (SymmetricPDFErr) then
            write (51,'(''   NColumn = '',I5)') NNuisance+1
     $        + DATASETBinningDimension(iset)

            write (51,
     $'(''   ColumnType = '',I1,''*"Bin","Theory",'',i3,''*"Error"'')')
     $       DATASETBinningDimension(iset), NNuisance

            write (51,'(''   ColumnName = '',200(''"'',A,''",''))'
     $           ,advance='no' )
     $           ( trim(DATASETBinNames(i,iset)),
     $           i=1,DATASETBinningDimension(iset) ), 'theory',
     $           ( trim(System(nsys+i)),i=1,NNuisance-1)
            write (51,'(A,''"'')')
     $           ( trim(System(nsys+i)),i=NNuisance,NNuisance)
            write (51,'(''   Percent = '',I3,''*True'')') NNuisance
         else
            write (51,'(''   NColumn = '',I5)') NNuisance*2+1
     $        + DATASETBinningDimension(iset)
            write (51,
     $'(''   ColumnType = '',I1,''*"Bin","Theory",'',i3,''*"Error"'')')
     $       DATASETBinningDimension(iset), NNuisance*2
            write (51,'(''   ColumnName = '',200(''"'',A,''",''))'
     $           ,advance='no' )
     $           ( trim(DATASETBinNames(i,iset)),
     $           i=1,DATASETBinningDimension(iset) ), 'theory',
     $           ( trim(System(nsys+i))//'+',
     $           trim(System(nsys+i))//'-',i=1,NNuisance-1)
            write (51,'(A,''",'',''"'',A,''"'')')
     $           ( trim(System(nsys+i))//'+',
     $           trim(System(nsys+i))//'-',i=NNuisance,NNuisance)
            write (51,'(''   Percent = '',I3,''*True'')') NNuisance*2
         endif

         write (51,'(''&End '')')

         do i = 1, NDATAPOINTS(iset)
            ipoint = Datasetidx(iset,i)
            if (SymmetricPDFErr) then
               write (51,'(200E12.4)')
     $  ( AbstractBins(j,ipoint),j=1,DATASETBinningDimension(iset)),
     $           theo_cent(ipoint),
     $  ( -Beta(j,ipoint)*100.0,               ! negative sign, since it is inverted in lhapdferrors.cc
     $           j=NSys+1
     $           ,NSys+NNuisance)
            else
               write (51,'(200E14.6)')
     $  ( AbstractBins(j,ipoint),j=1,DATASETBinningDimension(iset)),
     $           theo_cent(ipoint),
     $  ( -BetaAsym(j,1,ipoint)*100.0, -BetaAsym(j,2,ipoint)*100.,
     $           j=NSys+1
     $           ,NSys+NNuisance)
            endif
         enddo
         close (51)
      enddo

C---------------------------------------------------------

      end
