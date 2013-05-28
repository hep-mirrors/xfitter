
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

      include 'thresholds.inc'

      real*4 glun00h100,ubar00h100,uval00h100,dbar00h100
      real*4 dval00h100,cbar00h100,sbar00h100,bbar00h100

      double precision nf

      double precision pdf
      dimension pdf(-6:6)


      double precision pdfl,hf_get_alphas,alphas(0:160),q2valpdf(0:160)
     $     ,xvalpdf(0:160)
      dimension pdfl(-6:6)
      integer iq,jx,j

  ! Store how many PDFs are written out:
      integer NPdfs 
      parameter (NPdfs = 13)

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
         write (81,'(14(2x,A12))')
     $        ' x ',' g    ',' U    ',' D    ',' Ubar    ', ' Dbar    ',
     $        ' u_val    ', ' d_val    ', ' sea    ' ,' u_sea    ',
     $        ' d_sea    ', ' str    ',' chm    ',' bot    '

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

            write(81,810)
     +           x,gval,U,D,d_Ubar,d_Dbar,umin,dmin,sea,u_sea,d_sea,str,chm,bot
 810        format(14(2x,G12.6))
 811        format(I3,2x,23(2x,G12.6))


         enddo

         close(81)
      
 999  continue




cv store for LHAPDF
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
cv         grid(162+iq)=q2

         DO jx=0,160
            IF(Jx.LE.80)THEN
               X=10**(6D0/120D0*Jx-6D0)
            ELSE
               X=10**(2D0/80D0*(Jx-80)-201D-2)
            ENDIF
cv       grid(1+jx)=x
            call  hf_get_pdfs(x,q2,pdfl)
            write(76,666) PDFl(0), PDFl(2)-PDFl(-2), PDFl(1)-PDFl(-1),
     $           PDFl(-2), PDFl(-1),
     $           PDFl(3), PDFl(4), PDFl(5)


         enddo
      ENDDO

      close (76)

 666   format(8(2x,G12.6))




      return
      end

*   -------------------------------------------------
*   -------------------------------------------------


      SUBROUTINE WRITEFITTEDPOINTS

      implicit none
      
      include 'steering.inc'
      include 'ntot.inc'
      include 'datasets.inc'
      INCLUDE 'indata.inc'
      include 'systematics.inc'
      INCLUDE 'theo.inc'
      
      integer i,j,index,PlotVarColIdx,PreviousPlots
      double precision PlotVar

      integer GetBinIndex

      open(90,file=TRIM(OutDirName)//'/fittedresults.txt')
      write(90,*)ndatasets

      PreviousPlots = 0

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
     $        , ' pull   ', 'iset', 'iplot'
 17      format(1X,9(A11,1X),A4,A12)

         do j=1,NDATAPOINTS(i)
            index = DATASETIDX(i,j)

            PlotVarColIdx = GetBinIndex(i,TRIM(Gplotvarcol(i)))

            if(PlotVarColIdx.eq.0) then
               if(Gplotvarcol(i).eq.'undefined') then
                  call HF_Errlog(13021000,
     $                 'W:Plotting options not set for data set')
               else
                  call HF_Errlog(13012901,
     $                 'W:Plotting: Can not find one of the columns')
               endif
               PlotVar = 0.
            else 
               PlotVar = AbstractBins(PlotVarColIdx,index)
            endif

            write(90,'(1X,9(e11.5,1X),i4,i4,A1,E11.5)') 
     $              AbstractBins(1,index),
     $              AbstractBins(2,index),AbstractBins(3,index),
     &           DATEN(index),ALPHA_MOD(index),
     &           E_TOT(index)/100.*DATEN(index),THEO(index), THEO_MOD(index),
     &           (DATEN(index)-THEO_MOD(index))/ALPHA_MOD(index),
     &           DATASETNUMBER(i), JPLOT(index), '/',PlotVar

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

      subroutine get_lhapdferrors
C---------------------------------
      implicit none
C
      include 'ntot.inc'
      include 'couplings.inc'
      include 'steering.inc'
      include 'alphas.inc'
      include 'fcn.inc'
      include 'theo.inc'
      include 'indata.inc'
      include 'systematics.inc'

      integer iset
      integer nsets
      double precision chi2tot
C Function:
      double precision chi2data_theory
      double precision alphasPDF 
      character*4 c

C Store central, plus/minus theory predictions:
      double precision theo_cent(NTOT), 
     $     theo_minus(NTOT), theo_plus(NTOT)
      integer nsysLoc

C Some hack to store PDFs
      character*48 name
      character*48 base
      integer i,idx,j
      character*3 tag


C-----------------------------------------------------------
      nsets = nLHAPDF_Sets+1
C      nsets = 45
      nsysLoc = nsys
      print *,'Nsets=',nsets
      
      do iset=0,nLHAPDF_Sets 
         call InitPDF(iset)
         alphas = alphasPDF(Mz)
         if (iset.eq.0) then
            ! Initialise:
            chi2tot = chi2data_theory(1)      
C Store theory prediction:
            do i=1,Npoints
               THEO_CENT(i) = THEO(i)
            enddo
         else
            chi2tot = chi2data_theory(min(2,iset))
C Store theory predictions for  up/down variations:
            if ( mod(iset,2).eq.0) then
               do i=1,Npoints
                  THEO_PLUS(i) = THEO(i)
               enddo
            else
               do i=1,Npoints
                  THEO_MINUS(i) = THEO(i)
               enddo
            endif

C Also add "systematics" for these sources:
            if ( mod(iset,2).eq.0) then
               nsysLoc = nsysLoc + 1
               if ( nsysLoc - nsys .ge. 10 ) then
                  write (tag,'(''_'',i2)')  nsysLoc - nsys
               else
                  write (tag,'(''_0'',i1)')  nsysLoc - nsys
               endif

               SYSTEM(nsysLoc) = 'PDF_nuisance_param'//tag
               do i=1,Npoints
                  Beta(nsysLoc,i) = 
     $                 (THEO_PLUS(i)-THEO_MINUS(i))/2.0D0
     $                 / THEO_CENT(i)

                  if (Scale68) then
                     Beta(nsysLoc,i) = Beta(nsysLoc,i) / 1.64
                  endif
               enddo
            endif

         endif
         print '(''Got PDF set='',i5,'' chi2='',F10.1,'' ndf='',i5)',
     $        iset,chi2tot,ndfmini
         write(86,*) iset, ' ', chi2tot/ndfmini
         
         if (iset.lt.10) then
            write (c,'(''000'',I1)') iset
         else
            write (c,'(''00'',I2)') iset
         endif


         call WriteFittedPoints
         call my_system
     $ ('mv '//TRIM(OutDirName)//'/fittedresults.txt '
     & //TRIM(OutDirName)//'/fittedresults.txt_set'//c)

         if (iset.eq.0) then
            name = TRIM(OutDirName)//'/pdfs_q2val_'
         else
            i = (iset-1) / 2 + 1

            if (i.lt.10) then
               write (tag,'(''s0'',i1)') i
            elseif (i.lt.100) then
               write (tag,'(''s'',i2)') i
            endif

            base = TRIM(OutDirName)//'/pdfs_q2val_'//tag
            idx = index(base,' ')-1
            if ( mod(iset,2).eq.1) then
               name = base(1:idx)//'m_'
            else
               name = base(1:idx)//'p_'
            endif
         endif
         call store_pdfs(name)

      enddo

C Write some chi2 tests, vs central prediction or floating eighenvectors:
      Call initpdf(0)

      open (85,file=TRIM(OutDirName)//'/Results_00.txt'
     $     ,status='unknown')
      chi2tot = chi2data_theory(3)        
      close (85)

      do i=Nsys+1,NSysLoc
         n_syst_meas(i) = Npoints
         do j=1,Npoints
            syst_meas_idx(j,i) = j
         enddo
         SysScalingType(i) = isLinear ! scale theory
      enddo
      NSys = Nsysloc
      ResetCommonSyst = .true.
      
      print *,nsys
      open (85,file=TRIM(OutDirName)//'/Results_fitPDF.txt'
     $     ,status='unknown')
      chi2tot = chi2data_theory(3)        
      close (85)



      end

C--------------------------------------------------------------
      subroutine write_pars(ifcn3)
C-------------------------------------------------------------
C Extra output of PDF parameters
C-------------------------------------------------------------
      implicit none
      include 'fcn.inc'
      include 'endmini.inc'
      include 'steering.inc'
      integer ifcn3

      integer i
      double precision val,err,xlo,xhi
      integer ipar
      character*32 parname
      character*32 fname

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
      do i=1,mne
         call mnpout(i,parname,val,err,xlo,xhi,ipar)
         if (Trim(parname).ne.'undefined') then
            if (xlo.eq.0.and.xhi.eq.0) then
               write (71,72) i, Trim(parname), val,err
            else
               write (71,72) i, Trim(parname), val,err,xlo,xhi
            endif
         endif
      enddo
 72   format (I5,'   ','''',A,'''',4F12.6)
      close(71)

C-------------------------------------------------------------
      end

      Subroutine FindBestFCN3
C-------------------------------------------------------------------
C
C Find minuit train which has the best chi2 for the control sample.
C
C--------------------------------------------------------------------
      implicit none
      include 'endmini.inc'
      include 'steering.inc'
      integer i,iminCont, kflag
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

      !> Dump PDFs for this:
      call PDF_param_iteration(pkeep3(1,iminCont),2) !> Decode params.
C
C Fix some pars by sum-rules:
C
      kflag = 0
      call SumRules(kflag)
      call Evolution

C !> Ready to store: 
cv      open (76,file='output/lhapdf.block.txt',status='unknown')
cv      call store_pdfs('output/pdfs_q2val_')
C store the optimal values 

      open (76,file=TRIM(OutDirName)//'/opt_lhapdf.block.txt',
     &   status='unknown')
      call store_pdfs(TRIM(OutDirName)//'/opt_pdfs_q2val_')


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


      subroutine WriteCSforAverager
C----------------------------------------------------
C The subroutine writes the theoretical values of
C the cross sections according fo the input bins.
C The subroutine is used for integration with HERAverager
C----------------------------------------------------
C The functionality is constrained to H1 ZEUS data
C----------------------------------------------------
C P.Belov 19/11/2012
C----------------------------------------------------
      implicit none

      include 'steering.inc'
      include 'ntot.inc'
      include 'datasets.inc'
      INCLUDE 'indata.inc'
      include 'systematics.inc'
      INCLUDE 'theo.inc'

      integer i,j,index

      open(90,file='./'//TRIM(OutDirName)//'/heraverager.dat')
C      write(90,*)ndatasets




        write(90,*) '!* '
        write(90,*) '!* Swimming set from HERAFITTER for the HERAverager'
        write(90,*) '&Data'
        write(90,*) '  Name = ''Swimming'' '
        write(90,*) '  NData = ',NPOINTS
        write(90,*) '  NColumn = 4'
        write(90,*) '  ColumnType = 3*''Bin'',''Sigma'' '
        write(90,*) '  ColumnName = ''x'', ''Q2'', ''y '',' //
     & ' ''reduced x-section '' '
        write(90,*) '                                                    '
        write(90,*) '  IndexDataset = 666'
        write(90,*) '  Reaction = ''NC e+-p'' '
        write(90,*) '                                                    '
        write(90,*) '&END'

      do i=1,ndatasets

          do j=1,NDATAPOINTS(i)
             index = DATASETIDX(i,j)
             write(90,'(1X,9(e11.5,1X),i4)') 
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
