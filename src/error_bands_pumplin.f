
      subroutine Error_Bands_Pumplin

      implicit none

      include 'steering.inc'
      include 'pdfparam.inc'
      include 'endmini.inc'
      include 'alphas.inc'
      include 'thresholds.inc'
      include 'ntot.inc'
      include 'systematics.inc'
      include 'g_offset.inc'	

      integer shift_dir
      double precision a
      dimension a(MNE)

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

      integer mpar

      double precision shift
C Function
      double precision GetUmat
      ! double precision DecorVarShift

C---------------------------------------------------------------
      
C
C  Fix relation between internal and external params.
C
      mpar = 0



      do ind=1,MNE
         call mnpout(ind,parname,parval,parerr,parlolim,
     $        parhilim,iunint(ind))
         if (iunint(ind).gt.0) then
            write (6,*) 'Parameter',ind,' name=',parname
            write (6,*) 'Internal index=',iunint(ind)
            write (6,*) ' '
            mpar = mpar + 1
         endif
      enddo

      print *,mpar,' variable parameters'


      do ind=1,mpar
         do ind2=1,MNE
            if (iunint(ind2).eq.ind) then
               iexint(ind) = ind2
            endif
         enddo
      enddo

      print *,'Relation internal->external index:'

      do ind=1,mpar
         write (6,*) 'internal=',ind,' external=',iexint(ind)
      enddo


      npar = MNE !> npar runs over external parameters.


C
C Loop over de-correlated (diagonalised) errors:
C
      do j=1,mpar
  
         jext = iexint(j)

         base = TRIM(OutDirName)//'/pdfs_q2val_'//tag(j)
         idx = index(base,' ')-1
         
         base2 = TRIM(OutDirName)//'/pdfs_'//tag(j)
         idx2  = index(base2,' ')-1

         do shift_dir=-1,1,2
            if (shift_dir.eq.-1) then
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


C
C Shift variable paramters by the j-th de-correlated error:
C
            do i=1,npar
               a(i) = pkeep(i) 
               iint = iunint(i)
               if (iint.gt.0) then
                  if(doOffset) then
                    shift = shift_dir * DecorVarShift(iint, j)
                  else
                    shift = shift_dir * GetUmat(iint,j)
                  endif
                  a(i) = a(i) + shift
               endif
            enddo  ! i


C
C Decode "a". 2 stands for IFLag = 2, which is a normal iteration.
C
            call PDF_param_iteration(a,2)
            
C
C Fix some pars by sum-rules:
C
            kflag = 0
            call SumRules(kflag)
            call Evolution

C
C Write results out:
C
            open (76,file=name2,status='unknown')
            call store_pdfs(name)
            close (76)
            if(shift_dir.eq.-1) then
              call save_data_lhapdf6(j*2-1)
            else
              call save_data_lhapdf6(j*2)
            endif
            

         enddo  ! shift_dir

      enddo  ! j

      return
      end


