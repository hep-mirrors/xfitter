C----------------------------------------------------------------------
C> @file xfitter_gn_bands.f
C> @brief Output of the Pumplin eigenvector sets for the GaussNewton minimizer
C>
C> The output part of Error_Bands_Pumplin (src/error_bands_pumplin.f),
C> without the MINUIT calls: the module sets the parameters of each
C> member, these routines write the same files, LHAPDF members and
C> asymmetric theory uncertainties as for MINUIT doErrors: Pumplin.
C----------------------------------------------------------------------

C> @brief NTOT, the size of the theory arrays
      integer function gn_ntot()
      implicit none
#include "ntot.inc"
      gn_ntot = NTOT
      end

C----------------------------------------------------------------------
C> @brief One member: eigenvector j, direction idir (-1 or +1), at the
C>        parameters already set in gParameters
C> @param[inout] TheoVars theory of the members, NTOT x 2 x mpar
C----------------------------------------------------------------------
      subroutine gn_pumplin_member(j, idir, mpar, TheoVars)
      implicit none
#include "steering.inc"
#include "endmini.inc"
#include "ntot.inc"
#include "systematics.inc"
#include "fcn.inc"
#include "theo.inc"
      integer j, idir, mpar
      double precision TheoVars(NTOT,2,mpar)
      integer i, idx, idx2
      character*48 name, name2
      character*300 base, base2
      character tag(40)*3
      data (tag(i),i=1,40) /'s01','s02','s03','s04','s05',
     +     's06','s07','s08','s09','s10',
     +     's11','s12','s13','s14','s15',
     +     's16','s17','s18','s19','s20',
     +     's21','s22','s23','s24','s25',
     +     's26','s27','s28','s29','s30',
     +     's31','s32','s33','s34','s35',
     +     's36','s37','s38','s39','s40'/
      double precision chichi
      double precision chi2data_theory ! function
C---------------------------------------------------------------------
      base = TRIM(OutDirName)//'/pdfs_q2val_'//tag(j)
      idx = index(base,' ')-1
      base2 = TRIM(OutDirName)//'/pdfs_'//tag(j)
      idx2  = index(base2,' ')-1
      if (idir.eq.-1) then
         name  = base(1:idx)//'m_'
         name2 = base2(1:idx2)//'m.lhgrid'
      else
         name  = base(1:idx)//'p_'
         name2 = base2(1:idx2)//'p.lhgrid'
      endif

      ifcncount = ifcncount+1
      chichi = chi2data_theory(2)
      TheoVars(:,(idir+1)/2+1,j) = THEO

      open (76,file=name2,status='unknown')
      call store_pdfs(name)
      close (76)
      if (idir.eq.-1) then
         call save_data_lhapdf6(j*2-1)
         call error_band_action(j*2-1)
      else
         call save_data_lhapdf6(j*2)
         call error_band_action(j*2)
      endif
      end

C----------------------------------------------------------------------
C> @brief After the members: restore the theory of the final fcn call
C>        and write the asymmetric theory uncertainties
C----------------------------------------------------------------------
      subroutine gn_pumplin_finish(mpar, TheoVars)
      implicit none
#include "ntot.inc"
#include "systematics.inc"
#include "theo.inc"
      integer mpar
      double precision TheoVars(NTOT,2,mpar)
C---------------------------------------------------------------------
      Theo = TheoFCN3
      Theo_mod = TheoModFCN3
      ALPHA_Mod = ALphaModFCN3
      call GetTheoErrorsAsym(TheoVars,ntot,mpar)
      call writefittedpoints
      end
