c Adapted from LHAPDF uncertainties.f
      subroutine GetPDFUncType_HERAF(name,
     $     lMonteCarlo,lAsymhess,lSymmhess)
      implicit none
      logical lMonteCarlo,lAsymhess,lSymmhess
      character*30 name
      lMonteCarlo = .false.
      lAsymhess = .false.
      lSymmhess = .false.
      if ((name(:5).eq.'NNPDF')
     $     .or.(name(:7).eq.'Alekhin')
     $     .or.(name(:5).eq.'Botje')
     $     .or.(name(:5).eq.'Fermi')
     $     .or.(index(name,'_MC').gt.0)
     $     .or.(index(name,'MC_').gt.0)
     $     .or.(index(name,'-MC').gt.0)
     $     .or.(index(name,'MC-').gt.0)
     $     ) then               ! Monte Carlo PDF sets
         lMonteCarlo = .true.
      else if ((name(:4).eq.'A02M').or.(name(:4).eq.'a02m')
     $        .or.(name(:6).eq.'ABKM09').or.(name(:6).eq.'abkm09')
     $        .or.(name(:5).eq.'ABM11').or.(name(:5).eq.'abm11')
     $     .or.(index(name,'EIGSYM').gt.0)
     $        ) then            ! symmetric eigenvector PDF sets
         lSymmhess = .true.
      else                      ! default: assume asymmetric Hessian eigenvector PDF sets
         lAsymhess = .true.
      endif
      end subroutine GetPDFUncType_HERAF
