c Adapted from LHAPDF uncertainties.f
      subroutine GetPDFUncType_HERAF_lhapdf5(lMonteCarlo,lAsymhess,
     $        lSymmhess, name)
      implicit none
      character*(*) name
      logical lMonteCarlo,lAsymhess,lSymmhess
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
     $        .or.(name(:5).eq.'ABM11')
     $        .or.(name(:5).eq.'abm12')
     $        .or.(name(:6).eq.'ABMP15')
     $        .or.(name(:6).eq.'ABMP16')
     $     .or.(index(name,'EIGSYM').gt.0)
     $        ) then            ! symmetric eigenvector PDF sets
         lSymmhess = .true.
      else                      ! default: assume asymmetric Hessian eigenvector PDF sets
         lAsymhess = .true.
      endif
      end subroutine GetPDFUncType_HERAF_lhapdf5





      subroutine GetPDFUncType_HERAF_lhapdf6(lMonteCarlo,lAsymhess,
     $        lSymmhess,name)
      implicit none

      character*(*) name
      logical lMonteCarlo,lAsymhess,lSymmhess
      ! logical variables for lhapdf interface 
      logical lhapdf_mc, lhapdf_symmetric
      integer nset
      lMonteCarlo = .false.
      lAsymhess = .false.
      lSymmhess = .false.

#ifndef LHAPDF_ENABLED
      call hf_errlog(14081501, "S: Call to lhapdf function but"//
     $      "xFitter compiled without --enable-lhapdf switch")
#else

      call getnset(nset)
      call getpdfunctypem(nset, lhapdf_mc, lhapdf_symmetric)

              if(lhapdf_symmetric) then 
                      lMonteCarlo=.false.
                      lSymmhess=.true.
                      lAsymhess=.false.
              else
                      lMonteCarlo=.false.
                      lSymmhess=.false.
                      lAsymhess=.true.

              endif

              if(lhapdf_mc) then
                      lMonteCarlo=.true.
                      lSymmhess=.false.
                      lAsymhess=.false.
              endif
#endif
      

      end subroutine GetPDFUncType_HERAF_lhapdf6
      




      subroutine GetPDFUncType_HERAF(lMonteCarlo,lAsymhess,lSymmhess
     $     ,name)

              implicit none
#include "steering.inc"
              character*(*) name
              logical lMonteCarlo,lAsymhess,lSymmhess
              character version*32
              lMonteCarlo = .false.
              lAsymhess = .false.
              lSymmhess = .false.
              if(PDF_DECOMPOSITION.eq."LHAPDF") then
#ifndef LHAPDF_ENABLED
             call hf_errlog(26061547, "S: Call to lhapdf function but"//
     $      "xFitter compiled without --enable-lhapdf switch")
#else
                      call getlhapdfversion(version)
                      if(index(version, '5.').eq.1) then
                      call GetPDFUncType_HERAF_lhapdf5(lMonteCarlo,
     $                             lAsymhess, lSymmhess, name)
                      else if(index(version, '6.').eq.1) then
                      call GetPDFUncType_HERAF_lhapdf6(lMonteCarlo,
     $                             lAsymhess, lSymmhess, name)
                      else 
                      call hf_errlog(26061518, "S: lhapdf can not"//
     $                "determine error type")
                      endif
#endif
              else
                 if ( DoBandsSym ) then
                    lSymmhess=.true.
                 else
                    lAsymhess=.true.
                 endif
              endif
      end subroutine GetPDFUncType_HERAF
