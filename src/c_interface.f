      subroutine fill_c_common
#include "d506cm.inc"
#include "ntot.inc"
#include "c_interface.inc"
#include "steering.inc"
#include "couplings.inc"
#include "thresholds.inc"

        c_mz= mz

        c_LHAPDF6OutDir= LHAPDF6OutDir
        c_OutDirName= OutDirName
        c_lhapdfset=lhapdfset
        c_WriteAlphaSToMemberPDF= WriteAlphaSToMemberPDF

        do idx=1,NXGridMax
          c_xgrid(idx)= xgrid(idx)
        enddo

        c_nx= nxgrid
        c_read_xgrid= ReadXGrid

        c_itheory = itheory

        !c_dobands= dobands Broken since 2.2.0
        c_hf_mass(1)= hf_mass(1)
        c_hf_mass(2)= hf_mass(2)
        c_hf_mass(3)= hf_mass(3)
        c_i_fit_order= i_fit_order
        c_ipdfset= ipdfset
        c_lead= lead
        c_useGridLHAPDF5= useGridLHAPDF5
        c_writeLHAPDF6= writeLHAPDF6
        c_extrapdfs   = extrapdfs
        c_kmuc = kmuc
        c_kmub = kmub
        c_kmut = kmut

      end
