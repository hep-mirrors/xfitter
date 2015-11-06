      subroutine fill_c_common
#include "d506cm.inc"
#include "ntot.inc"
#include "c_interface.inc"
#include "steering.inc"
#include "couplings.inc"

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

        c_dobands= dobands
        c_hf_mass(1)= hf_mass(1)
        c_hf_mass(2)= hf_mass(2)
        c_hf_mass(3)= hf_mass(3)
        c_i_fit_order= i_fit_order
        c_ipdfset= ipdfset
        c_lead= lead
        c_useGridLHAPDF5= useGridLHAPDF5
        c_writeLHAPDF6= writeLHAPDF6
        c_extrapdfs   = extrapdfs
      end

      function get_nmembers()
#include "d506cm.inc"
#include "steering.inc"
              integer get_nmembers, nset, i
              get_nmembers=0

              if(PDF_DECOMPOSITION.eq."LHAPDF") then
#ifndef LHAPDF_ENABLED
            call hf_errlog(29061521, "S: Call to lhapdf function but"//
     $      "HERAfitter compiled without --enable-lhapdf switch")
#else
                      call getnset(nset)
                      call numberpdfm(nset, get_nmembers)
                      if(get_nmembers.gt.1) then
                              get_nmembers=get_nmembers+1
                      endif
#endif
              else
                      get_nmembers=1 !central value
                      
                      if(dobands) then
                      do i=1,MNE
                              if(nvarl(i).eq.1) then
                                      get_nmembers=get_nmembers+2
                              endif
                      enddo
                      endif
              endif
              end 
