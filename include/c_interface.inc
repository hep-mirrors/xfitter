      double precision xgrid(NXGridMax)
      integer NXgrid
      logical ReadXGrid
C
      common/ext_xgrid/xgrid,nxgrid,ReadXGrid

      !------------------------------------------------ c interface

      double precision c_Mz

      character*128 c_LHAPDF6OutDir
      character*128 c_OutDirName

      double precision c_xgrid(NXGridMax)
      integer c_nx
      logical c_read_xgrid

      integer c_nvarl(200) ! see d506cm.inc in minuit

      logical c_useGridLHAPDF5
      logical c_WriteAlphaSToMemberPDF
      logical c_lead
      integer c_ipdfset
      logical c_writeLHAPDF6
      integer c_i_fit_order
      real*4  c_hf_mass(3)
      logical c_dobands

      common/CCommonInterface/ c_mz, c_LHAPDF6OutDir, 
     $ c_OutDirName, c_xgrid, c_nx, c_read_xgrid,
     $  c_nvarl, c_dobands, c_hf_mass, c_i_fit_order, c_ipdfset, c_lead, 
     $  c_useGridLHAPDF5, c_writeLHAPDF6, c_WriteAlphaSToMemberPDF

     
