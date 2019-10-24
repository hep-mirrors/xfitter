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
