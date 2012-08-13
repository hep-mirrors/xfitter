C =========================================================
C === THE FILE CONTAINS PARAMETERS FOR THE DIPOLE MODEL ===
C === THE FILE IS CALLED BY DIPOLE MODEL SUBROUTINES    ===
C =========================================================
      
      double precision sig0,xlam,x0,xm,cBGK,eBGK
      
      common/pass1/ sig0,xlam,x0,xm,cBGK,eBGK
      
      
      double precision cm
      integer icharm
      
      common/chsteer/ cm,icharm

C
C     Namelist of dipole type variable
C
      character*32 TypeOfDipoleFit
      namelist/DipoleType/TypeOfDipoleFit
