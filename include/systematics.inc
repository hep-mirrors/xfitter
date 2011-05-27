      integer NSYS
* new2 jf : mandy test , to decorrelate all errors in H1Z data files
C--- 06/08/2010: NSYS was increased from 550 up to 559 for 'tab_Ep_460_575.dat'
C--- 13/08/2010: remark: NTOT should be 574 in case of 'tab.dat' for Dataset 24
C--- in order to enable systematics for Dataset 24
      parameter (NSYS = 600)

      character*5 SYSTEM
      dimension SYSTEM(NSYS)
      double precision BETA(NSYS,NTOT)
      double precision ALPHA(NTOT)
      double precision sysa(nsys,nsys)
      integer npoints,ndis
      common/systema/ALPHA,BETA,sysa,npoints,ndis,SYSTEM

      
