      Subroutine  GetIntegratedNCXsection(IDataSet)
C-----------------------------------------------------------------
C
C  Created by SG, 24/05/2011
C  Calculate x-integrated  cross section.
C
C-----------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'couplings.inc'
      include 'qcdnumhelper.inc'
      
      integer IDataSet

      integer idxQ2, idxX1, idxX2, i, j, idx

      integer NSplit
      parameter (NSplit=100)

      integer NQ2Max
      parameter (NQ2Max=60)

      integer Npt

      double precision charge,S,y,yplus,yminus
      double precision xmin,xmax,q2v,x(0:NSplit,NQ2Max)
     $     ,q2(0:NSplit,NQ2Max)
      double precision F2p(0:NSplit,NQ2Max),xF3p(0:NSplit,NQ2Max)
     $     ,FLp(0:NSplit,NQ2Max)
      double precision F2m(0:NSplit,NQ2Max),xF3m(0:NSplit,NQ2Max)
     $     ,FLm(0:NSplit,NQ2Max)
      double precision F2,xF3,FL
      double precision XSec(0:NSplit),Xsint


C
C Constants: should be somewhere else ...
C
      double precision alphaem,pi,convfac
      parameter (alphaem = 1./137.0D0)
      parameter (pi = 3.141592653589793D0)
      parameter (convfac = 3.893796D8)
      
      double precision alpha_run

      double precision eu,ed, e2u,e2d ,pz
      double precision ve,ae,au,ad,vu,vd,A_u,A_d,B_u,B_d
      parameter (eu =  2./3.)
      parameter (ed = -1./3.)
      parameter (e2u = eu*eu)
      parameter (e2d = ed*ed)      

C Functions:
      integer GetBinIndex
      integer GetInfoIndex
      double precision DSIMPS
      double precision AEMRUN
C-----------------------------------------------------------------


      ve = -0.5d0 + 2.*sin2thw
      ae = -0.5d0                  
      au = cau
      ad = cad
                  
      vu = au - (4.d0/3.d0)*sin2thw
      vd = ad + (2.d0/3.d0)*sin2thw
                  


                                  

      if (NDATAPOINTS(IDataSet).gt.NQ2Max) then
         print *,'ERROR IN GetIntegratedNCXsection'
         print *,'INCREASE NQ2Max to ',NDATAPOINTS(IDataSet)
         stop
      endif

C Get indexes for Q2, x1 and x2:
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX1 = GetBinIndex(IDataSet,'xmin')
      idxX2 = GetBinIndex(IDataSet,'xmax')

      if (idxX1.eq.0 .or. idxX2.eq.0 .or. idxQ2.eq.0) then
         Return
      endif


C prepare bins:
      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)

         xmin = AbstractBins(idxX1,idx)
         xmax = AbstractBins(idxX2,idx)
         q2v  = AbstractBins(idxQ2,idx)

         do j=0,NSplit
            q2(j,i) = q2v
            x(j,i)  = (xmax-xmin)/NSplit * j + xmin
         enddo

      enddo

      Npt = (NSplit+1)*NDATAPOINTS(IDataSet)


C QCDNUM, caclulate FL, F2 and xF3 for all bins:
      charge = DATASETInfo( GetInfoIndex(IDataSet,'e charge'), IDataSet)
      S = (DATASETInfo( GetInfoIndex(IDataSet,'sqrt(S)'), IDataSet))**2

      CALL ZMSTFUN(1,CNEP2F,X,Q2,FLp,Npt,0)
      CALL ZMSTFUN(2,CNEP2F,X,Q2,F2p,Npt,0)
      CALL ZMSTFUN(3,CNEP3F,X,Q2,XF3p,Npt,0)    
      CALL ZMSTFUN(1,CNEM2F,X,Q2,FLm,Npt,0)
      CALL ZMSTFUN(2,CNEM2F,X,Q2,F2m,Npt,0)
      CALL ZMSTFUN(3,CNEM3F,X,Q2,XF3m,Npt,0) 


      do i=1,NDATAPOINTS(IDataSet)

         PZ = 4.d0 * sin2thw * cos2thw * (1.+Mz**2/Q2(0,i))
         PZ = 1./Pz
         A_u = e2u - ve*PZ*2.*eu*vu +(ve**2 + ae**2)*PZ**2*(vu**2+au**2)
         A_d = e2d - ve*PZ*2.*ed*vd +(ve**2 + ae**2)*PZ**2*(vd**2+ad**2)
         B_u = -ae*PZ*2.*eu*au + 2.*ve*ae*(PZ**2)*2.*vu*au
         B_d = -ae*PZ*2.*ed*ad + 2.*ve*ae*(PZ**2)*2.*vd*ad

         alpha_run = AEMRUN(q2(0,i))

C Get x-sections:
         do j=0,NSplit
            y = q2(j,i)/(S*X(j,i))
            yplus  = 1+(1-y)**2
            yminus = 1-(1-y)**2


            XF3  = B_U*XF3p(j,i)  + B_D*XF3m(j,i)
            F2   = A_U*F2p(j,i)   + A_D*F2m(j,i)
            FL   = A_U*FLp(j,i)   + A_D*FLm(j,i)


            if (charge.gt.0) then
               XSec(j) = yplus*F2 - yminus*xF3 - y*y*FL
            else
               XSec(j) = yplus*F2 + yminus*xF3 - y*y*FL
            endif
            XSec(j) = XSec(j) * (2*pi*alpha_run**2)/ 
     $           (q2(j,i)**2*x(j,i))*convfac
C            print *,'hihi',q2(j,i),x(j,i),y,Xsec(j),F2
         enddo


C Integrate:
         xsint = DSIMPS(XSec(0),x(0,i),x(NSplit,i),NSplit)

C         print *,'hahaha',xsint,x(0,i),x(NSplit,i),S,XSec(5), A_U*F2p(0,i)   + A_D*F2m(0,i)
C     $        ,q2(0,i),x(0,i),y
         
         idx =  DATASETIDX(IDataSet,i)
         THEO(idx) =  xsint
      enddo



C-----------------------------------------------------------------
      end


      Subroutine GetReducedNCXsection(IDataSet)
C----------------------------------------------------------------
C
C  Created by SG, 25/05/2011
C  Start with zero mass implementation
C
C---------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'couplings.inc'
      include 'qcdnumhelper.inc'

      integer IDataSet

C
C Constants: should be somewhere else ...
C
      double precision alphaem,pi,convfac
      parameter (alphaem = 1./137.0D0)
      parameter (pi = 3.141592653589793D0)
      parameter (convfac = 3.893796D8)


      double precision eu,ed, e2u,e2d ,pz
      double precision ve,ae,au,ad,vu,vd,A_u,A_d,B_u,B_d
      parameter (eu =  2./3.)
      parameter (ed = -1./3.)
      parameter (e2u = eu*eu)
      parameter (e2d = ed*ed)      

      integer idxQ2, idxX, idxY, i,  idx
      
      integer NPmax
      parameter(NPmax=1000)

      double precision X(NPmax),Y(NPmax),Q2(NPmax)
      double precision yplus, yminus

      double precision F2p(NPmax),xF3p(NPmax),FLp(NPmax)
      double precision F2m(NPmax),xF3m(NPmax),FLm(NPmax)
      
      double precision F2,xF3,FL,XSec

      double precision Xv,Yv,Q2v

      double precision Charge, S

C Functions:
      integer GetBinIndex
      integer GetInfoIndex

C---------------------------------------------------------

      ve = -0.5d0 + 2.*sin2thw
      ae = -0.5d0                  
      au = cau
      ad = cad
                  
      vu = au - (4.d0/3.d0)*sin2thw
      vd = ad + (2.d0/3.d0)*sin2thw
 

      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN GetReducedNCXsection'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)
         stop
      endif

C Get indexes for Q2, x and y:
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX  = GetBinIndex(IDataSet,'x')
      idxY = GetBinIndex(IDataSet,'y')

      if (idxQ2.eq.0 .or. idxX.eq.0 .or. idxY.eq.0) then
         Return
      endif

C prepare bins:
      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)

         X(i)   = AbstractBins(idxX,idx)
         Y(i)   = AbstractBins(idxY,idx)
         Q2(i)  = AbstractBins(idxQ2,idx)
      enddo
     
C QCDNUM, caclulate FL, F2 and xF3 for all bins:
      charge = DATASETInfo( GetInfoIndex(IDataSet,'e charge'), IDataSet)
      S = (DATASETInfo( GetInfoIndex(IDataSet,'sqrt(S)'), IDataSet))**2

      CALL ZMSTFUN(1,CNEP2F,X,Q2,FLp,NDATAPOINTS(IDataSet),0)
      CALL ZMSTFUN(2,CNEP2F,X,Q2,F2p,NDATAPOINTS(IDataSet),0)
      CALL ZMSTFUN(3,CNEP3F,X,Q2,XF3p,NDATAPOINTS(IDataSet),0)    
      CALL ZMSTFUN(1,CNEM2F,X,Q2,FLm,NDATAPOINTS(IDataSet),0)
      CALL ZMSTFUN(2,CNEM2F,X,Q2,F2m,NDATAPOINTS(IDataSet),0)
      CALL ZMSTFUN(3,CNEM3F,X,Q2,XF3m,NDATAPOINTS(IDataSet),0) 

      do i=1,NDATAPOINTS(IDataSet)

         PZ = 4.d0 * sin2thw * cos2thw * (1.+Mz**2/Q2(i))
         PZ = 1./Pz
         A_u = e2u - ve*PZ*2.*eu*vu +(ve**2 + ae**2)*PZ**2*(vu**2+au**2)
         A_d = e2d - ve*PZ*2.*ed*vd +(ve**2 + ae**2)*PZ**2*(vd**2+ad**2)
         B_u = -ae*PZ*2.*eu*au + 2.*ve*ae*(PZ**2)*2.*vu*au
         B_d = -ae*PZ*2.*ed*ad + 2.*ve*ae*(PZ**2)*2.*vd*ad

C Get x-sections:
         yplus  = 1+(1-y(i))**2
         yminus = 1-(1-y(i))**2


         XF3  = B_U*XF3p(i)  + B_D*XF3m(i)
         F2   = A_U*F2p(i)   + A_D*F2m(i)
         FL   = A_U*FLp(i)   + A_D*FLm(i)
         

         if (charge.gt.0) then
            XSec = F2 - yminus/yplus*xF3 - y(i)*y(i)/yplus*FL
         else
            XSec = F2 + yminus/yplus*xF3 - y(i)*y(i)/yplus*FL
         endif

         idx =  DATASETIDX(IDataSet,i)
         THEO(idx) =  XSec

      enddo


      end
