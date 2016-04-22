
      BLOCK DATA DATEwk
C                       
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NSP = 2, NGN = 3, NBN = 4, NPOL = 3)
      PARAMETER (NFL = NSP * NGN, Ntt1=NSP*NBN*NSP, Ntt2=NFL*NBN*NFL)
      PARAMETER (NANG= NGN*(NGN-1)/2, NPHS=(NGN-1)*(NGN-2)/2)
 
      COMMON / EwkPAR / WMS, ZMS, SWG2, ALFE, ALFEW(NBN)
      COMMON / MIXPAR / CMX(NANG), DMX(NPHS)
      COMMON / EW1LCP / GLL(NSP,NBN,NSP), GLR(NSP,NBN,NSP), 
     >                  GLV(NSP,NBN,NSP), GLA(NSP,NBN,NSP)
      COMMON / EW1QCP / GQL(NSP,NBN,NSP), GQR(NSP,NBN,NSP), 
     >                  GQV(NSP,NBN,NSP), GQA(NSP,NBN,NSP)
      COMMON / EW2QCP / HQL(NFL,NBN,NFL), HQR(NFL,NBN,NFL), 
     >                  HQV(NFL,NBN,NFL), HQA(NFL,NBN,NFL)
 
ccpy 03/14/2012: According to arXiv:0911.0884v2,  SWG2_DIS = 0.2315
ccpy      DATA WMS, ZMS, SWG2, ALFE / 80.41, 91.187, 0.23124, 7.297353E-3 /
      DATA WMS, ZMS, SWG2, ALFE / 80.41, 91.187, 0.2315, 7.297353E-3 /
      Data (((GLL(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /
     >     (((GLR(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GLV(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GLA(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 / 
C
     >     (((GQL(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /
     >     (((GQR(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GQV(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GQA(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /   
C
     >     (((HQL(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
     >     (((HQR(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
     >     (((HQV(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
     >     (((HQA(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
 
C                        ****************************
      END

