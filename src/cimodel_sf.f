      subroutine cimodel_sf(p_in, etaull_in,etadll_in,etaulr_in,etadlr_in,
     $     etaurl_in,etadrl_in,etaurr_in,etadrr_in)



      implicit none

      include 'steering.inc'
      include 'pdfparam.inc'


      
      double precision p_in(*)
      double precision etaull_in,etadll_in,etaulr_in,etadlr_in
      double precision etaurl_in,etadrl_in,etaurr_in,etadrr_in


cv leptoquarks
      if (modelHP.eq.1) then
                        
         etauLL_in = p_in(74)/2
         etadLL_in= 0.d0
         
         etauLR_in=0.d0
         etadLR_in=0.d0
                       
         etauRL_in=0.d0
         etadRL_in=0.d0
                       
         etauRR_in=0.d0
         etadRR_in=0.d0
                       
      elseif (modelHP.eq.2) then
                       
         etauLL_in = 0.d0
         etadLL_in= 0.d0
                       
         etauLR_in=0.d0
         etadLR_in=0.d0
                       
         etauRL_in=0.d0
         etadRL_in=0.d0
         
         etauRR_in=p_in(74)/2
         etadRR_in=0.d0
      elseif (modelHP.eq.3) then
         
         etauLL_in = 0.d0
         etadLL_in= 0.d0
         
         etauLR_in=0.d0
         etadLR_in=0.d0
         
         etauRL_in=0.d0
         etadRL_in=0.d0
         
         etauRR_in=0.d0
         etadRR_in=p_in(74)/2
         
      elseif (modelHP.eq.4) then
         
         etauLL_in = 0.d0
         etadLL_in= 0.d0
         
         etauLR_in=-p_in(74)/2
         etadLR_in=0.d0
         
         etauRL_in=0.d0
         etadRL_in=0.d0
         
         etauRR_in=0.d0
         etadRR_in=0.d0
         
      elseif (modelHP.eq.5) then
         
         etauLL_in = 0.d0
         etadLL_in= 0.d0
         
         etauLR_in=0.d0
         etadLR_in=0.d0
         
         etauRL_in=-p_in(74)/2
         etadRL_in= etauRL_in
         
         etauRR_in=0.d0
         etadRR_in=0.d0
      elseif (modelHP.eq.6) then
         
         etauLL_in = 0.d0
         etadLL_in= 0.d0
         
         etauLR_in=0.d0
         etadLR_in=-p_in(74)/2
         
         etauRL_in=0.d0
         etadRL_in=0.d0
         
         etauRR_in=0.d0
         etadRR_in=0.d0
         
      elseif (modelHP.eq.7) then
         
         etauLL_in = p_in(74)/2
         etadLL_in=  2*etauLL_in 
         
         etauLR_in=0.d0
         etadLR_in=0.d0
         
         etauRL_in=0.d0
         etadRL_in=0.d0
         
         etauRR_in=0.d0
         etadRR_in=0.d0
         
      elseif (modelHP.eq.8) then
         
         etauLL_in = 0.d0
         etadLL_in= -p_in(74)
         
         etauLR_in=0.d0
         etadLR_in=0.d0
         
         etauRL_in=0.d0
         etadRL_in=0.d0
         
         etauRR_in=0.d0
         etadRR_in=0.d0
      elseif (modelHP.eq.9) then
         
         etauLL_in = 0.d0
         etadLL_in= 0.d0
         
         etauLR_in=0.d0
         etadLR_in=0.d0
         
         etauRL_in=0.d0
         etadRL_in=0.d0
         
         etauRR_in=0.d0
         etadRR_in=-p_in(74)
         
      elseif (modelHP.eq.10) then
         
         etauLL_in = 0.d0
         etadLL_in= 0.d0

         etauLR_in=0.d0
         etadLR_in=0.d0

         etauRL_in=0.d0
         etadRL_in=0.d0
         
         etauRR_in=-p_in(74)
         etadRR_in=0.d0
         
      elseif (modelHP.eq.11) then
         
         etauLL_in = 0.d0
         etadLL_in= 0.d0
         
         etauLR_in=0.d0
         etadLR_in=p_in(74)
         
         etauRL_in=0.d0
         etadRL_in=0.d0
         
         etauRR_in=0.d0
         etadRR_in=0.d0
         
      elseif (modelHP.eq.12) then
         
         etauLL_in = 0.d0
         etadLL_in= 0.d0
         
         etauLR_in=0.d0
         etadLR_in=0.d0
         
         etauRL_in=p_in(74)
         etadRL_in= etauRL_in
         
         etauRR_in=0.d0
         etadRR_in=0.d0
      elseif (modelHP.eq.13) then
         
         etauLL_in = 0.d0
         etadLL_in= 0.d0
         
         etauLR_in=p_in(74)
         etadLR_in=0.d0
         
         etauRL_in=0.d0
         etadRL_in= 0.d0
         
         etauRR_in=0.d0
         etadRR_in=0.d0
      elseif (modelHP.eq.14) then
         
         etauLL_in = -2*p_in(74)
         etadLL_in= -p_in(74)
         
         etauLR_in=0.d0
         etadLR_in=0.d0
         
         etauRL_in=0.d0
         etadRL_in= 0.d0
         
         etauRR_in=0.d0
         etadRR_in=0.d0
         
         
         
         
cv compositness
      elseif (modelHP.eq.15) then
         
cv                       print*,'cucu'
         etauLL_in =p_in(74)
         etadLL_in=etauLL_in 
ch                       print*,'cucu', etauLL_in, p_in(74)
         
         etauLR_in=0.d0
         etadLR_in=0.d0
         
         etauRL_in=0.d0
         etadRL_in=0.d0
         
         etauRR_in=0.d0
         etadRR_in=0.d0
         
         
         
      elseif (modelHP.eq.16) then
         
         etauLL_in = 0.d0
         etadLL_in= 0.d0 
         
         etauLR_in=p_in(74)
         etadLR_in=etauLR_in
         
         etauRL_in=0.d0
         etadRL_in=0.d0

         etauRR_in=0.d0
         etadRR_in=0.d0
                       


      elseif (modelHP.eq.17) then
                      
         etauLL_in = 0.d0
         etadLL_in= 0.d0 
                      
         etauLR_in=0.d0
         etadLR_in=0.d0

         etauRL_in=p_in(74)
         etadRL_in=etauRL_in

         etauRR_in=0.d0
         etadRR_in=0.d0
                       


      elseif (modelHP.eq.18) then

         etauLL_in = 0.d0
         etadLL_in= 0.d0 
                      
         etauLR_in=0.d0
         etadLR_in=0.d0


         etauRL_in=0.d0
         etadRL_in=0.d0

                        
         etauRR_in=p_in(74)
         etadRR_in=etauRR_in
                       



      elseif (modelHP.eq.19) then
                        
         etauLL_in = p_in(74)
         etadLL_in= etauLL_in
                      
         etauLR_in=etauLL_in
         etadLR_in=etauLL_in
                      

         etauRL_in=etauLL_in
         etadRL_in=etauLL_in

                        
         etauRR_in=etauLL_in
         etadRR_in=etauLL_in
                       

      elseif (modelHP.eq.20) then
                      
         etauLL_in = p_in(74)
         etadLL_in= etauLL_in
                      
         etauLR_in=-etauLL_in
         etadLR_in=-etauLL_in
                      

         etauRL_in=-etauLL_in
         etadRL_in=-etauLL_in

         
         etauRR_in=etauLL_in
         etadRR_in=etauLL_in
                       
      elseif (modelHP.eq.21) then
                      
         etauLL_in = p_in(74)
         etadLL_in= etauLL_in
                      
         etauLR_in=-etauLL_in
         etadLR_in=-etauLL_in


         etauRL_in= etauLL_in
         etadRL_in= etauLL_in
                      
                        
         etauRR_in=-etauLL_in
         etadRR_in=-etauLL_in
                       
      elseif (modelHP.eq.22) then
                      
                      etauLL_in = p_in(74)
                      etadLL_in= etauLL_in
                      
                      etauLR_in=-etauLL_in
                      etadLR_in=-etauLL_in


                      etauRL_in= 0.d0
                      etadRL_in= 0.d0

                        
                      etauRR_in=0.d0
                      etadRR_in=0.d0

                  elseif (modelHP.eq.23) then
                      
                      etauLL_in = p_in(74)
                      etadLL_in= etauLL_in
                      
                      etauLR_in=0.d0
                      etadLR_in=0.d0


                      etauRL_in= etauLL_in
                      etadRL_in= etauLL_in

                        
                      etauRR_in=0.d0
                      etadRR_in=0.d0

                 elseif (modelHP.eq.24) then
                      
                      etauLL_in = p_in(74)
                      etadLL_in= etauLL_in
                      
                      etauLR_in=0.d0
                      etadLR_in=0.d0


                      etauRL_in= 0.d0
                      etadRL_in= 0.d0

                        
                      etauRR_in=etauLL_in
                      etadRR_in=etauLL_in

                elseif (modelHP.eq.25) then
                      
                      etauLL_in = 0.d0
                      etadLL_in=  0.d0
                      
                      etauLR_in=p_in(74)
                      etadLR_in=etauLR_in


                      etauRL_in= etauLR_in
                      etadRL_in= etauLR_in

                        
                      etauRR_in=0.d0
                      etadRR_in=0.d0

               elseif (modelHP.eq.26) then
                      
                      etauLL_in = 0.d0
                      etadLL_in=  0.d0
                      
                      etauLR_in=p_in(74)
                      etadLR_in=etauLR_in


                      etauRL_in= 0.d0
                      etadRL_in= 0.d0

                        
                      etauRR_in=etauLR_in
                      etadRR_in=etauLR_in

                   elseif (modelHP.eq.27) then
                      
                      etauLL_in = 0.d0
                      etadLL_in=  0.d0
                      
                      etauLR_in= 0.d0
                      etadLR_in= 0.d0


                      etauRL_in= p_in(74)
                      etadRL_in= etauRL_in

                        
                      etauRR_in=-etauRL_in
                      etadRR_in=-etauRL_in

                   elseif (modelHP.eq.28) then
                      
                      etauLL_in = p_in(74)
                      etadLL_in=  0.d0
                      
                      etauLR_in= -etauLL_in 
                      etadLR_in= 0.d0


                      etauRL_in= 0.d0
                      etadRL_in= 0.d0

                        
                      etauRR_in=0.d0
                      etadRR_in=0.d0

                   elseif (modelHP.eq.29) then
                      
                      etauLL_in = p_in(74)
                      etadLL_in=  0.d0
                      
                      etauLR_in= 0.d0 
                      etadLR_in= 0.d0


                      etauRL_in= etauLL_in 
                      etadRL_in= 0.d0

                        
                      etauRR_in=0.d0
                      etadRR_in=0.d0
                   elseif (modelHP.eq.30) then
                      
                      etauLL_in = p_in(74)
                      etadLL_in=  0.d0
                      
                      etauLR_in= 0.d0 
                      etadLR_in= 0.d0


                      etauRL_in= 0.d0
                      etadRL_in= 0.d0

                        
                      etauRR_in=etauLL_in
                      etadRR_in=0.d0
                   elseif (modelHP.eq.31) then
                      
                      etauLL_in = 0.d0
                      etadLL_in=  0.d0
                      
                      etauLR_in= p_in(74)
                      etadLR_in= 0.d0


                      etauRL_in= etauLR_in
                      etadRL_in= 0.d0

                        
                      etauRR_in=0.d0
                      etadRR_in=0.d0
                   elseif (modelHP.eq.32) then
                      
                      etauLL_in = 0.d0
                      etadLL_in=  0.d0
                      
                      etauLR_in= p_in(74)
                      etadLR_in= 0.d0


                      etauRL_in= 0.d0
                      etadRL_in= 0.d0

                        
                      etauRR_in=etauLR_in
                      etadRR_in=0.d0
                   elseif (modelHP.eq.33) then
                      
                      etauLL_in = 0.d0
                      etadLL_in=  0.d0
                      
                      etauLR_in= 0.d0
                      etadLR_in= 0.d0


                      etauRL_in= p_in(74)
                      etadRL_in= 0.d0

                        
                      etauRR_in=-etauRL_in
                      etadRR_in=0.d0

                   endif


        return
      end

