      Subroutine qcdnumgetall (x, QQ, f)

*     return pdfs 


      Implicit double precision (a-h,o-z)

      double precision x, QQ
      double precision f(7)
      
      integer iset              ! iset --> unpolarised =1, polarised = 2 , fragmentation function = 3
      parameter (iset = 1)      ! custom = 4, external 5-9
     



      dimension pdf(-6:6)

      integer i


      call hf_get_pdfs(x,qq,pdf)
      
      do i=1,7         
         f(i) = pdf(i-1)
      enddo
c      write (*,1008)'x,qq,f', x,qq,(f(i),i=1,7)
 1008 FORMAT (A8,F7.4,1X,F6.1,1X,7(1PE11.4,1X))

      return 
      end

c ========================================================================

      double precision function qcdnumget (x, QQ, iparton)
      
cws      Implicit double precision (a-h,o-z)

      double  precision x, QQ
      integer iparton
cws      dimension pdf(-6:6)

      integer iset              ! iset --> unpolarised =1, polarised = 2 , fragmentation function = 3
      parameter (iset = 1)      ! custom = 4, external 5-9
     

      integer flag             
      parameter (flag = 1)      ! if set to zero  fpdfxq will return null when x and q2 are outside the 
                                ! grid boundaries

      qcdnumget = FVALXQ ( iset, iparton, x, QQ, flag )

      return
      end
