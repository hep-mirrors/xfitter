*     ---------------------------------------------------------  	 
*     Calculate chisquare:
*      - first get error matrix
*      - invert error matric to get errors
*      - calculate chisquare
*     ---------------------------------------------------------

      subroutine GetNewChisquare(flag_in,n0_in,fchi2_in,rsys_in,ersys_in,pchi2_in,fcorchi2_in)
      implicit none

      include 'ntot.inc'
      include 'steering.inc'
      include 'systematics.inc'
      
      integer n0_in, flag_in
      double precision fchi2_in, ERSYS_in(NSYSMax), RSYS_in(NSYSMax)
      double precision pchi2_in(nset), fcorchi2_in

      double precision chi2tmp
      integer i,jsys

c initialisation 
      do jsys=1,nsys
         rsys_in(jsys) = 0.d0
         ersys_in(jsys) = 0.d0
      enddo

      do i=1,nset
         pchi2_in(i)=0.d0
      enddo
      
      fchi2_in = 0.d0

      if(Chi2CorErr.eq.'Nuisance') then
         call calc_nuisance(rsys_in, ersys_in, fcorchi2_in)
c         fchi2_in = fcorchi2_in ! this is missing in the original code? why?
      endif 
      
      if (Chi2UncorErr.eq.'Simple') then
         call calc_simple_chi2(rsys_in, chi2tmp, pchi2_in)
         fchi2_in = chi2tmp + fchi2_in
      endif
      
      if((Chi2CorErr.eq.'Matrix').or.(Chi2UncorErr.eq.'Matrix')) then
         call GetCovChisquare(flag_in,n0_in, chi2tmp, pchi2_in)
         fchi2_in = chi2tmp + fchi2_in
      endif
            
      if(ICHI2.eq.41) then   ! Poissonian tail chi2 contribution
         call calc_poisson(n0_in, chi2tmp)
         fchi2_in = chi2tmp + fchi2_in
      endif

      return 
      end

      subroutine calc_nuisance(rsys_in, ersys_in, fcorchi2_in)
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'systematics.inc'
      
      double precision ERSYS_in(NSYSMax), RSYS_in(NSYSMax)
      double precision fcorchi2_in

      print *, 'not yet implemented'
      return 
      end
      
      subroutine calc_simple_chi2(rsys_in, chi2tmp, pchi2_in)
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'systematics.inc'
      
      double precision RSYS_in(NSYSMax)
      double precision chi2tmp, pchi2_in(nset)

      print *, 'not yet implemented'
      return 
      end

      subroutine calc_poisson(n0_in, chi2)
      implicit none
      include 'ntot.inc'
      include 'systematics.inc'
     
      double precision chi2
      integer n0_in, ipoint
      double precision stat, unc, const, fac, error

      chi2=0.d0
      do ipoint=1,n0_in
         call GetPointScaledErrors(ipoint,fac,stat,unc,const)
         error = dsqrt(stat**2+unc**2+const**2)
         chi2 = chi2 + 2.*log( error/alpha(ipoint)) 
      enddo

      return 
      end



      subroutine GetChisquare(flag_in,n0_in,fchi2_in,rsys_in,ersys_in,pchi2_in,fcorchi2_in)
      implicit none

      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'for_debug.inc'
      include 'systematics.inc'
      INCLUDE 'theo.inc'


      integer ir(nsysmax),n0_in 
      integer isys,ipoint,jpoint,ifail,flag_in
      double precision chisq,fchi2_in,chi2error,fchi2_error
      double precision d,t,error,errorunc,errorconst
      double precision errorsta, fac, facd, fcorchi2_in
      double precision d_i, t_i, error_i, d_j, error_j, t_j         
      integer i,j,jsys,h1iset
      double precision pchi2_in(nset)
      double precision EBSYS_in(NSYSMax),ERSYS_in(NSYSMax)
      double precision BSYS_in(NSYSMax),RSYS_in(NSYSMax)
      double precision sub

      double precision dNEvt, tNEvt

      double precision factor2, factor_1, factor_2

      integer npoisson, ngauss

C      integer getcachesize

*     ----------------------------------------------------------
*     Initialise
*     ----------------------------------------------------------
      npoisson = 0
      ngauss   = 0

      chisq=0.d0

      fchi2_error = 0.d0

      fchi2_in = 0.d0

 !> Also zero fit/control sample chi2s
      chi2_fit = 0.
      chi2_cont = 0.


      sub = 0.d0


      do jsys=1,nsys
         ir(nsys)=0.d0
         bsys_in(jsys) = 0.d0
         rsys_in(jsys) = 0.d0
         ebsys_in(jsys) = 0.d0
         ersys_in(jsys) = 0.d0
      enddo

      do i=1,nset
         pchi2_in(i)=0.d0
      enddo
*     -- for ICHI2=2, the matrix sysa is filled already in Systematics

      if (ICHI2.ne.2) then
         do isys=1,nsys
            do jsys=1,nsys
               sysa(isys,jsys) = 0.d0
               if (jsys.eq.isys) sysa(isys,jsys) = 1.d0
            enddo
         enddo
      endif


*     ----------------------------------------------------------
*     Pascaud-like chi2 plus error scaled modifications
*     ----------------------------------------------------------
      if (mod(ICHI2, 10).eq.1) then



*     ---------------------------------------------------------
*     now calculate the bsys(nsys) and the sysa matrix
*     ---------------------------------------------------------

            
         do ipoint=1,n0_in
               
            d = DATEN(ipoint)
            t = THEO(ipoint)
            error = ALPHA(ipoint)
                  

***   scale errors for chi2 calculation
***   in principle:  unc*(t/d),  sta*dsqrt(t/d)
            if (ICHI2.eq.11 .or. ICHI2.eq.41) then
               if ( (alpha(ipoint)/d.gt.Chi2MaxError)
     $              .or. (ControlFitSplit.and. .not. FitSample(ipoint))
     $              ) then
C     Turn off the point for the syst. errors shift estimation:
                  error = 1.D010
               else
***   mixed scaling - decompose - scale - recombine
                  errorunc = E_UNC(ipoint)*d/100.
                  errorconst = E_STA_CONST(ipoint)*d/100.

                  if (errorunc.gt.error) then
                     errorsta = 0.
                  else
                     errorsta = error**2-errorunc**2-errorconst**2
                     if (errorsta.gt.0) then
                        errorsta = sqrt(errorsta)
                     else
                        errorsta = 0.
                     endif
                  endif
                  if (t.gt.0) then
                     errorsta = errorsta*dsqrt(abs(t/d))
                     errorunc = errorunc*(abs(t/d))
                  endif
                  error = dsqrt(errorsta**2+errorunc**2+errorconst**2)
               endif

            else if (ICHI2.eq.21) then
***   linear scaling
               error = error*(abs(t/d))
            else if (ICHI2.eq.31) then
***   sqrt scaling
               error = error*dsqrt(abs(t/d))
            endif

            

            factor2 = 1.D0/error**2

            do isys = 1,nsys                                 

               if (beta(isys,ipoint).ne.0) then
                  if (SysAdditive(isys)) then
                     factor_1 = -d
                  else
                     factor_1 = t
                  endif

                  bsys_in(isys) = bsys_in(isys) 
     +                 + factor_1*(d-t)*BETA(isys,ipoint)*factor2
               
c                  ebsys_in(isys) = ebsys_in(isys)     !> ????
c     +                 + t * BETA(isys,ipoint)/error
               
                  do  jsys=isys,nsys
                     if (beta(jsys,ipoint).ne.0) then
                        if (SysAdditive(jsys)) then
                           factor_2 = -d
                        else
                           factor_2 = t
                        endif

                        sysa(isys,jsys) = sysa(isys,jsys)
     +                       + beta(isys,ipoint)*beta(jsys,ipoint)
     $                       *factor_1*factor_2*factor2

                     endif
                  enddo
               endif
            enddo

         enddo


         do isys=1,nsys
C            print '(5F10.2)',(sysa(isys,jsys),jsys=1,5)
            do jsys=isys+1,nsys
               sysa(jsys,isys) = sysa(isys,jsys)
            enddo
         enddo
*     ---------------------------------------------------------
*     inverse sysa and find the shifts
*     ---------------------------------------------------------


         
         if (nsys.gt.0) then
            if (flag_in.eq.3) then
               Call DEQINV(NSys,sysa,NSYSMAX,IR,IFAIL,1,bsys_in)
            else
c               CALL DPOSV('Upper',Nsys,1,sysa,NSYSMax,Bsys_in,NSYSmax
c     $              ,Ifail)
               Call DEQN(NSys,sysa,NSYSMAX,IR,IFAIL,1,bsys_in)
            endif
         endif

         do isys=1,nsys
            rsys_in(isys) = -bsys_in(isys)
         enddo

         if (DEBUG.and.flag_in.eq.1) then
            write(78,*)
            do isys=1,nsys
               write(78,*) 'isys rsys ',isys,rsys_in(isys)
            enddo
            write(78,*)
         endif




*     ---------------------------------------------------------
*     now calculate the chi2
*     ---------------------------------------------------------

         do ipoint=1,n0_in

            h1iset = JSET(ipoint)
            d = daten(ipoint)
            t = theo(ipoint)

*** Factor to scale theory/data:
            fac  = 1.d0
            facd = 1.d0
            do isys=1,nsys
               if (.not. SysAdditive(isys) ) then
                  fac  = fac  - rsys_in(isys)*beta(isys,ipoint)                
               else
                  facd = facd + rsys_in(isys)*beta(isys,ipoint)                
               endif
            enddo 


            error = alpha(ipoint)

***   scale errors for chi2 calculation - as above!
            if (ICHI2.eq.11 .or. ICHI2.eq.41) then
***   mixed scaling - decompose - scale - recombine
               errorunc = E_UNC(ipoint)*d/100.
               errorconst = E_STA_CONST(ipoint)*d/100.
               if (errorunc.gt.error) then
                  errorsta = 0.
               else
                  errorsta = error**2-errorunc**2-errorconst**2
                  if (errorsta.gt.0) then
                     errorsta = sqrt(errorsta)
                  else
                     errorsta = 0.
                  endif
               endif
               if (t.gt.0) then
                  if (iDH_MOD.ne.0) then
                     errorsta = errorsta*dsqrt(abs(t*fac/d))
                  else
                     errorsta = errorsta*dsqrt(abs(t/d))
                  endif
                  errorunc = errorunc*(abs(t/d))
               endif
               error = dsqrt(errorsta**2+errorunc**2+errorconst**2)
               !> Extra contribution due to 2xlog sigma term:
               chi2error =  2.*log( error/alpha(ipoint)) !> subtract un-modified error such that delta chi2=0 if errors are not modified.


            else if (ICHI2.eq.21) then
***   linear scaling
               error = error*(abs(t/d))
            else if (ICHI2.eq.31) then
***   sqrt scaling
               error = error*dsqrt(abs(t/d))
            endif

                                  
            if (DEBUG.and.flag_in.eq.1) then
               write(78,*) 'ipoint fac ',ipoint,fac
            endif


            t = t*fac
            d = d*facd

            THEO_MOD(ipoint)=t
            ALPHA_MOD(ipoint)=error

C
C Add Poisson log likelihood
C
            if (alpha(ipoint)/d.gt.Chi2MaxError) then
C use Poisson formula:
               dNEvt = (d*d)/(alpha(ipoint)*alpha(ipoint))
               tNEvt = (t*d)/(alpha(ipoint)*alpha(ipoint))
               chisq = 2* ( tNEvt - dNEvt - dNEvt*log(tNevt/dNEvt))
               
C               print *,dnevt,tnevt,chisq
               chi2error = 0.
               npoisson = npoisson+1
            else
               chisq = (d-t)**2/error**2
               ngauss = ngauss+1
            endif

            if (ControlFitSplit) then
               if (FitSample(ipoint)) then
                  chi2_fit  = chi2_fit + chisq
                  fchi2_in  = fchi2_in + chisq
               else
                  chi2_cont = chi2_cont + chisq
               endif
      
            else               
               fchi2_in = fchi2_in + chisq
            endif

            if ( ICHI2.eq.41) then
               fchi2_error = fchi2_error + chi2error
               fchi2_in = fchi2_in + chi2error
            endif

            if (flag_in.eq.3) then
               pchi2_in(h1iset) = pchi2_in(h1iset)+chisq
            endif

            if (flag_in.eq.1) then
               write(87,880) h1iset, 
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
            endif
            if (flag_in.eq.3) then
               write(88,880) h1iset,
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
 880           format(1x, i2, 2x, G12.6, 2x, G12.4, 2x, G12.6, 3(2x, G12.4))
            endif

*     -- errors on the shifts 

            if (flag_in.eq.3) then
               do isys=1,nsys
                  ersys_in(isys) = sqrt(sysa(isys,isys))
               enddo
            endif

         enddo

         do isys=1,nsys
            fchi2_in = fchi2_in + rsys_in(isys)**2
         enddo

         

c....print out the correlated chi2
         if (flag_in.eq.3) then
            fcorchi2_in=0.0d0
            do isys=1,nsys
               fcorchi2_in= fcorchi2_in +rsys_in(isys)**2
            enddo
         endif


c...........................


      elseif (ICHI2.eq.2) then  !  CTEQ-like chi2

         do ipoint=1,n0_in

            d_i = daten(ipoint)
            t_i = theo(ipoint)
            error_i = alpha(ipoint)

            chisq = (d_i-t_i)**2 / error_i**2
            fchi2_in = fchi2_in + chisq

            do jsys=1,nsys
               bsys_in(jsys) = bsys_in(jsys)
     +              + beta(jsys,ipoint)*(d_i-t_i)/error_i**2
            enddo

            if (flag_in.eq.3) then
               h1iset = JSET(ipoint)
               pchi2_in(h1iset) = pchi2_in(h1iset)+chisq
            endif


         enddo

*     
*     - Now calculate the term to subtract from chi2 (cf CTEQ)

         do i=1,nsys
            do j=1,nsys
               sub = sub + bsys_in(i) * sysa(i,j) * bsys_in(j)
            enddo
         enddo


         print *,'haha',fchi2_in,sub
         fchi2_In = fchi2_in - sub


*     -- and get the systematic shifts :

         do i=1,nsys
            RSYS_in(i) = 0.
            do j=1,nsys
               rsys_in(i) = rsys_in(i) + sysa(i,j) * bsys_in(j)
            enddo
         enddo

c....print out the correlated chi2
         if (flag_in.eq.3) then
            fcorchi2_in=0.0d0
            do isys=1,nsys
               fcorchi2_in= fcorchi2_in +rsys_in(isys)**2
            enddo
         endif
c...........................

      elseif (mod(ICHI2, 10).eq.3) then ! Offset
*     ---------------------------------------------------------
*     now calculate the chi2
*     ---------------------------------------------------------

         fcorchi2_in=0.0d0
         
         do ipoint=1,n0_in

            h1iset = JSET(ipoint)
            d = daten(ipoint)
            t = theo(ipoint)
            error = alpha(ipoint)

            if (DEBUG.and.flag_in.eq.1) then
               write(78,*) 'ipoint fac ',ipoint,fac
            endif


            THEO_MOD(ipoint)=t
            ALPHA_MOD(ipoint)=error

C
C Add Poisson log likelihood
C
            if (alpha(ipoint)/d.gt.Chi2MaxError) then
C use Poisson formula:
               dNEvt = (d*d)/(alpha(ipoint)*alpha(ipoint))
               tNEvt = (t*d)/(alpha(ipoint)*alpha(ipoint))
               chisq = 2* ( tNEvt - dNEvt - dNEvt*log(tNevt/dNEvt))
               
C               print *,dnevt,tnevt,chisq
               chi2error = 0.
               npoisson = npoisson+1
            else
               chisq = (d-t)**2/error**2
               ngauss = ngauss+1
            endif
            fchi2_in = fchi2_in + chisq
            if (flag_in.eq.3) then
              fcorchi2_in = fcorchi2_in + ((1-t/d)/E_TOT(ipoint)*1.d2)**2
            endif

            if (flag_in.eq.3) then
               pchi2_in(h1iset) = pchi2_in(h1iset)+chisq
            endif

            if (flag_in.eq.1) then
               write(87,880) h1iset, 
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
            endif
            if (flag_in.eq.3) then
               write(88,880) h1iset,
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
cws     880           format(1x, i2, 2x, G12.6, 2x, G12.4, 2x, G12.6, 3(2x, G12.4))
            endif

         enddo

      endif

      if (ichi2.eq.41) then
         print '(''Chi2 due to 2xlog sigma term'',F6.1)', fchi2_error
      endif

      if (npoisson.gt.0) then
         print '(''Use Poisson stats for '',I6,'', use Gauss stats for''
     $        ,I6,'' events'')',npoisson,ngauss
      endif

      return
      end
