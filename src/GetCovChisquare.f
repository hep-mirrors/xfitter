*     ---------------------------------------------------------  	 
*     Calculate chisquare with the help of covariance matrix:
*     ---------------------------------------------------------

      subroutine GetCovChisquare(flag_in,n0_in,fchi2_in,pchi2_in)

      implicit none

      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'for_debug.inc'
      include 'systematics.inc'
      include 'theo.inc'
      include 'covar.inc'

      integer n0_in, flag_in
      double precision pchi2_in(nset)
      double precision fchi2_in, chisq
      double precision stat1, unc1, const1
      double precision stat2, unc2, const2

      integer i,j,IFAIL, h1iset
      double precision cov_total(n0_in, n0_in), tmp(n0_in),d,t

*     ----------------------------------------------------------
*     Initialise
*     ----------------------------------------------------------

      chisq=0.d0
      fchi2_in=0.d0

      do i=1,nset
         pchi2_in(i)=0.d0
      enddo

      do i=1, n0_in
         do j=1, n0_in
            cov_total(i, j) = 0.d0
         enddo
      enddo

      do i=1,n0_in
         call GetPointScaledErrors(i,stat1,unc1,const1)
         do j=1,n0_in
            call GetPointScaledErrors(j,stat2,unc2,const2)
            
C     fill with statistical part            
            cov_total(i,j) = corr_stat(i,j) * stat1 * stat2
            
C     add correlated systematic part from fully correlated part
            cov_total(i,j) = cov_total(i,j) + (corr_syst(i,j) * THEO(i) * THEO(j))

C     add correlated systematic part from given covariance matrix
            cov_total(i,j) = cov_total(i,j) + cov(i,j)

C     add uncorrelated systematic part
            if(i.eq.j) then
               cov_total(i,j) = cov_total(i,j) + (unc1 * unc2)
            endif

         enddo
      enddo

      CALL DINV  (n0_in,cov_total,n0_in,tmp,IFAIL)

      if(ifail.eq.-1) then
         Call HF_ERRLOG(11040001,'S: Matrix inversion failed !')
      endif

      do i=1,n0_in
         do j=1,n0_in
            
            chisq = (DATEN(i)-THEO(i)) *
     $           cov_total(i,j) * 
     $           (DATEN(j)-THEO(j))
            
            fchi2_in = fchi2_in + chisq

            h1iset = JSET(i)
            pchi2_in(h1iset) = pchi2_in(h1iset)+chisq
         enddo
      enddo

      return
      end
