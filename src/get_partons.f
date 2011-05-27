      subroutine Get_Partons(x,q2,upv,dnv,uusea,ddsea,str,chm,bot,glu)
      include 'pdfparam.inc'
c      implicit double precision (a-h, o-z)
      
      double precision x,q2,upv,dnv,uusea,ddsea,str,chm,bot,glu,flav_number
      double precision pdf
      dimension pdf(-6:6)
      xflav = flav_number(q2)

	if (x.eq.1) goto 999

      call fpdfxq(1,x,q2,pdf ,0)

c      usea = 0.5* (uplus - upv + sing/xflav)
c      dsea = 0.5* (dplus - dnv + sing/xflav)
c      str = 0.5* (splus + sing/xflav)
c      chm = 0.5* (cplus + sing/xflav)
c      bot = 0.5* (bplus + sing/xflav)
      upv=pdf(2)-pdf(-2)
      dnv=pdf(1)-pdf(-1)
      glu=pdf(0)
      uusea=pdf(-2)!+pdf(-4)
      ddsea=pdf(-1)!+pdf(-3)
      str=pdf(-3)
      chm=pdf(-4)
      bot=pdf(-5)
 

     
      if(xflav.lt.5) bot=0.0d0
      if(xflav.lt.4) chm=0.0d0

c ===== PATCH: FORCE C AND B  .GT. ZERO
      if(bot.lt.0.0d0)  bot=0.0d0
      if(chm.lt.0.0d0)  chm=0.0d0

 999  continue
      return
      end
