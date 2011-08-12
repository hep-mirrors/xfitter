      subroutine minuit_ini

      character*72 minfile

      open ( 25, file='output/minuit.out.txt' )
      minfile='minuit.in.txt' 
      endif
      write(6,*) ' read minuit input params from file ',minfile
      open ( 24, file=minfile )
      open (  7, file='minuit.save.txt' )


      call mintio(24,25,7)


      return
      end
