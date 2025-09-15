C
C  Store inforamtion from minuit FCN 
C
      integer IFlagFCN          !> FCN minuit flag, set by minuit (1-init, 2-iteration, 3-finalisation)
      integer IfcnCount         !>  Count number of iterations
      integer nparFCN
      integer ndfMINI      

      Common/CFCN/IFlagFCN,nparFCN,ndfMINI,IfcnCount

C
C  Errors-on-errors
C

      double precision epsilon_value  ! Declare epsilon_value
      integer n_iterations
      Common/CEpsilon/epsilon_value,n_iterations 