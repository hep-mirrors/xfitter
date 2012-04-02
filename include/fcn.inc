C
C  Store inforamtion from minuit FCN 
C
      integer IFlagFCN          !> FCN minuit flag, set by minuit (1-init, 2-iteration, 3-finalisation)
      integer IfcnCount         !>  Count number of iterations
      integer nparFCN
      integer ndfMINI      

      Common/CFCN/IFlagFCN,nparFCN,ndfMINI,IfcnCount
