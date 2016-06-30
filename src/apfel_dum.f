c----------------------------------------------------------
c     this file contains dummy APFEL interface routines
c     which are called in case xFitter was not compiled
c     with --enable-apfel option
c----------------------------------------------------------
      subroutine print_apfel_error_message
      print *, '--------------------------------------------------'
      print *, 'APFEL: You have chosen to use APFEL but xFitter'
      print *, 'is not compiled with --enable-apfel option.'
      call exit(-10)
      return 
      end

      subroutine apfel_ini
      call print_apfel_error_message
      return 
      end

      subroutine alphaqcd
      call print_apfel_error_message
      end

      subroutine setpdfset
      call print_apfel_error_message
      end

      subroutine xpdfall
      call print_apfel_error_message
      end

      subroutine xpdfallphoton
      call print_apfel_error_message
      end

      subroutine setalphaqcdref
      call print_apfel_error_message
      end

      subroutine evolveapfel
      call print_apfel_error_message
      end

      subroutine initializeapfel_dis
      call print_apfel_error_message
      end

      subroutine initializeapfel
      call print_apfel_error_message
      end

      subroutine setmassscalereference
      call print_apfel_error_message
      end

      subroutine setpolarizationdis
      call print_apfel_error_message
      end

      subroutine setprojectiledis
      call print_apfel_error_message
      end

      subroutine setprocessdis
      call print_apfel_error_message
      end

      subroutine setpolemasses
      call print_apfel_error_message
      end

      subroutine setperturbativeorder
      call print_apfel_error_message
      end

      subroutine setnumberofgrids
      call print_apfel_error_message
      end

      subroutine setgridparameters
      call print_apfel_error_message
      end

      function f2total(x)
      double precision x,f2total
      call print_apfel_error_message
      end

      function fltotal(x)
      double precision x,fltotal
      call print_apfel_error_message
      end

      function f3total(x)
      double precision x,f3total
      call print_apfel_error_message
      end

      function f2charm(x)
      double precision x,f2charm
      call print_apfel_error_message
      end

      function flcharm(x)
      double precision x,flcharm
      call print_apfel_error_message
      end

      function f2bottom(x)
      double precision x,f2bottom
      call print_apfel_error_message
      end

      function flbottom(x)
      double precision x,flbottom
      call print_apfel_error_message
      end

      subroutine computestructurefunctionsapfel
      call print_apfel_error_message
      end

      subroutine setzmass
      call print_apfel_error_message
      end

      subroutine setwmass
      call print_apfel_error_message
      end

      subroutine setsin2thetaw
      call print_apfel_error_message
      end

      subroutine setgfermi
      call print_apfel_error_message
      end

      subroutine setckm
      call print_apfel_error_message
      end

      subroutine setmassscheme
      call print_apfel_error_message
      end

      subroutine setmsbarmasses
      call print_apfel_error_message
      end

      subroutine enablemassrunning
      call print_apfel_error_message
      end

      subroutine setpropagatorcorrection
      call print_apfel_error_message
      end

      subroutine setewcouplings
      call print_apfel_error_message
      end

      subroutine enabledynamicalscalevariations
      call print_apfel_error_message
      end

      subroutine setrenqratio
      call print_apfel_error_message
      end

      subroutine setfacqratio
      call print_apfel_error_message
      end

      subroutine setmaxflavourpdfs
      call print_apfel_error_message
      end

      subroutine setmaxflavouralpha
      call print_apfel_error_message
      end

      function HeavyQuarkMass(i,Q)
      integer i
      double precision Q,HeavyQuarkMass
      call print_apfel_error_message
      end
