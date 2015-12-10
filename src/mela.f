************************************************************************
*
*     Initialization routine for MELA. This subrotine writes the input
*     file to be fed to MELA for the evolution and the computation of
*     the sturecture functions.
*
*     Author: Valerio Bertone
*     Created: 08/09/2015
*
************************************************************************
      subroutine MELA_init()
*
      implicit none
*
#include "ntot.inc"
#include "steering.inc"
#include "couplings.inc"
*
      character*13 card
*
      double precision Q_ref,Alphas_ref
      double precision hf_get_alphas
*
      Q_ref      = mz
      Alphas_ref = hf_get_alphas(Q_ref*Q_ref)
*
      card = "MELA_init.dat"
      open(unit=7,file=card,status="unknown")
      write(7,'(a)') "#"
      write(7,'(a)') "# Configuration file. Comments start with #."
      write(7,'(a)') "#"
      write(7,'(a)') "# If any of the following parameters is missing or commented out,"
      write(7,'(a)') "# the default value will be used."
      write(7,'(a)') "#"
      write(7,'(a)') "#########################################"
      write(7,'(a)') "# Evolution parameters                  #"
      write(7,'(a)') "#########################################"
      write(7,'(a,i1,a)') "IPT      = ",I_FIT_ORDER-1,"                            # Perturbative order (0,1,2)"
      write(7,'(a)') "MODEV    = 'PTH'                        # Solution of the DGLAP equation ('ITE','TRN','PTH')"
      write(7,'(a)') "NS       = 'VFNS'                       # Mass scheme used for the LHA evoultion ('VFNS','FFNS')"
      write(7,'(a)') "NFMAX    = 6                            # Maximum number of flavours for the LHA evolution (3,4,5,6)"
      write(7,'(a)') "NFFN     = 3                            # Number of active flavours in the FFNS (3,4,5,6)"
      write(7,'(a)') "HQMASS   = 'POLE'                       # Heavy quark mass scheme ('POLE','MSBAR')"
      write(7,'(a)') "KRF      = 1                            # muR / muF"
      write(7,'(a)') "EVOL     = 'SPACE'                      # Evolution ('SPACE' = space-like (PDFs),'TIME' = time-like (FFs))"
      write(7,'(a)') "DISTF    = 'xFitter'                 # Input distributions"
      write(7,'(a)') "#########################################"
      write(7,'(a)') "# Coupling                              #"
      write(7,'(a)') "#########################################"
      write(7,'(a,f9.5,a)') "QREF     = ",Q_ref,",0                  # Reference scale for alpha_s [GeV]"
      write(7,'(a,f9.5,a)') "ASREF    = ",Alphas_ref,",0                  # alpha_s at the refence scale"
      write(7,'(a)') "#########################################"
      write(7,'(a)') "# Heavy quark masses                    #"
      write(7,'(a)') "#########################################"
      write(7,'(a,f9.5,a)') "MC       = ",mch,"                    # Charm threshold  [GeV]"
      write(7,'(a,f9.5,a)') "MB       = ",mbt,"                    # Bottom threshold [GeV]"
      write(7,'(a,f9.5,a)') "MT       = ",mtp,"                    # Top threshold    [GeV]"
      write(7,'(a)') "#########################################"
      write(7,'(a)') "# DIS                                   #"
      write(7,'(a)') "#########################################"
      write(7,'(a)') "PROC     = 'NC'                         # DIS process ('EM','NC')"
      write(7,'(a)') "KRENQ    = 1.0                          # muR / Q"
      write(7,'(a)') "KFACQ    = 1.0                          # muF / Q"
      write(7,'(a)') "#########################################"
      write(7,'(a)')
      close(7)
*
*     Read parameters of the evolution from the card
*
      call ReadParameters(card)
*
*     Initialize MELA
*
      call InitializeEvolution
*
      return
      end
