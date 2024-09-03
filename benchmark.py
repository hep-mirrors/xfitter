import os
import subprocess
import shutil

def get_steering(fname='steering.txt'):
  out = '''
*  Namelist to control input data
*

&InFiles
!  ! Input files:
 !   InputFileNames = 'examples/ploughshare/wplus-thexp.dat'
 
  InputFileNames = 
     !'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCem-thexp.dat',
     'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_920-thexp.dat',
     !'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_820-thexp.dat',
     !'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_575-thexp.dat',
     !'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_460-thexp.dat',
     !'test_data.dat'
&End

&InCorr
  ! Number of correlation (statistical, systematical or full) files
    NCorrFiles = 0
 
  ! Correlation files:
  !  CorrFileNames(1) = 'datafiles/hera/h1/jets/0904.3870/H1_NormInclJets_HighQ2_99-07___H1_NormInclJets_HighQ2_99-07.corr'
&End

&ReduceSyst
    ! even with tolerance =0 the following flag may speed up calculations
  do_reduce = .false.  ! turn-on to simplify/speedup chi2 calculation.
    ! tolerance = 0.0 for exact calculation, > 0.0 for improved speed.
  tolerance = 0.0
    ! depending on blas library, toggling the following flag may improve chi2 computation speed:
  useBlas = .false.
&End

&CovarToNuisance
   ! Global switch for using nuisance param representation for covariance mat.
  LConvertCovToNui = .False.

   ! Tolerance -- zero means exact transformation
  Tolerance = 0.0

   ! (Optional) -- try to subtract diagonal stat. uncertainties from total covariance when determining uncorrelated uncertainites
  LSubtractStat = .false.

   ! The following lines allow to adjust error scaling properties (default: :M)
  DataName     = 'CMS electon Asymmetry rapidity', 'CMS W muon asymmetry'
  DataSystType = ':A', ':A'
&End


*
* (Optional) List systematic sources, modify their scaling properties:
*
&Systematics
 !C      List sources, Results.txt file would list them first. Use the usual :A, :P, 
 !C      qualifiers to change the scalling properties
 !  ListOfSources = 'ATLAS_lumi2010', 'ATL_WZ2010_Source_13:A'
 !C      Modify the prior in chi2 definition (1.0 is default):
 !  PriorScaleName = 'ATLAS_lumi2010', 'ATL_WZ2010_Source_13'
 !  PriorScaleFactor = 0.0, 0.0 
&End


*
* Main steering cards
*
&xFitter 
 ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 !
 ! Chi2 definition. Following options are supported:
 !  
 ! -- Bias corrections for uncertainties --
 ! 'StatScale'    :  'Poisson',  'NoRescale' ( see also 'ExtraSystRescale' below )
 ! 'UncorSysScale':  'Poisson',  'Linear',  'NoRescale'
 ! 'CorSysScale'  :  'Linear',   'NoRescale'
 ! 
 ! -- Treatment of systematics in chi2 ---
 ! 'UncorChi2Type':  'Diagonal'  
 ! 'CorChi2Type'  :  'Hessian', 'Matrix', 'Offset'
 !
 ! -- Extra corrections ---
 !   are given as comma separated list for Chi2ExtraParam, they are off by default.
 !  'PoissonCorr'            : extra log correction accounting for changing uncertainties 
 !  'FirstIterationRescale' : re-scale uncertainties at the first iteration only 
 !  'ExtraSystRescale'      : additional re-scaling of stat. uncertainty to account for syst. shifts.

   CHI2SettingsName = 'StatScale', 'UncorSysScale', 'CorSysScale', 'UncorChi2Type', 'CorChi2Type'
   Chi2Settings     = 'Poisson'  , 'Linear',        'Linear'     , 'Diagonal'     , 'Hessian'
   Chi2ExtraParam = 'PoissonCorr'

 ! Debug flag
  LDEBUG     = False


 ! Quadratic approximation for asymmetric uncertainties
 ! AsymErrorsIterations = 10
&End


*
* Output steering cards
*
&Output 
  ! -- Q2 values at which the pdfs & errors are done (up to 20)
  Q2VAL = 1.9, 2., 3.0, 4.0, 5., 10., 100., 6464, 8317 
!  Q2VAL = 1.9, 4., 10., 100., 6464, 8317 

  ! How many x points to write (standard = 101)
  OUTNX = 101

  ! x-range of output (standard = 1E-4 1.0)
  OUTXRANGE = 1E-4, 0.9999
  ! Write out LHAPDF5 output
  ! WriteLHAPDF5 = true
&End



*
* Process dependent cuts
*
&Cuts

  !--------------------- NC ep  --------------------------

  ! Rule #1: Q2 cuts
   ProcessName(1)     = 'NC e+-p'
   Variable(1)        = 'Q2'
   CutValueMin(1)     = 3.8
   CutValueMax(1)     = 1000000.0

  ! Rule #2: x cuts
   ProcessName(2)     = 'NC e+-p'
   Variable(2)        = 'x'
   CutValueMin(2)     = 0.000001 
   CutValueMax(2)     = 1.0

  !---------------------  CC ep  ------------------

   ProcessName(3)     = 'CC e+-p'
   Variable(3)        = 'Q2'
   CutValueMin(3)     = 3.8
   CutValueMax(3)     = 1000000.0

   ProcessName(4)     = 'CC e+-p'
   Variable(4)        = 'x'
   CutValueMin(4)     = 0.000001 
   CutValueMax(4)     = 1.0

  !-------------------- DY pp  ----------------------

   ProcessName(5)     = 'CC pp'
   Variable(5)        = 'eta1'
   CutValueMin(5)     = -1.
   CutValueMax(5)     = 100.

  !------------------- Jets ---------------------------
   
   ProcessName(6)     = 'pp jets APPLGRID'
   Variable(6)        = 'pt1'
   CutValueMin(6)     = 20.
   CutValueMax(6)     = 1000000.

  !--------------------- Fixed target --------------------------

  ! Rule #7: Whad2 cut
   ProcessName(7)     = 'muon p'
   Variable(7)        = 'Whad2'
   CutValueMin(7)     = 15.   

  !--------------------- Fastnlo jets ----------------------

   ProcessName(8)     = 'FastNLO ep jets'
   Variable(8)        = 'kfac'
   CutValueMin(8)     = 0.0
   CutValueMax(8)     = 2.5

  !--------------------- NC ep charm ----------------

   ProcessName(9)     = 'NC e+-p charm'
   Variable(9)        = 'Q2'
   CutValueMin(9)     = 3.5 
   CutValueMax(9)     = 10000.0

   ProcessName(10)     = 'NC e+-p charm'
   Variable(10)        = 'x'
   CutValueMin(10)     = 0.000001 
   CutValueMax(10)     = 1.0

  !--------------------- by Dataset index ----------------

   ! applied to any (including unspecified) reaction by Dataset indices
   ! e.g. cut pTmax < 15.0 will be applied to datasets 996 and 996
   ! ProcessName must be 'DUMMY'
   !
   !ProcessName(11)     = 'DUMMY'
   !Dataset( 1,11)      = 995
   !Dataset( 1,11)      = 996
   !Variable(11)        = 'pTmax'
   !CutValueMin(11)     = 0.0
   !CutValueMax(11)     = 15.0
   
&End

*
* (Optional) MC errors steering cards
*
&MCErrors
  ! Activate MC method for error estimation if lRand = True
  lRAND   = False
  
  ! Use data (true, default) or theory (false) for the central values of the MC replica
  lRANDDATA = True

  ! MC method Seed
  ISeedMC = 123456 

  ! --- Choose what distribution for the random number generator 
  ! STATYPE (SYS_TYPE)  =   1  gauss
  ! STATYPE (SYS_TYPE)  =   2  uniform
  ! STATYPE (SYS_TYPE)  =   3  lognormal
  ! STATYPE (SYS_TYPE)  =   4  poisson (only for lRANDDATA = False !)
  !
  ! STATYPE = SYS_TYPE = 0 with lRANDDATA = False and lRAND   = True yeilds in substitution of data by theory (useful for sensitivity tests)
  !
  STATYPE =  1
  SYSTYPE =  1
&End
'''
  with open(fname, 'w') as fout:
    fout.write(out)

def get_constants(fname='constants.yaml'):
  out = f'''
# EW parameters
Mz : 91.1876
Mw : 80.385
Mh : 125.9

Wz : 2.4952
Ww : 2.085
Wh : 1d-3
Wtp: 2.0d0

gf : 1.16638e-5
convFac : 0.389379338e9
alphaem : 7.29735e-3

sin2thW :  0.23127
Vud : 0.97427
Vus : 0.2254
Vub : 0.00358
Vcd : 0.22520
Vcs : 0.97344 
Vcb : 0.04156
Vtd : 0.00872
Vts : 0.04076
Vtb : 0.999133

# lepton masses:
men : 1e-10
mel : 0.510998928e-3
mmn : 1e-10
mmo : 0.1056583715
mtn : 1e-10
mta : 1.77682

# light quark masses:
mup : 0.06983
mdn : 0.06983
mst : 0.150

# heavy quark masses:
mch : 1.43 
mbt : 4.50
mtp : 173.0

# QCD parameters
Order: {Order}
NFlavour: {NFlavour}
isFFNS: {isFFNS} # 0 (default): VFNS with number of flavours changing from 3 to NFlavour, 1: FFNS with number of flavours = NFlavour at any scale
'''
  with open(fname, 'w') as fout:
    fout.write(out)

def get_parameters(fname='parameters.yaml'):
  out = f'''
Minimizer: MINUIT # CERES
MINUIT:
  Commands: |
    set str 2
    call fcn 3
 #   migrad
 #   hesse
 #   call fcn 3
    
#  doErrors :  Hesse # None

CERES:
  offset: 2
  tolerance: 1e-5
  strategy: 0
  covariance: 1

Parameters:
  Ag   :  DEPENDENT
  Bg   : [ -0.061953, 0.27 ]
  Cg   : [  5.562367,  0.32 ]
  Agp  : [ 0.07311, 0.01 ]  # negative gluon ....
  Bgp  : [ -0.383100, 0.01 ]
  Cgp  : [ 25.0, 0.]  # fix C of negative gluon
  Auv  :  DEPENDENT
  Buv  : [ 0.810476, 0.016 ]
  Cuv  : [ 4.823512, 0.06 ]
  Duv  : [    0     ]
  Euv  : [ 9.921366, 0.8 ]
  Adv  :  DEPENDENT
  Bdv  : [ 1.029995, 0.06 ]
  Cdv  : [ 4.846279, 0.3 ]
  Aubar: [ 0.0, 0.0 ] # not used (Aubar=Adbar)
  Bubar: [ 0.0, 0.0  ] # not used (Bubar=Bdbar)
  Cubar: [ 7.059694, 0.8 ]
  Dubar: [ 1.548098, 1.0 ]
  Adbar: [ 0.1613, 0.01 ]
  Bdbar: [ -0.1273, 0.004  ]
  Cdbar: # another example of providing value, step etc.
    value: 9.586246
    step: 1.2345
    #min
    #max
    #pr_mean
    #pr_sigma
  ZERO : [ 0. ]          # zero
  fs   :   0.4   #no step means fixed
  DbarToS: "=fs/(1-fs)"

Parameterisations:
  par_uv:
    class: HERAPDF
    parameters: [Auv,Buv,Cuv,Duv,Euv]
  par_dv:
    class: HERAPDF
    parameters: [Adv,Bdv,Cdv]
  par_ubar:
    class: HERAPDF
    parameters: [Adbar,Bdbar,Cubar,Dubar]
  par_dbar:
    class: HERAPDF
    parameters: [Adbar,Bdbar,Cdbar]
  par_s: # s=fs/(1-fs) * Dbar
    class: Factor
    factor: DbarToS #name of parameter
    input: par_dbar
  par_g:
    class: NegativeGluon
    parameters: [Ag,Bg,Cg,ZERO,ZERO,Agp,Bgp,Cgp]
  #par_s:
    #class: HERAPDF
    #parameters: [As,Bs,Cs]
  #par_s:
    #class: Expression
    #expression: "Adbar*fs/(1-fs)*(x^Bdbar*(1-x)^Cdbar)"
  # Another example for Expression parameterisation
  #par_g:
    #class: Expression
    #expression: "Ag*(x^Bg*(1-x)^Cg-Agp*x^Bgp*(1-x)^Cgp)"

DefaultDecomposition: proton
Decompositions:
  proton:
    class: UvDvUbarDbarS
    xuv: par_uv
    xdv: par_dv
    xubar: par_ubar
    xdbar: par_dbar
    xs: par_s
    xg: par_g

DefaultEvolution: proton-{DefaultEvolution}

Evolutions:
  proton-APFELff:
    ? !include evolutions/APFEL.yaml
    qLimits : [1.0, 50000.0]
  proton-QCDNUM:
    ? !include evolutions/QCDNUM.yaml
    # The following allows QCDNUM to read PDFs from other evolutions:
    #EvolutionCopy: "proton-LHAPDF"
  proton-LHAPDF:
    class: LHAPDF
    set: "NNPDF30_nlo_as_0118"
    #set: "CT10nlo"
    member: 0
  proton-Hoppet:
    class: Hoppet

#  proton-APFEL:
#    ? !include evolutions/APFELxx.yaml
#    decomposition: proton
  antiproton:
    class: FlipCharge
    #input: proton-QCDNUM
    input: proton-LHAPDF
#  neutron:
#    class: FlipUD
#    input: proton-QCDNUM

Q0 : {Q0} # Initial scale =sqrt(1.9)

? !include constants.yaml

alphas : {alphas}

byReaction:
  # RT DIS scheme settings:
  RT_DISNC:
    ? !include reactions/RT_DISNC.yaml
  HOPPET_DISNC:
    ? !include reactions/HOPPET_DISNC.yaml
    # uncomment if defaultEvolution is not QCDNUM: RT_DISNC works with QCDNUM only, use EvolutionCopy
    #evolution: proton-QCDNUM
  # uncomment if defaultEvolution is not QCDNUM: RT_DISNC works with QCDNUM only, use EvolutionCopy
  #BaseDISCC:
  #  evolution: proton-QCDNUM
  # FONLL scheme settings:
  FONLL_DISNC:
    ? !include reactions/FONLL_DISNC.yaml
  FONLL_DISCC:
    ? !include reactions/FONLL_DISCC.yaml
  # FF ABM scheme settings:
  FFABM_DISNC:
    ? !include reactions/FFABM_DISNC.yaml
  FFABM_DISCC:
    ? !include reactions/FFABM_DISCC.yaml
  # AFB settings:
  AFB:
    ? !include reactions/AFB.yaml
  # APPLgrid settings:
  APPLgrid:
    ? !include reactions/APPLgrid.yaml
  # APPLgrid settings:
  # (optional) APFELgrid module settings:
  #  ? !include reactions/APFELgrid.yaml
  # (optional) Fractal module settings:
  Fractal_DISNC:
    ? !include reactions/Fractal_DISNC.yaml
#  DYTurbo:
#    ? !include reactions/DYTurbo.yaml

#byDataset: #Here one can redefine some parameters for specific datasets
#  #Parameter definitions here have the highest priority: they override both "byReaction" and "TermInfo"
#  "HERA1+2 NCep 920":
#    epolarity: 2

# Specify HF scheme used for DIS NC processes:
hf_scheme_DISNC :
 defaultValue : '{hf_scheme_DISNC}'
# defaultValue : 'RT_DISNC'        # global specification
  #defaultValue : 'BaseDISNC'       # global specification
#  defaultValue : 'FONLL_DISNC'     # global specification
#  defaultValue : 'FFABM_DISNC'
#  'HERA1+2 NCep 920' : 'BaseDISNC' # datafile specific (based on name)
#  1 : BaseDISNC
#  'HERA1+2 NCep 920' : 'Fractal_DISNC'  # Fractal model. Add parameters file if you want to try it (see above)

# Specify HF scheme used for DIS CC processes:
hf_scheme_DISCC :
  defaultValue : '{hf_scheme_DISCC}'       # global specification
#  defaultValue : 'FONLL_DISCC'     # global specification
#  defaultValue : 'FFABM_DISCC'     # global specification

#
# Profiler allows to add variations of parameters and PDF eigenvectors as additional nuisance parameters
#
Profiler:
  Parameters:
    alphas: [ 0.118, 0.119, 0.117 ]  # central, up, (down) variation. If down is not given, uses symmetrizes Up variation 
  #Evolutions:
  #  proton-LHAPDF:
  #    sets:    [CT10]
  #    members: [[0,1,end]]
  Status: "Off"                 # "Off" to turn off profiler
  WriteTheo: "Off"              # Can be "Off", "On" or "Asymmetric" (to store asymmetric variations)
  getChi2: "Off"                # determine and report chi2 for each variation
  enableExternalProfiler: "Off" # enable creation of additional files, needed for xfitter draw

OutputDirectory: "output" #Can be omitted, default is "output"
#OutputDirectory: "ApfelFF" #Can be omitted, default is "output"


#
# Possible levels to stop program execution:
#  1 - will stop on warnings
#  2 - will stop on errors (default)
#  3 - will stop on severe errors
#  4 - will stop on fatal
#  5 - will not stop on any error

MaxErrAllowed: 2
'''
  with open(fname, 'w') as fout:
    fout.write(out)

def run_cmd(cmd):
  print(cmd)
  subprocess.call(cmd.split())

if __name__ == '__main__':
  evolutions = ['APFELff']#, 'QCDNUM', 'Hoppet']

  Orders = ['LO']
  isFFNSs = [0]
  NFlavours = [5]
  
  alphas = 0.118
  Q0 = 1.378404875209
  hf_scheme_DISNC = 'HOPPET_DISNC'
  hf_scheme_DISCC = 'BaseDISCC'
  basedirname = 'bench'
  xfitter = '/home/zenaiev/soft/xfitter-hoppet/bin/xfitter'
  xfitterdraw = '/home/zenaiev/soft/xfitter-hoppet/bin/xfitter-draw'
  datafiles = '/home/zenaiev/soft/xfitter-datafiles'

  cwd = os.getcwd()
  for Order in Orders:
    for isFFNS in isFFNSs:
      for NFlavour in NFlavours:
        for DefaultEvolution in evolutions:
          dirname = f'{basedirname}/run_Order{Order}_isFFNS{isFFNS}_NFlavour{NFlavour}/{DefaultEvolution}'
          os.chdir(cwd)
          if os.path.exists(dirname):
            shutil.rmtree(dirname)
          os.makedirs(dirname)
          os.chdir(dirname)
          get_steering()
          get_constants()
          get_parameters()
          os.symlink(datafiles, './datafiles')
          os.symlink(xfitter, './xfitter')
          run_cmd(f'./xfitter')
          #os.symlink(xfitter, './xfitter')
          #run_cmd(f'./{xfitterdraw}')
