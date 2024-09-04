#!/usr/bin/env python
 
import os
import subprocess
import shutil
import sys

def get_steering(fname='steering.txt'):
  out = '''
*  Namelist to control input data
*

&InFiles
!  ! Input files:
    InputFileNames = '''
  for data in InputFileNames:
    out += f"\n      '{data}',"
  out += f'''
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
  !Q2VAL = 1.9, 2., 3.0, 4.0, 5., 10., 100., 6464, 8317 
  !Q2VAL = 1.9, 2., 10., 8317 
  Q2VAL = 1.9, 8317 
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
  proton-{DefaultEvolution}:'''
  if DefaultEvolution in evolutions_with_yaml_file:
    out += f'''
    ? !include evolutions/{DefaultEvolution}.yaml'''
  else:
    out += f'''
    class: {DefaultEvolution}'''
  out += f'''

Q0 : {Q0} # Initial scale =sqrt(1.9)

? !include constants.yaml

alphas : {alphas}
'''
  if hf_scheme_DISNC in reactions_with_yaml_file:
    out += f'''
byReaction:
  {hf_scheme_DISNC}:
    ? !include reactions/{hf_scheme_DISNC}.yaml
'''
  if hf_scheme_DISCC in reactions_with_yaml_file:
    out += f'''
byReaction:
  {hf_scheme_DISCC}:
    ? !include reactions/{hf_scheme_DISCC}.yaml
'''
  out += f'''
# Specify HF scheme used for DIS NC processes:
hf_scheme_DISNC :
 defaultValue : '{hf_scheme_DISNC}'
# Specify HF scheme used for DIS CC processes:
hf_scheme_DISCC :
  defaultValue : '{hf_scheme_DISCC}'       # global specification

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
  logfile = cmd.split()[0] + '.log'
  with open(logfile, 'w') as fout:
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in proc.stdout:
      sys.stdout.write(line.decode('utf-8'))
      fout.write(line.decode('utf-8'))
    for line in proc.stderr:
      sys.stderr.write(line.decode('utf-8'))
      fout.write(line.decode('utf-8'))
    proc.wait()
  if proc.returncode != 0:
    sys.exit(1)

def recreate_dir(dirname, cd=False, cwd=None):
  if cwd is not None:
    os.chdir(cwd)
  if os.path.exists(dirname):
    shutil.rmtree(dirname)
  os.makedirs(dirname)
  if cd:
    os.chdir(dirname)

def benchmark_run(suffix=None):
  #hf_scheme_DISNC = {'QCDNUM': 'BaseDISNC', 'APFEL': 'FONLL_DISNC', 'APFELxx': 'N3LO_DISNC', 'Hoppet': 'HOPPET_DISNC'}[DefaultEvolution]
  if suffix is None:
    suffix = f'DefaultEvolution{DefaultEvolution}_hf_scheme_DISNC{hf_scheme_DISNC}_hf_scheme_DISCC{hf_scheme_DISCC}'
  recreate_dir(f'{dirname}/{suffix}', cd=True, cwd=startdir)
  get_steering()
  get_constants()
  get_parameters()
  os.symlink(datafiles, './datafiles')
  os.symlink(xfitter, './xfitter')
  if DefaultEvolution == 'APFEL' and Order == 'NNNLO':
    os.remove('./xfitter')
    os.symlink(xfitter_n3lo, './xfitter')
  run_cmd(f'./xfitter')
  return suffix

def benchmark_plot(outputs, labels, extraopts):
  if len(outputs) == 0: return
  os.chdir(startdir)
  os.chdir(dirname)
  os.symlink(xfitterdraw, './xfitter-draw')
  run_cmd(f'./xfitter-draw {xfitterdraw_opts} {extraopts} --outdir plots ' + ' '.join(f'{outputs[ioutput]}/output:{labels[ioutput]}' for ioutput in range(len(outputs))))
  print(f'All plots stored in {os.getcwd()}/plots/plots.pdf')
  os.chdir('plots')
  if len(outputs) > 1:
    for q2 in ['1.9', '8317']:
      run_cmd(f'pdfunite q2_{q2}_pdf_uv.pdf q2_{q2}_pdf_dv.pdf q2_{q2}_pdf_g.pdf q2_{q2}_pdf_Sea.pdf q2_{q2}_pdf_uv_ratio.pdf q2_{q2}_pdf_dv_ratio.pdf q2_{q2}_pdf_g_ratio.pdf q2_{q2}_pdf_Sea_ratio.pdf pdfs_{q2}.pdf')
      run_cmd(f'pdfjam --nup 4x2 pdfs_{q2}.pdf --outfile pdfs_compact_{q2}.pdf --papersize {{24cm,12cm}}')
      print(f'Compact PDF plots stored in {os.getcwd()}/pdfs_compact_{q2}.pdf')


if __name__ == '__main__':
  fname_coms = 'bench_coms.txt'
  #run_cmd('ls bin ')
  #aaa
  reactions_with_yaml_file = ['FFABM_DISCC', 'FONLL_DISCC', 'HOPPET_DISNC', 'FFABM_DISNC', 'FONLL_DISNC', 'RT_DISNC']
  evolutions_with_yaml_file = ['APFELxx', 'APFEL', 'Hoppet', 'QCDNUM']

  #evolutions = ['APFELxx', 'APFEL', 'QCDNUM', 'Hoppet']
  #evolutions = ['QCDNUM', 'Hoppet']
  #evolutions = ['QCDNUM']
  #evolutions = ['APFELff', 'QCDNUM', 'Hoppet']
  #evolutions = ['APFEL']
  #evolutions = ['APFELff', 'APFEL']
  #evolutions = ['APFEL', 'Hoppet']
  DefaultEvolutions = []

  hf_scheme_DISNCs = ['BaseDISNC', 'HOPPET_DISNC']

  #Orders = ['LO']
  #Orders = ['LO', 'NLO', 'NNLO']
  Orders = ['NNLO']
  #Orders = ['NNNLO']
  isFFNSs = [0]
  NFlavours = [5]
  #isFFNSs = [1]
  #NFlavours = [3]
  #isFFNSs = [0,1]
  #NFlavours = [3,4,5]
  
  alphas = 0.118
  Q0 = 1.378404875209
  #hf_scheme_DISNC = 'HOPPET_DISNC'
  #hf_scheme_DISNC = 'BaseDISNC'
  hf_scheme_DISCC = 'BaseDISCC'
  basedirname = 'runs_bench'
  xfitter = '/home/zenaiev/soft/xfitter-hoppet/bin/xfitter'
  xfitter_n3lo = '/home/zenaiev/soft/xfitter-hoppet/bin/xfitter_n3lo'
  xfitterdraw = '/home/zenaiev/soft/xfitter-hoppet/bin/xfitter-draw'
  datafiles = '/home/zenaiev/soft/xfitter-datafiles'
  xfitterdraw_opts = '--no-logo'
  #xfitterdraw_opts += ' --q2all'
  xfitterdraw_opts += ' --splitplots-pdf'
  #xfitterdraw_opts += ' --ratiorange 0.99:1.01'
  #xfitterdraw_opts += ' --ratiorange 0.999:1.001'
  xfitterdraw_opts += ' --no-shifts'
  #xfitterdraw_opts += ' --no-tables'

  startdir = os.getcwd()
  for Order in Orders:
    FONLLVariant = 'B' if Order == 'NLO' else 'C'
    for isFFNS in isFFNSs:
      for NFlavour in NFlavours:
        if isFFNS == 0 and NFlavour != 5: continue
        outputs = []
        # loop over evolutions
        InputFileNames = None
        hf_scheme_DISNC = None
        dirname = f'{basedirname}/Order{Order}_isFFNS{isFFNS}_NFlavour{NFlavour}/evolutions'
        recreate_dir(dirname, cd=True, cwd=startdir)
        for DefaultEvolution in DefaultEvolutions:
          outputs.append(benchmark_run(suffix=DefaultEvolutions))
        benchmark_plot(outputs, DefaultEvolutions, extraopts=' --q2all --ratiorange 0.99:1.01 --no-tables')
        # loop over reactions
        DefaultEvolution = 'QCDNUM'
        InputFileNames = ['datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_920-thexp.dat']
        outputs = []
        dirname = f'{basedirname}/Order{Order}_isFFNS{isFFNS}_NFlavour{NFlavour}/reactions'
        recreate_dir(dirname, cd=True, cwd=startdir)
        for hf_scheme_DISNC in hf_scheme_DISNCs:
          outputs.append(benchmark_run(suffix=hf_scheme_DISNC))
        benchmark_plot(outputs, hf_scheme_DISNCs, extraopts='--no-pdfs --only-theory')
