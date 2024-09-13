#!/usr/bin/env python
 
import os
import subprocess
import shutil
import sys
import argparse
import glob
import time
import numpy as np

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

def enable_evolution(evolution):
  out = f'''
  proton-{evolution}:'''
  if evolution in evolutions_with_yaml_file:
    out += f'''
    ? !include evolutions/{evolution}.yaml'''
  else:
    out += f'''
    class: {evolution}'''
  return out

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

Evolutions:'''
  out += enable_evolution(DefaultEvolution)
  if extraEvolutionLines is not None:
    out += extraEvolutionLines
  out += f'''

Q0 : {Q0} # Initial scale =sqrt(1.9)

? !include constants.yaml

alphas : {alphas}
'''
  if hf_scheme_DISNC in reactions_with_yaml_file:
    out += f'''
byReaction:
  {hf_scheme_DISNC}:
    ? !include reactions/{hf_scheme_DISNC}.yaml'''
  if hf_scheme_DISCC in reactions_with_yaml_file:
    out += f'''
  {hf_scheme_DISCC}:
    ? !include reactions/{hf_scheme_DISCC}.yaml
'''
    if extraReactionLines is not None:
      out += extraReactionLines
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

def make_datafile_dis(fname, current, Q2s, xs, sqrts, charge):
  charge_title = {1: '+', -1: '-'}[charge]
  out = f'''* PSEUDODATA for DIS SF benchmark
&Data 
  Name = "PSEUDODATA" 
  Reaction = "{current} e+-p"
  TermName = 'R'
  TermSource = 'use:hf_scheme_DIS{current}'
  TermInfo = 'type=sigred:flav=incl:echarge={charge}:epolarity=0'
  TheorExpr = 'R'
  NData = {sum(len(xs[v]) for v in Q2s)}
  NColumn =   5
  ColumnType = 3*"Bin","Sigma", 1*"Error"
  ColumnName = "Q2","x","y","Sigma", "stat"
  Percent = 1*true
&End 
&PlotDesc
   PlotN = {len(Q2s)}
   PlotDefColumn = 'Q2'
   PlotDefValue = {','.join([str(v-0.001) for v in Q2s] + [str(Q2s[-1]+0.001)])}
   PlotVarColumn = \'x\''''
  for iq2 in range(len(Q2s)):
    out += f'''
   PlotOptions({iq2+1})  = 'Experiment:PSEUDO @ExtraLabel:e^{{{charge_title}}}p #rightarrow e^{{{charge_title}}}X ({current}) Q^{{2}} = {str(Q2s[iq2])} GeV^{{2}} #sqrt{{s}} = {sqrts:.0f} GeV @XTitle: x @YTitle: #sigma_{{red}}  @Title: @Xlog\''''
  out += f'''
&End'''
  out += f'''
*{"Q2":>12s}{"x":>12s}{"y":>12s}{"Sigma":>12s}{"stat":>12s}'''
  for q2 in Q2s:
    for x in xs[q2]:
      y = sqrts**2/(q2*x)
      sigma = 1.0
      stat = 50.0
      out += f'''
 {q2:12.4e}{x:12.4e}{y:12.4e}{sigma:12.4e}{stat:12.4e}'''
  with open(fname, 'w') as fout:
    fout.write(out)

def run_cmd(cmd):
  start = time.time()
  print(f'running command: {cmd} [{os.getcwd()}] ... ', end = '', flush = True)
  logfile = cmd.split()[0] + '.log'
  with open(logfile, 'w') as fout:
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in proc.stdout:
      fout.write(line.decode('utf-8'))
      if args.verbose:
        sys.stdout.write(line.decode('utf-8'))
    for line in proc.stderr:
      fout.write(line.decode('utf-8'))
      if args.verbose:
        sys.stderr.write(line.decode('utf-8'))
    proc.wait()
  print(f'[took {time.time()-start:.2f}s]', end = '', flush = True)
  if proc.returncode != 0:
    print(f' FAILED')
    sys.exit(1)
  with open(logfile) as fout:
    if any('******** HF_STOP forced program termination ********' in l for l in open(logfile).readlines()):
      print(f' FAILED')
      sys.exit(1)
  print(f'')


def recreate_dir(dirname, cd=False, cwd=None):
  if cwd is not None:
    os.chdir(cwd)
  if os.path.exists(dirname):
    shutil.rmtree(dirname)
  os.makedirs(dirname)
  if cd:
    os.chdir(dirname)

def benchmark_run(label=None):
  if label is None:
    label = f'DefaultEvolution{DefaultEvolution}_hf_scheme_DISNC{hf_scheme_DISNC}_hf_scheme_DISCC{hf_scheme_DISCC}'
  recreate_dir(f'{dirname}/{label}', cd=True, cwd=startdir)
  get_steering()
  get_constants()
  get_parameters()
  os.symlink(datafiles, './datafiles')
  os.symlink(xfitter, './xfitter')
  #if DefaultEvolution == 'APFELxx' and Order == 'NNNLO' or hf_scheme_DISNC == 'N3LO_DISNC' or hf_scheme_DISCC == 'N3LO_DISCC':
  #  os.remove('./xfitter')
  #  os.symlink(xfitter_n3lo, './xfitter')
  run_cmd(f'./xfitter')
  return label

def calc_max_dif(thpreds, thpreds_ref, eps=None):
  maxdif = 0
  for ith in range(len(thpreds)):
    if thpreds_ref[ith] != 0.:
      reldif = abs(thpreds[ith] / thpreds_ref[ith] - 1.)
    elif thpreds[ith] == 0.:
      reldif = 0.
    else:
      reldif = thpreds[ith]
    maxthpred = max(abs(thpreds[ith]), abs(thpreds_ref[ith]))
    # ignore differences if numbers are smaller than eps
    if eps is None or maxthpred > eps:
      dif = reldif
    else:
      dif = 0
      #dif = maxthpred
    maxdif = max(maxdif, dif)
  return maxdif

def read_pdfs(fname):
  thpreds = []
  for fname in sorted(glob.glob(f'{fname}/output/pdfs_q2val_*.txt')):
    with open(fname) as f:
      for l in f.readlines():
        if len(l.split()) == 16 and l.split()[0] != 'x':
          #thpreds += [float(v) for v in l.split()]
          # check PDFs only in a stable range of x
          x = float(l.split()[0])
          if x > 1e-4 and x < 0.1:
            thpreds += [float(v) for v in l.split()[1:]]
  return thpreds

def read_thpreds(fname):
  with open(fname) as f:
    thpreds = [float(l.split()[6]) for l in f.readlines() if len(l.split()) == 13]
  return thpreds

def make_plots(outputs, extraopts):
  os.symlink(xfitterdraw, './xfitter-draw')
  run_cmd(f'./xfitter-draw {xfitterdraw_opts} {extraopts} --outdir plots ' + ' '.join(f'{outputs[ioutput]}/output:{outputs[ioutput]}' for ioutput in range(len(outputs))))
  print(f'All plots stored in {os.getcwd()}/plots/plots.pdf')
  # make compact PDF plots for Q2 = 1.9, 8317 GeV2
  if '--no-pdfs' not in f'{xfitterdraw_opts} {extraopts}':
    if all(shutil.which(com) for com in ['pdfunite', 'pdfjam']):
      os.chdir('plots')
      if len(outputs) > 1:
        for q2 in ['1.9', '8317']:
          run_cmd(f'pdfunite q2_{q2}_pdf_uv.pdf q2_{q2}_pdf_dv.pdf q2_{q2}_pdf_g.pdf q2_{q2}_pdf_Sea.pdf q2_{q2}_pdf_uv_ratio.pdf q2_{q2}_pdf_dv_ratio.pdf q2_{q2}_pdf_g_ratio.pdf q2_{q2}_pdf_Sea_ratio.pdf pdfs_{q2}.pdf')
          run_cmd(f'pdfjam --nup 4x2 pdfs_{q2}.pdf --outfile pdfs_compact_{q2}.pdf --papersize {{30cm,15cm}}')
          print(f'Compact PDF plots stored in {os.getcwd()}/pdfs_compact_{q2}.pdf')
    else:
      print(f'skipping production of compact PDF plots because commands "pdfunite" or "pdfjam" are not available')

def benchmark_results(outputs, extraopts):
  if len(outputs) == 0: return
  os.chdir(startdir + '/' + dirname)
  # calculate max. differences
  if len(outputs) >= 2:
    # max. difference in PDFs
    thpreds_ref = read_pdfs(outputs[0])
    maxdiffs = [0] * len(outputs)
    for ioutput in range(1, len(outputs)):
      thpreds = read_pdfs(outputs[ioutput])
      maxdiffs[ioutput] = calc_max_dif(thpreds, thpreds_ref, eps = 1e-2)
    print(f'max. relative difference in PDFs ' + ', '.join(f'{outputs[0]} vs {outputs[ioutput]} = {maxdiffs[ioutput]:.1e}' for ioutput in range(1, len(outputs))))
    # max. difference in theory predictions
    thpreds_ref = read_thpreds(f'{outputs[0]}/output/fittedresults.txt')
    if len(thpreds_ref) > 0:
      maxdiffs = [0] * len(outputs)
      for ioutput in range(1, len(outputs)):
        thpreds = read_thpreds(f'{outputs[ioutput]}/output/fittedresults.txt')
        maxdiffs[ioutput] = calc_max_dif(thpreds, thpreds_ref)
      print(f'max. relative difference in theory predictions ' + ', '.join(f'{outputs[0]} vs {outputs[ioutput]} = {maxdiffs[ioutput]:.1e}' for ioutput in range(1, len(outputs))))
  # make plots
  if args.plot:
    make_plots(outputs, extraopts)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Benchamark xFitter evolutions and/or reactions')
  parser.add_argument('-v', '--verbose', action='store_true', help='Verbose (print output of commands)')
  parser.add_argument('-p', '--plot', action='store_true', help='Produce plots')
  args = parser.parse_args()

  reactions_with_yaml_file = ['FFABM_DISNC', 'FFABM_DISCC', 'FONLL_DISNC', 'FONLL_DISCC', 'HOPPET_DISNC', 'HOPPET_DISCC', 'RT_DISNC', 'N3LO_DISNC', 'N3LO_DISCC']
  evolutions_with_yaml_file = ['APFELxx', 'APFEL', 'HOPPET', 'QCDNUM']

  #DefaultEvolutions = ['QCDNUM']
  #DefaultEvolutions = ['APFELxx', 'APFEL', 'QCDNUM', 'HOPPET']
  #DefaultEvolutions = ['APFELxx', 'HOPPET']
  DefaultEvolutions = []

  # hf_scheme_DISNCs is a list where every entry is [label, hf_scheme_DISNC, extraEvolutionLines, extraReactionLines]
  # label: directory name and plotting label
  # hf_scheme_DISNC: name for parameters.yaml
  # extraEvolutionLines: string to be appended to node "Evolutions" (None if nothing to append)
  # extraReactionLines: string to be appended to node "<scheme>:" (None if nothing to append)
  hf_scheme_DISNCs = [
    #['HOPPET_N3LO', 'HOPPET_DISNC', None, 'Order_HOPPET_Evolution: NNLO'],
    #['APFELxx_ZMVFNS_N3LO', 'N3LO_DISNC', '\n  proton-APFELxx:\n    ? !include evolutions/APFELxx.yaml\n', '    massive: 0\nOrder_HOPPET_Evolution: NNLO'],
    #['HOPPET_N2LO', 'HOPPET_DISNC', None, 'Order: NNLO'],
    #['APFELxx_ZMVFNS_N2LO', 'N3LO_DISNC', '\n  proton-APFELxx:\n    ? !include evolutions/APFELxx.yaml\n', '    massive: 0\nOrder: NNLO'],
    ['HOPPET', 'HOPPET_DISNC', None, None],
    ['FONLL_ZMVFNS', 'FONLL_DISNC', '\n  proton-APFEL:\n    ? !include evolutions/APFEL.yaml\n    FONLLVariant: {FONLLVariant}\n    MassScheme: \'ZM-VFNS\'', None],
    ['Base', 'BaseDISNC', None, None],
    ##['FFABM', 'FFABM_DISNC', None, None],
    ##['FONLL', 'FONLL_DISNC', '\n  proton-APFEL:\n    ? !include evolutions/APFEL.yaml', None],
  ]
  #hf_scheme_DISNCs = []

  #Orders = ['LO']
  Orders = ['NLO']
  #Orders = ['NNLO']
  #Orders = ['NNNLO']
  #Orders = ['LO', 'NLO', 'NNLO']
  isFFNSs = [0]
  NFlavours = [5]
  #isFFNSs = [1]
  #NFlavours = [3]
  #isFFNSs = [0,1]
  #NFlavours = [3,4,5]
  
  alphas = 0.118
  Q0 = 1.378404875209
  basedirname = 'runs_bench'
  xfitter = '/home/zenaiev/soft/xfitter-hoppet/bin/xfitter'
  #xfitter = '/home/zenaiev/soft/xfitter-n3lo/bin/xfitter'
  #xfitter_n3lo = '/home/zenaiev/soft/xfitter-n3lo/bin/xfitter' # currently this is in a separate branch
  xfitterdraw = '/home/zenaiev/soft/xfitter-hoppet/bin/xfitter-draw'
  datafiles = '/home/zenaiev/soft/xfitter-datafiles'
  xfitterdraw_opts = '--no-logo'
  xfitterdraw_opts += ' --splitplots-pdf'
  xfitterdraw_opts += ' --no-shifts'
  #xfitterdraw_opts += ' --no-tables'

  startdir = os.getcwd()
  for Order in Orders:
    FONLLVariant = 'B' if Order == 'NLO' else 'C' # needed for APFELff, otherwise it complains
    for isFFNS in isFFNSs:
      for NFlavour in NFlavours:
        if isFFNS == 0 and NFlavour != 5: continue # skip VFNS with nf<5 
        
        # benchmark evolutions
        InputFileNames = []
        hf_scheme_DISNC = None
        hf_scheme_DISCC = None
        dirname = f'{basedirname}/Order{Order}_isFFNS{isFFNS}_NFlavour{NFlavour}/evolutions'
        recreate_dir(dirname, cd=True, cwd=startdir)
        outputs = []
        for DefaultEvolution in DefaultEvolutions:
          outputs.append(benchmark_run(label=DefaultEvolution))
        benchmark_results(outputs, extraopts=' --q2all --ratiorange 0.99:1.01 --no-tables')
        
        # benchmark reactions
        DefaultEvolution = 'QCDNUM'
        #DefaultEvolution = 'APFELxx'
        if Order == 'NNNLO':
          DefaultEvolution = 'HOPPET'
        xs_NC = {
          5: np.logspace(np.log10(5e-5), np.log10(0.65), 50),
          50: np.logspace(np.log10(5e-4), np.log10(0.65), 50),
          500: np.logspace(np.log10(5e-3), np.log10(0.65), 50),
          30000: np.logspace(np.log10(0.3), np.log10(0.75), 50),
        }
        xs_CC = {
          5: np.logspace(np.log10(5e-5), np.log10(0.15), 50),
          50: np.logspace(np.log10(5e-4), np.log10(0.15), 50),
          500: np.logspace(np.log10(5e-3), np.log10(0.15), 50),
          30000: np.logspace(np.log10(0.3), np.log10(0.75), 50),
        }
        dirname = f'{basedirname}/Order{Order}_isFFNS{isFFNS}_NFlavour{NFlavour}/reactions'
        recreate_dir(dirname, cd=True, cwd=startdir)
        make_datafile_dis('NCep.dat', 'NC', list(xs_NC.keys()), xs_NC, 318., +1)
        make_datafile_dis('NCem.dat', 'NC', list(xs_NC.keys()), xs_NC, 318., -1)
        make_datafile_dis('CCep.dat', 'CC', list(xs_CC.keys()), xs_CC, 318., +1)
        make_datafile_dis('CCem.dat', 'CC', list(xs_CC.keys()), xs_CC, 318., -1)
        InputFileNames = [
          '../NCep.dat',
          '../NCem.dat',
          '../CCep.dat',
          '../CCem.dat',
          #'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCem-thexp.dat',
          #'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_920-thexp.dat',
          #'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_820-thexp.dat',
          #'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_575-thexp.dat',
          #'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_460-thexp.dat',
        ]
        outputs = []
        for entry in hf_scheme_DISNCs:
          hf_scheme_DISNC = entry[1]
          hf_scheme_DISCC = hf_scheme_DISNC.replace('NC', 'CC')
          extraEvolutionLines = entry[2]
          if extraEvolutionLines is not None:
            extraEvolutionLines = extraEvolutionLines.format(FONLLVariant=FONLLVariant)
          extraReactionLines = entry[3]
          outputs.append(benchmark_run(label=entry[0]))
        #benchmark_results(outputs, extraopts='--no-pdfs --only-theory')
        benchmark_results(outputs, extraopts='--no-pdfs --only-theory --no-tables')
