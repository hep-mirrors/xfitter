#!/bin/sh
# the next line restarts using tclsh \
if [ -z $DISPLAY ]; then exec tclsh "$0" "$@"; else exec wish "$0" -- "$@"; fi

catch {source $tcl_rcFileName}
source [file join [file dirname [info script]] job_farm_fns.tcl]

#-------------------------------------------------
#  Retrieve HERAFitter jobs from a farm
#-------------------------------------------------
set Version 1.0

set _TEST_JS 0
set Verbose 1

# --- redefine proc TestJobResults
# ===========================
proc TestJobResults {args} {
  set outzip [lindex $args 0]
  if {[catch {exec unzip -l $outzip params_????.txt} ans]} {
    # Say "*** ERROR in TestJobResults:" err
    return 1
  }
  return 0
}

# ###############################################

MainStartUp

Say "|-----  Using farm $FarmName  -----|" info

set Args [GetOptions {
  lbl "Last"
  peek 0
  ?run 0
  jid  0
} Opts $argv]
# Say $Args
# Say [array get Opts]

set Label $Opts(lbl)
set nWaitSec $Opts(peek)
if {$nWaitSec > 0 && $nWaitSec < $MinWaitSec} {
  Say "Peek interval increased to $MinWaitSec seconds." warnb
  set nWaitSec $MinWaitSec
}

if {$Opts(jid)} {
  JobGetResults $Opts(jid)
  exit 0
}

while {[set nRunning [CollectResults]]} {
  Say "==> Number of unfinished jobs:  $nRunning\n" warn
  if {$nWaitSec < $MinWaitSec} break
  after [expr $nWaitSec*1000]
}

if {$nRunning} {Stop 0 "DONE"}

if {$Opts(run) || [Yes "Results collected.\nPerform final calculations?"]} {
  Say "Starting final run ..." info
  # make_ControlCards 2 0
  # exec $MainExe > FitPDF.log 2> FitPDF.err
  # file rename -force -- $CCorg $CCfname
  run_local 2 0
  Stop 0 "DONE"
}

exit 0
