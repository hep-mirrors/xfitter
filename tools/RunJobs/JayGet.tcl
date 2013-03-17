catch {source $tcl_rcFileName}
source [file join [file dirname [info script]] job_farm_fns.tcl]

#-------------------------------------------------
#  Retrieve HERAFitter jobs from a farm
#-------------------------------------------------

set Version 1.2
set _TEST_JS 0
# set Verbose 1

# --- redefine proc TestJobResults
# ===========================
proc TestJobResults {args} {
  set outzip [lindex $args 0]
  if {[catch {exec unzip -l $outzip params_*.txt} ans]} {
    return 1
  }
  return 0
}

# ###############################################

MainStartUp
Initialize 0

# --- Options. Each line contains: name def_value type [comment]
set options_ {
  {lbl  $DefaultLabel   any "Label"}
  {peek 0 int "check for results every 'peek' minutes"}
  {jid 0 int "job id to retrieve; 0 means ALL"}
  {verbose 1 int "Verbosity level"}
  {run 0 bool "perform final run if possible"}
}
set options_ [subst $options_]
set argv [opts::GetOptsX $options_ Opts $argv 1]
# ...................................................

Initialize

set Label $Opts(lbl)
if {$Label == ""} {set Label $DefaultLabel}
set Verbose $Opts(verbose)

ShowHeader

set nWait $Opts(peek)
if {$nWait > 0 && $nWait < $MinWait} {
  Say "Peek interval increased to $MinWait minutes." warn
  set nWait $MinWait
}

if {$Opts(jid)} {
  JobGetResults $Opts(jid)
  exit 0
}

while {[set nRunning [CollectResults]]} {
  Say "Number of unfinished jobs:  $nRunning"
  if {$nWait < $MinWait} break
  Say "->->->  will check again in $nWait minutes..."
  after [expr int($nWait*60000)]
}

# if {$nRunning} {Stop 0 "DONE"}
if {$nRunning} {exit}
if {[check_results]} {
  Say "Still not enough data to finalize Offset calculation" info
  exit 0
}

if {$Opts(run) || [Yes "Results collected." "Perform final calculations?"]} {
  Say "Starting final run ..." info
  run_local 2 0
  # Stop 0 "DONE"
}

exit 0
