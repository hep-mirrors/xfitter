#!/bin/sh
# the next line restarts using tclsh or wish \
if [ -z $DISPLAY ]; then exec tclsh "$0" "$@"; else exec wish "$0" -- "$@"; fi

catch {source $tcl_rcFileName}
source [file join [file dirname [info script]] job_farm_fns.tcl]
# source [file join [file dirname [info script]] cfg.tcl]
# source [file join [file dirname [info script]] $FarmName.tcl]

#-------------------------------------------------
#  Submit jobs for Offset calculations
#-------------------------------------------------

set Version 1.0
set _TEST_JS 0

# ###############################################

set Args [GetOptions {
  lbl "Last"
  i ""
  x ""
  verbose 1
} Opts $argv]
# puts $Args
# puts [array get Opts]

set Label $Opts(lbl)
set Verbose $Opts(verbose)

MainStartUp

Say "|-----  Using farm $FarmName  -----|" info

set JobList [file join $job_wk_dir $Label.jobs]
if [file exists $JobList] {
  Stop 1 "There are some not retrieved jobs!\nRemove '$JobList' to ignore previous submissions
  or submit with a new label (use option -lbl <Label>)"
}

set res [file join $output_dir Results_0.txt]
if {![file exists $res]} {
  Stop 1 "*** ERROR: Run the central fit first."
}
if {![regexp -lineanchor {^ Systematic shifts\s+(\d+)$} [ReadFile $res] m nCorSys]} {
  Stop 1 "Cannot read nCorSys from $res"
}
Say "nCorSys = $nCorSys"
set muList {}
for {set mu -$nCorSys} {$mu <= $nCorSys} {incr mu} {
  if {$mu == 0} continue
  if {[lsearch $Opts(x) $mu] >= 0} continue
  if {($Opts(i) != "") && ([lsearch $Opts(i) $mu] < 0)} continue
  lappend muList $mu
}
if {$Verbose > 1} {Say "Sending jobs for CorSysIndex = $muList"}
file delete input.zip
set nts [llength $muList]
foreach mu $muList {
  if {[JobSubmit $mu]} break
  incr nts -1
}


Stop 0 "DONE"

