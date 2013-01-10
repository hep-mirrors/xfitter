#!/bin/sh
# the next line restarts using tclsh \
if [ -z $DISPLAY ]; then exec tclsh "$0" "$@"; else exec wish "$0" -- "$@"; fi

catch {source $tcl_rcFileName}
source [file join [file dirname [info script]] job_farm_fns.tcl]
# source [file join [file dirname [info script]] zarah-new.tcl]

#-------------------------------------------------
#  Run job locally
#-------------------------------------------------

set Version 1.0

set Verbose 1
set _TEST_JS 0

# ----------------------------------------------------------------

# ###############################################

MainStartUp

set Args [GetOptions {
  ind ""
  use 0
} Opts $argv]
# puts $Args
# puts [array get Opts]

set mu [string trim $Opts(ind)]
if {![string is integer $mu]} {Stop 1 "Illegal 'ind' value: $mu"}
run_local $Opts(use) $mu

Stop 0 "DONE"
