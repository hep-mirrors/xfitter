catch {source $tcl_rcFileName}
source [file join [file dirname [info script]] job_farm_fns.tcl]

#-------------------------------------------------
#  Run job locally
#-------------------------------------------------

set Version 1.3

set Verbose 1
set _TEST_JS 0

# ----------------------------------------------------------------

# ###############################################

MainStartUp
Initialize 0

# Pause [PathRelTo [pwd] $MainCodeDir 1]

# --- Options. Each line contains: name def_value type [comment]
set options_ {
  {ind  ""  int "CS index"}
  {use  0   int "Use previous results (0,1,2)"}
}
set argv [opts::GetOptsX $options_ Opts $argv ]
# ...................................................

Initialize

set mu [string trim $Opts(ind)]
run_local $Opts(use) $mu

Stop 0 "JayRun finished"
