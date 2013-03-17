catch {source $tcl_rcFileName}
source [file join [file dirname [info script]] job_farm_fns.tcl]

#-------------------------------------------------
#  Submit jobs for Offset calculations
#-------------------------------------------------

set Version 1.2
set _TEST_JS 0

# ###############################################

MainStartUp
Initialize 0

# --- Options. Each line contains: name def_value type [comment]
set options_ {
  {lbl  $DefaultLabel   any "Label"}
  {i "" any "include indices e.g. -i 1,3:6; default is ALL except 0"}
  {x "" any "exclude indices e.g. -x -5:-1"}
  {verbose 1 int "Verbosity level"}
  {quiet 0 bool "No confirmation of selected indices"}
}
set options_ [subst $options_]
set argv [opts::GetOptsX $options_ Opts $argv 1]
# ...................................................

Initialize

set Label $Opts(lbl)
if {$Label == ""} {set Label $DefaultLabel}
set Verbose $Opts(verbose)

ShowHeader

set JobList [file join $job_wk_dir $Label.jobs]
if [file exists $JobList] {
  Stop 1 "There are some not retrieved jobs!" "Remove '$JobList' to ignore previous submissions or submit with a new label (use option -lbl <Label>)" 0
}

if {$Opts(i) == 0} {
  set nCorSys 0; # required for the job list
  set muList [list 0 ]
} {
  if {[GetNCS]} {Stop 1 "Bye."}
  Say "nCorSys = $nCorSys"
  
  set muList $Opts(i)
  if {$muList == ""} {set muList ":"}
  if {[catch {set muList [ExpandNumList $muList $nCorSys -$nCorSys]} ans]} {Stop 1 $ans}
  set muList [ListRemove $muList 0]
  if {[catch {set xList [ExpandNumList $Opts(x) $nCorSys -$nCorSys]} ans]} {Stop 1 $ans}
  set muList [ListRemove $muList $xList]
}

# Pause [CompressNumList $muList]
if {![Yes "Ready to send jobs" "for CorSysIndex = [CompressNumList $muList]" yes "?"  $isMSWin]} {exit}
# if {$Verbose > 1} {Say "Sending jobs for CorSysIndex = $muList"}
# file delete input.zip
set nts [llength $muList]
set nTot $nts
foreach mu $muList {
  if {[JobSubmit $mu]} break
  incr nts -1
}

if {$nts} {Stop 1 "Bye!"} {Stop 0 "$nTot jobs submitted"}

