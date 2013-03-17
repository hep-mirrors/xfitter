#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

catch {source $tcl_rcFileName}
# source [file join [file dirname [info script]] utils.tcl]
source [file join [file dirname [info script]] job_farm_fns.tcl]

# -------------------------------------------------
#  Configure Jay -- HERAFitter batch utilities
# -------------------------------------------------

set Version 1.3
set _TEST_JS 0
# set Verbose 1

# ###############################################

MainStartUp
Initialize 0

set thisDir [file normalize [file dirname [info script]]]
set tail ".cfg.tcl"
set clist [glob -dir $thisDir -tails *$tail]
foreach f [lsort $clist] {
  # set cname [file root [file root $f]]
  # lappend FarmNames [string range $f 0 end-[string length $tail]]
  append FarmNames "\n  [string range $f 0 end-[string length $tail]]"
  # Say $cname
}

# --- Options. Each line contains: name def_value type [comment]
  # {quiet  0   bool "configure in the current folder only"}
set options_ {
  {<help_text>  "Prepare for running Offset mode fits locally and optionally on a farm.\nUsage: [file rootname [file tail [info script]]] \[options\] \[farm_name\]
Available farm_names:$FarmNames"}
  {local  0   bool "configure in the current folder only"}
}
set argv [opts::GetOptsX [subst $options_] Opts $argv 1]
set argc [llength $argv]
# ...................................................



# Pause "$thisDir\n[pwd]"
if {$Opts(local) && [string eq $thisDir [pwd]]} {Stop 1 "JayConfig -local\ncannot be run in the installation directory" "" $isMSWin}
set EnvOK [info exists env(HERAFITTER_SYS)]
if $EnvOK {set MainCodeDir $env(HERAFITTER_SYS)} {
  set MainCodeDir [file normalize [file join $thisDir ../..]]
  # Say "It is recommended to set env. var.\n   HERAFITTER_SYS=installation folder of the HERA Fitter\n" warn
  # Say "I assume it is '$MainCodeDir'"
  # if {![cYes "I assume that the HERA Fitter installation directory is:\n$MainCodeDir\n"]} {exit}
}

set UserHFdir "~/.HERAFitter"
set UserHFdir [file normalize $UserHFdir]
file mkdir $UserHFdir
if {$Opts(local)} {set dest_dir [pwd]} {set dest_dir $UserHFdir}

set sty err

# Say [PathRelTo [pwd] $MainCodeDir 1]
if {[file pathtype [PathRelTo $dest_dir $MainCodeDir 1]] != "absolute"} {
  # if {![Yes "Current path '[pwd]'\nis within the HERAFitter installation tree." "Are you sure to continue?" no "?" $isMSWin]} {exit}
  if {![Yes "You have chosen to save current configuration in the folder\n$dest_dir\nwhich is within the HERAFitter installation tree." "Are you sure to continue?" no "?" $isMSWin]} {exit}
}
set trg_cfg [file join $dest_dir "jay.cfg.tcl"]
# Pause $trg_cfg
if [file exists $trg_cfg] {
  if {![Yes "Configuration file\n$trg_cfg\nalready exists." "Overwrite?" yes "?" $isMSWin]} {exit}
}
  
if {$argc} {
  set cname [lindex $argv 0]
  set src_cfg [file join $thisDir $cname$tail]
  
  if {![file exists $src_cfg]} {
    Stop 1 "Unknown farm name '$cname'" "Available farm names:$FarmNames" $isMSWin
  }
  
  file copy -force $src_cfg $trg_cfg
  if {!$EnvOK} {
    # set rp [PathRelTo $MainCodeDir [pwd]]
    # WriteFile "set MainCodeDir $rp\n\n[ReadFile $trg_cfg utf-8]" $trg_cfg utf-8 lf
    WriteFile "set MainCodeDir $MainCodeDir\n\n[ReadFile $trg_cfg utf-8]" $trg_cfg utf-8 lf
  }
  
  # Stop 0 "Configuration saved to $trg_cfg\nYou can edit this file to make further adjustments." "" $isMSWin
  # set msg "Configuration saved to\n$trg_cfg"
  set msg ""
  set msgd "You can edit this file to make further adjustments."
} {
  file delete $trg_cfg
  if {!$EnvOK} {
    WriteFile "set MainCodeDir $MainCodeDir\n" $trg_cfg utf-8 lf
  }
  
  # if {$Farm(Type) == ""} {
    # set msg "No farm configured"
  # } {
    # set msg "Your farm: $Farm(Type) $Farm(Name)"
  # }
  set msg "No farm selected"
  set msgd "Run:\n  JayConfig <farm_name>    to select a new farm\n  JayConfig -?             for help"
}

set rList {JayRun JayGet JaySub}
# --- make links in the install. dir.
set here [pwd]
cd $thisDir
# Pause [pwd]
foreach f $rList {
  file delete $f
  file link $f ".jay_go"
}
cd $here
  
# --- make scripts in $HOME/bin
set dir "~/bin"
file mkdir $dir
foreach f $rList {
  set p [file join $dir $f]
  if $isMSWin {
    append p ".bat"
    file delete $p
    WriteFile "@start [file join $::thisDir $f.tcl] %*" $p
    # file attributes $p -permissions "+x"
  } {
    file delete $p
    WriteFile "#!/bin/sh\nexec [file join $::thisDir $f] \$@" $p
    file attributes $p -permissions "+x"
  }
}

Stop 0 "Configuration saved to\n$trg_cfg\n$msg" $msgd $isMSWin
