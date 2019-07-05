#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@" 

catch {source $tcl_rcFileName}
# package require wsutils

if {[info commands lassign] == ""} {
  proc lassign {Vals args} {
    foreach k $args v $Vals {upvar $k x; set x $v}
  }
}

#===========================
proc ReadTextFile {fn {Enc ""} {SplitByLines 0}} {
  # utf-8 identity
  #if {![file exists $fn]} {Warn "No file '$fn'"}
  if {![file exists $fn]} {return}
  set ch [open $fn]
  if {$Enc != ""} {fconfigure $ch -encoding $Enc}
  set xtx [string trim [read $ch]]
  close $ch
  if $SplitByLines {set xtx [split $xtx \n]}
  return $xtx
}

# ============================
proc Stop {rc msg} {
  if $rc {
    set sty error
    set tit Error
  } {
    set sty info
    set tit Done
  }
  puts $msg
  if {$::isMSWin && [info exists ::tk_version]} {
    tk_messageBox -icon $sty -title $tit -message $msg -type ok
  }
  exit $rc
}

# =============================
proc IgnoreSec {tit} {
  foreach ignpat [list "Stat." "Syst." "Covariance"] {
    if [string match "$ignpat *" $tit ] {return 1}
  }
  return 0
}

# =============================
proc compare {fn0 fn} {
  global Verbose
  # set fn0 [file normalize $fn0]
  set tx0 [ReadTextFile $fn0 utf-8 1]
  if {$tx0 == ""} {Stop 1 "Missing or empty file\n'$fn0'"}
  set tx [ReadTextFile $fn utf-8 1]
  if {$tx == ""} {Stop 1 "Missing or empty file\n'$fn'"}

  puts "\nComparing numerical values in"
  puts $fn
  puts "to"
  puts $fn0
  puts ""

  set fmt "%.3g"
  set asym 0
  set dmin 0.0
  set dmax 0.0
  set N 0
  set N0 0
  set Sec 0
  set nBad 0
  set DeltaMax 1.5e-3
  set Relative 1
  lappend tx0 ""
  foreach t0 $tx0 t1 $tx {
    incr N
    set t0 [string trim $t0]
    if {![string is double [lindex $t0 0]]} {
      # --- start of Section
      set Ignore [IgnoreSec $t0]
      incr Sec
      set dmin 0.0
      set dmax 0.0
      if $Ignore continue
      puts $t0
      continue
    }
    if {$t0 == ""} {
      # --- end of Section
      set BAD [expr {$N0 && ($dmax >= $DeltaMax)}]
      if {$Verbose || $BAD} {
        if {$N0} {
          puts -nonewline "Max. relative difference encountered at line $N0, item $k0:  "
          if {$asym} {puts [format "$fmt,  $fmt" $dmin $dmax]} {puts [format $fmt $dmax]}
          if {$dmax < $DeltaMax} {puts "--- OK\n"} {incr nBad; puts "*** TOO BAD\n"}
        } {puts "--- No differences found\n"}
      }
      set DeltaMax 0.01
      set Relative 0
      continue
    }
    if $Ignore continue
    set t1 [string trim $t1]
    set k 0
    foreach v0 $t0 v1 $t1 {
      incr k
      if $Relative {set d [expr double($v1)/$v0 -1]} {set d [expr double($v1) - $v0]}
      if {!$asym} {set d [expr abs($d)]}
      if {$d > $dmax} {set dmax $d; set N0 $N; set k0 $k}
      if {$asym && $d < $dmin} {set dmin $d}
      # if {abs($d) > 2e-2} {Pause "$N:$k\n$v1\n$v0\n$d"}
    }
  }

  if {$nBad} {set msg "Some bad results found"} {set msg "All OK"}
  return $msg
}

# =============================
proc FixFname {fn} {
  global resName
  # --- [file isdirectory $fn] is true also when $fn is a link to a directory
  if {![file exists $fn]} {Stop 1 "'$fn'\ndoes not exist!"}
  if {![file isdirectory $fn]} {return $fn}
  set fn1 [file join $fn $resName]
  if {[file exists $fn1]} {return $fn1}
  set fn1 [file join $fn output $resName]
  if {[file exists $fn1]} {return $fn1}
  Stop 1 "Cannot locate $resName\nin $fn"
}

# #########################################

set isMSWin [string match "windows" $tcl_platform(platform)]
if {$::isMSWin && [info exists ::tk_version]} {
  console show; update
  wm withdraw .
}

set Verbose 0

lassign $argv fn0 fn
# puts "'$fn0'   '$fn'" ; exit
set resName "offset.save.txt"

if {$fn == ""} {set fn output}
# --- fn is relative to the current dir
# set fn [file normalize $fn]
# if [file isdirectory $fn] {set fn [file join $fn $resName]}
set fn [FixFname $fn]

if {$fn0 == ""} {set fn0 "output"}
# --- fn0 is relative to this script dir
set fn0 [file join [file dirname [info script]] $fn0]
# if [file isdirectory $fn0] {set fn0 [file join $fn0 $resName]}
set fn0 [FixFname $fn0]

Stop 0 [compare $fn0 $fn]
