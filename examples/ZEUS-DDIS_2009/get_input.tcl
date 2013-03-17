#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@" 

catch {source $tcl_rcFileName}
lappend auto_path [file join [file dirname [info script]] ../../common/TkTcl]
package require wsutils
# source [file join [file dirname [info script]] ../../common/TkTcl/wsutils.tcl]

# #########################################

# set isMSWin [string match "windows" $tcl_platform(platform)]
# if {$::isMSWin && [info exists ::tk_version]} {
  # console show; update
  # wm withdraw .
# }
MainStartUp

set Verbose 0

set Name DIFFRACTION
set srcdir [file join [file dirname [info script]] ../../input_steering]
set fl0 [glob -nocomplain -dir $srcdir -tails *.$Name]
set flist {}
set len [string length ".$Name"]
foreach f $fl0 {
  lappend flist [string range $f 0 end-$len]
}
if {![Yes "Setup name: $Name" "Going to copy\n    [join $flist "\n    "]\nto the current dir." yes "?" $isMSWin]} {exit}
# puts $flist
foreach fin $fl0 f $flist {file copy -force [file join $srcdir $fin] $f}
# Stop 0 [join $flist \n]
exit
