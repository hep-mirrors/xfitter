set isGUI [info exists ::tk_version]

#============================
proc find_file {fnam {pathList ""}} {
  if {$pathList == ""} {set pathList [list "." [file dirname [info script]]]}
  foreach path $pathList {
    set fn [file join $path $fnam]
    if {[file exists $fn]} {return $fn}
  }
  return -code error $fnam
}

#============================
proc source_if_exists {fn} {
if {[file exists $fn]} {uplevel "source $fn"}
}

#============================
proc ReadFile {fnam} {
  set ch [open $fnam]
  set tx [string trim [read $ch]]
  close $ch
  #return [split $tx \n]
  return $tx
}

#============================
proc Say {txt {style normal}} {
  puts $txt; update
}

#=========================================================
proc ConsAsk {msg qstr} {
  set qstr [string tolower $qstr]
  set ans "A"
  # Say $msg
  while {[string first $ans $qstr] < 0} {
    Say "$msg \[$qstr\]"
    if {[gets stdin ans] < 0} {return}
    if {$ans == ""} {return [string index $qstr 0]}
    set ans [string tolower [string index $ans 0]]
  }
  return $ans
}

#=========================================================
proc Yes {msg {def yes} {tit "?"}} {
  if [info exists ::tk_version] {
    set p [focus]
    if {$p == ""} {set p .}
    return [string eq yes [tk_messageBox -parent $p -icon question -title $tit -message $msg -type yesno -default $def]]
  } {
    return [expr {[ConsAsk "$msg" "YN"] == "y"}]
  }
}
 
#=========================================================
proc Pause {msg {tit "Pause..."} {CancProc {exit}}} {
  if [info exists ::tk_version] {
    set p [focus]
    if {$p == ""} {set p .}
    if {[tk_messageBox -parent $p -title $tit -message $msg -type okcancel] != "cancel"} {return 0}
  } {
    if {[ConsAsk "$msg\nContinue?" "YN"] == "y"} {return 0}
  }
  return [uplevel 1 $CancProc]
}

#=========================================================
proc Stop {rc msg {tit "Stop!"}} {
  if $rc {set sty error} {set sty info}
  if [info exists ::tk_version] {
    set p [focus]
    if {$p == ""} {set p .}
    tk_messageBox -parent $p -icon $sty -title $tit -message $msg -type ok
  } {
    Say $msg $sty
  }
  exit $rc
}

#==============================================
proc GetOptions {optdefs oparn alist {esch -}} {
  upvar $oparn oarr
  # Pause $alist
  array unset oarr
  # array set oarr $optdefs
  foreach {k v} $optdefs {
    set isBool($k) 0
    if {[string match {\?*} $k]} {
      set k [string range $k 1 end]
      set isBool($k) 1
      set v 0
    }
    set oarr($k) $v
  }
  # Pause [array get oarr]
  set Keys [array names oarr]
  set nparsed 0
  set narg [llength $alist]
  # for {set i 0} {$i < $narg} {incr i} {}
  while {$nparsed < $narg} {
    set k [lindex $alist $nparsed]
    #  Pause "$k = $v"
    if {[string eq "$esch$esch" $k]} {incr nparsed; break}
    if {![string match "$esch*" $k]} break
    set k [string range $k 1 end]
    set key [array names oarr $k*]
    if {[llength $key] != 1} {
      Stop 1 "Unresolved option: $k\nMust be one of: $Keys"
    }
    incr nparsed
    if {$isBool($key)} {
      set oarr($key) 1
    } {
      set oarr($key) [lindex $alist $nparsed]
      incr nparsed
    }
  }
  # Pause "nparsed = $nparsed"
  return [lrange $alist $nparsed end]
}

#==============================================
proc MainStartUp {} { uplevel {
  set isMSWin [string match "windows" $tcl_platform(platform)]
  set SystemConsole [catch	{console show; update}] 
  if $::isGUI {wm withdraw . }
}}

