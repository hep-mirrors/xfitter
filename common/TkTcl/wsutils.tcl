package provide wsutils 2.1.1

# set isGUI [info exists ::tk_version]

if {[info commands lassign] == ""} {
  proc lassign {Vals args} {
    foreach k $args v $Vals {upvar $k x; set x $v}
  }
}

#============================
proc link_here {fn} {
  # MainCodeDir
  # set fn datafiles
  if {[file exists $fn]} return
  set f [file join $::MainCodeDir $fn]
  # Say "Data = $f"
  if {![file exists $f]} {return -code error "Cannot locate $fn"}
  file link $fn $f
}

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
proc ReadFile {fnam {enc ""}} {
  # utf-8 identity
	# if {![file exists $fnam]} {return}
  set ch [open $fnam]
	if {$enc != ""} {fconfigure $ch -encoding $enc}
  set tx [string trim [read $ch]]
  close $ch
  #return [split $tx \n]
  return $tx
}

# ================================================
proc WriteFile {tx fn {enc ""} {eol ""}} {
	# eol: "lf" for UNIX, "crlf" for DOS
  set ch [open $fn w]
  if {$enc != ""} {fconfigure $ch -encoding $enc}
  if {$eol != ""} {fconfigure $ch -translation $eol}
  puts $ch "$tx"
  close $ch
}

#============================
proc Say_n {txt {style normal}} {
  switch -glob $style {
    err* {set pre "*** ERROR: "}
    warn* {set pre "*** WARNING: "}
    info* {set pre "--> "}
    default {set pre ""}
  }
  puts -nonewline "$pre$txt"; flush stdout; update
}

#============================
proc Say {txt {style normal}} {Say_n "$txt\n" $style}

#=========================================================
proc ConsAsk {msg qstr} {
  set qstr [string tolower $qstr]
  set ans "A"
  # Say $msg
  while {[string first $ans $qstr] < 0} {
    Say_n "$msg \[$qstr\]: "
    if {[gets stdin ans] < 0} {return}
    if {$ans == ""} {return [string index $qstr 0]}
    set ans [string tolower [string index $ans 0]]
  }
  return $ans
}

#=========================================================
proc Yes {msg {det ""} {def yes} {tit "?"} {w 1}} {
  if {$w && [info exists ::tk_version]} {
    set p [focus]
    if {$p == ""} {set p .}
    if {[catch {set ans [tk_messageBox -parent $p -icon question -title $tit -message $msg -detail $det -type yesno -default $def ]}]} {
      set ans [tk_messageBox -parent $p -icon question -title $tit -message "$msg\n$det" -type yesno -default $def ]
    }
    update
    return [string eq yes $ans]
  } {
    if {$det != ""} {append msg "\n$det"}
    return [expr {[ConsAsk "$msg" "YN"] == "y"}]
  }
}
 
#=========================================================
# proc cYes {msg} {
  # return [expr {[ConsAsk "$msg" "YN"] == "y"}]
# }
 
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
proc Stop {rc msg {det ""} {w 1}} {
  if $rc {
    set sty error
    set tit Error
  } {
    set sty info
    set tit Done
  }
  
  if {$w && [info exists ::tk_version]} {
    Say $msg $sty
    if {$det != ""} {Say $det}
    set p [focus]
    if {$p == ""} {set p .}
    if {[catch {tk_messageBox -parent $p -icon $sty -title $tit -message $msg -detail $det -type ok}]} {
      tk_messageBox -parent $p -icon $sty -title $tit -message "$msg\n$det" -type ok
    }
  } {
    # if {$det != ""} {append msg "\n$det"}
    Say $msg $sty
    Say $det
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
      Stop 1 "Unresolved option: $k" "Must be one of: $Keys"
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

#=================================================
proc TempFileName {{subdir ""}} {
  global env
  if [info exists env(TEMP)] {
    set tdir $env(TEMP)
  } elseif [info exists env(TMP)] {
    set tdir $env(TMP)
  } else {
    set tdir [pwd]
  }
  #puts <$tdir>
  set i 0
  while {[file exist [set tn [file join $tdir $subdir _tmp_$i]]]} {incr i}
  return $tn
}

#==============================================
proc MainStartUp {{isGUI 0}} { 
  uplevel #0 {
    set isMSWin [string match "windows" $tcl_platform(platform)]
    set SystemConsole [catch	{console show; update}] 
  }
  if {!$isGUI && [info exists ::tk_version]} {wm withdraw . }
}

# =============================================
proc ExpandNumList {lst {n_max 16} {n_min 1}} {
  set L {}
  foreach n [split $lst ,] {
    set rng [split $n :]
    # puts $rng
    if {[llength $rng] == 1} {
      set n0 $n
      set n1 $n
    } {
      lassign $rng n0 n1
      if {$n0 == ""} {set n0 $n_min}
      if {$n1 == ""} {set n1 $n_max}
      # puts "$n0 - $n1"
      if {$n0 > $n1} {return -code error "Illegal range '$n'"}
    }
    if {$n0 < $n_min || $n1 > $n_max} {return -code error "'$n' is outside \[$n_min,$n_max\]"}
    for {set i $n0} {$i <= $n1} {incr i} {lappend L $i}
  }
  return [lsort -integer -unique $L]
}

# =============================================
proc ListRemove {lst rlst} {
  set L {}
  foreach a $lst {
    if {[lsearch $rlst $a] < 0} {lappend L $a}
  }
  return $L
}

# =============================================
proc CompressNumList {lst} {
  if {![llength $lst]} return
  set lst [lsort -integer -unique $lst]
  set n0 [lindex $lst 0]
  set p $n0
  foreach n $lst {
    if {$n > $p+1} {
      # if {$p > $n0} {lappend L $n0:$p} {lappend L $n0}
      if {$p == $n0} {lappend L $n0} elseif {$p == $n0+1} {lappend L $n0,$p} else {lappend L $n0:$p}
      set n0 $n
    }
    set p $n
  }
  if {$p > $n0} {lappend L $n0:$p} {lappend L $n0}
  return [join $L ,]
}

# ============================================
proc PathRelTo {path base {strict 0}} {
  # --- path and base must be absolute; if not, [pwd] is prepended
  # --- 'file normalize' makes abs. path, converts \\ to / and removes trailing slashes:
  # --- file normalize m:/ == M:/
  # --- file normalize / == <CUR_DRIVE>:/
  # puts "--- '$path' rel to '$base'"
  if {$base == ""} {return -code error "empty base"}
  if {$path == ""} {return -code error "empty path"}
  set path [file normalize $path]
  set base [file normalize $base]
  if {[string equal $path $base]} {return "."}
  if {$strict} {
    if {![string match "*/" $base]} {append base "/"}
    if {![string match $base* $path]} {return $path}
    return [string range $path [string length $base] end]
  }
  set Lpath [file split $path]
  set Lbase [file split $base]
  set nc 0
  foreach p $Lpath b $Lbase {
    # puts "  $p $b"; update
    if {![string equal $p $b]} break
    incr nc
  }
  if {!$nc} {return $path}
  set rpath [string repeat "../" [expr [llength $Lbase] - $nc]]
  # Pause "rpath = $rpath"
  return [eval file join \$rpath [lrange $Lpath $nc end]]
}
