package provide cl_opts 1.2
package require wsutils

namespace eval opts {

  variable Options
  variable Type
  variable Comment
  variable HelpText
  
  #==============================================
  proc help_text {} {
    variable Options
    variable Type
    variable Comment
    variable HelpText
    if {$HelpText != ""} {set tx "$HelpText\n"}
    append tx "Available options \(name \[default\] <type> comment\):\n"
    append tx "  -?, -help, --help    show this text\n"
    foreach k [lsort [array names Options]] {
      append tx "  -$k \[$Options($k)\] <$Type($k)>   $Comment($k)\n"
    }
    return $tx
  }

  #==============================================
  proc show_options {{emsg ""}} {
    # Stop 0 [help_text]
    if {$emsg != ""} {
      Stop 1 $emsg "Use: [info script] -?\nfor help" $::isMSWin
      # Say "$emsg" err
      # Say "Use: [info script] -?\nfor help"
      # exit 1
    }
    Stop 0 [help_text] "" $::isMSWin
    # if {$emsg != ""} {Stop 1 "$emsg\nUse -? or --help"}
    # if {$emsg != ""} {Stop 1 "$emsg"}
    # if {$emsg != ""} {exit 1}
    # exit
  }

  #==============================================
  proc GetOptsX {optdefs oparn alist {short_bool 0}} {
    variable Options
    variable Type
    variable Comment
    variable HelpText
    
    array set Msg {
      badopt {"Illegal option: -$k"}
      noval {"No value for option -$key"}
      badval {"Illegal value for option -$key: $v\nMust be of type $Type($key)"}
      ambopt {"Ambiguous option: -$k\nmatches: $key"}
    }
    upvar $oparn oarr
    #	Pause $alist
    array unset Options
    array unset Type
    array unset Comment
    array unset oarr
    # --- optdefs = list of lists
    # --- {name defval type [comment]}
    # --- reserved options: help ?
    
    # set Opts {
      # {dash       0            bool }
      # {dashscl    1            double "dash-length scale" }
      # {fscl       0.5          real "font-size scale" }
      # {size       {17cm,26cm}  any }
      # {count       6  int }
      # {font       1  any= "" {Helvetica "Times-Roman"}}
    # }

    # set optdefs [linsert $optdefs 0 {quiet 0 bool}]
    set HelpText ""
    foreach opt $optdefs {
      lassign $opt k v t c
      if {[string eq {<help_text>} $k]} {
        set HelpText $v
      } {
        set Options($k) $v
        set oarr($k) $v
        set Type($k) $t
        set Comment($k) $c
      }
    }
    set nparsed 0
    set nArgs [llength $alist]
    set RealArgs {}
    for {set nparsed 0} {$nparsed < $nArgs} {incr nparsed} {
      set k [lindex $alist $nparsed]
      if {[string eq {--} $k]} {
        # --- stop parsing options
        incr nparsed
        break
      }
      if {![string match {-*} $k]} {
        lappend RealArgs $k
        continue
      }
      set k [string range $k 1 end]
      if {$k == "help" || $k == "-help" || $k == "?"} {
        show_options
        continue
      }
      set key [array names oarr $k]
      if {$key == ""} {set key [array names oarr $k*]}
      switch [llength $key] {
        0 {eval show_options $Msg(badopt)}
        1 {
          if {$short_bool && [string match bool* $Type($key)]} {
            set oarr($key) 1
          } {
            incr nparsed
            if {$nparsed >= $nArgs} {eval show_options $Msg(noval)}
            set v [lindex $alist $nparsed]
            if {![string eq $Type($key) any]} {
              if {![string is $Type($key) $v]} {eval show_options $Msg(badval)}
            }
            if {[string match bool* $Type($key)]} {
              set oarr($key) [string is true $v]
            } else {
              set oarr($key) $v
            }
          }
        }
        default {eval show_options $Msg(ambopt)}
      }
      #	Pause "nparsed = $nparsed"
    }
    return [concat $RealArgs [lrange $alist $nparsed end]]
  }

  # namespace export *
  # namespace export GetOptsX
}
