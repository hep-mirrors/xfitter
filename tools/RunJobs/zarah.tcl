# --- ZARAH commands used:
# zarah-jobsub
# zarah-jobq
# zarah-jobget
# ---------------------------------------------------------------

# ===================================
proc farm_submit {driver args} {
  if {[catch {eval exec zarah-jobsub -q M -s \$driver $args 2>@1} ans ]} {
    return -code error $ans
  }
  # --- parse output
  if {[regexp {Job\s+(\d+)\s+submitted\s*} $ans m jno]} {
    scan $jno %d jno
    return $jno
  }
  return -code error $ans
}

# ===================================
proc farm_query {what} {
# --- zarah-jobq ALWAYS returns status 1
# --- work around this bug:
  set tmp "_ooo_[clock seconds].tmp"
  catch {exec zarah-jobq -$what > $tmp}
  after 500
  set ch [open $tmp]
  set ans [string trim [read $ch]]
  close $ch
  file del $tmp

  set qList {}
  foreach qi [lrange [split [string trim $ans] \n ] 2 end] {
    scan [lindex $qi 0] %d j
    lappend qList $j
  }
  return $qList
}

# ===================================
proc farm_query_all {} {return [farm_query a]}

# ===================================
proc farm_query_running {} {return [farm_query r]}

# ===================================
proc farm_getresults {jobid args} {
  # --- get std. output and files specified in args
  # --- if args is empty then get all files
  set cmd "exec zarah-jobget --yes -n -j \$jobid"
  if {$args != ""} {append cmd " stdout $args"} {append cmd " -a"}
  if {[catch {eval $cmd 2>@1} ans ]} {
    return -code error $ans
  }
  return $ans
}

