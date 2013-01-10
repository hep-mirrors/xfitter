# --- NQS commands used:
# jobsub
# jobq
# jobls
# jobget
# ---------------------------------------------------------------

# ===================================
proc farm_cmd {cmd} {
  # --- return farm command with optional host and user names
  set fc "exec $cmd"
  if {$::FarmHost != ""} {append fc " -h $::FarmHost"}
  if {$::FarmUser != ""} {append fc " -u $::FarmUser"}
  return $fc
}

# ===================================
proc farm_submit {driver args} {
  # --- Submit job to the farm
  # --- return job identifier needed to retrieve the results
  if {[catch {eval [farm_cmd jobsub] -q M \$driver -f $args 2>@1} ans ]} {
    return -code error $ans
  }
  # --- parse output
  if {[regexp {Job\s+(\d+)\s+submitted.\s*} $ans m jno]} {
    scan $jno %d jno
    return $jno
  }
  return -code error $ans
}

# ===================================
proc farm_query_all {} {
  set jl [string trim [eval [farm_cmd jobls]]]
  set jl [split $jl \n]
  set jall {}
  foreach rr $jl {
    scan [lindex $rr 0] %d jid
    lappend jall $jid
  }
  return $jall
}

# ===================================
proc farm_query_running {} {
  set jrun [string trim [eval [farm_cmd jobq]]]
  set jrun [lrange [split $jrun \n ] 3 end]
  set Running {}
  if {[string match "No unfinished *" [lindex $jrun 0]]} {return }
  foreach rr $jrun {
    if {[string match {\**} $rr]} {break}
    scan [lindex $rr 1] %d rid
    lappend Running $rid
  }
  return $Running
}

# ===================================
proc farm_getresults {jobid args} {
  # --- get std. output and files specified in args
  # --- if args is empty then get all files
  set cmd [farm_cmd jobget]
  if {[catch {eval $cmd -s $jobid > stdout} ans]} {
    return -code error $ans
  }
  append cmd " -n $jobid"
  if {$args != ""} {append cmd " -f $args"}
  if {[catch {eval $cmd 2>@1} ans ]} {
    return -code error $ans
  }
  return $ans
}

