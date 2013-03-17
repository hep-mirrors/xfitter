# --- NQS commands used:
# jobsub
# jobq
# jobls
# jobget
# ---------------------------------------------------------------

# if {$Farm(User) == ""} {set Farm(User) $tcl_platform(user)}
if {$Farm(Driver) == ""} {set Farm(Driver) "jnqs.sh"}

# ===================================
proc farm_cmd {cmd} {
  # --- return farm command with optional host and user names
  set fc "exec $cmd"
  if {$::Farm(Host) != ""} {append fc " -h $::Farm(Host)"}
  if {$::Farm(User) != ""} {append fc " -u $::Farm(User)"}
  return $fc
}

# ===================================
proc farm_submit {driver args} {
  # --- Submit job to the farm
  # --- return job identifier needed to retrieve the results
  global Farm Label
  set cmd [farm_cmd jobsub]
  if {$Farm(Queue) != ""} {append cmd " -q $Farm(Queue)"}
  append cmd " -j $Label"
  if {[catch {eval $cmd \$driver -f $args 2>@1} ans ]} {
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
proc farm_getresults {jrec args} {
  # --- get std. output and files specified in args
  # --- jrec = job record = list: nCS iCS jid tmpdir ...
  # --- if args is empty then get all files
  set jid [lindex $jrec 2]
  set cmd [farm_cmd jobget]
  if {[catch {eval $cmd -s $jid > stdout} ans]} {
    return -code error $ans
  }
  append cmd " -n $jid"
  if {$args != ""} {append cmd " -f $args"}
  if {[catch {eval $cmd 2>@1} ans ]} {
    return -code error $ans
  }
  return $ans
}

