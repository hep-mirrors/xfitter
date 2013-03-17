# Portable Batch System
# uses local file storage
# --- PBS commands used:
# qsub

# qstat
  # JobID Username Queue Jobname SessID  ReqNodes ReqTasks ReqMemory ReqTime Stat Time
  # 813965.h1farmer0.des zohrab   SL512H3G task.dg568  15136     1  --    --    --  R 00:23
  # 26311704.batch.g     plgwues     plgrid-t test03              --    --     --     --  01:00 Q   --
  # SessID becomes a number once the job is running (Stat == 'R')
    # the job identifier assigned by PBS. 
    # the job owner. 
    # The queue in which the job currently resides. 
    # The job name given by the submitter. 
    # The session id (if the job is running). 
    # The number of nodes requested by the job. 
    # The number of cpus or tasks requested by the job. 
    # The amount of memory requested by the job. 
    # Either the cpu time, if specified, or wall time requested by the job, (hh:mm). 
    # The job's current state. 
    # The amount of cpu time or wall time used by the job (hh:mm). 
    
# ---------------------------------------------------------------

set LFSok 1

if {$Farm(User) == ""} {set Farm(User) $tcl_platform(user)}
if {$Farm(Driver) == ""} {set Farm(Driver) "jpbs.sh"}

# ===================================
proc farm_cmd {cmd} {
  # --- return farm command with optional host and user names
  set fc "exec $cmd"
  # if {$::Farm(Host) != ""} {append fc " -h $::Farm(Host)"}
  if {$::Farm(User) != ""} {append fc " -u $::Farm(User)"}
  return $fc
}

# ===================================
proc farm_submit {driver args} {
  # --- Submit job to the farm
  # --- return job identifier needed to retrieve the results
  global Farm Label
  set cmd [farm_cmd qsub]
  if {$Farm(Queue) != ""} {append cmd " -q $Farm(Queue)"}
  if {$Farm(TimeLimit) != ""} {append cmd " -l walltime=$Farm(TimeLimit)"}
  if {$Farm(WorkDir) != ""} {append cmd " -v RunTimeDir=$Farm(WorkDir)"}
  append cmd " -o _job_stdout -e _job_stderr -N $Label"
  # append cmd " -k eo" ; # nie wiem gdzie to zostaje
  if {[catch {eval $cmd \$driver 2>@1} ans ]} {
    return -code error $ans
  }
  # --- parse output
  set ans [string trim $ans]
  if {[string match "Args-check:*" $ans]} {set ans [lindex [split $ans \n] 1]}
  if {[regexp {\A(\d+)\.} $ans m jno]} {
    scan $jno %d jno
    return $jno
  }
  return -code error $ans
}

# ===================================
proc farm_query_running {} {
  global Verbose
  # --- returns list of unfinished jobs
  set Keys {JobID Username Queue Jobname SessID  ReqNodes ReqTasks ReqMemory ReqTime Stat Time}
  set jrun [string trim [eval [farm_cmd qstat]]]
  if {$jrun == ""} {return}
  if {$Verbose > 1} {Say $jrun}
  # set jrun [lrange [split $jrun \n ] 3 end]
  set Running {}
  # if {[string match "No unfinished *" [lindex $jrun 0]]} {return }
  foreach rr [split $jrun \n ] {
    if {$Verbose > 2} {Say "rr='$rr'"}
    set rr [string trim $rr]
    regsub -all {\s+} $rr "\t" rr
    set rdat [split $rr \t]
    if {$Verbose > 2} {Say "rdat=$rdat"}
    if {[llength $rdat] != [llength $Keys]} continue
    foreach k $Keys v $rdat {set $k $v}
    if {![regexp {\A(\d+)\..+} $JobID mm jno]} continue
    # if {[string match {\**} $rr]} {break}
    # scan [lindex $rr 1] %d rid
    lappend Running $jno
  }
  if {$Verbose > 1} {Say $Running}
  return $Running
}

# ===================================
proc farm_getresults {jrec args} {
  # --- get std. output and files specified in args
  # --- jrec = job record = list: nCS iCS jid tmpdir ...
  # --- if args is empty then get all files
  if {$::Verbose > 2} {Say $jrec}
  lassign $jrec nCS iCS jid rdir
  Say "CS index = $iCS"
  # set jid [lindex $jrec 2]
  # set rdir [lindex $jrec 3]
  # --- rdir contains all job results
  set jobdir [pwd]
  cd $rdir
  # --- wait for job output files
  set ntry 3
  while {$ntry} {
    glob -nocomplain _job_std*
    set fc 0
    foreach f {_job_stderr _job_stdout} { if {[file exists $f]} {incr fc} }
    if {$fc == 2} break
    incr ntry -1
    if {!$ntry} break
    after 3000
  }
  if {$fc < 2} {return -code error "_job_stderr and _job_stdout not available (yet)"}
  
  foreach f [glob -nocomplain *.?$jid] {file rename $f $jobdir}
  foreach f [glob -nocomplain _job_*] {file rename $f $jobdir}
  foreach f $args {file rename $f $jobdir}
  cd $jobdir
  # file delete -force $rdir
  return
}

