# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  # \date 2012--2013
  # \copyright Creative Commons license CC-BY-NC 3.0
# _____________________________________________________________

# -------------------------------------------------
#   Routines to run HERAFitter jobs on a farm
# -------------------------------------------------
set Version 1.0

source [file join [file dirname [info script]] utils.tcl]

set CCfname "steering.txt"
set FarmUser ""
set FarmHost ""

source [file join [file dirname [info script]] jay.cfg.tcl]
# source [find_file jay.cfg.tcl]
source_if_exists jay.cfg.tcl
source [file join [file dirname [info script]] $FarmName.tcl]

#======================================
proc get_datfiles {cc_fn} {
  set dfList {}
  foreach tx [split [string trim [ReadFile $cc_fn]] \n ] {
    set tx [string trim $tx]
    if [string match {[\*!]*} $tx] continue
    # -- InputFileNames(1) = 'datafiles/hera/ZEUS_LRG_98-00.dat'
    # -- CorrFileNames(1) = 'datafiles/hera/H1_NormInclJets_HighQ2_99-07___H1_NormInclJets_HighQ2_99-07.corr'
    if {![regexp {FileNames\s*\(\s*\d+\s*\)\s*=\s*(.+)\Z} $tx m dn]} continue
    set dn [string trim $dn]
    lappend dfList [string trim $dn ' ]
  }
  return $dfList
}

#======================================
proc job_sub1 {} {
  # --- Submit job to the farm
  # --- return job identifier needed to retrieve the results
  # ---        or "" on error
  set drv [find_file $::Driver]
  if {[catch {farm_submit $drv $::MainExe input.zip} ans ]} {
    Say "*** ERROR submitting job:\n$ans" err
    return
  }
  return $ans
}

#======================================
proc Suffix {CSind} {
  if {!$CSind} {return "0"}
  set c [expr {$CSind > 0 ? "p" : "m"}]
  set mu [expr abs($CSind)]
  return "[format %03d $mu]$c"
}

#======================================
proc JobSubmit {CSind} {
	global output_dir job_wk_dir Label Verbose
	global CCfname CCorg
  Say "CorSysIndex = $CSind"

  set sv [file join $output_dir minuit.save_0.txt]
  if {![file exists $sv]} {
    Say "*** ERROR: File $sv not found.\nRun the central fit first."
    return 1
  }
  make_ControlCards 1 $CSind
  set flist [get_datfiles $CCfname]
  lappend flist $CCfname ewparam.txt minuit.in.txt $sv
  if {$::QuickInput && [file exists input.zip]} {
    set ans [eval exec zip -u input.zip $flist]
  } {
    file delete input.zip
    set ans [eval exec zip input.zip $flist]
  }
  file rename -force -- $CCorg $CCfname
  if {$Verbose > 1} {Say $ans}
  set jno [job_sub1]

  file mkdir $job_wk_dir
  if {$jno != ""} {
    Say "    $jno submitted."
    set JobList [file join $job_wk_dir $Label.jobs]
  	set ch [open $JobList a]
  	# puts $ch "[format %010d $jno] $logfn"
  	puts $ch "$jno"
  	close $ch
    file delete [file join $output_dir "params_[Suffix $CSind].txt"]
    return 0
  }
  Say "*** submission failed."
  return 1
}

# ===========================
proc TestJobResults {args} {return 0}

#===================================
proc GenJobDir {jid} {
	set dn [clock format [clock seconds] -format %Y-%m-%d]_$jid
	return $dn
}

#===================================
proc CollectResults {} {
  # --- retrieve completed jobs and remove them from the list
  # --- return # of still running jobs
  
	global output_dir job_wk_dir Label
  
  set nSent 0
  set JobList [file join $job_wk_dir $Label.jobs]
  if {[file exists $JobList]} {
    set SentJobs [split [ReadFile $JobList] \n ]
    set nSent [llength $SentJobs]
  }
  if {!$nSent} {
    Say "No jobs submitted or all jobs already retrieved." warn
    return 0
  }
  
  Say "Connecting to farm..."
  set AllJobs [farm_query_all]
	set nAll [llength $AllJobs]
  if {!$nAll} {
    Say "  No jobs on farm!" warn
    return 0
  }
  
	set AllRunningJobs [farm_query_running]
	# set nRun [llength $AllRunningJobs]
  
	set FinishedJobs {}
	set RunningJobs {}
  foreach j $SentJobs {
    if {[lsearch $AllJobs $j] < 0} continue
    if {[lsearch $AllRunningJobs $j] < 0} {lappend FinishedJobs $j} {lappend RunningJobs $j}
  }
  
  if {$FinishedJobs == ""} {Say "  No completed jobs to retrieve." warn}
  foreach j $FinishedJobs {
    set rc [JobGetResults $j]
    if {$rc > 0} {
      Say "  *** ERROR\n  See FitPDF.err and FitPDF.log files.\nYou may also need to look at the job output files packed in 'output.zip'." err
    } elseif {$rc < 0} {
      Say "  will try to get it later."
      lappend RunningJobs $j
    }
  }
  
  # --- update JobList
  if {$RunningJobs == ""} {
    file del $JobList
    return 0
  } else {
    set och [open $JobList w]
    puts $och [join $RunningJobs \n]
    close $och
    return [llength $RunningJobs]
  }
}

#======================================
proc job_get1 {jid} {
  # --- Get job results from the farm
  # --- return 0 on success
  if {[catch {farm_getresults $jid $::CCfname ewparam.txt minuit.in.txt FitPDF.log FitPDF.err output.zip} ans ]} {
    Say "*** ERROR: [string trim $ans]" err
    return 1
  }
  return 0
}

#===================================
proc JobGetResults {jid} {
  Say "Retrieving $jid ..." info
  global output_dir
  file mkdir $output_dir
  set here [pwd]
  cd $output_dir
  set jobdir [GenJobDir $jid]
  file mkdir $jobdir
  cd $jobdir
  set ntry 3
  
  while {[set rc [job_get1 $jid]] && $ntry} {
    after 1000
    incr ntry -1
    Say "  trying once more..."
  }
  cd ..
  set rc [expr -($rc)]
  if {!$rc} {
    set rc [TestJobResults $jobdir/output.zip]
    if {!$rc} {
      Say [exec unzip -o $jobdir/output.zip -x fittedresults.txt pulls.first.txt pulls.last.txt parsout_* 2>@1]
      Say "... done."
    }
  }
  cd $here
  return $rc
}


#======================================
proc make_ControlCards {UsePrev {CorSysNo ""} } {
	global CCfname CCorg
  set CCorg "org_$CCfname"
  if {![file exists $CCorg]} {file copy $CCfname $CCorg}
  set och [open $CCfname w]
  puts $och "&CSOffset"
  if {$CorSysNo != ""} {puts $och "  CorSysIndex = $CorSysNo"}
  puts $och "  UsePrevFit = $UsePrev"
  puts $och "&End"

	puts $och [ReadFile $CCorg]
  close $och
}

# ====================================
proc run_local {useprev {mu ""}} {
  global CCorg CCfname MainExe output_dir
  set msg "Running $MainExe for"
  if {$mu != ""} {append msg " CorSysIndex = $mu,"}
  append msg " UsePrev = $useprev ..."
  Say $msg  info
  file mkdir $output_dir
  make_ControlCards $useprev $mu
  set ofn [file join $output_dir "FitPDF"]
  if {$mu != ""} {append ofn "_[Suffix $mu]"}
  set rc [catch {exec $MainExe > $ofn.log 2> $ofn.err} ans]
  if {$rc} {
      Say "*** ERROR running '$MainExe':\n$ans" err
  } {
    if {[file exists $ofn.err] && ([file size $ofn.err] == 0)} {file delete $ofn.err} {
      Say "*** ERROR running '$MainExe':\nsee '$ofn.err'" err
    }
  }
  file rename -force -- $CCorg $CCfname
}

