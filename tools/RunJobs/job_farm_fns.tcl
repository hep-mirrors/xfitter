# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  # \date 2012--2013
  # \copyright Creative Commons license CC-BY-NC 3.0
# _____________________________________________________________

# -------------------------------------------------
#   Routines to run HERAFitter jobs on a farm
# -------------------------------------------------
set Version 1.3

# source [file join [file dirname [info script]] utils.tcl]
lappend auto_path [file join [file dirname [info script]] ../../common/TkTcl]
package require cl_opts

# ----------------
#     Defaults
# ----------------

# ===============================
proc Initialize {{lvl 2}} {
  set ::InitMode $lvl
  uplevel #0 {
    foreach key {Type Name User Host Queue TimeLimit WorkDir Driver} {set Farm($key) ""}
    set QuickInput 0 ; # --- create new input.zip for each job, otherwise update only

    source [file join [file dirname [info script]] jay_def_cfg.tcl]
    # Say "MainCodeDir = $MainCodeDir"
    set MainExe [file join $MainCodeDir $MainExeRel]

    set UserHFdir "~/.HERAFitter"
    source_if_exists [file join $UserHFdir jay.cfg.tcl]
    source_if_exists jay.cfg.tcl

    if {$InitMode} {
      if {[file pathtype [PathRelTo [pwd] $MainCodeDir 1]] != "absolute"} {
        if {![Yes "Current working directory\nis within the HERAFitter installation tree." "Do you really want to run here?" no "?" $isMSWin]} {exit}
      }
      link_here datafiles
      if {$InitMode > 1} {
        if {$Farm(Type) == ""} {Stop 1 "Farm type not set.\nRun JayConfig"}
        source [file join [file dirname [info script]] $Farm(Type).tcl]

        foreach k [array names Farm] {
          set Farm($k) [subst [regsub -all {<(.+?)>} $Farm($k) "$\{\\1\}"]]
        }
        set Farm(WorkDir) [string trimright $Farm(WorkDir) /]
      }
    }
  }
}

# ===============================
proc ShowHeader {} {
  global Farm
  Say "------------------------------------------------------------"
  Say "Farm  = $Farm(Type) $Farm(Name)" info
  Say "Queue = $Farm(Queue)" info
  Say "------------------------------------------------------------\n"
}

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
proc Suffix {CSind} {
  if {!$CSind} {return "0"}
  set c [expr {$CSind > 0 ? "p" : "m"}]
  set mu [expr abs($CSind)]
  return "[format %03d $mu]$c"
}

#======================================
# --- fills global nCorSys
# --- returns RC
proc GetNCS {} {
  global nCorSys
  set res [file join $::output_dir Results_0.txt]
  if {![file exists $res]} {
    Say "File '$res' not found.\nRun the central fit first." err
    return 1
  }
  if {![regexp -lineanchor {^ Systematic shifts\s+(\d+)$} [ReadFile $res] m nCorSys]} {
    Say "Cannot read nCorSys from '$res'" errb
    return 1
  }
  return 0
}

#======================================
proc job_sub1 {} {
  # --- Submit job to the farm
  # --- return job identifier needed to retrieve the results
  # ---        or "" on error
  set drv [find_file $::Farm(Driver)]
  if {[catch {farm_submit $drv $::MainExe input.zip} ans ]} {
    # Say "Failed to submit job:\n$ans" err
    Say "Submission failed:\n$ans" err
    return
  }
  return $ans
}

#======================================
proc pack_input {flist {optflist ""} } {
  set rc 0
  foreach f $flist {
    if {![file exists $f]} {
      Say "File $f not found." errb
      incr rc
    }
  }
  if {$rc} {return 1}
  
  foreach f $optflist {
    if {[file exists $f]} {lappend flist $f}
  }
  
  if {$::QuickInput && [file exists input.zip]} {
    # set ans [eval exec zip -u input.zip $flist]
    set rc [catch {eval exec zip -u input.zip $flist} ans]
  } {
    file delete input.zip
    # set ans [eval exec zip input.zip $flist]
    set rc [catch {eval exec zip input.zip $flist} ans]
  }
  if {$rc} {Say $ans errb} {  if {$::Verbose > 1} {Say $ans} }
  return $rc
}

#======================================
proc JobSubmit {CSind} {
	global Farm output_dir job_wk_dir Label Verbose
	global CCfname CCorg
  # foreach {k v} [array get ::Farm] {puts "$k=$v"}
  # exit
  Say "CorSysIndex = $CSind"

  make_ControlCards 1 $CSind
  set flist [get_datfiles $CCfname]
  lappend flist $CCfname ewparam.txt minuit.in.txt
  foreach f {minuit.save_0.txt params_0.txt statcov_0.txt MI_saved_0.txt} {
    lappend flist0 [file join $output_dir $f]
  }
  if {$CSind} {
    set flist [concat $flist $flist0]
    set flist0 {}
  }
  if [pack_input $flist $flist0] {return 1}
  
  file rename -force -- $CCorg $CCfname
  set tmpdir ""
  if {$::LFSok} {
    set tmpdir [TempFileName [file norm $job_wk_dir]]
    file mkdir $tmpdir
    file rename input.zip $tmpdir
    file copy -force $::MainExe $tmpdir
    set drv [find_file $Farm(Driver)]
    file copy -force $drv $tmpdir
    set here [pwd]
    cd $tmpdir
  }
  set jno [job_sub1]
  # set jno 123456
  if {$::LFSok} {cd $here}

  if {$jno == ""} {return 1}
  Say "    $jno submitted." info
  
  file mkdir $job_wk_dir
  set JobList [file join $job_wk_dir $Label.jobs]
  set ch [open $JobList a]
  # puts $ch "[format %010d $jno] $logfn"
  puts $ch [join [list $::nCorSys $CSind $jno $tmpdir] \t ]
  close $ch
  # === ???
  file delete [file join $output_dir "params_[Suffix $CSind].txt"]
  return 0
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
  
	global output_dir job_wk_dir Label LFSok
  
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
  
  if {!$LFSok} {
    Say "Connecting to farm..."
    set AllJobs [farm_query_all]
    set nAll [llength $AllJobs]
    if {!$nAll} {
      Say "  No jobs on farm!" warn
      return 0
    }
  }
  
	set AllRunningJobs [farm_query_running]
	# set nRun [llength $AllRunningJobs]
  
	set FinishedJobs {}
	set RunningJobs {}
  foreach sj $SentJobs {
    set j [lindex [split $sj \t] 2]
    if {!$LFSok && ([lsearch $AllJobs $j] < 0)} continue
    if {[lsearch $AllRunningJobs $j] < 0} {lappend FinishedJobs $sj} {lappend RunningJobs $sj}
  }
  
  if {$FinishedJobs == ""} {Say "  No completed jobs to retrieve."}
  foreach sj $FinishedJobs {
    # lassign [split $sj] j tmpd
    # set rc [JobGetResults $j $tmpd]
    set rc [JobGetResults [split $sj \t]]
    if {$rc > 0} {
      Say "Sthg. wrong or missing in the job results..." err
      Say "See FitPDF.err and FitPDF.log files.\nYou may also need to look at the job output files packed in 'output.zip'."
    } elseif {$rc < 0} {
      Say "  will try to get it later."
      lappend RunningJobs $sj
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
proc job_get1 {jrec} {
  # --- Get job results from the farm
  # --- jrec = job record = list: nCS iCS jid tmpdir ...
  # --- return 0 on success
  if {[catch {farm_getresults $jrec FitPDF.err FitPDF.log $::CCfname ewparam.txt minuit.in.txt output.zip} ans ]} {
    Say "[string trim $ans]" err
    return 1
  }
  return 0
}

#===================================
proc JobGetResults {jrec} {
  # --- jrec = job record = list: nCS iCS jid tmpdir ...
  lassign $jrec nCS iCS jid tempdir
  # set jid [lindex $jrec 2]
  Say "Retrieving $jid ..." info
  global output_dir
  file mkdir $output_dir
  set here [pwd]
  cd $output_dir
  set jobdir [GenJobDir $jid]
  file mkdir $jobdir
  cd $jobdir
  set ntry 3
  
  while {$ntry && [set rc [job_get1 $jrec]]} {
    after 6000
    incr ntry -1
    Say "  trying once more..."
  }
  cd ..
  if {$::Verbose > 2} {Say "JobGetResults: $rc"}
  set rc [expr {-$rc}]
  if {!$rc} {
    set rc [TestJobResults [file join $jobdir output.zip]]
    if {!$rc} {
      set ans [exec unzip -o $jobdir/output.zip -x fittedresults.txt pulls.first.txt pulls.last.txt parsout_* 2>@1]
      if {$::Verbose > 1} {Say $ans}
      if {$tempdir != ""} {file delete -force $tempdir}
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
  if {[file exists $CCorg]} {
    if {![Yes "File '$CCorg' exists\nas if the last run stopped unexpectedly!" "Use '$CCorg' as the primary input?" yes "Ooops..."]} {exit}
  } else {file copy $CCfname $CCorg}
  set och [open $CCfname w]
  puts $och "&CSOffset"
  if {$CorSysNo != ""} {puts $och "  CorSysIndex = $CorSysNo"}
  puts $och "  UsePrevFit = $UsePrev"
  puts $och "&End"
	puts $och [ReadFile $CCorg]
  close $och
}

# ====================================
proc check_results {{nCS 0}} {
  # set here [pwd]
  global output_dir nCorSys
  if {!$nCS} {
    if {[GetNCS]} {return 1}
    set nCS $nCorSys
  }
  if {![file exists [file join $output_dir statcov_0.txt]]} {return 1}
  for {set ics -$nCS} {$ics <= $nCS} {incr ics} {
    set suf "_[Suffix $ics].txt"
    foreach f {MI_saved params Results} {
      if {![file exists [file join $output_dir "$f$suf"]]} {return 1}
    }
  }
  return 0
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
      Say "error running '$MainExe':\n$ans" err
  } {
    if {[file exists $ofn.err] && ([file size $ofn.err] == 0)} {file delete $ofn.err} {
      Say "error running '$MainExe':\nsee '$ofn.err'" err
    }
  }
  file rename -force -- $CCorg $CCfname
}

