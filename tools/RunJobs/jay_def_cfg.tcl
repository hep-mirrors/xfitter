# -------------------------------------------------
#     Default configuration
# -------------------------------------------------

# --- Can Farm access the local file system?
# --- PBS yes, NQS no
set LFSok 0

if [info exists env(HERAFITTER_SYS)] {
  set MainCodeDir $env(HERAFITTER_SYS)
} {
  set MainCodeDir [file normalize [file join [file dirname [info script]] ../..]]
}
set MainExeRel bin/FitPDF ;  # --- relative to $MainCodeDir

# --- paths relative to the current RUN folder:
set job_wk_dir "jay_work"
# --- Warning! The following names are currently hard-coded in the HERAFitter
set output_dir "output"
set CCfname "steering.txt"

set DefaultLabel HFit
set MinWait 1
