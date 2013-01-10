#-------------------------------------------------
#  Current configuration
#-------------------------------------------------

# --- paths relative to the current RUN folder
set output_dir "output";  # --- Warning! Currently 'output' is hard-coded in the HERAFitter
set job_wk_dir "jay_work"
set MainExe bin/FitPDF

# --- driving script, Driver is first searched in the current folder
# --- and next in the folder of this script
set Driver jay.sh

set QuickInput 0 ; # --- create new input.zip for each job, otherwise update only
# set FarmName "zarah"
set FarmName "NQS"
# set FarmHost zenithsub
set MinWaitSec 8
