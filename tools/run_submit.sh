#!/bin/sh

working_directory="/scratch/jtomasz/HERAFITTER/h1fitter/RUN"
JobSubmit="jobsub -h zenithsub  -qM eaze_run.sh -f steering.txt FitPDF ewparam.txt minuit.in.txt datafiles/hera/ZEUS_LRG_98-00.dat datafiles/hera/ZEUS_LPS_98-00.dat"
JobStatus="jobq -h zenithsub $name"
JobGet="jobget"

time=50

./run_fitter.sh 0 0 
#cat >  $working_directory/minuit.ini.txt <<EOF_CARDS
#return 
#EOF_CARDS

./FitPDF $working_directory/stering.txt

NrParameters=`head -1 $working_directory/output/params_0.txt`
echo "NrParameters " $NrParameters

\rm -r $working_directory/minuit.ini.txt
#create minuit.ini.txt
head -61 $working_directory/output/minuit.save_0.txt  > $working_directory/minuit.ini.txt
cat >> $working_directory/minuit.ini.txt <<EOF_CARDS1

set eps 1.0e-11
*set print 3
call fcn 3
migrad 200000
*hesse
set print 3
save
return
EOF_CARDS1

i=1
\rm  $working_directory/ListOfJobs
until [ $i == $NrParameters ]    ; do 
     name=steering.txt
     ./run_fitter.sh $i 0
     $JobSubmit > $working_directory/file_tmp
     gawk -F" " '{ print $2 }' $working_directory/file_tmp >> $working_directory/ListOfJobs
     i=`expr $i + 1`
done

i=0
NrNegParameters=`expr 0 - $NrParameters`

until [ $i =  $NrNegParameters ]   ; do 
    i=`expr $i - 1`
     name=steering.txt
     ./run_fitter.sh $i 0
     $JobSubmit > $working_directory/file_tmp
     gawk -F" " '{ print $2 }' $working_directory/file_tmp >> $working_directory/ListOfJobs
done


# take jobs 

for file in $working_directory/ListOfJobs ; do
  echo `basename  $working_directory/ListOfJob ` >> file_list
done


\rm  $working_directory/output_jobq

sleep 1 
$JobStatus > $working_directory/output_jobq
NrJobs="`wc -l $working_directory/output_jobq |  gawk -F" " '{ print $1 }'`"


# 4 number of output line in comman jobq if all files are ready

until [ $NrJobs == 4 ] ; do 
    $JobStatus > $working_directory/output_jobq
    NrJobs="`wc -l $working_directory/output_jobq |  gawk -F" " '{ print $1 }'`"
    sleep $time 
done 

NoRunJobs="`head -1 $working_directory/output/params_0.txt`"
echo "NoRunJobs" $NoRunJobs

until [ "$NoRunJobs" == "No unfinished job found" ] ; do
$JobStatus > $working_directory/output_jobq
NoRunJobs="`tail -1 $working_directory/output_jobq`"
sleep $time 
done

file_list=" `awk '{print $1}' $working_directory/ListOfJobs`"

i=0
for name in $file_list ; do
    i=`expr $i + 1`
    echo $name $i
    cd  $working_directory/output
    $JobGet $name
    
done 


cd $working_directory

# collect information from all jobs
$working_directory/run_fitter.sh 0 2
$working_directory/FitPDF steering.txt
