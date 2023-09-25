#First, perform profiling for all study cases.
for dir in "stat" "stat+syst"; do
    cd ${dir}
    sh run.sh
    cd ..
done
