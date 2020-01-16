#The format of input files was changed while migrating to cmake (commit 084a0ca01d51bd26dbd10677bc80b629e080e007)
#This script migrates parameters.yaml and steering.txt from before that merge to after that merge
#Run it from the directory where the parameters.yaml and steering.txt are
sed -i "s/\(\? \!include .*\)\/parameters.yaml/\1.yaml/;s/UvDvubardbars/UvDvUbarDbarS/" parameters.yaml
sed -i "/\(TheoryType\|RunningMode\|HF_SCHEME\|LUseAPPLgridCKM\|WriteLHAPDF6\) = /d" steering.txt
