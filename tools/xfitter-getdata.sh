#!/bin/bash
#########################################################################
# script to download data files from hepforge stored in xFitter package #
# author: RP, OZ, March 2017                                            #
#########################################################################
#
if [ $# -eq 0 ]
  then
      echo "no argument supplied, please type ./xfitter-getdata.sh -h for help"
      exit 0
fi    

if [[ "$1" = "help" ]] || [[ "$1" = "--help" ]] || [[ "$1" = "-h" ]] || [[ "$1" = "h" ]] ; then
    echo "------------------------------------------------------------"
    echo "xfitter-getdata.sh script allows to download public data"
    echo "sets from the xFitter package (from hepforge)"  
    echo " "  
    echo "For the list of available data sets please check:"
    echo "  http://xfitter.hepforge.org/data.html "
    echo " "  
    echo "or input_steering/steering.txt.ALLdata in xFitter package "
    echo " "  
    echo "To download all available data sets do: " 
    echo " ./xfitter-getdata.sh ALL "
    echo " "  
    echo "To download a specific data set do: "
    echo " ./xfitter-getdata.sh arXivNumber(reportNumber) "
    echo " "  
    echo "this will download the indicated dataset with the directory " 
    echo "structure ready to be used in xFitter "
    echo "(dirCollider/dirExperiment/dirReactionType/dirArxivNumber)"
    echo " "  
    echo "Note, that related theory or correlation files (if any) "
    echo "will be downloaded together with the specified data file "
    echo " "  
    echo "===> to see available data sets to download type: "
    echo " ./xfitter-getdata.sh --print (or -p) "
    echo "------------------------------------------------------------"
    exit 0
fi

if [[ "$1" = "-p" ]] || [[ "$1" = "--print" ]] ; then
    echo "-------------------------------------------------------------------------------------"
    echo -e "\033[34m available data sets in xFitter (to download use arXivNumber or reportNumber) \e[0m"
    printf "%15s%20s%20s%30s\n" "Collider" "Experiment" "Reaction" "arXivNumber or reportNumber";
    # read html page with datafiles to tmp file
    lynx --dump http://xfitter.hepforge.org/data.html > o
    # remove all lines above References AND remove README lines 
    cat o | sed  -e '1,/References/d' -e '/README/d' > o1
    # remove all lines below References 
    cat o | sed  -e '/References/q' > o2
    # loop over arxiv numbers and retrieve additional inforamtion
    # print only text after last slash AND remove tar.gz at the end of string
    for num in `awk '{gsub(".tar.gz", " ", $0); print $NF}' FS=/ o1`; do
      #echo $num;
      cat o2 | grep $num | sed -e 's/\[\([^]]*\)\]//g' | awk '{printf "%15s%20s%20s%30s\n", $2, $3, $4, $5}'
    done
    # rm not needed file 
    rm o o1 o2
    echo " "  
    echo -e "\033[92m example: ./xfitter-getdata.sh 1412.2862 \e[0m"
    echo "------------------------------------------------------------"
    exit 0
fi


# first dowload all data if requested:
if [[ $1 = "ALL" ]] ; then
    if wget -q --spider www.hepforge.org/archive/xfitter/data_files/ALL.tar.gz > /dev/null; then
       echo " ----------------------------------------------------------- "
       echo "ALL.tar.gz file found on hepforge, will proceed with the download "
       echo ""
       wget www.hepforge.org/archive/xfitter/data_files/ALL.tar.gz
       #option -k does not allow to overwrite existing diretory/files
       tar -k -zxvf ALL.tar.gz
       if [ $? -ne 0 ]; then
            echo " ------------------------------------------------------ "
            echo -e "\033[92m datafiles directory exists, do you want to overwrite it? (y/n) \e[0m"
            read answer
            # of user wants to overwrite existing files, let him/her
            if echo "$answer" | grep -iq "^y" ;then
               echo " overwriting  existing datafiles ..."
               tar -zxvf ALL.tar.gz
            else 
                echo " exiting... (no data were extracted) "
            fi
       fi
       # ask if user wants to remove tar.gz file (needed?)       
       echo -e "\033[31m do you want to remove ALL.tar.gz? \e[0m"
       read answer
       if echo "$answer" | grep -iq "^y" ;then
          echo " removing ALL.tar.gz ..."
          rm ALL.tar.gz
       fi
       echo ""
       echo -e "\033[92m  ------ done. xFitter support: xfitter-help@desy.de ----------------  \e[0m"
    else 
        echo "ALL.tar.gz file does not exist, please send email to xfitter-help@desy.de"
        exit
    fi
# now dowload separate files rif equested:
elif [[ $1 != "ALL" ]] ; then
    for arg in "$@"; do
        if wget -q --spider www.hepforge.org/archive/xfitter/data_files/$arg.tar.gz > /dev/null; then
           echo " ----------------------------------------------------------- "
           echo "$arg found on hepforge, will proceed with the download "
           echo ""
           wget www.hepforge.org/archive/xfitter/data_files/$arg.tar.gz
           #option -k does not allow to overwrite existing diretory/files
           tar -k -zxvf $arg.tar.gz 
           if [ $? -ne 0 ]; then
                echo " ------------------------------------------------------ "
                echo -e "\033[92m File you trying to extract already exists, do you want to overwrite it? (y/n) \e[0m"
                read answer
                # of user wants to overwrite existing files, let him/her
                if echo "$answer" | grep -iq "^y" ;then
                   echo " overwriting  existing file $arg ..."
                   tar -zxvf $arg.tar.gz
                else 
                    echo " exiting... (no data were extracted) "
                fi
           fi
           # ask if user wants to remove tar.gz file (needed?)       
           echo -e "\033[92m do you want to remove $arg.tar.gz? (y/n) \e[0m"
           read answer
           if echo "$answer" | grep -iq "^y" ;then
              echo " removing $arg.tar.gz ..."
              rm $arg.tar.gz
           fi
           echo ""
           echo -e "\033[92m  ------ done. xFitter support: xfitter-help@desy.de ----------------  \e[0m"
        # if data set doesn't exist...   
        else 
            echo -e "\033[31m $arg does not exist on http://xfitter.hepforge.org/data.html, please check!  \e[0m"
            exit
        fi
    done
fi
