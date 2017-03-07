#!/bin/bash
#########################################################################
# script to download data files from hepforge stored in xFitter package #
# author: RP, March 2017                                                #
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
    echo "------------------------------------------------------------"
    echo "available data sets in xFitter to download (arXivNumber or reportNumber) "
    # read html page with datafiles to tmp file
    lynx --dump http://xfitter.hepforge.org/data.html > o
    # remove all lines above References AND remove README lines 
    sed  -i -e '1,/References/d' -e '/README/d' o
    # print only text after last slash AND remove tar.gz at the end of string
    awk '{gsub(".tar.gz", " ", $0); print $NF}' FS=/ o
    # rm not needed file 
    rm o
    echo " "  
    echo "example: ./xfitter-getdata.sh 1412.2862 "
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
            echo " datafiles diretory exists, do you want to overwrite it? (y/n) "
            read answer
            # of user wants to overwrite existing files, let him/her
            if echo "$answer" | grep -iq "^y" ;then
               echo " overwriting  existing datafiles ..."
               tar -zxvf ALL.tar.gz
            else 
                echo " exiting... (no data were extracted) "
                return 0
            fi
       fi
       # ask if user wants to remove tar.gz file (needed?)       
       echo " do you want to remove ALL.tar.gz? (y/n) "
       read answer
       if echo "$answer" | grep -iq "^y" ;then
          echo " removing ALL.tar.gz ..."
          rm ALL.tar.gz
       else 
           return 0
       fi
       echo ""
       echo " ------ done. xFitter support: xfitter-help@desy.de ---------------- "
    else 
        echo "ALL.tar.gz file does not exist, please send email to xfitter-help@desy.de"
        return
    fi
# now dowload separate files rif equested:
elif [[ $1 != "ALL" ]] ; then
    if wget -q --spider www.hepforge.org/archive/xfitter/data_files/$1.tar.gz > /dev/null; then
       echo " ----------------------------------------------------------- "
       echo "$1 found on hepforge, will proceed with the download "
       echo ""
       wget www.hepforge.org/archive/xfitter/data_files/$1.tar.gz
       #option -k does not allow to overwrite existing diretory/files
       tar -k -zxvf $1.tar.gz 
       if [ $? -ne 0 ]; then
            echo " ------------------------------------------------------ "
            echo " File you trying to extract already exists, do you want to overwrite it? (y/n) "
            read answer
            # of user wants to overwrite existing files, let him/her
            if echo "$answer" | grep -iq "^y" ;then
               echo " overwriting  existing file $1 ..."
               tar -zxvf $1.tar.gz
            else 
                echo " exiting... (no data were extracted) "
                return 0
            fi
       fi
       # ask if user wants to remove tar.gz file (needed?)       
       echo " do you want to remove $1.tar.gz? (y/n) "
       read answer
       if echo "$answer" | grep -iq "^y" ;then
          echo " removing $1.tar.gz ..."
          rm $1.tar.gz
       else 
           return 0
       fi
       echo ""
       echo " ------ done. xFitter support: xfitter-help@desy.de ---------------- "
    # if data set doesn't exist...   
    else 
        echo "$1 does not exist on http://xfitter.hepforge.org/data.html, please check!"       
        return
    fi
fi


