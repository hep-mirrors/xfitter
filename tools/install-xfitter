#!/bin/bash
#####################################################################
## Configuration ####################################################

## Programs versions
lhapdfver=6.1.6
hoppetver=1.2.0
applgridver=1.5.15
qcdnumver=17-01-14
apfelver=3.0.0
melaver=2.0.1
apfelgridver=1.0.1

## Optional manual configurations
#MANUALCONF=1

## Option 1: Use local environment
## These settings assume that root is already in your PATH and boost lib is installed
## For versions older than 1.1.0, need to specify the location of CERNLIB
#MODE=local
#export CERN_ROOT=/usr/lib

## Option 2: Setup environment from CERN cvmfs
## These settings include compiler, cernlib and root
#MODE=cern
#gccv=4.6
#os=slc6
#arch=x86_64
#rootversion=5.34.18

## End of Configuration ####################################################
#####################################################################
#Check for dummies
if [[ $0 == bash || $0 = csh ]]
then
    echo "Please don't source me, I am an executable!"
    echo "Give me attributes with:"
    echo "chmod +x install-xfitter"
    echo "and run me with:"
    echo "./install-xfitter"
    return 2
fi

if [[ -z $1 ]]
then
    echo
    echo "usage:"
    echo "$0 <version|deps>"
    echo
    echo "available versions:"
    vers=`git ls-remote --tags https://gitlab.cern.ch/fitters/xfitter.git | sed 's|/| |g; s|\^| |' | awk '{print $4}' | uniq`
    echo "$vers"
    echo "master"
    echo
    echo "to reinstall only dependences, run:"
    echo "$0 deps"
    echo
    exit
fi
mode=$1
shift 1

#in deps mode, read version from the version file
if [[ $mode != "deps" ]]
    then
    version=$mode
else
    if [[ ! -e version ]]
	then
	echo
	echo "could not find file \"version\""
	echo "cannot determine current version"
	echo "run first:"
	echo "$0 <version>"
	echo
	exit
    fi
    version=`cat version`
    echo "reinstalling dependencies for xFitter version $version"
fi
# strip the xfitter version
stripversion=`echo $version |sed "s/\.//g"`


#check that requested version exists
if [[ $version != "master" ]]
    then
    exist=0
    for ver in ``
    do
	if [[ $version == $ver ]]
	then
	    exist=1
	fi
    done

    vers=`git ls-remote --tags https://gitlab.cern.ch/fitters/xfitter.git | sed 's|/| |g; s|\^| |' | awk '{print $4}' | uniq`
 
    for ver in $vers
    do
	if [[ $version == $ver ]]
	then
	    exist=1
	fi
    done

    if [[ $exist == 0 ]]
    then
	echo
	echo "version $version not found, available versions:"
	echo ""
	echo "$vers"
	echo "master"
	echo
	exit
    fi
fi

if [[ $mode != "deps" && -e xfitter-${version} &&  -e herafitter-${version} ]]
then
    echo
    echo "xfitter-${version} already exists, remove it first"
    echo "To reinstall only dependences, run:"
    echo "$0 deps"
    echo
    exit
fi
if [[ $mode == "deps" && ! -e xfitter-${version}  && ! -e herafitter-${version}  ]]
then
    echo
    echo "xfitter-${version} or herafitter-${version}   does not exist, install it first with:"
    echo "$0 $version"
    echo
    exit
fi

#automatically detect system:
if [[ -z $MANUALCONF ]]
then
    which sw_vers >& /dev/null
    if [[ $? == 0 ]]
    then
	echo "Detected Mac OS X system"
	MODE=local
    else
	SYS=$(echo `lsb_release -i |cut -d: -f2`)
	ver=$(echo `lsb_release -r |cut -d: -f2`)
	if [[ $SYS == Scientific* && $ver == 6.* ]]
	then
	    echo "Detected SL6 Linux distribution"
	    MODE=cern
	    gccv=4.8
	    os=slc6
	    arch=x86_64
	    rootversion=5.34.36
	    boostver=1.53.0
	    pyth=2.7
        elif [[ $SYS == CentOS* && $ver == 7.* ]]
        then
            echo "Detected CentOS7 Linux distribution"
            MODE=cern
            gccv=4.8
            os=centos7
            arch=x86_64
            rootversion=6.08.04
            boostver=1.53.0
            pyth=2.7
	elif [[ $SYS == Scientific* && $ver == 5.* ]]
	then
	    echo "Detected SL5 Linux distribution"
	    MODE=cern
	    gccv=4.3
	    os=slc5
	    arch=x86_64
	    rootversion=5.34.00
	    boostver=1.48.0
	    python=2.6.5
	    pyth=2.6
	elif [[ $SYS == "Ubuntu" ]]
	then
	    echo "Detected Ubuntu distribution"
	    MODE=local
	else
	    echo "Sorry, I don't recognize your system:"
	    echo "$SYS $ver"
	    echo "I will assume you have root installed in your system,"
            echo "gcc version >= 4.3, python, boost libraries, and wget"
	    echo "If this doesn't work, and you have /cvmfs/sft.cern.ch mounted"
	    echo "edit me (I am $0) and try to setup appropriate settings"
	    echo "in the section: manual configuration"
	    echo
	    MODE="local"
	fi
    fi
fi
if [[ $MODE == "cern" ]]
    then
#    if [[ ! -e /afs/cern.ch ]]
    if [[ ! -e /cvmfs/sft.cern.ch ]]
	then
	echo
	echo "/cvmfs/sft.cern.ch not mounted, forcing local MODE"
	echo "Fasten you seat belt"
	echo "I hope you have root, gcc >= 4.3, python and boost libraries"
	echo "all installed in your system"
	echo
	MODE="local"
    fi
fi

if [[ $MODE == "cern" ]]
then
    compiler=`echo gcc${gccv} | sed "s/\.//"`
#    . /afs/cern.ch/sw/lcg/contrib/gcc/${gccv}/${arch}-${os}/setup.sh
    . /cvmfs/sft.cern.ch/lcg/contrib/gcc/${gccv}/${arch}-${os}/setup.sh
#    . /afs/cern.ch/sw/lcg/app/releases/ROOT/${rootversion}/${arch}-${os}-${compiler}-opt/root/bin/thisroot.sh
    . /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/${rootversion}/${arch}-${os}-${compiler}-opt/root/bin/thisroot.sh
    if [[ $os == slc5 ]]
    then
        echo "LEGACY SL5 ! using afs"
	PYTHONBIN=/afs/cern.ch/sw/lcg/external/Python/${python}/${arch}-${os}-${compiler}-opt/bin
	PATH=$PYTHONBIN:$PATH
	export BOOST=--with-boost=/afs/cern.ch/sw/lcg/external/Boost/${boostver}_python${pyth}/${arch}-${os}-${compiler}-opt
    fi
    if [[ $os == slc6 ]]
    then
	export BOOST=--with-boost=/cvmfs/sft.cern.ch/lcg/external/Boost/${boostver}_python${pyth}/${arch}-${os}-${compiler}-opt
    fi
fi

#check some basic dependendencies before starting the installation
which git >& /dev/null
if [[ $? != 0 ]]
then
    echo "Error, git not found"
    exit
fi

which root >& /dev/null
if [[ $? != 0 ]]
then
    echo "Error, root not found"
    exit
fi

which wget >& /dev/null
if [[ $? == 0 ]]
then
    http=wget
else
    which curl >& /dev/null
    if [[ $? == 0 ]]
    then
	http=curl
    else
	echo "Error, wget or curl not found"
	exit
    fi
fi
    
#directory:
CURRENTDIR=`pwd`

#clean up
rm version setup.sh compile quickstart.readme.txt >& /dev/null
rm install.log >& /dev/null



# keep for debugging
installDeps=1


if [[ $installDeps == 0 ]]
then
   echo "Skip installation of dependences"
else    
#Make all dependencies
    rm -rf deps >& /dev/null
    mkdir deps
    cd deps
#lhapdf:
    echo "Installing LHAPDF $lhapdfver..."
    if (( `echo $lhapdfver |cut -d. -f1` >= 6 ))
    then
	lhapdf="LHAPDF"
	withboost=$BOOST
    else
	lhapdf="lhapdf"
    fi

    if [[ $http == "curl" ]]
    then
	curl https://www.hepforge.org/archive/lhapdf/${lhapdf}-${lhapdfver}.tar.gz > ${lhapdf}-${lhapdfver}.tar.gz 2>> $CURRENTDIR/install.log
    else
	wget https://www.hepforge.org/archive/lhapdf/${lhapdf}-${lhapdfver}.tar.gz >> $CURRENTDIR/install.log 2>&1
    fi
    tar xfz ${lhapdf}-${lhapdfver}.tar.gz  >> $CURRENTDIR/install.log 2>&1
    cd ${lhapdf}-${lhapdfver} 
    ./configure --prefix=$CURRENTDIR/deps/lhapdf $withboost >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    make -j 9 install >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    cd - >& /dev/null

 #hoppet:
    echo "Installing HOPPET $hoppetver..."
    if [[ $http == "curl" ]]
    then
	curl http://hoppet.hepforge.org/downloads/hoppet-${hoppetver}.tgz > hoppet-${hoppetver}.tgz 2>> $CURRENTDIR/install.log
    else
	wget http://hoppet.hepforge.org/downloads/hoppet-${hoppetver}.tgz >> $CURRENTDIR/install.log 2>&1
    fi
    tar xfz hoppet-${hoppetver}.tgz  >> $CURRENTDIR/install.log  2>&1
    cd hoppet-${hoppetver}
    ./configure --prefix=$CURRENTDIR/deps/hoppet  >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    make -j 9 install  >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    cd - >& /dev/null

 # setup paths for applgrid:
    export PATH=$CURRENTDIR/deps/hoppet/bin/:$PATH
    export PATH=$CURRENTDIR/deps/lhapdf/bin/:$PATH

 #applgrid:
    echo "Installing APPLGRID $applgridver..."
    if [[ $http == "curl" ]]
    then
	curl https://www.hepforge.org/archive/applgrid/applgrid-$applgridver.tgz > applgrid-$applgridver.tgz  2>> $CURRENTDIR/install.log
    else
	wget https://www.hepforge.org/archive/applgrid/applgrid-$applgridver.tgz >> $CURRENTDIR/install.log 2>&1
    fi
    tar xfz applgrid-$applgridver.tgz  >> $CURRENTDIR/install.log  2>&1
    cd applgrid-$applgridver
    ./configure --prefix=$CURRENTDIR/deps/applgrid  >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    make   >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    make -j 9 install  >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    cd - >& /dev/null

    export PATH=$CURRENTDIR/deps/applgrid/bin/:$PATH

 #apfel
    echo "Installing APFEL $apfelver..."
    if [[ $http == "curl" ]]
    then
	curl https://github.com/scarrazza/apfel/archive/${apfelver}.tar.gz > ${apfelver}.tar.gz 2 >> $CURRENTDIR/install.log
    else
	wget https://github.com/scarrazza/apfel/archive/${apfelver}.tar.gz >> $CURRENTDIR/install.log 2>&1
    fi
    mv ${apfelver}.tar.gz apfel-${apfelver}.tar.gz
    tar xfvz apfel-${apfelver}.tar.gz >> $CURRENTDIR/install.log 2>&1
    cd apfel-${apfelver}
    ./configure --prefix=$CURRENTDIR/deps/apfel --disable-lhapdf >> $CURRENTDIR/install.log  2>&1

    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    make -j 9 install  >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    cd - >& /dev/null
 # setup paths for apfel:
    export PATH=$CURRENTDIR/deps/apfel/bin/:$PATH

 #apfelgrid
    lhapdf get NNPDF30_nlo_as_0118
    echo "Installing APFELgrid $apfelgridver..."
    if [[ $http == "curl" ]]
    then
	curl https://github.com/nhartland/APFELgrid/archive/v${apfelgridver}.tar.gz > v${apfelgridver}.tar.gz 2 >> $CURRENTDIR/install.log
    else
	wget https://github.com/nhartland/APFELgrid/archive/v${apfelgridver}.tar.gz >> $CURRENTDIR/install.log 2>&1
    fi
    mv v${apfelgridver}.tar.gz APFELgrid-${apfelgridver}.tar.gz
    tar xfvz APFELgrid-${apfelgridver}.tar.gz >> $CURRENTDIR/install.log 2>&1
    cd APFELgrid-${apfelgridver}
    ./setup.sh  >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 1 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    cd - >& /dev/null

 #mela
    echo "Installing MELA $melaver..."

    if [[ $http == "curl" ]]
    then
	curl https://github.com/vbertone/MELA/archive/${melaver}.tar.gz > ${melaver}.tar.gz 2 >> $CURRENTDIR/install.log
    else
	wget https://github.com/vbertone/MELA/archive/${melaver}.tar.gz >> $CURRENTDIR/install.log 2>&1
    fi
    mv ${melaver}.tar.gz MELA-${melaver}.tar.gz
    tar xfvz MELA-${melaver}.tar.gz >> $CURRENTDIR/install.log 2>&1
    cd MELA-${melaver}
    ./configure --prefix=$CURRENTDIR/deps/mela  >> $CURRENTDIR/install.log  2>&1

    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi  
    make -j 9 install  >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    cd - >& /dev/null
 # setup paths for mela:
    export PATH=$CURRENTDIR/deps/mela/bin/:$PATH

 #qcdnum
    echo "Installing QCDNUM $qcdnumver..."
    qcdnumstripver=`echo $qcdnumver |sed "s/-//g"`
    if [[ $http == "curl" ]]
    then
	curl http://www.nikhef.nl/user/h24/qcdnum-files/download/qcdnum${qcdnumstripver}.tar.gz > qcdnum${qcdnumstripver}.tar.gz 2>> $CURRENTDIR/install.log
    else
	wget http://www.nikhef.nl/user/h24/qcdnum-files/download/qcdnum${qcdnumstripver}.tar.gz >> $CURRENTDIR/install.log 2>&1
    fi
    tar xfz qcdnum${qcdnumstripver}.tar.gz  >> $CURRENTDIR/install.log  2>&1
    cd qcdnum-${qcdnumver}

    ./configure --prefix=$CURRENTDIR/deps/qcdnum  >> $CURRENTDIR/install.log  2>&1
    export PATH=$CURRENTDIR/deps/qcdnum/bin/:$PATH
    
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi 
    make -j 9 install  >> $CURRENTDIR/install.log  2>&1
    if [[ $? != 0 ]]
    then
	echo "Error, check install.log for details"
	exit
    fi
    cd - >& /dev/null
fi
cd $CURRENTDIR


 #xfitter
if [[ $mode != "deps" ]]
then
    echo "Installing xFitter $version..."

    
    if [[ $version == "master" ]]
    then
	git clone https://gitlab.cern.ch/fitters/xfitter.git >> $CURRENTDIR/install.log  2>&1 
	mv xfitter xfitter-master
    else
	if [[ $http == "curl" ]]
	then
#	    curl  https://gitlab.cern.ch/fitters/xfitter/repository/archive.tar.gz\?ref=$version > xfitter-$version.tar.gz
	    curl https://www.hepforge.org/archive/xfitter/tar_files/xfitter-$version.tar.gz > xfitter-$version.tar.gz
	else
#	    wget  https://gitlab.cern.ch/fitters/xfitter/repository/archive.tar.gz\?ref=$version >> $CURRENTDIR/install.log 2>&1
#	    mv archive.tar.gz\?ref=$version xfitter-$version.tar.gz
	    wget https://www.hepforge.org/archive/xfitter/tar_files/xfitter-$version.tar.gz  >> $CURRENTDIR/install.log 2>&1
	fi
       # unpack nicely:
	rm -fr xfitter-${version} 
	mkdir xfitter-${version} ; tar xfz xfitter-${version}.tar.gz -C xfitter-${version} --strip-components 1
    fi
else
    make -C xfitter-${version} clean >> $CURRENTDIR/install.log  2>&1
fi


#make a setup run enviroment script
echo $version > version
echo "export CURRENTDIR=`pwd`" > setup.sh
echo "export version=\`cat version\`" >> setup.sh
echo "export PATH=$CURRENTDIR/xfitter-$version/bin:\$PATH" >> setup.sh
echo "export PATH=$CURRENTDIR/deps/hoppet/bin:\$PATH" >> setup.sh
echo "export PATH=$CURRENTDIR/deps/applgrid/bin:\$PATH" >> setup.sh
echo "export PATH=$CURRENTDIR/deps/lhapdf/bin:\$PATH" >> setup.sh
echo "export PATH=$CURRENTDIR/deps/apfel/bin:\$PATH" >> setup.sh
echo "export PATH=$CURRENTDIR/deps/mela/bin:\$PATH" >> setup.sh
echo export LD_LIBRARY_PATH=\$CURRENTDIR/deps/hoppet/lib/:\$LD_LIBRARY_PATH   >> setup.sh
echo export LD_LIBRARY_PATH=\$CURRENTDIR/deps/lhapdf/lib/:\$LD_LIBRARY_PATH   >> setup.sh
echo export LD_LIBRARY_PATH=\$CURRENTDIR/deps/applgrid/lib/:\$LD_LIBRARY_PATH >> setup.sh
echo export LD_LIBRARY_PATH=\$CURRENTDIR/deps/apfel/lib/:\$LD_LIBRARY_PATH >> setup.sh
echo export LD_LIBRARY_PATH=\$CURRENTDIR/deps/mela/lib/:\$LD_LIBRARY_PATH >> setup.sh
echo export LD_LIBRARY_PATH=\$CURRENTDIR/deps/qcdnum/lib/:\$LD_LIBRARY_PATH >> setup.sh
echo "export PATH=$CURRENTDIR/deps/qcdnum/bin:\$PATH" >> setup.sh


if [[ $MODE == "cern" ]]
then
#    echo . /afs/cern.ch/sw/lcg/contrib/gcc/${gccv}/${arch}-${os}/setup.sh            >> setup.sh
    echo . /cvmfs/sft.cern.ch/lcg/contrib/gcc/${gccv}/${arch}-${os}/setup.sh            >> setup.sh
#    . /afs/cern.ch/sw/lcg/app/releases/ROOT/${rootversion}/${arch}-${os}-${compiler}-opt/root/bin/thisroot.sh
#    echo cd /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/${rootversion}/${arch}-${os}-${compiler}-opt/root >> setup.sh

    echo "cd /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/${rootversion}/${arch}-${os}-${compiler}-opt/root/" >> setup.sh
    echo ". /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/${rootversion}/${arch}-${os}-${compiler}-opt/root/bin/thisroot.sh ">> setup.sh
    echo "cd -" >>setup.sh
fi

#make a compilation script
echo source ./setup.sh > compile
echo export PATH=\$CURRENTDIR/deps/hoppet/bin/:\$PATH       >> compile
echo export PATH=\$CURRENTDIR/deps/lhapdf/bin/:\$PATH       >> compile
echo export PATH=\$CURRENTDIR/deps/applgrid/bin/:\$PATH     >> compile
echo export PATH=\$CURRENTDIR/deps/apfel/bin/:\$PATH     >> compile
echo export PATH=\$CURRENTDIR/deps/mela/bin/:\$PATH     >> compile


echo cd xfitter-\$version                                        >> compile
echo autoreconf --install                                 >> compile

echo ./configure --enable-applgrid --enable-lhapdf --enable-apfel --enable-mela  --enable-apfelgrid     >> compile

echo make -j 9 install                                    >> compile

chmod +x compile
echo "Compiling xFitter $version..."

./compile >> $CURRENTDIR/install.log  2>&1
if [[ $? != 0 ]]
then
    echo "Error, check install.log for details"
    exit
fi

# run test
source ./setup.sh
cd xfitter-${version}

bin/xfitter >> $CURRENTDIR/install.log  2>&1



if [[ $? != 0 ]]
then
    echo "Error in testing xfitter executable, check install.log for details"
    exit
fi
cd - >& /dev/null
echo "xFitter installation successful!"
echo "Check install.log file for details"
echo
# setup a run dir
if [[ ! -e run ]]
then
    mkdir -p run
    cp  xfitter-${version}/steering.txt \
	xfitter-${version}/minuit.in.txt \
	xfitter-${version}/ewparam.txt \
	run
    rsync -a --exclude=".*" xfitter-${version}/datafiles run/
    rsync -a --exclude=".*" xfitter-${version}/theoryfiles run/
else
    echo "\"run\" directory already exists, I won't touch it"
    echo
fi

#quick start guide
if [[ ! -e quickstart-readme.txt ]]
then
    echo "for a quick start do:" >> quickstart-readme.txt
    echo >> quickstart-readme.txt
    echo "source setup.sh #setup environment" >> quickstart-readme.txt
    echo "cd run          #enter your run directory" >> quickstart-readme.txt
    echo "xfitter          #run your fit" >> quickstart-readme.txt
    echo >> quickstart-readme.txt
    echo "The fit options are controlled by the 3 files steering.txt, minuit.in.txt, and ewparam.txt" >> quickstart-readme.txt
    echo "Data and theory files are in datafiles directory" >> quickstart-readme.txt
    echo >> quickstart-readme.txt
    echo "After the fit, display you results with:" >> quickstart-readme.txt
    echo "xfitter-draw output/" >> quickstart-readme.txt
    echo "A file output/plots.pdf will be created" >> quickstart-readme.txt
    echo "You can find further documentation in xfitter-$version/README and xfitter-$version/doc/tex/manual/" >> quickstart-readme.txt
    echo "-----------------------------------------------------------"
    cat quickstart-readme.txt
    echo "-----------------------------------------------------------"
    echo 
    echo "To read again these instructions, see quickstart-readme.txt"
    echo "Have fun!"
fi
