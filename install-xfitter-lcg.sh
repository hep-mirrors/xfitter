#!/bin/bash

#####################################################################
# 
#
# Installation and setup script for xFitter using the LCG environment
# from CERN through /cvmfs
#
#
#####################################################################



# ----------------------------------------------------------------- #
#                USER INPUT STARTS HERE                             #
# ----------------------------------------------------------------- #
#  Please specify versions
# ----------------------------------------------------------------- #
xfitterbranch=fastNLO-v2.6        # default: main [or master?]
yamlver=0.2.5                     # default: 0.2.5
qcdnumver=18-00-00                # default: 18-00-00
applgridver=1.6.36                # default: 1.6.36
apfelxxver=                       # default: 4.8.0
pineapplver="0.6.0-alpha.17"      # default: 0.6.0-alpha.17
dyturbover=                       # default: 1.4.2
# ----------------------------------------------------------------- #
#                      END OF USER INPUT                            # 
# ----------------------------------------------------------------- #



#####################################################################
#  setup the install script and environment
#####################################################################
# --- test environmental variables
#      - PLATFORM
#      - LCG_VERSION
if [[ ! -f /cvmfs/sft.cern.ch/lcg/views/$LCG_VERSION/$PLATFORM/setup.sh ]]; then
    echo " +-----------------------------------------------------------------------------------------------+"
    echo " | Error! Cannot find LCG setup script. Probably variables PLATFORM and/or LCG_VERSION are not set!"
    echo " |        Please set them, and run this install script again."
    echo " |        See https://lcginfo.cern.ch for further details and options."
    echo " |        As example:"
    echo " +-----------------------------------------------------------------------------------------------+"
    echo "export PLATFORM=x86_64-el9-gcc14-opt"
    echo "export LCG_VERSION=LCG_107a"
    echo " +-----------------------------------------------------------------------------------------------+"
    exit 1
fi


#####################################################################
# --- prepare  setup script
#####################################################################
export CURRENTDIR=$PWD
INSTALLDIR=$CURRENTDIR/deps

# --- Write setup file
cat << EOF >setup-lcgenv-$PLATFORM-$LCG_VERSION.sh
#!/bin/bash  
PLATFORM=$PLATFORM
LCG_VERSION=$LCG_VERSION
source /cvmfs/sft.cern.ch/lcg/views/$LCG_VERSION/$PLATFORM/setup.sh
export PATH=$INSTALLDIR/bin:\$PATH
export LD_LIBRARY_PATH=$INSTALLDIR/lib:\$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$INSTALLDIR/lib/pkgconfig:\$PKG_CONFIG_PATH
export yaml_DIR=$INSTALLDIR
gcc --version | head -1
which gcc
echo "Info. Using:"
echo "      PLATFORM    = $PLATFORM"
echo "      LCG_VERSION = $LCG_VERSION"
echo "      gcc         = \$(which gcc)"
echo "      version     = \$(gcc --version | head -1)"
EOF

chmod +x setup-lcgenv-$PLATFORM-$LCG_VERSION.sh

#####################################################################
# --- Print info message
#####################################################################
echo " +-----------------------------------------------------------------------------------------------+"
echo " | Info. Using:"
echo " |       PLATFORM    = $PLATFORM"
echo " |       LCG_VERSION = $LCG_VERSION"
echo " +-----------------------------------------------------------------------------------------------+"
echo " | "
echo " | Call 'source setup-lcgenv-$PLATFORM-$LCG_VERSION.sh' when you come back."
echo " | "
echo " +-----------------------------------------------------------------------------------------------+"
echo " source setup-lcgenv-$PLATFORM-$LCG_VERSION.sh"
echo " +-----------------------------------------------------------------------------------------------+"


# running setup
source setup-lcgenv-$PLATFORM-$LCG_VERSION.sh

# --- Test setup
if [[ "$CC" != "$(which gcc)" ]]; then
    echo "Error! CC and gcc path are not consistent! Please check your setup!"
    echo $CC
    echo $(which gcc)
    echo "$CC"
    echo "$(which gcc)"
    exit 1
fi
if [[ "$COMPILER_PATH/bin/gcc" != $(which gcc) ]]; then
    echo "Error! COMPILER_PATH and gcc path are not consistent! Please check your setup!"
    exit 1
fi


#####################################################################
# printout
#####################################################################
echo " +-----------------------------------------------------------------------------------------------+"
echo " | Package versions (if unset, then package will not be installed/provided): "
echo " | "
echo " | xFitter:    $xfitterbranch"
echo " | yaml:       $yamlver"
echo " | QCDNUM:     $qcdnumver"
echo " | APPLgrid:   $applgridver"
echo " | Apfel++:    $apfelxxver"
echo " | Pineappl:   $pineapplver"
echo " | DYTurbo:    $dyturbover"
echo " | "
echo " | Available packages from lcg"
echo " | gcc:        $(gcc --version | head -1)"
echo " | ROOT:       $(root-config --version)"
echo " | GSL:        $(gsl-config --version)"
echo " | LHAPDF:     $(lhapdf-config --version)"
echo " | Hoppet:     $(hoppet-config 2>&1 | head -1)"
echo " | Apfel:      $(apfel-config --version)"
echo " | yaml-cpp:   $(pkg-config --variable=libdir yaml-cpp)"
echo " | "
echo " | Info:"
echo " | Skipping Hoppet, since master version is needed, but its installation is not available in this skipt (do it!) " # -DCMAKE_DISABLE_FIND_PACKAGE_HOPPET=TRUE
echo " | Skipping DYTurbo due to compilation/linking issue in xFitter" # -DCMAKE_DISABLE_FIND_PACKAGE_DYTurbo=TRUE
echo " | Install directory:  $INSTALLDIR"
echo " | "
echo " +-----------------------------------------------------------------------------------------------+"


read -p "Do you want to continue? [y/n]: " choice
case "$choice" in
  [yY]) echo -e "Continuing...\n";;
  [nN]) echo "Exiting..."; exit 1;;
  *) echo "Invalid input. Please enter y or n."; exit 1;;
esac

mkdir -p $INSTALLDIR

### debugs:
# xfitterbranch= #fastNLO-v2.6        # default: main [or master?]
# yamlver= #0.2.5                     # default: 0.2.5
# qcdnumver= #18-00-00                # default: 18-00-00
# applgridver= #1.6.36                # default: 1.6.36
# apfelxxver=                       # default: 4.8.0
# pineapplver= #"0.6.0-alpha.17"      # default: 0.6.0-alpha.17
# 


cd $CURRENTDIR
#####################################################################
# installing applgrid
#####################################################################
if [[ ! -n $applgridver ]]; then
    echo "Skipping installation of APPLgrid."
else
    # ----------------------------------------------------------------- #
    echo "Installing APPLGRID $applgridver..."
    # ----------------------------------------------------------------- #
    cd $CURRENTDIR
    APPLGRID_URL=https://applgrid.hepforge.org/downloads/applgrid-"$applgridver".tgz
    APPLlog=$CURRENTDIR/install_applgrid.log
    # --- download and upack
    wget --no-check-certificate $APPLGRID_URL >> $APPLlog 2>&1  || { echo "Error. Fetching of APPLgrid failed. Check $APPLlog for details"; exit 1; }
    tar xfz applgrid-$applgridver.tgz  >> $APPLlog  2>&1   || { echo "Error. Unpacking of APPLgrid failed. Check $APPLlog for details"; exit 1; }
    # --- configure, compile and install
    cd applgrid-$applgridver
    # need to supply c++1X flag from root explicitly
    root_c_flags=`root-config --cflags`
    ./configure CXXFLAGS="${root_c_flags}" --prefix=$INSTALLDIR  >> $APPLlog  2>&1  || { echo "Error. Configure of APPLgrid failed. Check $APPLlog for details"; exit 1; }
    make -j 8  >> $APPLlog  2>&1  || { echo "Error. Compilation of APPLgrid failed. Check $APPLlog for details"; exit 1; }
    make install  >> $APPLlog  2>&1
    # --- return
    echo -e "APPLgrid installed successfully. See $APPLlog for installation details.\n"
    cd $CURRENTDIR
    # ----------------------------------------------------------------- #
fi 

cd $CURRENTDIR
#####################################################################
# installing pineappl
#####################################################################
if [[ ! -n $pineapplver ]]; then
    echo "Skipping installation of Pineappl."
else
    # ----------------------------------------------------------------- #
    echo "Installing PineAPPL $pineapplver..."
    # ----------------------------------------------------------------- #
    arch=x86_64
    rusttag="unknown-linux-gnu"
    pinelog=$CURRENTDIR/install_pineappl.log
    # ----------------------------------------------------------------- #
    wget https://github.com/NNPDF/pineappl/releases/download/v1.0.0/pineappl_capi-${arch}-${rusttag}.tar.gz  >> $pinelog 2>&1  || \
        { echo "Error. Fetching of Pineappl failed. Check $pinelog for details"; exit 1; }
    tar xzvf pineappl_capi-${arch}-${rusttag}.tar.gz -C $INSTALLDIR >> $pinelog 2>&1 || { echo "Error. Unpacking of Pineappl failed. Check $pinelog for details"; exit 1; }
    # --- return
    echo -e "Pineappl installed successfully. See $pinelog for installation details.\n"
    cd $CURRENTDIR
    # ----------------------------------------------------------------- #
fi


cd $CURRENTDIR
#####################################################################
# installing yaml
#####################################################################
if [[ ! -n $yamlver ]]; then
    echo "Skipping installation of yaml."
else
    # ----------------------------------------------------------------- #
    echo "Get yaml version $yamlver"
    # ----------------------------------------------------------------- #
    yamllog=install_yaml.log
    # --- download and upack
    wget http://pyyaml.org/download/libyaml/yaml-${yamlver}.tar.gz >> $yamllog  2>&1  || { echo "Error. Fetching of YAML failed. Check $yamllog for details"; exit 1; }
    tar xzvf yaml-${yamlver}.tar.gz >> $yamllog  2>&1
    cd yaml-${yamlver}  >> $yamllog
    # --- configure, compile and install ...
    ./configure --prefix=$INSTALLDIR >> $yamllog 2>&1 || { echo "Error. Configure of YAML failed. Check $yamllog for details"; exit 1; }
    make >> $yamllog 2>&1   || { echo "Error. Compilation of YAML failed. Check $yamllog for details"; exit 1; }
    make install >> $yamllog 2>&1
    # --- return
    echo -e "Yaml installed successfully. See $yamllog for installation details.\n"
    cd $CURRENTDIR
    # ----------------------------------------------------------------- #
fi



cd $CURRENTDIR
#####################################################################
# installing QCDNUM
#####################################################################
if [[ ! -n $qcdnumver ]]; then
    echo "Skipping installation of QCDNUM."
else
    # ----------------------------------------------------------------- #
    echo "Installing QCDNUM $qcdnumver..."
    # ----------------------------------------------------------------- #
    qcdnumstripver=`echo $qcdnumver |sed "s/-//g"`
    qcdnumlog=install_qcdnum.log
    # --- download and unpack
    wget https://gitlab.cern.ch/fitters/xfitter/-/raw/master/tools/qcdnum${qcdnumstripver}.tar.gz >> $qcdnumlog 2>&1
    tar xzvf qcdnum${qcdnumstripver}.tar.gz  >> $qcdnumlog  2>&1
    cd qcdnum-${qcdnumver}
    # --- configure
    ./configure --prefix=$INSTALLDIR  >> $qcdnumlog  2>&1  || { echo "Error. Configure of QCDNUM failed. Check $qcdnumlog for details"; exit 1; }
    # --- compile
    make -j 8 install  >> $qcdnumlog  2>&1  || { echo "Error. Compilation of QCDNUM failed. Check $qcdnumlog for details"; exit 1; }
    # --- return
    echo -e "QCDNUM installed successfully. See $qcdnumlog for installation details.\n"
    cd $CURRENTDIR
    # ----------------------------------------------------------------- #
fi



#####################################################################
# installing apfelxx
#####################################################################
if [[ ! -n $apfelxxver ]]; then
    echo "Skipping installation of Apfel++."
else
    # ----------------------------------------------------------------- #
    echo "Installing APFELxx $apfelxxver..."
    # ----------------------------------------------------------------- #
    apfelxxlog=$CURRENTDIR/install_apfelxx.log
    # ---- fetch Apfel++
    wget --no-check-certificate https://github.com/vbertone/apfelxx/archive/refs/tags/${apfelxxver}.tar.gz  >> $apfelxxlog 2>&1  || \
        { echo "Error. Fetching of APFELxx failed. Check $apfelxxlog for details"; exit 1; }
    mv ${apfelxxver}.tar.gz apfelxx-${apfelxxver}.tar.gz
    tar xfvz apfelxx-${apfelxxver}.tar.gz >> $apfelxxlog 2>&1 || { echo "Error. Unpacking of APFELxx failed. Check $apfelxxlog for details."; exit 1; }
    cd apfelxx-${apfelxxver}
    # ---- build Apfel++
    mkdir build; cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALLDIR >> $apfelxxlog  2>&1  || { echo "Error. Build of APFELxx failed. Check $apfelxxlog for details."; exit 1; }
    make -j 8 >> $apfelxxlog  2>&1  || { echo "Error. Build of APFELxx failed. Check $apfelxxlog for details."; exit 1; }
    # specific fix for rel 4.8.0
    #    echo "Info:  'error: 'transform' is not a member of 'std'' is related to a missing include of #include <algorithm> in ./src/kernel/matrix.cc."
    #    echo "       fix it, cd to apfelxx-4.8.0/build, and type 'make' and 'make install' again."
    make -j 8 install >> $apfelxxlog  2>&1
    # --- return
    echo -e "APFELxx installed successfully. See $apfelxxlog for installation details.\n"
    cd $CURRENTDIR
fi


cd $CURRENTDIR
#####################################################################
# installing DYTurbo
#####################################################################
if [[ ! -n $dyturbover ]]; then
    echo "Skipping installation of DYTurbo."
else
    # ----------------------------------------------------------------- #
    echo "Installing DYTurbo $dyturbover..."
    # ----------------------------------------------------------------- #
    dyturbolog=$CURRENTDIR/install_dyturbo.log
    wget --no-check-certificate https://dyturbo.hepforge.org/downloads/dyturbo-${dyturbover}.tar.gz >> $dyturbolog 2>&1 || \
        { echo "Error. Fetching of DYTurbo failed. Check $dyturbolog for details"; exit 1; }
    tar xfz dyturbo-${dyturbover}.tar.gz  >> $CURRENTDIR/install.log  2>&1 || { echo "Error. Unpacking of DYTurbo failed. Check $dyturbolog for details"; exit 1; }
    cd dyturbo-${dyturbover}
    #export CXXFLAGS="${CXXFLAGS} -std=c++0x"
    ./configure --enable-Ofast --disable-tlsopt --prefix=$INSTALLDIR >> $dyturbolog  2>&1 || { echo "Error. Configure of DYTurbo failed. Check $dyturbolog for details"; exit 1; }
    make -j 8 install  >> $dyturbolog  2>&1 || { echo "Error. Compilation of DYTurbo failed. Check $dyturbolog for details"; exit 1; }
    # --- return
    echo -e "DYTurbo installed successfully. See $dyturbolog for installation details.\n"
    cd $CURRENTDIR
fi



cd $CURRENTDIR
#####################################################################
# installing xfitter
#####################################################################
if [[ ! -n $xfitterbranch ]]; then
    echo "Skipping installation of xFitter."
else
    cd $CURRENTDIR
    git clone --branch $xfitterbranch --single-branch https://gitlab.cern.ch/fitters/xfitter.git xfitter_git
    mkdir xfitter-build; cd xfitter-build
    cmake ../xfitter_git  -DCMAKE_INSTALL_PREFIX=$INSTALLDIR -DCMAKE_DISABLE_FIND_PACKAGE_Ceres=TRUE -DCMAKE_DISABLE_FIND_PACKAGE_Ploughshare=TRUE -DCMAKE_DISABLE_FIND_PACKAGE_APFELgrid=TRUE -DCMAKE_DISABLE_FIND_PACKAGE_Hathor=TRUE -DCMAKE_DISABLE_FIND_PACKAGE_HOPPET=TRUE -DCMAKE_DISABLE_FIND_PACKAGE_DYTurbo=TRUE  || { echo "Error. CMake of xFitter failed. See above's log for details"; exit 1; }
    make -j8   || { echo "Error. Compilatoin of xFitter failed. See above's log for details"; exit 1; }
    make install
fi



cd $CURRENTDIR
#####################################################################
#  prepare a run directory
#####################################################################
if [[ ! -e run ]]; then
    mkdir -p run
    cp  xfitter_git/steering.txt     run
    cp  xfitter_git/parameters.yaml  run
    cp  xfitter_git/constants.yaml   run
    rsync -a --exclude=".*" xfitter_git/datafiles run/.
else
    echo "\"run\" directory already exists, I won't touch it"
    echo ""
fi


#####################################################################
#  quick start guide
#####################################################################
if [[ ! -e quickstart-readme.txt ]]
then
    echo "for a quick start do:" >> quickstart-readme.txt
    echo >> quickstart-readme.txt
    echo "source setup-lcgenv-$PLATFORM-$LCG_VERSION.sh  #setup environment" >> quickstart-readme.txt
    echo "cd run          #enter your run directory" >> quickstart-readme.txt
    echo "xfitter          #run your fit" >> quickstart-readme.txt
    echo >> quickstart-readme.txt
    echo ""  >> quickstart-readme.txt
    echo "Installation directory of all packages: $INSTALLDIR"  >> quickstart-readme.txt
    echo "Path of executables:                    $INSTALLDIR/bin"  >> quickstart-readme.txt
    echo ""  >> quickstart-readme.txt
    echo "The fit options are controlled by the 3 files steering.txt, parameters.yaml, and constants.yaml" >> quickstart-readme.txt
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

