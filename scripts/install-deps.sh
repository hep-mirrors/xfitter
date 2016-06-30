#!/bin/bash
#####################################################################
# Installer scripts for xFitter dependencies 

## Configuration ####################################################

## Programs versions
lhapdfver=6.1.5
hoppetver=1.1.5
applgridver=1.4.70
qcdnumver=17-01-12
apfelver=2.7.0
melaver=2.0.1

#. scripts/setup.sh

if [ -e $_LOCALINST ] ; then
  echo "Dependencies taken from cache, skip installation"
  exit
fi

#export BOOST="--with-boost=/afs/cern.ch/sw/lcg/external/Boost/1.48.0_python2.6/x86_64-slc6-gcc46-opt"

#Make all dependencies
rm -rf $_LOCALINST >& /dev/null
_LOG="$_CURRENTDIR/install.log"
rm -rf $_LOG
mkdir $_LOCALINST && cd $_LOCALINST
pwd

#lhapdf:

echo "Installing LHAPDF $lhapdfver..."
if (( `echo $lhapdfver |cut -d. -f1` >= 6 ))
then
    lhapdf="LHAPDF"
#    withboost=$BOOST
else
    lhapdf="lhapdf"
fi
wget https://www.hepforge.org/archive/lhapdf/${lhapdf}-${lhapdfver}.tar.gz >> $_LOG 2>&1
tar xzf ${lhapdf}-${lhapdfver}.tar.gz  >> $_LOG 2>&1
cd ${lhapdf}-${lhapdfver} 
./configure --prefix=$_LOCALINST >> $_LOG  2>&1
_check_result $? "LHAPDF configure" $_LOG
make -j 9 install >> $_LOG  2>&1
_check_result $? "LHAPDF compile" $_LOG
cd - >& /dev/null

#hoppet:
echo "Installing HOPPET $hoppetver..."
wget http://hoppet.hepforge.org/downloads/hoppet-${hoppetver}.tgz >> $_LOG 2>&1
tar xzf hoppet-${hoppetver}.tgz  >> $_LOG  2>&1
cd hoppet-${hoppetver}
./configure --prefix=$_LOCALINST  >> $_LOG  2>&1
_check_result $? "HOPPET configure" $_LOG
make -j 9 install  >> $_LOG  2>&1
_check_result $? "HOPPET compile" $_LOG
cd - >& /dev/null

#applgrid:
echo "Installing APPLGRID $applgridver..."
wget https://www.hepforge.org/archive/applgrid/applgrid-$applgridver.tgz >> $_LOG 2>&1
tar xzf applgrid-$applgridver.tgz  >> $_LOG  2>&1
cd applgrid-$applgridver
# adjust gfortran library
sed -e '/^FRTLIB/s/ =.*$/ = /' -e '/^FRTLLIB/s/ -L.*)//' -i ./src/Makefile.in
./configure --prefix=$_LOCALINST  >> $_LOG  2>&1
_check_result $? "APPLgrid configure" $_LOG
make install >> $_LOG  2>&1
_check_result $? "APPLgrid compile" $_LOG
cd - >& /dev/null

#apfel
echo "Installing APFEL $apfelver..."
#git clone https://github.com/scarrazza/apfel.git >>$_LOG 2>&1
wget https://github.com/scarrazza/apfel/archive/${apfelver}.tar.gz >>$_LOG 2>&1
tar xzf ${apfelver}.tar.gz >>$_LOG 2>&1
cd  apfel-${apfelver}
#git checkout tags/${apfelver} >>$_LOG 2>&1
./configure --prefix=$_LOCALINST  >> $_LOG  2>&1
_check_result $? "APFEL configure" $_LOG
make -j 9 install  >> $_LOG  2>&1
_check_result $? "APFEL compile" $_LOG
cd - >& /dev/null

#mela
echo "Installing MELA $melaver..."
wget https://github.com/vbertone/MELA/archive/${melaver}.tar.gz >>$_LOG 2>&1
#git clone https://github.com/vbertone/MELA.git >>$_LOG 2>&1
tar xzf ${melaver}.tar.gz  >>$_LOG 2>&1
#mv MELA mela-${melaver}
cd MELA-${melaver}
#git checkout tags/${melaver} >>$_LOG 2>&1
#echo "Installing MELA $melaver..."
#melastripver=`echo $melaver |sed "s/-//g"`
#if [[ $http == "curl" ]]
#then
#curl http://apfel.hepforge.org/downloads/mela-${melastripver}.tar.gz > mela-${melastripver}.tar.gz 2>> $_LOG
#else
#wget http://apfel.hepforge.org/downloads/mela-${melastripver}.tar.gz >> $_LOG 2>&1
#fi
#tar xfz mela-${melastripver}.tar.gz  >> $_LOG  2>&1
#cd  mela-${melaver}
./configure --prefix=$_LOCALINST  >> $_LOG  2>&1
_check_result $? "MELA configure" $_LOG
make -j 9 install  >> $_LOG  2>&1
_check_result $? "MELA compile" $_LOG
cd - >& /dev/null

#qcdnum
echo "Installing QCDNUM $qcdnumver..."
qcdnumstripver=`echo $qcdnumver |sed "s/-//g"`
wget http://www.nikhef.nl/user/h24/qcdnum-files/download/qcdnum${qcdnumstripver}.tar.gz >> $_LOG 2>&1
tar xzf qcdnum${qcdnumstripver}.tar.gz  >> $_LOG  2>&1
cd qcdnum-${qcdnumver}
./configure --prefix=$_LOCALINST  >> $_LOG  2>&1
_check_result $? "QCDNUM configure" $_LOG
make -j 9 install  >> $_LOG  2>&1
_check_result $? "QCDNUM compile" $_LOG

cd - >& /dev/null

