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

CURRENTDIR=$(pwd)
LOCALINST="$CURRENTDIR/deps"
#export BOOST="--with-boost=/afs/cern.ch/sw/lcg/external/Boost/1.48.0_python2.6/x86_64-slc6-gcc46-opt"

#. /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6-gcc46-opt/setup.sh

function _check_result () {
  if [[ $1 != 0 ]]; then
    echo "$2 error, check install.log for details"
    exit
  fi
}

#Make all dependencies
rm -rf $LOCALINST >& /dev/null
rm $CURRENTDIR/install.log
mkdir $LOCALINST && cd $LOCALINST

#lhapdf:

echo "Installing LHAPDF $lhapdfver..."
if (( `echo $lhapdfver |cut -d. -f1` >= 6 ))
then
    lhapdf="LHAPDF"
    withboost=$BOOST
else
    lhapdf="lhapdf"
fi
wget https://www.hepforge.org/archive/lhapdf/${lhapdf}-${lhapdfver}.tar.gz >> $CURRENTDIR/install.log 2>&1
tar xfz ${lhapdf}-${lhapdfver}.tar.gz  >> $CURRENTDIR/install.log 2>&1
cd ${lhapdf}-${lhapdfver} 
./configure --prefix=$LOCALINST >> $CURRENTDIR/install.log  2>&1
_check_result $? "LHAPDF configure"
make -j 9 install >> $CURRENTDIR/install.log  2>&1
_check_result $? "LHAPDF compile"
cd - >& /dev/null

#hoppet:
echo "Installing HOPPET $hoppetver..."
wget http://hoppet.hepforge.org/downloads/hoppet-${hoppetver}.tgz >> $CURRENTDIR/install.log 2>&1
tar xfz hoppet-${hoppetver}.tgz  >> $CURRENTDIR/install.log  2>&1
cd hoppet-${hoppetver}
./configure --prefix=$LOCALINST  >> $CURRENTDIR/install.log  2>&1
_check_result $? "HOPPET configure"
make -j 9 install  >> $CURRENTDIR/install.log  2>&1
_check_result $? "HOPPET compile"
cd - >& /dev/null

#applgrid:
echo "Installing APPLGRID $applgridver..."
wget https://www.hepforge.org/archive/applgrid/applgrid-$applgridver.tgz >> $CURRENTDIR/install.log 2>&1
tar xfz applgrid-$applgridver.tgz  >> $CURRENTDIR/install.log  2>&1
cd applgrid-$applgridver
./configure --prefix=$LOCALINST  >> $CURRENTDIR/install.log  2>&1
_check_result $? "APPLgrid configure"
make install >> $CURRENTDIR/install.log  2>&1
_check_result $? "APPLgrid compile"
cd - >& /dev/null

#apfel
echo "Installing APFEL $apfelver..."
git clone https://github.com/scarrazza/apfel.git >>$CURRENTDIR/install.log 2>&1
mv apfel apfel-${apfelver}
cd  apfel-${apfelver}
git checkout tags/${apfelver} >>$CURRENTDIR/install.log 2>&1
./configure --prefix=$LOCALINST  >> $CURRENTDIR/install.log  2>&1
_check_result $? "APFEL configure"
make -j 9 install  >> $CURRENTDIR/install.log  2>&1
_check_result $? "APFEL compile"
cd - >& /dev/null

#mela
echo "Installing MELA $melaver..."
#wget https://github.com/vbertone/MELA/archive/2.0.1.tar.gz >>$CURRENTDIR/install.log 2>&1
git clone https://github.com/vbertone/MELA.git >>$CURRENTDIR/install.log 2>&1
#tar xzf 2.0.1.tar.gz  >>$CURRENTDIR/install.log 2>&1
mv MELA mela-${melaver}
cd mela-${melaver}
git checkout tags/${melaver} >>$CURRENTDIR/install.log 2>&1
#echo "Installing MELA $melaver..."
#melastripver=`echo $melaver |sed "s/-//g"`
#if [[ $http == "curl" ]]
#then
#curl http://apfel.hepforge.org/downloads/mela-${melastripver}.tar.gz > mela-${melastripver}.tar.gz 2>> $CURRENTDIR/install.log
#else
#wget http://apfel.hepforge.org/downloads/mela-${melastripver}.tar.gz >> $CURRENTDIR/install.log 2>&1
#fi
#tar xfz mela-${melastripver}.tar.gz  >> $CURRENTDIR/install.log  2>&1
#cd  mela-${melaver}
./configure --prefix=$LOCALINST  >> $CURRENTDIR/install.log  2>&1
_check_result $? "MELA configure"
make -j 9 install  >> $CURRENTDIR/install.log  2>&1
_check_result $? "MELA compile"
cd - >& /dev/null

#qcdnum
echo "Installing QCDNUM $qcdnumver..."
qcdnumstripver=`echo $qcdnumver |sed "s/-//g"`
wget http://www.nikhef.nl/user/h24/qcdnum-files/download/qcdnum${qcdnumstripver}.tar.gz >> $CURRENTDIR/install.log 2>&1
tar xfz qcdnum${qcdnumstripver}.tar.gz  >> $CURRENTDIR/install.log  2>&1
cd qcdnum-${qcdnumver}
./configure --prefix=$LOCALINST  >> $CURRENTDIR/install.log  2>&1
_check_result $? "QCDNUM configure"
make -j 9 install  >> $CURRENTDIR/install.log  2>&1
_check_result $? "QCDNUM compile"

cd - >& /dev/null

export PATH="$LOCALINST/bin/:$PATH"
export LD_LIBRARY_PATH="$LOCALINST/lib/:$LD_LIBRARY_PATH"

