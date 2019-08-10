cmake=$(which cmake)
SOURCE_DIR=$(dirname $(readlink -e $0)) #absolute path to directory of this script
BUILD_DIR=$SOURCE_DIR/build
INSTALL_DIR=$SOURCE_DIR
#INSTALL_DIR=$SOURCE_DIR
if [ "$1" == "reset" ];then
  if [ -d $BUILD_DIR ];then
		rm $(<$BUILD_DIR/install_manifest.txt)
		rm -r $BUILD_DIR
		rmdir --ignore-fail-on-non-empty $INSTALL_DIR
	fi
else
  mkdir -p $BUILD_DIR
  cd $BUILD_DIR
  if [ ! -f Makefile ];then
    CXXFLAGS=-I~/opt/usr/include $cmake $SOURCE_DIR -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR
  fi
  #make -j$(nproc) install
  make
  if [ "$1" == "run" ];then LD_LIBRARY_PATH=$INSTALL_DIR/lib/:$LD_LIBRARY_PATH $INSTALL_DIR/bin/xfitter;fi
fi
