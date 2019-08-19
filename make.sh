#Wrapper script for cmake

#./make.sh clean     - delete all build files
#./make.sh uninstall - delete all installed files
#./make.sh build     - configure and compile
#./make.sh           - same (configure and compile)
#./make.sh install   - configure, compile, and install
#./make.sh run       - configure, compile, install, and run

#Default install is in source:
#executable               goes to ./bin
#main library             goes to ./lib
#dynamically loaded modules go to ./lib/xfitter

CMAKE_OPTIONS="-DCMAKE_BUILD_TYPE=Debug"

#Uncommect to disable some some of the optional packages
#CMAKE_OPTIONS=$CMAKE_OPTIONS" -DCMAKE_DISABLE_FIND_PACKAGE_APFEL=TRUE"
#CMAKE_OPTIONS=$CMAKE_OPTIONS" -DCMAKE_DISABLE_FIND_PACKAGE_LHAPDF=TRUE"

SOURCE_DIR=$(dirname $(readlink -e $0)) #absolute path to directory of this script
BUILD_DIR=$SOURCE_DIR/build
INSTALL_DIR=$SOURCE_DIR

if [ "$1" == "clean" ];then
  #Delete the build directory
  if [ -d $BUILD_DIR ];then
		rm -r $BUILD_DIR
		rmdir --ignore-fail-on-non-empty $INSTALL_DIR
	fi
elif [ "$1" == "uninstall" ];then
  #Delete the installed exeuctable and libraries
  if [ $INSTALL_DIR == $SOURCE_DIR ];then #if xFitter is installed in-source
    #then only remove lib and bin
    rm -r lib bin share
  else
    #else remove the whole install dir
    rm -r $INSTALL_DIR
  fi
elif [ "$1" == "install" ] || [ "$1" == "run" ] || [ -z "$1" ];then
  #make sure build directory exists
  mkdir -p $BUILD_DIR
  #cd to build directory and invoke cmake theere
  cd $BUILD_DIR
  if [ ! -f Makefile ];then
    cmake $CMAKE_OPTIONS $SOURCE_DIR -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR
  fi
  if [ "$1" == "install" ] || [ "$1" == "run" ];then
    make -j$(nproc) install
  else
    make -j$(nproc)
  fi
  if [ "$?" -eq 0 ] && [ "$1" == "run" ];then
    cd $SOURCE_DIR #cd to where steering files are
    LD_LIBRARY_PATH=$INSTALL_DIR/lib/:$LD_LIBRARY_PATH $INSTALL_DIR/bin/xfitter
  fi
else
  echo "Unknown command \"$1\""
fi
