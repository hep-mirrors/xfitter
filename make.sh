#!/bin/bash

#Wrapper script for cmake

#./make.sh clean     - delete all build files
#./make.sh uninstall - delete all installed files
#./make.sh build     - configure and compile
#./make.sh           - same (configure and compile)
#./make.sh install   - configure, compile, and install
#./make.sh run       - configure, compile, install, and run
#./make.sh reconfigure - configure from scratch (can be shortened up to ./make.sh rec)

#Default install is in source:
#executable               goes to ./bin
#main library             goes to ./lib
#dynamically loaded modules go to ./lib/xfitter

#CMAKE_FLAGS=$CMAKE_FLAGS" -DCMAKE_BUILD_TYPE=Release"
CMAKE_FLAGS=$CMAKE_FLAGS" -DCMAKE_BUILD_TYPE=Debug"

#Uncommect to disable some some of the optional packages
#CMAKE_FLAGS=$CMAKE_FLAGS" -DCMAKE_DISABLE_FIND_PACKAGE_APFEL=TRUE"
#CMAKE_FLAGS=$CMAKE_FLAGS" -DCMAKE_DISABLE_FIND_PACKAGE_APFELxx=TRUE"
#CMAKE_FLAGS=$CMAKE_FLAGS" -DCMAKE_DISABLE_FIND_PACKAGE_Ceres=TRUE"

SOURCE_DIR=`pwd` #absolute path to directory of this script
BUILD_DIR=$SOURCE_DIR/build
INSTALL_DIR=$SOURCE_DIR

cmd=$1

if [ "$(echo "$cmd"|cut -b1-3)" == "rec" ];then #reconfigure
  cmd=reconfigure
fi

if [ "$cmd" == "clean" ];then
  #Delete the build directory
  if [ -d $BUILD_DIR ];then
		rm -r $BUILD_DIR
		rmdir --ignore-fail-on-non-empty $INSTALL_DIR
	fi
elif [ "$cmd" == "uninstall" ];then
  #Delete the installed exeuctable and libraries
  if [ $INSTALL_DIR == $SOURCE_DIR ];then #if xFitter is installed in-source
    #then only remove lib and bin
    rm -r lib bin share
  else
    #else remove the whole install dir
    rm -r $INSTALL_DIR
  fi
elif [ "$cmd" == "reconfigure" ] || [ "$cmd" == "install" ] || [ "$cmd" == "run" ] || [ "$cmd" == "build" ] || [ -z "$cmd" ];then
  #make sure build directory exists
  mkdir -p $BUILD_DIR
  #cd to build directory and invoke cmake theere
  cd $BUILD_DIR
  if [ "$cmd" == "reconfigure" ] && [ -e CMakeCache.txt ];then
    rm CMakeCache.txt
  fi
  if [ ! -f Makefile ] || [ ! -f CMakeCache.txt ];then
    cmake $CMAKE_FLAGS  $SOURCE_DIR -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR || exit
  fi
  if [ "$cmd" == "reconfigure" ];then
    exit 0
  fi
  if [ "$cmd" == "install" ] || [ "$cmd" == "run" ];then
    make VERBOSE=yes -j$(nproc) install || exit
  else
    make VERBOSE=yes -j$(nproc) || exit
  fi
  if [ "$cmd" == "run" ];then
    cd $SOURCE_DIR #cd to where steering files are
    $INSTALL_DIR/bin/xfitter || exit
  fi
else
  echo "Unknown command \"$1\""
  exit 10
fi
