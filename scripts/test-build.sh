#!/bin/bash

# This script tests building the xFitter code in different configurations

print_usage() {
    echo "USAGE: ./test-build.sh [options]"
    echo "where [options] is a list of configure modules, automatically preceeded by -enable "
    echo "Example ./test-build.sh applgrid lhapdf"
    echo "will run ./configure --enable-applgrid --enable-lhapdf"
    echo "Empty list build default configuration"
    echo "To disable an option, use: ./test-build.sh root=no"
}

####
if (grep -i help <<< $@); then
    print_usage
    return
fi

_DATE=$(date +%F.%H%I)
_LOG="buildlog_${_DATE}.log"
echo ${_DATE} ${_LOG}
touch ${_LOG}

# combine config string
conf_args=""
for arg in $@; do
  conf_args=$conf_args" --enable-"$arg
done

# autotools function
function run_autotools () {
  echo_separator "AUTOTOOLS TEST"
  autoreconf --force --install >> ${_LOG} 2>&1
  _check_result $? "Autotools test" ${_LOG}
}

# configure function
function run_configure () {
  echo_separator "CONFIGURE TEST"
  ./configure --prefix=$(pwd) $@ >> ${_LOG} 2>&1
  _check_result $? "Configure test" ${_LOG}
}

# make compile function
function run_make () {
  echo_separator "MAKE TEST"
  make $@ >> ${_LOG} 2>&1
  _check_result $? "Make test" ${_LOG}
}

# make install function
function run_make_install () {
  echo_separator "MAKE INSTALL TEST"
  make install >> ${_LOG} 2>&1
  _check_result $? "Make install test" ${_LOG}
}

# basic config
#if [[ $1 == "basic" ]]; then
  make maintainer-clean
  run_autotools
  echo "Running ./configure $conf_args"
  run_configure $conf_args
  run_make
  run_make_install
  if [[ ${_RESULT} -eq 0 ]]; then
    echo "Basic compile test OK."
  else
    echo "Basic compile test FAIL."
  fi
#fi

exit
