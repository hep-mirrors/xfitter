#!/bin/bash

# This script tests building the xFitter code in different configurations

_LOG="buildlog_${_DATE}.log"
touch ${_LOG}

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

if [[ $1 -eq "basic" ]]; then
  run_autotools
  run_configure
  run_make
  run_make_install
  if [[ ${_RESULT} -eq 0 ]]; then
    echo "Basic compile test OK."
  else
    echo "Basic compile test FAIL."
  fi
fi
