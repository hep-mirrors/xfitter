#!/bin/bash

# This script tests building the xFitter code in different configurations

_SUCCESS=0
_FAIL=1
_RESULT=${_SUCCESS}
_DATE=$(date +%F_%H%I)
_LOG="testlog_${_DATE}.log"
_ERR="testlog_${_DATE}.err"
touch ${_LOG} ${_ERR}

# echo separator function
function echo_separator () {
  echo "${1}..."
  printf "=================================================================\n \
========================= $1 ========================\n \
=================================================================\n" >> ${_LOG}
  printf "=================================================================\n \
========================= $1 ========================\n \
=================================================================\n" >> ${_ERR}
}

# autotools function
function run_autotools () {
  echo_separator "AUTOTOOLS TEST"
  autoreconf --force --install 1>> ${_LOG} 2>> ${_ERR}
  [[ $? -eq 0 ]] && _RESULT=${_SUCCESS} || _RESULT=${_FAIL}

  if [[ ${_RESULT} -eq 1 ]]; then
    echo "Autotools test fail. See error log $ERR"
    exit 1
  else 
    echo " - OK"
  fi
}

# configure function
function run_configure () {
  echo_separator "CONFIGURE TEST"
  ./configure --prefix=$(pwd) $@ 1>> ${_LOG} 2>> ${_ERR}
  [[ $? -eq 0 ]] && _RESULT=${_SUCCESS} || _RESULT=${_FAIL}

  if [[ ${_RESULT} -eq 1 ]]; then
    echo "configure test fail. See error log $ERR"
    exit 1
  else 
    echo " - OK"
  fi
}

# make compile function
function run_make () {
  echo_separator "MAKE TEST"
  make $@ 1>> ${_LOG} 2>> ${_ERR}
  [[ $? -eq 0 ]] && _RESULT=${_SUCCESS} || _RESULT=${_FAIL}

  if [[ ${_RESULT} -eq 1 ]]; then
    echo "Make compile test fail. See error log $ERR"
    exit 1
  else 
    echo " - OK"
  fi
}

# make install function
function run_make_install () {
  echo_separator "MAKE INSTALL TEST"
  make install 1>> ${_LOG} 2>> ${_ERR}
  [[ $? -eq 0 ]] && _RESULT=${_SUCCESS} || _RESULT=${_FAIL}

  if [[ ${_RESULT} -eq 1 ]]; then
    echo "Make install test fail. See error log $_ERR"
    exit 1
  else 
    echo " - OK"
  fi
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
