
# set environment
export _CURRENTDIR=$(pwd)
export _LOCALINST="$_CURRENTDIR/deps"

#export $CURRENTDIR $LOCALINST

export PATH="$_LOCALINST/bin:$PATH"
export LD_LIBRARY_PATH="$_LOCALINST/lib:$LD_LIBRARY_PATH"

# set variables 
_SUCCESS=0
_FAIL=1
_RESULT=${_SUCCESS}
_DATE=$(date +%F_%H%I)

export _SUCCESS _FAIL _RESULT _DATE

_check_result () {

  [[ $1 -eq 0 ]] && _RESULT=${_SUCCESS} || _RESULT=${_FAIL}

  if [[ ${_RESULT} -eq 1 ]]; then
    echo "$2 fail. See log $3"
    cat $3
    exit 1
  else 
    echo "$2 - OK"
  fi
}

export -f _check_result

# echo separator function
echo_separator () {
  echo "${1}..."
  printf "=================================================================\n \
========================= $1 ========================\n \
=================================================================\n" >> ${_LOG}
}

export -f echo_separator
