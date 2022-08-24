#!/bin/bash
set -e

# Fix issue with singularity not honoring WORKDIR in docker
if [ -d /run ]; then cd /run; fi

# Setup xFitter
source ${XFITTER_INSTALL_DIR}/setup.sh

# Allow people to add LHAPDF data files
if [[ -d /pdfdata ]]; then
    if [[ -f /pdfdata/lhapdf.conf && -f /pdfdata/pdfsets.index ]]; then
        echo "Found LHAPDF data in /pdfdata, telling LHAPDF to use this directory."
        echo "Note: lhapdf-config --datadir will not be updated to this path. Don't panic."
        export LHAPDF_DATA_PATH=/pdfdata
    else
        echo "Invalid PDF folder found at /pdfdata. Please check your bindings."
        sleep 3
    fi
fi

# Add a symlink to the /data directory to prevent the need of modifying steering files
createdSymlink=0
if [[ -d /data && ! -d /run/datafiles ]]; then
    echo "Linking /data to /run/datafiles."
    ln -s /data /run/datafiles && createdSymlink=1
fi

# if the first argument starts with - we assume it is a xfitter-draw argument
if [[ $* == -* ]]; then
    xfitter
    xfitter-draw "$@"
else
    # Otherwise we run the command given
    # Note: the default command is 'xfitter && xfitter-draw output' from the dockerfile 
    eval $@
fi

# Clean up symlink if we made one
if [[ createdSymlink -eq 1 ]]; then
    unlink /run/datafiles
fi
