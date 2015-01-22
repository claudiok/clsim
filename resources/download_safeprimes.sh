#!/bin/bash
set -e

## DO NOT COPY THIS SCRIPT INTO YOUR PROJECT
## Ask on the dataclass e-mail list how to do this

# helper script to automatically download and extract safeprimes_base32.gz
# (it takes a long time to build this file..)

BASE_URL="http://code.icecube.wisc.edu/tools/clsim"
MD5_SUM_EXPECTED="911bfc8f19be55a76c5d0ac30c773d0a"

# require exactly one argument
if [ $# -ne 1 ]; then
    echo "Usage: `basename $0` {resources_dir}"
    exit 1
fi
RESOURCES_DIR=$1

# if the file exists, exit immediately
if [ -f $RESOURCES_DIR/safeprimes_base32.gz ]; then
    # if possible check the md5 checksum to see if the file is ok
    # (if not, just assume it is ok and exit)
    
    FOUND_MD5=0
    UNAME=`uname`
    md5=md5
    if [ x"$UNAME" = "xDarwin" ]; then
	md5=/sbin/md5
    fi
    command -v $md5 >/dev/null && FOUND_MD5=1 || FOUND_MD5=0
    if [ $FOUND_MD5 -eq 1 ]; then
        MD5_SUM_FILE=`$md5 -q $RESOURCES_DIR/safeprimes_base32.gz`
    else
        # exit if there is no md5 or md5sum tool available
        command -v md5sum >/dev/null || { exit 0; }
        md5=`md5sum $RESOURCES_DIR/safeprimes_base32.gz`
        MD5_SUM_FILE="${md5%% *}" # remove the first space and everything after it
    fi
    
    if [ "$MD5_SUM_FILE" != "$MD5_SUM_EXPECTED" ]; then
        echo "file $RESOURCES_DIR/safeprimes_base32.gz does exist, but checksum is not okay! removing and downloading again."
        echo "     found: $MD5_SUM_FILE"
        echo "  expected: $MD5_SUM_EXPECTED"
        
        rm $RESOURCES_DIR/safeprimes_base32.gz
    else
        # checksum is OK!
        exit 0
    fi
fi 


echo "trying to download $RESOURCES_DIR/safeprimes_base32.gz..."

set +e
command -v curl >/dev/null && USE_CURL=1 || USE_CURL=0
if [ $USE_CURL -eq 0 ]; then
    # no curl, maybe wget?
    command -v wget >/dev/null || { echo >&2 "  this script requires either curl or wget to be installed on your system"; exit 2; }
    # USE_CURL=0 imples wget
fi
set -e

FILENAME="safeprimes_base32.gz"

# download
echo "  downloading $BASE_URL/$FILENAME ..."
if [ $USE_CURL -eq 1 ]; then
    curl -L -o $RESOURCES_DIR/$FILENAME $BASE_URL/$FILENAME
else
    wget -O $RESOURCES_DIR/$FILENAME $BASE_URL/$FILENAME
fi

echo "  $RESOURCES_DIR/safeprimes_base32.gz downloaded."
