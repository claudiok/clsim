#!/bin/sh
set -e

# helper script to automatically download and extract safeprimes_base32.txt
# (it takes a long time to build this file..)

BASE_URL="http://www.icecube.wisc.edu/~ckopper"
MD5_SUM_EXPECTED="295ae93631b7ff627ae42c1a4b2a7d75"

# require exactly one argument
if [ $# -ne 1 ]; then
	echo "Usage: `basename $0` {resources_dir}"
	exit 1
fi
RESOURCES_DIR=$1

# if the file exists, exit immediately
if [ -f $RESOURCES_DIR/safeprimes_base32.txt ]; then
	# if possible check the md5 checksum to see if the file is ok
	# (if not, just assume it is ok and exit)
	
	hash md5 &> /dev/null
	if [ $? -eq 1 ]; then
		# no md5 tool available
		exit 0
	fi
	
	MD5_SUM_FILE=`md5 -q $RESOURCES_DIR/safeprimes_base32.txt`
	
	if [ "$MD5_SUM_FILE" != "$MD5_SUM_EXPECTED" ]; then
		echo "file $RESOURCES_DIR/safeprimes_base32.txt does exist, but checksum is not okay! removing and downloading again."
		echo "     found: $MD5_SUM_FILE"
		echo "  expected: $MD5_SUM_EXPECTED"
		
		rm $RESOURCES_DIR/safeprimes_base32.txt
	else
		# checksum is OK!
		exit 0
	fi
fi 


echo "trying to download $RESOURCES_DIR/safeprimes_base32.txt..."

set +e
hash curl &> /dev/null
set -e
if [ $? -eq 1 ]; then
	# no curl, maybe wget?
	
	hash wget &> /dev/null
	if [ $? -eq 1 ]; then
		# nothing found
	    echo >&2 "  this script requires either curl or wget to be installed on your system"
		exit 2
	else
		USE_CURL=0 # implies: use wget
	fi
else
	USE_CURL=1
fi

set +e
hash unxz &> /dev/null
set -e
if [ $? -eq 1 ]; then
	# xz not found, use .gz file
	# (just assume that gunzip is available,
	# if it is not, this system won't be
	# usable for IceTray anyway.)

	USE_XZ=0
	FILENAME="safeprimes_base32.txt.gz"
else
	USE_XZ=1
	FILENAME="safeprimes_base32.txt.xz"
fi

# download
echo "  downloading $BASE_URL/$FILENAME ..."
if [ "$USE_CURL" == "1" ]; then
	curl -o $RESOURCES_DIR/$FILENAME $BASE_URL/$FILENAME
else
	wget -O $RESOURCES_DIR/$FILENAME $BASE_URL/$FILENAME
fi

# and unzip
echo "  extracting $BASE_URL/$FILENAME ..."
if [ "$USE_XZ" == "1" ]; then
	unxz $RESOURCES_DIR/$FILENAME
else
	gunzip $RESOURCES_DIR/$FILENAME
fi

echo "  $RESOURCES_DIR/safeprimes_base32.txt downloaded."
