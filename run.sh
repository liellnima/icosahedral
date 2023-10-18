#!/bin/bash

# First argument: which file should be converted
# Second argument: NI parameter (distance between hexagon grid cells in km)
# Third argument: type of remapping (default: remapcon2)
# Fourth argument: how many days should be meaned together
# Fifth argument: indicate if you have one datafile with all the monthly data
NI=${2:-24}
REMAP=${3:-"remapcon"}
DAYS=${4:-7}
MON=${5:-"false"}

ROOT=${1/.nc}

if [[ $MON == "true"]]; then
	cdo splityear $1 $ROOT
	DAYS=1
fi


echo "Remap $1 with NI=$NI and remapping function $REMAP..."

if [[ -f $1 ]]; then
	./remap2icosahedral.sh $1 $NI $REMAP $DAYS
elif [[ -d $1 ]]; then
	for file in $1/*.nc; do
		file_path=`readlink -f $file`
		./remap2icosahedral.sh $file_path $NI $REMAP $DAYS
	done
else
	echo "$1 is not valid"
	exit 1
fi

./getlonlat.sh

echo "Finished."
