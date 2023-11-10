#!/bin/bash

# First argument: which file should be converted
# Second argument: indicate if you have one datafile with all the monthly data
# Third argument: how the climate variable is named in the data file(s)
# Fourth argument: how the climate variable needs to be renamed for GRIB2 conventions

# Fifth argument: NI parameter (distance between hexagon grid cells in km)
# Sixth argument: type of remapping (default: remapcon2)
# Seventh argument: how many days should be meaned together

MON=${2:-"nan"}
VARNC=${3:-"slp"}
VARGRIB=${4:-"msl"}
NI=${5:-24}
REMAP=${6:-"remapcon"}
DAYS=${7:-7}

FILES=$1
ROOT=${1/.nc}

echo $ROOT
echo $FILES

if [[ $MON == "mon" ]]; then
	cdo splityear $1 $ROOT
	FILES="$(dirname "$1")"
	DAYS=1
	rm $1
fi

echo "Remap $FILES with NI=$NI and remapping function $REMAP..."

if [[ -f $FILES ]]; then
	./remap2icosahedral.sh $FILES $NI $REMAP $DAYS $VARNC $VARGRIB
elif [[ -d $FILES ]]; then
	for file in $FILES/*.nc; do
		file_path=`readlink -f $file`
		./remap2icosahedral.sh $file_path $NI $REMAP $DAYS $VARNC $VARGRIB
	done
else
	echo "$FILES is not valid"
	exit 1
fi

# if needed, create vertex to lonlat mapping:
# (as long as 6250 steps -> use the one that is already uploaded)
# ./getlonlat.sh [SPECIFIC GRIB FILE] 

echo "Finished."
