#!/bin/bash

# First argument: which file should be converted
# Second argument: NI parameter (distance between hexagon grid cells in km) 
# Third argument: type of remapping (default: remapcon2)
NI=${2:-24}
REMAP=${3:-"remapcon"}
DAYS=${4:-7}

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
	
echo "Finished."


