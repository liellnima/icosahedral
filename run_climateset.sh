#!/bin/bash

# script to transform climateset data into icosahedral grid. assumes annual data.
# stores the new data in a seperate folder given by user

# First argument: folder with all subdirs, where each dir is a seperate year. assumes annual data.
# Second argument: folder where the new grib files should be stored
# Third argument: how the climate variable is named in the data file(s)
# Fourth argument: how the climate variable needs to be renamed for GRIB2 conventions

# Fifth argument: NI parameter (distance between hexagon grid cells in km)
# Sixth argument: type of remapping (default: remapcon2)
INPUT_DIR=$1
STORE_DIR=$2
VARNC=${3:-"tas"}
VARGRIB=${4:-"t"}
NI=${5:-24}
REMAP=${6:-"remapcon"}

echo "Remap files in $INPUT_DIR with NI=$NI and remapping function $REMAP..."

echo $INPUT_DIR
echo $STORE_DIR

# make new dir if it does not exist yet
# create subdir icosahedral and subdir longlat
mkdir -p $STORE_DIR/icosahedral
mkdir -p $STORE_DIR/longlat

# iterate through all dirs & files in INPUT_DIR
for dir in $INPUT_DIR/*/; do

  if [[ -d $dir ]]; then

    # create year as subdir in STORE_DIR
    dir=${dir%*/}
    year=${dir##*/}
    mkdir -p $STORE_DIR/icosahedral/$year
    mkdir -p $STORE_DIR/longlat/$year

    # run run.sh on each file
    for file in $dir/*.nc; do
  		file_path=`readlink -f $file`
  		./remap2icosahedral.sh $file_path $NI $REMAP 1 $VARNC $VARGRIB "grib"
  	done

    # Move grib files: "_icosahedral.grib" & "_longlat.grib"
    mv $dir/*_icosahedral.grib $STORE_DIR/icosahedral/$year
    mv $dir/*_longlat.grib $STORE_DIR/longlat/$year

  else
  	echo "$dir is not valid"
  	exit 1
  fi

done

# if needed, create vertex to lonlat mapping:
# (as long as 6250 steps -> use the one that is already uploaded)
# ./getlonlat.sh [SPECIFIC GRIB FILE]

echo "... finished. Stored in $STORE_DIR."
