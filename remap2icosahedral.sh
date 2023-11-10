#!/bin/bash

# First argument: which file should be converted
# Second argument: NI parameter (distance between hexagon grid cells in km)
# Third argument: type of remapping (default: remapcon2)
# Fourth argument: how many days should be meaned together
# Fifth argument: how the climate variable is named in the data file(s)
# Sixth argument: how the climate variable needs to be renamed for GRIB2 conventions
# Seventh argument: if one of the intermediate files should be kept.
# 		"nan" to keep nothing (default)
#			"nc" to keep the renamed nc file (vars following grib convention)
#			"grib" to keep the pure grib files with the original data but new var naming

NI=${2:-24}
REMAP=${3:-"remapcon"}
DAYS=${4:-1}
ROOT=${1/.nc}
VARNC=${5:-"slp"}
VARGRIB=${6:-"msl"}
KEEP=${7:-"nan"}

# convert netcdf file to grib file

# rename variable for mean sea level pressure (GRIB convention)
ncrename -h -O -v $VARNC,$VARGRIB $1 $ROOT"_renamed.nc"
# drop time_bnds variable (GRIB convention)
cdo select,name=$VARGRIB $ROOT"_renamed.nc" $ROOT"_updated.nc"
# convert to GRIB2
cdo -f grb2 copy $ROOT"_updated.nc" $ROOT"_longlat.grib"

# remap to icosahedral-hexagon map
cdo $REMAP,gme$NI $ROOT"_longlat.grib" $ROOT"_icosahedral.grib"

# delete unnecessary files
rm $ROOT"_renamed.nc"
if [[ ! $KEEP == "nc" ]]; then
	rm $ROOT"_updated.nc"
fi
if [[ ! $KEEP == "grib" ]]; then
	rm $ROOT"_longlat.grib"
fi

# average over given number of days
if [[ $DAYS > 1 ]]; then
	cdo timselmean,$DAYS $ROOT"_icosahedral.grib" $ROOT".grib"
	rm $ROOT"_icosahedral.grib"
fi
