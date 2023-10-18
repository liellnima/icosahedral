#!/bin/bash

# First argument: which file should be converted
# Second argument: NI parameter (distance between hexagon grid cells in km)
# Third argument: type of remapping (default: remapcon2)
# Fourth argument: how many days should be meaned together
# Fifth argument: how the climate variable is named in the data file(s)
# Sixth argument: how the climate variable needs to be renamed for GRIB2 conventions
NI=${2:-24}
REMAP=${3:-"remapcon"}
DAYS=${4:-1}
ROOT=${1/.nc}
VARNC=${5:"slp"}
VARGRIB=${6:"msl"}

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
rm $ROOT"_updated.nc"
rm $ROOT"_longlat.grib"

# average over given number of days
if [[ $DAYS > 1 ]]; then
	cdo timselmean,$DAYS $ROOT"_icosahedral.grib" $ROOT".grib"
	rm $ROOT"_icosahedral.grib"
fi
