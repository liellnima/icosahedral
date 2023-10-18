#!/bin/bash

# First argument: which file should be converted
# Second argument: NI parameter (distance between hexagon grid cells in km)
# Third argument: type of remapping (default: remapcon2)
# Fourth argument: how many days should be meaned together
NI=${2:-24}
REMAP=${3:-"remapcon"}
DAYS=${4:-1}
ROOT=${1/.nc}

# convert netcdf file to grib file

# rename variable for mean sea level pressure (GRIB convention)
ncrename -h -O -v slp,msl $1 $ROOT"_renamed.nc"
# drop time_bnds variable (GRIB convention)
cdo select,name=msl $ROOT"_renamed.nc" $ROOT"_updated.nc"
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
