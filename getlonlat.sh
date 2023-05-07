#!/bin/bash

# use this script only after having created the icosahedral files
# the script assumes that all the data files have the same resolution and grid
# has to be applied separately if you are operating with different grid sizes / resolutions

DATA_GME=${1:-"icosahedral_weekly/slp.1948.grib"}

# drop time-dimension and work with single timestep
cdo seltimestep,1 $DATA_GME mappings/gme_grid_example.grib

# create the vertex - lonlat mapping
cdo -s outputkey,xind,lon,lat mappings/gme_grid_example.grib > mappings/vertex_lonlat_mapping.txt

# remove grid examples
rm mappings/*.grib

echo "created mapping between GME vertices and longitude latitude values "
