# Remapper: Longitude/Latitude to Icosahedral-Hexagonal Grid

Bashscript to convert mean sea level pressure files downloaded from [NOAA](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html).
The longitude / latitude netCDF4 files are remapped to an icosahedral-hexagonal grid (GRIB2 files).

If desired the data can be also averaged over a certain amount of days (e.g. weekly average).

Check out `load.py` for an example of how to load the data in python.

# Usage

```
./run.sh $INPUT $NI $REMAP_FUNC $MEAN_DAYS
```

`$INPUT` can be both a file (e.g. `slp.1948.nc`) or a directory (e.g. `raw`) that contains a bunch of netcdf files.

`$NI` is the distance between the hexagonal grid cells in kilometers. The default is `NI=25` (25 km), leading to 6760 grid cells.

`$REMAP_FUNC` is the function that should be used for remapping. Default is `remapcon` (conservative first order remapping). The following options are possible: "remapbil", "remapbic", "remapnn", "remapdis", "remapcon", "remapcon2".

`$MEAN_DAYS` is the number of days over which the data should be averaged. Default is `$MEAN_DAYS=7`, i.e. a weekly average is calculated. Put `1` if you do not want that any average is calculated.

# Installation and Setup
## Command line tools

The following command line tools must be installed:

- **cdo**: Climate Data Operators version 2.0.4 (https://mpimet.mpg.de/cdo)
- **ncrename**: netCDF Renamer version 5.0.6 (https://linux.die.net/man/1/ncrename)

## Data
Download the data (e.g. mean level sea pressure) from [NOAA](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html) and put it in the directory `raw/`.

## Optional: Python packages
OPTIONAL: In case you want to run `load.py`:

1. Create environment: `conda create --name icosahedral python=3.9`
2. Activate environment: `conda activate icosahedral`
3. Install requirements: `pip install python_requirements.txt`
4. Run code snippet: `python load.py`

# Note
The attributes of the GRIB2 files will automatically say that the variable in question is "surface pressure" (sp). This is not true - this is mean sea level pressure (msl). The short-name "msl" is GRIB2-conform, but for some reason the variable name is not accepted and instead read as surface pressure. Please refer to the NOAA website linked above or the raw data files to retrieve the correct attributes.

