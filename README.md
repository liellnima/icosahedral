# Remapper: Longitude/Latitude to Icosahedral-Hexagonal Grid

Bashscript to convert mean sea level pressure files downloaded from [NOAA](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html).
The longitude / latitude netCDF4 files are remapped to an icosahedral-hexagonal grid (GRIB2 files).

If desired the data can be also averaged over a certain amount of days (e.g. weekly average).

Check out `load.py` for an example of how to load the data in python.

# Usage

```
./run.sh $INPUT $MON $VARNC $VARGRIB $NI $REMAP_FUNC $MEAN_DAYS
```

`$INPUT` can be both a file (e.g. `slp.1948.nc`) or a directory (e.g. `raw`) that contains a bunch of netcdf files.

`$MON` is a string flag `mon` indicating if you have monthly data summed up into one file. Default is `$MON=nan`. If set to `mon`, the data is split up into annual files. `$INPUT` has to be one specific file in this case. Example command: `./run.sh /raw/air.mon.mean.nc mon`.

`$VARNC` is the name of the variable you downloaded. The name is indicated in the file names. The default is `slp`.

`$VARGRIB` is how the variable should be renamed in order to be GRIB2 conform. See also section 'Rename Variables' to find out how you need to rename your original variable.

`$NI` is the distance between the hexagonal grid cells in kilometers. The default is `$NI=25` (25 km), leading to 6760 grid cells.

`$REMAP_FUNC` is the function that should be used for remapping. Default is `remapcon` (conservative first order remapping). The following options are possible: "remapbil", "remapbic", "remapnn", "remapdis", "remapcon", "remapcon2".

`$MEAN_DAYS` is the number of days over which the data should be averaged. Default is `$MEAN_DAYS=7`, i.e. a weekly average is calculated. Put `1` if you do not want that any average is calculated (or drop this argument).



The script automatically creates a mapping from the hexgonal grid cells to longitude/latitude values. The resulting mapping can be found in the folder `mappings`. This assumes that all data files had the same resolution and are operating on the exact same hexagonal grid.

# Installation and Setup
## Command line tools

The following command line tools must be installed:

- **cdo**: Climate Data Operators version 2.0.4 (https://mpimet.mpg.de/cdo)
- **ncrename**: netCDF Renamer version 5.0.6 (https://linux.die.net/man/1/ncrename)

## Data
Download the data (e.g. mean level sea pressure) from [NOAA](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html) and put it in the directory `raw/`.

### Download many files
If you want to download several files, e.g. all annual files of Air Temperature 2m above surface, you can first move into the appropiate directory (`raw/`) and use wget here:

```
wget -r -nd --no-parent -A 'air.2m.gauss.*.nc' https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Dailies/surface_gauss/
```

### Split file
If you have one file that collects monthly data for a range of years (1948 - 2023), you can run the following command to get a seperate file for each year:

```
cdo splityear $VARIABLE.mon.mean.nc $VARIABLE.mon.mean.
```
With `$VARIABLE` being for example 'air'.

## Rename variables
If you are using different variables from the NOAA NCEP reanlysis data, you need to adapt the script a little bit to your files and new variable names. Specifically, you will need to rename your variable (e.g. 'temp') to a GRIB2 conform variable naming (e.g. 't'). To find out how you have to rename your variables, check out the [ECMWF GRIB Parameter Database](https://codes.ecmwf.int/grib/param-db/).

## Optional: Python packages
OPTIONAL: In case you want to run `load.py`:

1. Create environment: `conda create --name icosahedral python=3.9`
2. Activate environment: `conda activate icosahedral`
3. Install requirements: `pip install python_requirements.txt`
4. Run code snippet: `python load.py`

In case you run into problems with `eccodes`, please try:
`conda install -c conda-forge eccodes`


# Note
The attributes of the GRIB2 files will automatically say that the variable in question is "surface pressure" (sp). This is not true - this is mean sea level pressure (msl). The short-name "msl" is GRIB2-conform, but for some reason the variable name is not accepted and instead read as surface pressure. Please refer to the NOAA website linked above or the raw data files to retrieve the correct attributes.

The file `grid_gme.cc` was uploaded for reference in case someone needs to look into the orginal code of how the data is projected to the grid. The file is not used within the code and stems from the cdo package.
