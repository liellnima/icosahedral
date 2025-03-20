# Remapper: Longitude/Latitude to Icosahedral-Hexagonal Grid

Bashscript to convert mean sea level pressure files downloaded from [NOAA](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html).
The longitude / latitude netCDF4 files are remapped to an icosahedral-hexagonal grid (GRIB2 files).

If desired the data can be also averaged over a certain amount of days (e.g. weekly average).

Check out `load.py` for an example of how to load the data in python.

# Installation and Setup
## Command line tools

The following command line tools must be installed:

- **cdo**: Climate Data Operators version 2.0.4 (https://mpimet.mpg.de/cdo)
- **ncrename**: netCDF Renamer version 5.0.6 (https://linux.die.net/man/1/ncrename)

I recommend installting cdo not from source, but by using conda and install cdo in your environment, see e.g. [this conda forge cdo installation](https://anaconda.org/conda-forge/cdo).

`conda install -c conda-forge eccodes`

## Data
Download the data (e.g. mean level sea pressure) from [NOAA](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html) and put it in the directory `raw/`.

### Download many files
If you want to download several files, e.g. all annual files of Air Temperature 2m above surface, you can first move into the appropriate directory (`raw/`) and use wget here:

```
wget -r -nd --no-parent -A 'air.2m.gauss.*.nc' https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Dailies/surface_gauss/
```

## Optional: Python packages
OPTIONAL: In case you want to run `load.py`:

1. Create environment: `conda create --name icosahedral python=3.9`
2. Activate environment: `conda activate icosahedral`
3. Install requirements: `pip install python_requirements.txt`
4. Run code snippet: `python load.py`

In case you run into problems with `eccodes`, please try:
`conda install -c conda-forge eccodes`


# Usage

```
./run.sh $INPUT $MON $VARNC $VARGRIB $NI $REMAP_FUNC $MEAN_DAYS
```

`$INPUT` can be both a file (e.g. `slp.1948.nc`) or a directory (e.g. `raw`) that contains a bunch of netcdf files.

`$MON` is a string flag `mon` indicating if you have monthly data summed up into one file. Default is `$MON=nan`. If set to `mon`, the data is split up into annual files. `$INPUT` has to be one specific file in this case. Example command: `./run.sh /raw/air.mon.mean.nc mon`.

`$VARNC` is the name of the variable you downloaded. The name is indicated in the file names. The default is `slp`.

`$VARGRIB` is how the variable should be renamed in order to be GRIB2 conform. See also section 'Rename Variables' to find out how you need to rename your original variable.

`$NI` is used to set the distance between the hexagonal grid cells (in km). For more information, refer to section "Setting the NI parameter". 


For further explanations on NI and the GME, refer e.g. to [Majewski et al](https://www.ecmwf.int/sites/default/files/elibrary/2000/10942-global-icosahedral-hexagonal-grid-point-model-gme-operational-version-and-high-resolution.pdf).

`$REMAP_FUNC` is the function that should be used for remapping. Default is `remapcon` (conservative first order remapping). The following options are possible: "remapbil", "remapbic", "remapnn", "remapdis", "remapcon", "remapcon2".

`$MEAN_DAYS` is the number of days over which the data should be averaged. Default is `$MEAN_DAYS=7`, i.e. a weekly average is calculated. Put `1` if you do not want that any average is calculated (or drop this argument).



The script automatically creates a mapping from the hexgonal grid cells to longitude/latitude values. The resulting mapping can be found in the folder `mappings`. This assumes that all data files had the same resolution and are operating on the exact same hexagonal grid.

## Split file
If you have one file that collects monthly data for a range of years (1948 - 2023), you can run the following command to get a seperate file for each year:

```
cdo splityear $VARIABLE.mon.mean.nc $VARIABLE.mon.mean.
```
With `$VARIABLE` being for example 'air'.

## Rename variables
If you are using different variables from the NOAA NCEP reanlysis data, you need to adapt the script a little bit to your files and new variable names. Specifically, you will need to rename your variable (e.g. 'temp') to a GRIB2 conform variable naming (e.g. 't'). To find out how you have to rename your variables, check out the [ECMWF GRIB Parameter Database](https://codes.ecmwf.int/grib/param-db/).


## Setting the NI parameter

NI "specifies the number of intervals on a main triangle side". See to the picture below for a better understanding. When building an icosahedral grid, we start with the "original triangular grid", i.e. we just split up the globe into 20 equal triangles. This is NI=1, i.e. the number of intervals on a major triangle edge is 1. The length of a major triangle edge for the Earth is ~7054 km. This triangle can now be split up into a finer and finer grid. Depending on what resolution you want to have, you can decide how many times the original grid should be split up.


Here, the **default** is `$NI=25`, leading to a spacing of ~280km between the grid cells, i.e. this is appropriate for datasets with a 2.5 degree resolution. With an NI=25, we get `N=6760` icosahedral grid cells. This follows the formula `N = 10 * (NI)^2`. 


If you have a dataset with a different resolution, you need to adapt your NI, if the icosahedral grid is expected to reflect the same resolution. You can **calculate your own NI** parameter via: `NI = 7054 km / $YOUR_RESOLUTION km`. The number 7054 refers to the length of a side of the "original icosahedral triangle" (the ones that are split over and over again). "Your resolution" is the resolution of your original data files. For example, for ERA5 data, with a resolution of 0.25 degree, i.e. around 27.83 km, we would expect `NI = 7054 km / 27.83 km`, so ~253.


You can **calculate your resolution**: Simply divide the circumference of the Earth (40075 km) by the number of longitude grid cells you have. For ERA5 data, you have 1440 longitude grid cells, i.e. your resolution is: `RES = 40075 km / 1440 = 27.83 km`. I.e. your resolution (at the equator) is 27.83 km.


Oftentimes you will run into **missing values** in your output. In that case adjust your NI to the next power of 2 (e.g. if you have NI=253, use 256 instead), or make you sure you can divide it by 16 or 24. If you are reading this and can explain me why that helps, please email me or open an issue!


**Step-by-step Guide**:
1. Look up the number of grid cells along the longitude (e.g. (144 x 73) grid --> 144 is the number you are looking for):
``cdo sinfo $YOUR_FILE``
2. Calculate the resolution at the circumference of Earth
``RES = 40075 / $YOUR_LON_GRID_CELLS``
3. Calculate the appropriate NI to maintain the same resolution
``NI = 7054 / RES``
4. Take the integer of NI (round down)
5. Do your remapping 
``cdo remapcon,gme$NI $INPUT_FILE $OUTPUT_FILE``
6. Check whether the remapping worked or if you got missing values (Miss : 0) with
``cdo -s -infov $YOUR_REMAPPED_FILE``
7. If you have missing values, try an NI that can be divided by 16 or 24.



# Notes
The attributes of the GRIB2 files will automatically say that the variable in question is "surface pressure" (sp). This is not true - this is mean sea level pressure (msl). The short-name "msl" is GRIB2-conform, but for some reason the variable name is not accepted and instead read as surface pressure. Please refer to the NOAA website linked above or the raw data files to retrieve the correct attributes.

The file `grid_gme.cc` was uploaded for reference in case someone needs to look into the orginal code of how the data is projected to the grid. The file is not used within the code and stems from the cdo package.
