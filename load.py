import xarray as xr
from pathlib import Path

def load(file_path: Path) -> xr.Dataset:
	""" Load a GRIB2 file and return an xarray object.
	Args:
		file_path (Path): Path to the file that should be loaded.
	Returns:
		xr.Dataset: The loaded dataset
	"""
	return xr.load_dataset(file_path, engine="cfgrib")
	
def print_infos(ds: xr.Dataset):
	""" Shows how to access the data on the temporal, spatial and variable axis.
	"""
	print("\nGeneral information:")
	print(ds)
	
	print("\nDimensions:")
	print(ds.sizes)
	
	print("\nAll weeks of sea level pressure for the first grid cell:")
	print(ds["sp"][0:52, 0])
	
	print("\nThe first week for all grid cells:")
	print(ds["sp"][0, 0:6759])
	
	print("\nThe last week for the last grid cell:")
	print(ds["sp"][52, 6000:6759])
	
def has_nans(ds: xr.Dataset) -> bool:
	""" Returns True if dataset contains any nans, else returns False.
	Args:
		ds (xr.Dataset): the dataset to analyse
	Returns:
		bool: if dataset contains nans or not
	"""
	return ds["sp"].isnull().any()
	
	
def debug():
	""" Internal usage only
	"""
	# 1. Does averaging over the weeks create nans? -> checked it, existed before weekly average
	# 2. Does the original data contain any nans? -> nope, this has no nans...
	file_path = Path("./raw/slp.1948.nc")
	ds = xr.load_dataset(file_path, engine="netcdf4")
	contains_nan_list = ds["slp"].isnull()
	print(contains_nan_list.any())
	

def main():
	""" Running examples of loading the data. Making some sanity checks.
	"""
	file_path = Path("./icosahedral_weekly/slp.1948.grib")
	
	ds = load(file_path)
	print(has_nans(ds))
	
	print_infos(ds)
	


if __name__ == "__main__":
	main()
