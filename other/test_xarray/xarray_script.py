import numpy as np
import pandas as pd
import xarray as xr
import time
import pickle


def my_fun():
	
	# NetCDF file
	nc_file = "CMIP6_CanESM5_Omon_piControl_r1i1p1f1_zos_6000-6199.nc"
	
	# Open the netCDF file
	with xr.open_dataset(nc_file, chunks={"lat":30, "lon":30}) as ds:
	#with xr.open_dataset(nc_file) as ds:
		
		# Calculate the "average time" 
		time_idx = xr.DataArray(np.arange(len(ds.coords["time"]))+1.0, dims="time")
		
		# Calculate the slope over at each gridpoint
		num = (ds.zos - ds.zos.mean(dim="time")) * (time_idx - time_idx.mean(dim="time"))
		denom = (time_idx - time_idx.mean(dim="time")) ** 2
		slope = num.sum(dim="time") / denom.sum(dim="time")
		
		# Make an annual average ZOS value
		zos_annual = ds.zos.groupby("time.year").mean()
		
		# n years
		baseyear = 2015
		endyear = baseyear + len(zos_annual)
		projyears = xr.DataArray(np.arange(endyear - baseyear), dims="time")
		
		# Drift values
		drift = xr.Dataset({"drift": (("lat", "lon", "time"), np.around(slope * projyears * 1000, decimals=1))},
			coords={"lat": ds.coords["lat"], "lon": ds.coords["lon"], "time": projyears+baseyear})
		
		
		# Write the output file
		drift.to_netcdf("test_netcdf_file.nc", encoding={"drift": {"dtype": "i2", "scale_factor": 0.1, "zlib": True, "complevel":9}})
		#outfile = open("test_pickle.pkl", 'wb')
		#pickle.dump(drift.compute(), outfile, protocol=-1)
		#outfile.close()
	
	return(None)

if __name__ == "__main__":
	
	my_fun()
	
	exit()
	