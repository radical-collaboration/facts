import numpy as np
import pandas as pd
import xarray as xr
from mpl_toolkits.basemap import Basemap
import glob
import os
import shutil
import re
import sys
import cartopy
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import fnmatch
#
import argparse
import fnmatch
from netCDF4 import Dataset
from pathlib import Path
import dask.array as da
import time
#
PD=os.getcwd(); #PD


# Function block

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# create folder with JnB dir.
def FOcr8(folder_name,subfolder_name=None,subsubfolder_name=None):
    # folder_path, subfolder_path , subsubfolder_path = FOcr8("000_full_sample_workflows","wf_1e","ssp585")
    # I/P: folder_name = '000_full_sample_components'
    folder_path = os.path.join(PD, folder_name)
    if os.path.exists(folder_path): shutil.rmtree(folder_path)
    os.mkdir(folder_path)
    # SUB Folder
    if subfolder_name is not None:
        subfolder_path = os.path.join(folder_path, subfolder_name)
        if os.path.exists(subfolder_path): shutil.rmtree(subfolder_path)
        os.mkdir(subfolder_path)
        # Sub Sub Folder
        if subsubfolder_name is not None:
            subsubfolder_path = os.path.join(subfolder_path, subsubfolder_name)
            if os.path.exists(subsubfolder_path): shutil.rmtree(subsubfolder_path)
            os.mkdir(subsubfolder_path)
            #
            return folder_path, subfolder_path, subsubfolder_path   
            #
        return folder_path, subfolder_path
    else:
        return folder_path
# ^^^


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def FOcr8_1(path, folder_name):
    folder_path = os.path.join(path, folder_name)
    if os.path.exists(folder_path): shutil.rmtree(folder_path)
    os.mkdir(folder_path)
    # Return the folder path
    return folder_path
# ^^^


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# deletes all files in a folder that match a particular pattern, except those that contain both pattern and pattern2:
def delete_files_with_pattern(folder, pattern, pattern2):
    for filename in os.listdir(folder):
        if fnmatch.fnmatch(filename, pattern):
            if fnmatch.fnmatch(filename, pattern2):
                continue  # Skip files matching the exclusion pattern
            file_path = os.path.join(folder, filename)
            os.remove(file_path)
            # print(f"Deleted file: {file_path}")
# ^^^

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Copy all files to the FOLDER that match a pattern
def cp(source_dir,destination_dir,pattern):
    # source_dir = expF
    # destination_dir = folder_path
    source_file_pattern = os.path.join(source_dir, pattern)
    matching_files = glob.glob(source_file_pattern)
    for file_path in matching_files:
        shutil.copy2(file_path, destination_dir)
# ^^^   
  

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Copy all files in a folder to another folder
def cp_dir2dir(srcDIR,dstnDIR):
    file_list = os.listdir(srcDIR)
    for file_name in file_list:
        source_file = os.path.join(srcDIR, file_name)
        destination_file = os.path.join(dstnDIR, file_name)
        shutil.copy2(source_file, destination_file)
# ^^^   
 
 
 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def fileNAME0(patH,name):
    folder_path = patH
    search_term = name   # replace with the word you want to search for
    file_pattern = f"{folder_path}/*{search_term}*"  # create a pattern to match files containing the search term
    matching_files = glob.glob(file_pattern)
    if len(matching_files)>1: 
        raise ValueError("There are 2 files with same keyword")
    fnme = os.path.basename(matching_files[0])
    return fnme
# ^^^


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def fileNAME(patH,name,name1):
    folder_path = patH
    search_term = name   # replace with the word you want to search for
    search_term1 = name1
    file_pattern = f"{folder_path}/*{search_term}*{search_term1}*"  # create a pattern to match files containing the search term
    matching_files = glob.glob(file_pattern)
    if len(matching_files)>1: 
        raise ValueError("There are 2 files with same keyword")
    fnme = os.path.basename(matching_files[0])
    return fnme
# ^^^
  


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find a file with these 3 keywds (will fail when you have AIS as different cmponents.)
def fileNAME1(folder_path,search_term, search_term1, search_term2):
    matching_files = []
    for file_path in glob.glob(f"{folder_path}/*"):
        file_name = os.path.basename(file_path)
        if all(term in file_name for term in [search_term, search_term1, search_term2]):
            matching_files.append(file_path)
    if len(matching_files) > 1:
        raise ValueError("There are 2 or more files with the same keyword")
    fnme = os.path.basename(matching_files[0])
    return fnme
# ^^^


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def ConvertSamplesToDist(infile, outfile):

	# Initialize valid variable names and units
	x_varname = None
	x_units = None
	x_scalefactor = None
	valid_varnames = ["sea_level_change", "sea_level_change_rate"]
	valid_varunits = {"sea_level_change": "mm", "sea_level_change_rate": "mm per year"}
	valid_scalefactors = {"sea_level_change": 1.0, "sea_level_change_rate": 0.1}

	# Quantiles to extract from the data
	q = np.unique(np.append(np.round(np.linspace(0,1,101), 3), (0.001, 0.005, 0.01, 0.05, 0.167, 0.5, 0.833, 0.95, 0.99, 0.995, 0.999)))

	# Missing value
	nc_missing_value = np.iinfo(np.int16).min

	# Open the file
	with xr.open_dataset(infile) as nc:

		# Check which variable we're dealing with
		for this_varname in valid_varnames:
			if this_varname in nc.keys():
				x_varname = this_varname
				x_scalefactor = valid_scalefactors[this_varname]
				x_units = valid_varunits[this_varname]
				break
		if x_varname is None:
			print("No valid variable name exists in {}".format(infile))
			return(1)


		# If the variable exists, extract and convert

		x_lats = nc['lat']
		x_lons = nc['lon']
		x_years = nc['years']
		x_ids = nc['locations'].values

		# Extract the attributes for this data
		x_attrs = nc.attrs

		# Loop over the chunks of samples
		nsites = len(x_ids)
			
		site_idx = np.arange(nsites)
		# Extract the data and make the quantile
		x = np.nanquantile(nc[x_varname].isel({"locations":site_idx}), q, axis=0)

		# Create the output array
		x_tempout = xr.Dataset({x_varname: (("quantiles", "years", "locations"), x.data, {"units":x_units, "missing_value":nc_missing_value}),
								"lat": (("locations"), x_lats[site_idx].data),
								"lon": (("locations"), x_lons[site_idx].data)},
			coords={"years": x_years.data, "locations": x_ids[site_idx].data, "quantiles": q.data}, attrs=x_attrs)

		# Write theoutput file
		x_tempout.to_netcdf(outfile, encoding={x_varname: {"scale_factor":x_scalefactor, "dtype": "i2", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})
			

	return(None)
# ^^^


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def LoadInfiles(infiles, years):

	# Initialize the return variables
	localsl_q = []
	file_list = []
	varname = None
	varunit = None
	ids = None
	lats = None
	lons = None
	qvar = None

	# Valid variables over which to generate pboxes
	valid_varnames = ["sea_level_change", "sea_level_change_rate"]
	valid_varunits = {"sea_level_change": "mm", "sea_level_change_rate": "mm per year"}
	valid_varscale = {"sea_level_change": 1.0, "sea_level_change_rate": 0.1}

	# Initialize the first file flag
	first_file = True

	# Loop over the input files
	for infile in infiles:

		# Open the input file for reading
		nc = Dataset(infile, "r")

		# Find which years in the projections match the requested years
		ncyears = nc.variables['years'][:]
		_, _, years_idx = np.intersect1d(years, ncyears, return_indices=True)

		# If this is the first file...
		# 1) determine which variable we're dealing with
		# 2) extract the non-pbox information that gets passed to the final pbox component file
		if first_file:

			# Test if a valid variable exists in this file. If so, set varname and
			# varunit and break loop.
			for this_varname in valid_varnames:
				if this_varname in nc.variables.keys():
					varname = this_varname
					varunit = valid_varunits[varname]
					varscale = valid_varscale[varname]
					break

			# If we couldn't find any valid variable name in the first file, exit with error
			if varname is None:
				raise Exception("No valid variable name exists in {}".format(infile))

			# Variable was found. Continue by collecting the rest of the non-pbox information
			ids = nc.variables['locations'][:]
			lats = nc.variables['lat'][:]
			lons = nc.variables['lon'][:]
			qvar = nc.variables['quantiles'][:]

			# Set the first file flag to false
			first_file = False

		# Extract the projection quantiles for this file
		localsl_q.append(nc.variables[varname][:,years_idx,:])

		# Close the netcdf file
		nc.close()

	# Convert everything into numpy arrays
	localsl_q = np.array(localsl_q)
	ids = np.array(ids)
	lats = np.array(lats)
	lons = np.array(lons)
	#qvar = np.array(qvar)
	qvar = np.round(np.array(qvar), 3)

	# Return everything
	return(localsl_q, varname, varunit, varscale, ids, lats, lons, qvar)


def main(infiles, outfile, pyear_start=2020, pyear_end=2100, pyear_step=10):

	# Read in all the infiles and generate a large data array
	years = np.arange(pyear_start, pyear_end+1, pyear_step)
	component_data, varname, varunit, varscale, ids, lats, lons, qvar = LoadInfiles(infiles, years)

	# Make the infile string for the output netcdf file
	infile_string = ", ".join([os.path.basename(x) for x in infiles])

	# Which indices are the median, above the median, and below the median
	median_idx = np.flatnonzero(qvar == 0.5)
	above_idx = np.arange(median_idx + 1, len(qvar))
	below_idx = np.arange(median_idx)

	# Initialize the output pbox array
	pbox = np.full(component_data.shape[1:], np.nan)

	# Use the mean of the medians of the components as the median for the pbox
	pbox[median_idx,:,:] = np.mean(component_data[:,median_idx,:,:], axis=0)

	# Use the minimum for below the median and maximum for above
	pbox[below_idx,:,:] = np.amin(component_data[:,below_idx,:,:], axis=0)
	pbox[above_idx,:,:] = np.amax(component_data[:,above_idx,:,:], axis=0)

	# Write the output netcdf file
	rootgrp = Dataset(outfile, "w", format="NETCDF4")

	# Define Dimensions
	site_dim = rootgrp.createDimension("locations", pbox.shape[2])
	year_dim = rootgrp.createDimension("years", pbox.shape[1])
	q_dim = rootgrp.createDimension("quantiles", pbox.shape[0])

	# Populate dimension variables
	lat_var = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var = rootgrp.createVariable("lon", "f4", ("locations",))
	id_var = rootgrp.createVariable("locations", "i4", ("locations",))
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	q_var = rootgrp.createVariable("quantiles", "f8", ("quantiles",))

	# Create a data variable
	localslq = rootgrp.createVariable(varname, "i2", ("quantiles", "years", "locations"), zlib=True, complevel=4)
	localslq.missing_value = np.iinfo(np.int16).min
	localslq.scale_factor = varscale

	# Assign attributes
	rootgrp.description = "Pbox generated from a FACTS sea-level change projection workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "Generated with files: {}".format(infile_string)
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees East"
	localslq.units = varunit

	# Put the data into the netcdf variables
	lat_var[:] = lats
	lon_var[:] = lons
	id_var[:] = ids
	year_var[:] = years
	q_var[:] = qvar
	localslq[:,:,:] = pbox

	# Close the netcdf
	rootgrp.close()


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CL fun
def GetScenarios(dir):

	# Get the scenario names from the available scenario directories
	# Ignore any hidden directories (i.e. .DS_Store)
	pb1e_scenarios = [x for x in os.listdir(os.path.join(dir, "pb_1e")) if not re.search(r"^\.", x)]
	pb1f_scenarios = [x for x in os.listdir(os.path.join(dir, "pb_1f")) if not re.search(r"^\.", x)]
	pb2e_scenarios = [x for x in os.listdir(os.path.join(dir, "pb_2e")) if not re.search(r"^\.", x)]
	pb2f_scenarios = [x for x in os.listdir(os.path.join(dir, "pb_2f")) if not re.search(r"^\.", x)]

	# Find the overlapping scenarios
	med_scenarios = list(set(pb1e_scenarios) & set(pb1f_scenarios))
	low_scenarios = list(set(pb2e_scenarios) & set(pb2f_scenarios))

	# Return the overlapping scenarios
	med_scenarios.sort()
	low_scenarios.sort()
	return(med_scenarios, low_scenarios)



def GetFiles(dir):

	# There should be a file for glaciers, landwaterstorage, oceandynamics,
	# AIS, GIS, and total...and for regional projections verticallandmotion.
	file_keys = ["glaciers", "landwaterstorage", "sterodynamics", "AIS", "GIS", "total", "verticallandmotion"]

	# Initialize list of matched file keys
	match_files = {}

	# Loop over the keys and find the associated files
	for this_key in file_keys:

		# Locate this file in the directory
		pattern = "*{}*.nc".format(this_key)
		this_file = fnmatch.filter(os.listdir(dir), pattern)

		# There should be only one match
		if len(this_file) == 1:

			# Append the match
			match_files[this_key] = os.path.join(dir, this_file[0])

		elif len(this_file) > 1:
			raise Exception("More than one file matched in {} for key {}".format(dir, this_key))

		else:

			match_files[this_key] = None

	# Return the dictionary of files
	return(match_files)



# def MakeConfidenceFile(infile_e=None, infile_f=None, f_years=np.arange(2020,2301,10), outfile=None, is_rates=False):
#-pk def MakeConfidenceFile(infile_e=None, infile_f=None, f_years=np.arange(2020,2151,10), outfile=None, is_rates=False):
def MakeConfidenceFile(infile_e=None, infile_f=None, f_years=np.arange(2020,2101,10), outfile=None, is_rates=False):

	# If both infile_e and infile_f are None, then there's no data for this component key.
	# Return and let the code move onto the next component key
	if infile_f is None and infile_e is None:
		return(1)

	# Variable names and attributes
	if is_rates:
		varname = "sea_level_change_rate"
		varscale = 0.1
	else:
		varname = "sea_level_change"
		varscale = 1.0

	# Open and subset the f file
	with xr.open_dataset(infile_f) as nc_f:
	#with xr.open_dataset(infile_f) as nc_f:
		nc_out = nc_f.sel(years=f_years)

	# Add the f file to the source list
	source_files = [infile_f]

	# If there's an e file, overlap it with the f file
	if infile_e is not None:
		with xr.open_dataset(infile_e) as nc_e:
		#with xr.open_dataset(infile_e) as nc_e:
			nc_out = nc_e.combine_first(nc_f.sel(years=f_years))

		# Append the e file to the source file list
		source_files.append(infile_e)

	# Define the missing value for the netCDF files
	nc_missing_value = np.iinfo(np.int16).min

	# Attributes for the output file
	nc_attrs = {"description": "Combined confidence output file for AR6 sea-level change projections",
			"history": "Created " + time.ctime(time.time()),
			"source": "Files Combined: {}".format(",".join(source_files))}

	# Put the attributes onto the output file
	nc_out.attrs = nc_attrs

	# Write the output file
	nc_out.to_netcdf(outfile, encoding={varname: {"scale_factor": varscale, "dtype": "i2", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	# Done
	return(None)

# ==> main loop
def GenerateConfidenceFiles(pboxdir, outdir):

	# Are we working with values or rates?
	is_rates = True if re.search(r"rates", pboxdir) is not None else False

	# Get the overlapping scenarios for each confidence level
	med_scenarios, low_scenarios = GetScenarios(pboxdir)

	# If these are rate pboxes...
	if is_rates:

		# Medium-confidence: pb_1f through 2150
		# Low-confidence: pb_2f through 2300

		# Loop over the medium scenarios
		for this_scenario in med_scenarios:

			# Get the list of files for this scenario
			pb1f_infiles = GetFiles(os.path.join(pboxdir, "pb_1f", this_scenario))

			# Loop over the available components
			for this_key in pb1f_infiles.keys():

				# Define the output file name
				outpath = Path(os.path.join(outdir, "medium_confidence", this_scenario))
				Path.mkdir(outpath, parents=True, exist_ok=True)
				outfile = os.path.join(outpath, "{}_{}_medium_confidence_rates.nc".format(this_key, this_scenario))

				# Make the output file
				#-pk  MakeConfidenceFile(infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2151,10), outfile=outfile, is_rates=is_rates)
				MakeConfidenceFile(infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2101,10), outfile=outfile, is_rates=is_rates)

		# Loop over the low scenarios
		for this_scenario in low_scenarios:

			# Get the list of files for this scenario
			pb1f_infiles = GetFiles(os.path.join(pboxdir, "pb_2f", this_scenario))

			# Loop over the available components
			for this_key in pb1f_infiles.keys():

				# Define the output file name
				outpath = Path(os.path.join(outdir, "low_confidence", this_scenario))
				Path.mkdir(outpath, parents=True, exist_ok=True)
				outfile = os.path.join(outpath, "{}_{}_low_confidence_rates.nc".format(this_key, this_scenario))

				# Make the output file
				# MakeConfidenceFile(infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2301,10), outfile=outfile, is_rates=is_rates, chunksize=chunksize)
				#-pk MakeConfidenceFile(infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2151,10), outfile=outfile, is_rates=is_rates)
				MakeConfidenceFile(infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2101,10), outfile=outfile, is_rates=is_rates)

	# These are value files...
	else:

		# Medium-confidence: pb_1e through 2100, pb_1f through 2150
		# Low-confidence: pb_2e through 2100, pb_2f through 2300

		# Loop over the medium scenarios
		for this_scenario in med_scenarios:

			# Get the list of files for this scenario
			pb1e_infiles = GetFiles(os.path.join(pboxdir, "pb_1e", this_scenario))
			pb1f_infiles = GetFiles(os.path.join(pboxdir, "pb_1f", this_scenario))

			# Loop over the available components
			for this_key in pb1e_infiles.keys():

				# Define the output file name
				outpath = Path(os.path.join(outdir, "medium_confidence", this_scenario))
				Path.mkdir(outpath, parents=True, exist_ok=True)
				outfile = os.path.join(outpath, "{}_{}_medium_confidence_values.nc".format(this_key, this_scenario))

				# Make the output file
				#-pk MakeConfidenceFile(infile_e=pb1e_infiles[this_key], infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2151,10), outfile=outfile)
				MakeConfidenceFile(infile_e=pb1e_infiles[this_key], infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2101,10), outfile=outfile)

		# Loop over the low scenarios
		for this_scenario in low_scenarios:

			# Get the list of files for this scenario
			pb1e_infiles = GetFiles(os.path.join(pboxdir, "pb_2e", this_scenario))
			pb1f_infiles = GetFiles(os.path.join(pboxdir, "pb_2f", this_scenario))

			# Loop over the available components
			for this_key in pb1e_infiles.keys():

				# Define the output file name
				outpath = Path(os.path.join(outdir, "low_confidence", this_scenario))
				Path.mkdir(outpath, parents=True, exist_ok=True)
				outfile = os.path.join(outpath, "{}_{}_low_confidence_values.nc".format(this_key, this_scenario))

				# Make the output file
				# MakeConfidenceFile(infile_e=pb1e_infiles[this_key], infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2301,10), outfile=outfile, chunksize=chunksize)
				#-pk MakeConfidenceFile(infile_e=pb1e_infiles[this_key], infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2151,10), outfile=outfile)
				MakeConfidenceFile(infile_e=pb1e_infiles[this_key], infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2101,10), outfile=outfile)

	# Done
	return(None)
# ^^^



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
region="global"
#
total = {
    "wf_1e": [f'total.workflow.wf1e.{region}.nc'],
    "wf_1f": [f'total.workflow.wf1f.{region}.nc'],
    "wf_2e": [f'total.workflow.wf2e.{region}.nc'],
    "wf_2f": [f'total.workflow.wf2f.{region}.nc'],
    "wf_3e": [f'total.workflow.wf3e.{region}.nc'],
    "wf_3f": [f'total.workflow.wf3f.{region}.nc'],
    "wf_4":  [f'total.workflow.wf4.{region}.nc']
}


common_files = [
    f'lws.ssp.landwaterstorage_{region}sl.nc',
    f'ocean.tlm.sterodynamics_{region}sl.nc',
    f'k14vlm.kopp14.verticallandmotion_localsl.nc'
]

component = {
    "wf_1e": [
        f'emuGrIS.emulandice.GrIS_{region}sl.nc',
        f'emuAIS.emulandice.AIS_{region}sl.nc',
        f'emuglaciers.emulandice.glaciers_{region}sl.nc',
        *common_files   
    ],
    "wf_1f": [
        f'GrIS1f.FittedISMIP.GrIS_GIS_{region}sl.nc',
        f'ar5AIS.ipccar5.icesheets_AIS_{region}sl.nc',
        f'ar5glaciers.ipccar5.glaciers_{region}sl.nc',
        *common_files
    ],
    "wf_2e": [
        f'emuGrIS.emulandice.GrIS_{region}sl.nc',
        f'larmip.larmip.AIS_{region}sl.nc',
        f'emuglaciers.emulandice.glaciers_{region}sl.nc',
        *common_files
    ],
    "wf_2f": [
        f'GrIS1f.FittedISMIP.GrIS_GIS_{region}sl.nc', 
        f'larmip.larmip.AIS_{region}sl.nc', 
        f'ar5glaciers.ipccar5.glaciers_{region}sl.nc', 
        *common_files
    ],
    "wf_3e": [
        f'emuGrIS.emulandice.GrIS_{region}sl.nc', 
        f'deconto21.deconto21.AIS_AIS_{region}sl.nc', 
        f'emuglaciers.emulandice.glaciers_{region}sl.nc', 
        *common_files
    ],
    "wf_3f": [
        f'GrIS1f.FittedISMIP.GrIS_GIS_{region}sl.nc', 
        f'deconto21.deconto21.AIS_AIS_{region}sl.nc', 
        f'ar5glaciers.ipccar5.glaciers_{region}sl.nc', 
        *common_files
    ],
    "wf_4": [
        f'bamber19.bamber19.icesheets_GIS_{region}sl.nc', 
        f'bamber19.bamber19.icesheets_AIS_{region}sl.nc', 
        f'ar5glaciers.ipccar5.glaciers_{region}sl.nc', 
        *common_files
    ],
}
# ^^^


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pBox = {
    "pb_1e": [
        f'emuAIS.emulandice.AIS_{region}sl.nc',
        f'larmip.larmip.AIS_{region}sl.nc',
        f'emuGrIS.emulandice.GrIS_{region}sl.nc',
        f'emuglaciers.emulandice.glaciers_{region}sl.nc',
        *common_files   
    ],
    "pb_1f": [
        f'ar5AIS.ipccar5.icesheets_AIS_{region}sl.nc',
        f'larmip.larmip.AIS_{region}sl.nc',
        f'GrIS1f.FittedISMIP.GrIS_GIS_{region}sl.nc',
        f'ar5glaciers.ipccar5.glaciers_{region}sl.nc',
        *common_files
    ],
    "pb_2e": [
        f'emuAIS.emulandice.AIS_{region}sl.nc',
        f'larmip.larmip.AIS_{region}sl.nc',
        f'deconto21.deconto21.AIS_AIS_{region}sl.nc',
        f'bamber19.bamber19.icesheets_AIS_{region}sl.nc',
        #
        f'emuGrIS.emulandice.GrIS_{region}sl.nc',
        f'bamber19.bamber19.icesheets_GIS_{region}sl.nc',
        #
        f'emuglaciers.emulandice.glaciers_{region}sl.nc',
        f'ar5glaciers.ipccar5.glaciers_{region}sl.nc',
        *common_files
    ],
    "wf_2f": [
        f'ar5AIS.ipccar5.icesheets_AIS_{region}sl.nc', 
        f'larmip.larmip.AIS_{region}sl.nc',
        f'deconto21.deconto21.AIS_AIS_{region}sl.nc',
        f'bamber19.bamber19.icesheets_AIS_{region}sl.nc',
        #
        f'GrIS1f.FittedISMIP.GrIS_GIS_{region}sl.nc', 
        f'bamber19.bamber19.icesheets_GIS_{region}sl.nc',
        #
        f'emuglaciers.emulandice.glaciers_{region}sl.nc',
        f'ar5glaciers.ipccar5.glaciers_{region}sl.nc',
        *common_files
    ],
}
# ^^^