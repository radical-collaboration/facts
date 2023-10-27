import os
import re
import sys
import time
import glob
import shutil
import fnmatch
from pathlib import Path
#
import numpy as np
import pandas as pd
import xarray as xr
import dask.array as da
from netCDF4 import Dataset
#
import argparse
#
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.basemap import Basemap
import cartopy
#
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
common_files = ['lws.ssp.landwaterstorage', 'ocean.tlm.sterodynamics']
local_files  =  'k14vlm.kopp14.verticallandmotion'
#
WF_file_patterns = {
    "wf1e": ['emuGrIS.emulandice.GrIS', 'emuAIS.emulandice.AIS', 'emuglaciers.emulandice.glaciers'],
    "wf2e": ['emuGrIS.emulandice.GrIS', 'larmip.larmip.AIS', 'emuglaciers.emulandice.glaciers'],
    "wf3e": ['emuGrIS.emulandice.GrIS', 'deconto21.deconto21.AIS_AIS', 'emuglaciers.emulandice.glaciers'],
    #
    "wf1f": ['GrIS1f.FittedISMIP.GrIS_GIS', 'ar5AIS.ipccar5.icesheets_AIS', 'ar5glaciers.ipccar5.glaciers'],
    "wf2f": ['GrIS1f.FittedISMIP.GrIS_GIS', 'larmip.larmip.AIS', 'ar5glaciers.ipccar5.glaciers'],
    "wf3f": ['GrIS1f.FittedISMIP.GrIS_GIS', 'deconto21.deconto21.AIS_AIS', 'ar5glaciers.ipccar5.glaciers'],
    "wf4": ['bamber19.bamber19.icesheets_GIS', 'bamber19.bamber19.icesheets_AIS', 'ar5glaciers.ipccar5.glaciers']
}
PB_file_patterns = {
    "pb_1e": ['emuAIS.emulandice.AIS', 'larmip.larmip.AIS', 'emuGrIS.emulandice.GrIS', 'emuglaciers.emulandice.glaciers'],
    "pb_2e": ['emuAIS.emulandice.AIS', 'larmip.larmip.AIS', 'deconto21.deconto21.AIS_AIS', 'bamber19.bamber19.icesheets_AIS',
              'emuGrIS.emulandice.GrIS', 'bamber19.bamber19.icesheets_GIS', 'emuglaciers.emulandice.glaciers', 'ar5glaciers.ipccar5.glaciers'],
    #
    "pb_1f": ['ar5AIS.ipccar5.icesheets_AIS', 'larmip.larmip.AIS', 'GrIS1f.FittedISMIP.GrIS_GIS','ar5glaciers.ipccar5.glaciers'],
	"pb_2f": ['ar5AIS.ipccar5.icesheets_AIS','larmip.larmip.AIS','deconto21.deconto21.AIS_AIS','bamber19.bamber19.icesheets_AIS',
			  'GrIS1f.FittedISMIP.GrIS_GIS','bamber19.bamber19.icesheets_GIS','emuglaciers.emulandice.glaciers', 'ar5glaciers.ipccar5.glaciers']
}
# ^^^
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Generic Function block
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
def create_directory(dir_name: str, action: str = "") -> str:
    full_path = os.path.join(os.getcwd(), dir_name)
    if action == "remove":
        shutil.rmtree(full_path, ignore_errors=True)
    os.makedirs(full_path, exist_ok=True)
    return full_path
#
def delete_files_with_pattern(folder_path, pattern, exclusion_pattern):
    folder_path = Path(folder_path) #Convert folder path string to a Path object
    for file_path in folder_path.glob(pattern):
        if not file_path.match(exclusion_pattern):
            file_path.unlink()
#
def copy_filename_with_pattern(source_dir,destination_dir,pattern):
    source_file_pattern = os.path.join(source_dir, pattern)
    matching_files = glob.glob(source_file_pattern)
    for file_path in matching_files:
        shutil.copy2(file_path, destination_dir)   
#
def copy_all_files_from(srcDIR,dstnDIR):
    file_list = os.listdir(srcDIR)
    for file_name in file_list:
        source_file = os.path.join(srcDIR, file_name)
        destination_file = os.path.join(dstnDIR, file_name)
        shutil.copy2(source_file, destination_file)
# 
def find_filename_with_pattern(path, search_term):
    path_obj = Path(path)
    matching_files = list(path_obj.glob(f"*{search_term}*"))
    if len(matching_files) > 1:
        raise ValueError("There are multiple files with the same keyword.")
    elif not matching_files:
        return None
    else:
        return matching_files[0].name
# ^^^

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## save the Files and folders.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def generate_pbox_nc_files_(expFolder,region,ssps,FOLDER1):

    def copy_and_convert(source_file, dest_folder, region):
        #shutil.copy2(source_file, dest_folder) #copy the og(sample) files.
        out_file_path = os.path.join(dest_folder, os.path.basename(source_file).split(".nc")[0] + "_quantiles.nc")
        Samples_to_Quantiles(source_file, out_file_path)
    
    path = create_directory(FOLDER1,"remove")
    
    skip = {
        'wf3e': ['ssp119', 'ssp370'],
        'wf3f': ['ssp119', 'ssp370'],
        'wf4':  ['ssp119', 'ssp245', 'ssp370']
    }
    
    for workflow, file_patterns in WF_file_patterns.items():
        for ssp in ssps:
            if ssp in skip.get(workflow, []):
                continue
            
            ssp_path = os.path.join(path, workflow, ssp)
            os.makedirs(ssp_path, exist_ok=True)
            
            # Copy Component files & convert to Quantiles.
            for file in file_patterns:
                component_file_path = f'{expFolder}/coupling.{ssp}/output/coupling.{ssp}.{file}_{region}sl.nc'
                copy_and_convert(component_file_path, ssp_path, region)
            
            # Copy Common files & convert to Quantiles.
            for common_file_item in common_files:
                common_file_path = f'{expFolder}/coupling.{ssp}/output/coupling.{ssp}.{common_file_item}_{region}sl.nc'
                copy_and_convert(common_file_path, ssp_path, region)
            
            # Copy Total files & convert to Quantiles.
            total_file_path = f'{expFolder}/coupling.{ssp}/output/coupling.{ssp}.total.workflow.{workflow}.{region}.nc'
            copy_and_convert(total_file_path, ssp_path, region)
            
            # Copy VLM file & convert to Quantiles.  
            if region == 'local':
                vlm_file_path = f'{expFolder}/coupling.{ssp}/output/coupling.{ssp}.{local_files}.{region}sl.nc'
                copy_and_convert(vlm_file_path, ssp_path, region)




# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Convert Samples to Quantiles.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Samples_to_Quantiles(in_file, out_file):
	#
	valid_variables = {
		"sea_level_change": {"units": "mm", "scale_factor": 1.0},
		"sea_level_change_rate": {"units": "mm per year", "scale_factor": 0.1}
	}
	quantiles = np.unique(np.append(np.round(np.linspace(0,1,101), 3), (0.001, 0.005, 0.01, 0.05, 0.167, 0.5, 0.833, 0.95, 0.99, 0.995, 0.999)))
	#
	missing_value = np.iinfo(np.int16).min
	#
	# Open the file
	with xr.open_dataset(in_file) as dataset:
		variable_name = next((var for var in valid_variables if var in dataset.variables), None)
		if not variable_name:
			print(f"No valid variable name exists in {in_file}")  # Changed input_file to in_file
			return 1
		#
		variable_data = {
			"units": valid_variables[variable_name]["units"],
			"missing_value": missing_value
		}
		scale_factor = valid_variables[variable_name]["scale_factor"]
		#
		lats = dataset['lat']
		lons = dataset['lon']
		years = dataset['years']
		location_ids = dataset['locations'].values
		quantile_values = np.nanquantile(dataset[variable_name], quantiles, axis=0)
		#
		output_data = xr.Dataset({
			variable_name: (("quantiles", "years", "locations"), quantile_values, variable_data),
			"lat": (("locations"), lats.data), "lon": (("locations"), lons.data)
		}, coords={"years": years.data, "locations": location_ids.data, "quantiles": quantiles}, attrs=dataset.attrs)
		#
		encoding = {
			variable_name: {
				"scale_factor": scale_factor,
				"dtype": "i2",
				"zlib": True,
				"complevel": 4,
				"_FillValue": missing_value
			}
		}
		output_data.to_netcdf(out_file, encoding=encoding)
	#
	return None
#^^^


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Generate P-box files.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
def load_infiles(infiles, years):
    valid_varnames = ["sea_level_change", "sea_level_change_rate"]
    valid_varunits = {"sea_level_change": "mm", "sea_level_change_rate": "mm per year"}
    valid_varscale = {"sea_level_change": 1.0, "sea_level_change_rate": 0.1}
	#
    localsl_q = []
    ds = xr.open_dataset(infiles[0], engine='netcdf4')
	#
    varname = next((v for v in valid_varnames if v in ds.variables), None)
    if not varname:
        raise ValueError(f"No valid variable name exists in {infiles[0]}")
	#
    varunit = valid_varunits[varname]
    varscale = valid_varscale[varname]
    ids = ds['locations'].values
    lats = ds['lat'].values
    lons = ds['lon'].values
    qvar = np.round(ds['quantiles'].values, 3)
	#
    for infile in infiles:
        with xr.open_dataset(infile, engine='netcdf4') as ds:
            localsl_q.append(ds[varname].sel(years=years).values)
	#
    return np.array(localsl_q), varname, varunit, varscale, ids, lats, lons, qvar
#
# ...............................................................................................
def generate_pbox(infiles, outfile, pyear_start, pyear_end, pyear_step):
    years = np.arange(pyear_start, pyear_end+1, pyear_step)
    component_data, varname, varunit, varscale, ids, lats, lons, qvar = load_infiles(infiles, years)
	#
    median_idx = np.flatnonzero(qvar == 0.5)
    above_idx = np.arange(median_idx + 1, len(qvar))
    below_idx = np.arange(median_idx)
	#
    pbox = np.full(component_data.shape[1:], np.nan)
    pbox[median_idx,:,:] = np.mean(component_data[:,median_idx,:,:], axis=0)
    pbox[below_idx,:,:] = np.amin(component_data[:,below_idx,:,:], axis=0)
    pbox[above_idx,:,:] = np.amax(component_data[:,above_idx,:,:], axis=0)
	#
    dataset = xr.Dataset({
        varname: (["quantiles", "years", "locations"], pbox, {
            "units": varunit,
            "missing_value": np.iinfo(np.int16).min,
            "scale_factor": varscale
        }),
        "lat": (["locations"], lats),
        "lon": (["locations"], lons),
        "locations": (["locations"], ids),
        "years": (["years"], years),
        "quantiles": (["quantiles"], qvar)
    }, attrs={
        "description": "Pbox generated from a FACTS sea-level change projection workflow",
        "history": f"Created {time.ctime(time.time())}",
        "source": f"Generated with files: {', '.join([os.path.basename(x) for x in infiles])}"
    })

    dataset.to_netcdf(outfile)
#
# ...............................................................................................
#
def create_Pbox(ssps,specify_pbox_files,FOLDER1):
    #
    cp_file_with_pattern    = lambda ssp_val, wf, pattern: copy_filename_with_pattern(f"{os.getcwd()}/1_workflow/{wf}/{ssp_val}", ssp_path, f"*{pattern}*")
    copy_all_wf_files       = lambda p, s, d: copy_all_files_from(os.path.join(os.getcwd(), f'1_workflow/wf{p.split("_")[1]}/{s}'), d)
    rm_file_with_pattern    = lambda ssp_path, sl_file, pbox_file: delete_files_with_pattern(ssp_path, sl_file, pbox_file)
    #    
    def create_pb_files(ssp_path, patterns, outfile_prefix, pbox, ssp_val):
        infiles = [f'{ssp_path}/{find_filename_with_pattern(ssp_path, pattern)}' for pattern in patterns]
        outfile = f"{ssp_path}/{outfile_prefix}-pb{pbox.split('_')[1]}-{ssp_val}_globalsl.nc"
        generate_pbox(infiles, outfile, pyear_start=2020, pyear_end=2100, pyear_step=10)
        #
    #
    path = create_directory(FOLDER1,"remove")
    #
    for pbox in PB_file_patterns:
        create_directory(f'{os.path.join(path, pbox)}')
        for ssp_val in ssps:
            # ..............................................................................
            if pbox in ['pb_2e', 'pb_2f'] and ssp_val in['ssp119','ssp245','ssp370']: continue
            # ..............................................................................
            ssp_path = create_directory(f'{os.path.join(path, pbox,ssp_val)}')
            #
            # ---> loop over pbox
            if pbox in specify_pbox_files:
                # Copy all wf files (ie. pb1e <-- wf1e)
                copy_all_wf_files(pbox, ssp_val, ssp_path)
                #
                pbox_files = specify_pbox_files[pbox]
                #
                # Copy additional needed files.
                for indv_wf, pattern in pbox_files['files_to_cp']:
                    cp_file_with_pattern(ssp_val, indv_wf, pattern)   
                #
                # Create P-box files and save.
                create_pb_files(ssp_path, pbox_files['AIS_pb'], "icesheets-AIS", pbox, ssp_val)
                create_pb_files(ssp_path, pbox_files['total_pb'], "total-workflow", pbox, ssp_val)
                if pbox not in ['pb_1e','pb_1f']:
                    create_pb_files(ssp_path, pbox_files['GIS_pb'], "icesheets-GIS", pbox, ssp_val)
                    if pbox in ['pb_2e']: create_pb_files(ssp_path, pbox_files['glacier_pb'], "glaciers", pbox, ssp_val)
                #
                # Remove OG files used to create final P-box files
                for sl_file, pbox_file in pbox_files['rm_pattern_xcpt']:
                    rm_file_with_pattern(ssp_path, sl_file, pbox_file)
                #
        print(f'... created {pbox} folder')
#^^^

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Generate Confidence level files.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_scenarios(directory: str) -> tuple:
	#
    scenario_types = list(PB_file_patterns.keys())
    scenarios = {stype: [x for x in os.listdir(os.path.join(directory, stype)) if not x.startswith('.')] for stype in scenario_types}
	#
    med_scenarios = sorted(list(set(scenarios["pb_1e"]) & set(scenarios["pb_1f"])))
    low_scenarios = sorted(list(set(scenarios["pb_2e"]) & set(scenarios["pb_2f"])))
	#
    return med_scenarios, low_scenarios

# ............................................................
def make_confidence_file(infile_e=None, infile_f=None, f_years=np.arange(2020, 2101, 10), outfile=None, is_rates=False):
    if infile_f is None and infile_e is None:
        return 1
    varname = "sea_level_change_rate" if is_rates else "sea_level_change"
    varscale = 0.1 if is_rates else 1.0
	#
    with xr.open_dataset(infile_f, engine='netcdf4') as nc_f:
        nc_out = nc_f.sel(years=f_years)
	#
    source_files = [infile_f]
    if infile_e is not None:
        with xr.open_dataset(infile_e, engine='netcdf4') as nc_e:
            nc_out = nc_e.combine_first(nc_f.sel(years=f_years))
        source_files.append(infile_e)
	#
    nc_missing_value = np.iinfo(np.int16).min
    nc_attrs = {
        "description": "Combined confidence output file for AR6 sea-level change projections",
        "history": f"Created {time.ctime(time.time())}",
        "source": f"Files Combined: {','.join(source_files)}"
    }
	#
    nc_out.attrs = nc_attrs
    nc_out.to_netcdf(outfile, encoding={
        varname: {
            "scale_factor": varscale,
            "dtype": "i2",
            "zlib": True,
            "complevel": 4,
            "_FillValue": nc_missing_value
        }
    })
# ............................................................
def generate_confidence_files(pboxdir: str, outdir: str):
    #
    is_rates = "rates" in pboxdir
    med_scenarios, low_scenarios = get_scenarios(pboxdir)
	#
    confidence_map = {"medium_confidence": med_scenarios,
                      "low_confidence": low_scenarios}
    sl_files = ["glaciers", "landwaterstorage", "sterodynamics", "AIS", "GIS", "total", "verticallandmotion"]
	#
    for conf_level, scenarios in confidence_map.items():
        for scenario in scenarios:
            scenario_dir = os.path.join(pboxdir, "pb_1f" if "medium" in conf_level else "pb_2f", scenario)
            #       
            files = {sl: f'{scenario_dir}/{find_filename_with_pattern(scenario_dir, sl)}' for sl in sl_files}
            #
            for key, file_path in files.items():
                if file_path is None or 'None' in file_path:
                    continue  
                #
                outpath = Path(outdir, conf_level, scenario); outpath.mkdir(parents=True, exist_ok=True)
                filename = f"{key}_{scenario}_{conf_level}_{'rates' if is_rates else 'values'}.nc"
                outfile = outpath / filename
                #
                make_confidence_file(infile_f=file_path, f_years=np.arange(2020, 2101, 10), outfile=str(outfile), is_rates=is_rates)

# ^^^


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Data frame for Workflows.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def workflow_dataframe():
    df_reshaped = pd.DataFrame.from_dict(WF_file_patterns, orient='index').transpose().melt(var_name="Workflow", value_name="Pattern").dropna()

    # Grouping by the 'Workflow' column and aggregating the patterns into separate columns
    df_grouped = df_reshaped.groupby("Workflow")['Pattern'].apply(list).reset_index()
    df_grouped_patterns = pd.DataFrame(df_grouped['Pattern'].to_list())
    df_final = pd.concat([df_grouped['Workflow'], df_grouped_patterns], axis=1)

    # Renaming the columns
    unique_patterns = df_reshaped['Pattern'].unique()
    column_labels = [pattern.split('.')[-1] for pattern in unique_patterns]
    new_column_names = ['Workflow'] + column_labels[:df_final.shape[1]-1]
    df_final.columns = new_column_names
    return df_final