import numpy as np
import os
import sys
from netCDF4 import Dataset
import cftime

''' IncludeCMIP6ZOSModels.py

This script parses through a directory of models and loads annual mean 'zos' data from each model.
A directory structure of 'variable'>'Model' is expected.
Historical and SSP files are expected for each model, if not the model is excluded from the ensemble.

Parameters:
model_dir       = Directory of model output. Each model is a subdirectory within this one.
varname        = Name of the variables of interest
years           = Years of interest.
scenario    = SSP of interest

Return:
model_list  = Vector of model names that are to be included (nmodels)
ZOS     = Sea-level change (years, nmodels)

'''

def angd(lat0, lon0, lat, lon):

	# Convert the input from degrees to radians
	(lat0, lon0) = np.radians((lat0, lon0))
	(lat, lon) = np.radians((lat, lon))

	# Calculate the angle between the vectors
	temp = np.arctan2(np.sqrt((np.cos(lat)*np.sin(lon-lon0))**2 + \
	(np.cos(lat0)*np.sin(lat) - np.sin(lat0)*np.cos(lat) * np.cos(lon-lon0))**2),\
	(np.sin(lat0)*np.sin(lat) + np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0)))

	# Convert the results from radians to degrees and return
	return(np.degrees(temp))

def CalcWeights(qlat, qlon, lat, lon, idwrad, idwpow, idwmin):

	# Calculate the distance between the query lat/lon and the grid points
	ad = angd(lat, lon, qlat, qlon)
	ad[np.nonzero(ad <= idwmin)] = idwmin

	# Find the indices that are below the IDW threshold
	idw_idx = np.nonzero(ad <= idwrad)

	# Calculate the weights
	idw_weights = ad[idw_idx]**-idwpow

	return(idw_idx, idw_weights)

def IDW(val, idw_weights, idw_idx):

	# Apply the weights to the data
	num = val[idw_idx[0], idw_idx[1],:] * idw_weights[:,np.newaxis]

	# Initialize the output
	qval = []

	# Loop over the times
	for i in np.arange(num.shape[1]):

		# Identify which data points are valid for this year
		good_idx = np.flatnonzero(~np.isnan(num[:,i]))

		# Normalize num by the sum of valid weights at this time
		if(len(good_idx) > 0):
			qval.append(np.sum(num[good_idx,i]) / np.sum(idw_weights[good_idx]))
		else:
			qval.append(np.nan)

	# Convert the output to a numpy array
	qval = np.array(qval)

	# Apply the weights to the values
	#qval = np.nansum(val[idw_idx[0], idw_idx[1],:] * idw_weights[:,np.newaxis], axis=0) / np.nansum(idw_weights)
	#qval = np.sum(val[idw_idx[0], idw_idx[1],:] * idw_weights[:,np.newaxis], axis=0) / np.sum(idw_weights)

	return(qval)

#-----------------------------------------------------------------------------------------

def IncludeCMIP6ZOSModels(model_dir, varname, years, include_models, include_scenarios, focus_sites_lats, focus_sites_lons):

	# Initialize the model list and data matrix
	model_list = []
	scenario_list = []
	init_zos = True

	# Initialize variables to hold the IDW weights and indices
	all_idx = []
	all_weights = []
	ZOS = []

	# Initialize IDW parameters
	idw_rad = 3.5
	#idw_rad = 5.0
	idw_pow = 3.0
	idw_min = 0.005

	# Loop through available models in model_dir
	for i in np.arange(len(include_models)):

		# Try this model/scenario pair
		model = include_models[i]
		scenario = include_scenarios[i]

		# Skip if the folder/file found in this directory is hidden or directory is a file with a . extension
		if '.' in model:
			continue

		# Skip if this model is not available
		if model not in os.listdir(model_dir):
			continue

		# initialize lists to store data
		runtype_data = {'historical':[],scenario:[]}
		runtype_datayrs = {'historical':[],scenario:[]}

		incorporate = True # incorporate model or not

		# Read in historical and ssp data
		for runtype in ('historical',scenario):

			#start of filename for runtype currently processed
			filename_id = varname + '_Omon_' + model + '_' + runtype

			# find the historical or ssp file you want to read in for this model (exact filename depends on the experiment years)
			filename=[]
			for files_forModel in os.listdir(os.path.join(model_dir,model)): # loop through files in model folder
				if files_forModel[0:len(filename_id)] == filename_id:
					filename = files_forModel # assign filename

			if not filename: #if the right filename cannot be found:
				incorporate=False

			if(incorporate):

				# Open the netCDF file
				nc_fid = Dataset(os.path.join(model_dir,model,filename), 'r')

				# If this is the first model, collect the model lats and lons
				# Note, this assumes all models have been put on the same grid
				if(init_zos and runtype == 'historical'):

					# Raw model lats/lons
					model_lats = nc_fid.variables['lat'][:]
					model_lons = nc_fid.variables['lon'][:]
					n_model_lats = len(model_lats)
					n_model_lons = len(model_lons)

					# Reshape for use in IDW functions
					model_lats = np.tile(model_lats, (n_model_lons, 1))
					model_lons = np.tile(model_lons, (n_model_lats, 1)).T

					# Calculate the weights for all the sites
					for i in np.arange(len(focus_sites_lats)):

						# Calculate the weights
						(idx, weights) = CalcWeights(focus_sites_lats[i], focus_sites_lons[i], model_lats, model_lons, idw_rad, idw_pow, idw_min)
						all_idx.append(idx)
						all_weights.append(weights)

					# Done with initialization
					init_zos = False

				# read out the data
				datatime = nc_fid.variables['time']
				dat = nc_fid.variables[varname][:]
				dat = np.array(dat.filled(np.nan))

				# rearrange per year and compute average along year axis
				dat = np.mean( np.reshape(dat, (int(dat.shape[0]/12),12,dat.shape[1],dat.shape[2])), axis=1 )

				# Calculate the years
				nctime = cftime.num2date(datatime, datatime.units, datatime.calendar)
				datayrs = [int(x.strftime("%Y")) for x in nctime]
				#start_time = datetime.datetime(year=1850, month=6, day=1, hour=0, minute=0, second=0)
				#datayrs = [(start_time+datetime.timedelta(x)).year for x in datatime.data[0::12]]

				#store into dict for each cmip6 runtype
				runtype_data[runtype] = np.array(dat).T
				runtype_datayrs[runtype] = np.array(datayrs[::12])
				#runtype_data[runtype] = np.ma.array(dat).T
				#runtype_datayrs[runtype] = np.ma.array(datayrs[::12])

		if(incorporate):
			# check for overlap of historical and scenario datasets using data years
			overlap=np.isin(runtype_datayrs['historical'],runtype_datayrs[scenario],invert=True) # if historical years are in scenario
			runtype_datayrs['historical'] = runtype_datayrs['historical'][overlap] #remove these from the historical arrays
			runtype_data['historical'] = runtype_data['historical'][:,:,overlap]

			#concatenate historical and scenario years into full series
			fullyrs = np.concatenate((runtype_datayrs['historical'],runtype_datayrs[scenario]))
			fulldata = np.concatenate((runtype_data['historical'],runtype_data[scenario]), axis=2)

			# Put the ZOS data onto the requested years
			reduced_data = np.apply_along_axis(lambda fp, xp: np.interp(years, xp, fp, left=np.nan, right=np.nan), axis=2, arr=fulldata, xp=fullyrs)

			# Calculate the zos values for all sites from this model
			model_zos = map(lambda idx, w: IDW(reduced_data, w, idx), all_idx, all_weights)

			# Add this model to the overall data structure
			ZOS.append(np.array(list(model_zos)))

			# Append the model to the model list
			model_list.append(model)
			scenario_list.append(scenario)

	# Convert ZOS to a numpy array and reshape (years, models, sites)
	ZOS = np.array(ZOS)
	ZOS = np.transpose(ZOS, axes=(2,0,1))

	# Return variables
	return(model_list, scenario_list, ZOS)



if __name__ == '__main__':

	import pickle

	modeldir = "./data/cmip6/zos"
	varname = "zos"
	datayears = np.arange(1861,2300)
	include_models = os.listdir(modeldir)
	#include_models = ["IPSL-CM6A-LR", "CanESM5", "GISS-E2-1-G", "CNRM-ESM2-1", "MRI-ESM2-0"]
	include_scenarios = np.repeat("ssp119", len(include_models))

	#site_lats = [40.70, 40.70, 40.70, 40.70, 40.70, 40.70, 40.70]
	#site_lons = [x-0 for x in [-78.01, -77.01, -76.01, -75.01, -74.01, -73.01, -72.01]]

	#site_lats = [-17,-17,-17,-17,-17,-17,-17]
	#site_lons = [x+0 for x in [130, 131, 132, 133, 134, 135, 136]]

	#site_lats = [-17,-17,-17]
	#site_lons = [126,127,128]

	#site_lats = [-17,-17,-17,-17]
	#site_lons = [96,97,98,99]

	site_lats = [-7,-7,-7]
	site_lons = [-108,-107,-106]


	(modellist, scenariolist, ZOS) = IncludeCMIP6ZOSModels(modeldir, varname, datayears, include_models, include_scenarios, site_lats, site_lons)
	#output = {"ZOS": ZOS, "modellist": modellist}
	#outfile = open("includeCMIP6zos_temp.pkl", 'wb')
	#pickle.dump(output, outfile)
	#outfile.close()

	#with open("includeCMIP6zos_temp.pkl", 'rb') as f:
	#	pkl_file = pickle.load(f)
	#	ZOS = pkl_file["ZOS"]
	#	modellist = pkl_file["modellist"]


	print(modellist)

	#year_idx = np.flatnonzero(datayears == 2099)
	#baseyear_idx = np.flatnonzero(datayears == 2005)

	#for i in np.arange(len(modellist)):
		#temp_string = ",".join([str(x) for x in ZOS[year_idx,i,:]-ZOS[baseyear_idx,i,:]])
		#print("{}: {}".format(modellist[i], temp_string))

	#print("Mean")
	#print(np.mean(ZOS[year_idx,:,:]-ZOS[baseyear_idx,:,:], axis=1))
	#print("SD")
	#print(np.std(ZOS[year_idx,:,:]-ZOS[baseyear_idx,:,:], axis=1))



	sys.exit()
