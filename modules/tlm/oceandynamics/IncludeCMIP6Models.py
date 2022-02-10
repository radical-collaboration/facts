import numpy as np
import os
import sys
import cftime
from netCDF4 import Dataset

''' IncludeCMIP6Models.py

This script parses through a directory of models and loads annual mean 'zostoga' data from each model.
A directory structure of 'variable'>'Model' is expected.
PiControl, Historical and SSP files are expected for each model, if not the model is excluded from the ensemble.

Parameters:
model_dir		  = Directory of model output. Each model is a subdirectory within this one.
varname			  = Name of the variables of interest
years			  = Years of interest
include_models	  = List of models to attempt to include
include_scenario  = List of scenarios to attempt to include

Return:
model_list	= Vector of model names that are to be included (nmodels)
ZOSTOGA		= Global average thermosteric sea-level change (years, nmodels)
CONTROL_ZOSTOGA = Global average thermosteric sea-level change of piControl run (years,nmodels)

Note: The original code (Matlab version from Kopp14) includes both ZOSTOGA and ZOSGA,
though it's not apparent to me why these are being treated as the same quantity. ZOSTOGA
is the total thermosteric sea-level change while ZOSGA is total sea-level change including
mass flux from land-ice.

'''
def IncludeCMIP6Models(model_dir, varname, years, include_models, include_scenarios):

	# Initialize the model list and data matrix
	model_list = []
	scenario_list = []
	init_zostoga = True

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

		# Read in control, historical and ssp data
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
				# read out data
				nc_fid = Dataset(os.path.join(model_dir,model,filename), 'r')
				#datayrs = nc_fid.variables['year'][:]
				nctime = nc_fid.variables["time"]
				temp_nctime = cftime.num2date(nctime, nctime.units, nctime.calendar)
				datayrs = np.array([int(x.strftime("%Y")) for x in temp_nctime])
				dat = nc_fid.variables[varname][:]

				# if monthly means, convert to annual
				if datayrs[0] == datayrs[1]:
					# rearrange per year and compute average along year axis
					dat = np.mean( np.reshape(dat,(int(len(dat)/12), 12)), axis=1 )
					datayrs = datayrs[0::12]

				#store into dict for each cmip6 runtype
				runtype_data[runtype] = np.array(dat)
				runtype_datayrs[runtype] = np.array(datayrs)

		if(incorporate):
			# check for overlap of historical and scenario datasets using data years
			overlap=np.isin(runtype_datayrs['historical'],runtype_datayrs[scenario],invert=True) # if historical years are in scenario
			runtype_datayrs['historical'] = runtype_datayrs['historical'][overlap] #remove these from the historical arrays
			runtype_data['historical'] = runtype_data['historical'][overlap]

			#concatenate historical and scenario years into full series
			fullyrs = np.concatenate((runtype_datayrs['historical'],runtype_datayrs[scenario]))
			fulldata = np.concatenate((runtype_data['historical'],runtype_data[scenario]))

			#interpolate to requested years, add nans where no data available
			data_to_append = np.interp(years, fullyrs, fulldata, left=np.nan, right=np.nan)

			# If this model produces nan for ZOSTOGA or CONTROL_ZOSTOGA, remove the model
			if(np.all(np.isnan(data_to_append))):
				continue

			if(init_zostoga): # if first model
				ZOSTOGA = data_to_append
				init_zostoga = False
			else:
				ZOSTOGA = np.vstack((ZOSTOGA, data_to_append)) #stack arrays for the different models

			model_list.append(model) #list of model names
			scenario_list.append(scenario) #list of scenarios

	return(model_list, scenario_list, np.transpose(ZOSTOGA))


	''' #from previous version, Greg's code for quality-flagging pre-processed CMIP5:

		# If the change is too small, ignore this
		good_inds = np.nonzero(~np.isnan(dat[:,1]))[0]
		totalchange = np.sum(np.diff(dat[good_inds[0]:good_inds[-1], 1]))/(dat[good_inds[-1],0] - dat[good_inds[0],0])
		if(totalchange < 2e-4):
			incorporate=False

		# If there's no data near year 2000, ignore this
		if(np.min(np.abs(dat[good_inds,0]-2000)) > 2):
				incorporate=False
	'''

if __name__ == '__main__':

	modeldir = "./data/cmip6/zostoga"
	varname = "zostoga"
	datayears = np.arange(1861,2100)

	models = os.listdir(modeldir)
	scenarios = np.repeat("ssp119", len(models))

	(modellist, scenariolist, CZOSTOGA) = IncludeCMIP6Models(modeldir, varname, datayears, models, scenarios)

	print(modellist)


	exit()
