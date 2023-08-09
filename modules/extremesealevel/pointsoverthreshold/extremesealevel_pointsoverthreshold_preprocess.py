import os
import numpy as np
import pandas as pd
from datetime import datetime as dt
from netCDF4 import Dataset
import re
import pickle
import argparse

# Do not warn about chained assignments
pd.options.mode.chained_assignment = None  # default='warn'

''' extremesealevel_preprocess_pointsoverthreshold.py

Runs the pre-processing stage for the default extreme sea-level workflow.

Parameters:
total_localsl_file = Total localized sea-level projection file. Site lats/lons are taken from this file and mapped to the GESLA database [default="./total-workflow_localsl.nc"]
minDays = Minimum number of days in a valid year [default=250]
minYears = Minimum number of years available [default=20]
match_lim = Radius around requested locations to find a matching tide gauge in GESLA database [default=0.1]
center_year = This year +/- 9 years for centering data must be available [default=2000]
pctPot = Percentile for Point Over Threshold analysis [default=95]
gpd_pot_threshold = Percentile for GPD analysis [default=99.7]
cluster_lim = Maximum number of hours that define a cluster for extreme events [default=72]
target_years = Space-delimited list of projection years of interest (e.g. 2050 2100) [default=2100]
gesla_dir = Directory containing GESLA database [default="./gesla_data"]
pipeline_id = Unique identifier for this instance of the module [default="1234test"]

'''

def readmeta(filename):
	
	with open(filename,encoding = 'raw_unicode_escape') as myfile:
		head = [next(myfile) for x in range(6)]
	station_name = head[1][12:-1].strip()	 
	station_lat = float(head[4][11:-1].strip())
	station_lon = float(head[5][12:-1].strip()) 
	
	return(station_name, station_lat, station_lon)



def extract_gesla_locations(gesladir):
	
	# Get a list of the gesla database files
	geslafiles = os.listdir(gesladir)
	
	# Initialize the station variables
	station_names = []
	station_lats = []
	station_lons = []
	station_filenames = []
	
	# Loop over the gesla files
	for this_file in geslafiles:
		
		# Extract the station header information
		this_name, this_lat, this_lon = readmeta(os.path.join(gesladir, this_file))
		
		# Append this information to appropriate lists
		station_names.append(this_name)	 
		station_lats.append(float(this_lat))
		station_lons.append(float(this_lon))
		station_filenames.append(this_file)
	
	return(station_names, station_lats, station_lons, station_filenames)



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



def mindist(qlat, qlon, lats, lons, limit=0.1):
	
	# Calculate the angular distance
	dist = angd(lats, lons, qlat, qlon)
	
	# If the minimum distance is beyond the limit, print a warning and return None
	if(np.amin(dist) > limit):
		print("Warning: No match for lat={0} lon={1}. Minimum distance to nearest tide gauge is {2:.4f} (limit={3}). Returning \"None\"".format(qlat, qlon, np.amin(dist), limit))
		return(None)
	
	else:
		
		# Perform an indirect sort of the distances
		sort_idx = np.argsort(dist)
		
		# Find which ones fall within the radius limit
		min_idx = sort_idx[np.flatnonzero(dist[sort_idx] <= limit)]
		
		return(min_idx)


def load_tg(geslafile, minDays, minYears, center_year, pctPot, gpd_pot_threshold, cluster_lim):
	
	# Initialize data constraint flags
	include_file = False
	pass_nyears = False
	pass_centeryears = False
	
	#read meta data
	with open(geslafile,encoding = 'raw_unicode_escape') as myfile:
		head = [next(myfile) for x in range(6)]
	station_name = head[1][12:-1].strip()    
	station_lat = float(head[4][11:-1].strip())
	station_lon = float(head[5][12:-1].strip())
    
	#read raw data
	rawdata=np.loadtxt(geslafile, skiprows=32, dtype={'names':('date','hour','val','qf','ua'),'formats': ('S10', 'S8', 'f4', 'i4', 'i4')}, encoding = 'raw_unicode_escape') 
	
	#remove invalid or missing data (see also GESLA2 definitions)
	rawdata = rawdata[rawdata['val'] >-9.9] # skip remaining missing data
	rawdata = rawdata[rawdata['qf'] == 1] # data with quality-flag 1 is correct
	
	#store date and time in datetime object with hourly resolution	  
	time_hourly = [dt(year=int(mydate[0:4]),month=int(mydate[5:7]),day=int(mydate[8:10]),
					hour=int(myhour[0:2])) for mydate,myhour in zip(rawdata['date'],rawdata['hour'])]
	
	#compute hourly means as means of all time entries starting with same hour
	df = pd.DataFrame({'mytime': time_hourly, 'height': rawdata['val']})
	data_hourly = df.groupby('mytime',as_index=False).mean()
	
	#unique years in datetimes
	unique_years=np.array( list( {(i.year) for i in data_hourly['mytime']})) 
	
	#calculate annual means
	annual_means = data_hourly.groupby(data_hourly.mytime.dt.year).mean()
	
	#count number of observations in each year
	obs_pyear = data_hourly.groupby(data_hourly.mytime.dt.year,as_index=False).count()
	
	#select only years with enough obs
	goodYears_ind = obs_pyear>minDays*24 #nad hoc, no inferences made about how these hours are distributed in time
	good_years = unique_years[goodYears_ind.mytime.values]
	nYears = len(good_years)
	
	# Make sure there are 19 good year surrounding year 2000
	center_years = np.arange(center_year-9, center_year+10)
	if np.sum(np.isin(center_years, good_years)) >= 16:
		pass_centeryears = True
	
	# Does this file pass the minimum data check
	if nYears >= minYears:
		pass_nyears = True
	
	# Does everything pass?
	if not (pass_centeryears and pass_nyears):
		raise Exception("File does not pass data constraint")
	
	# loop over unique years with enough observations and compute anomalies wrt annual mean
	good_data=[]
	for myyear in good_years:
		yearly_anom = data_hourly[data_hourly['mytime'].dt.year == myyear ]
		yearly_anom.height -= annual_means.loc[myyear].height #subtract annual means
		
		good_data.append(yearly_anom)
		
	good_data = pd.concat(good_data)
	good_data.reset_index(drop=True,inplace=True)
	
	# compute Peak-over-Threshold (POT) height and fetch observed heights exceeding POT
	pot_threshold = np.percentile(good_data.height.values,pctPot)
	
	extremes = good_data[good_data['height'] > pot_threshold]
	extremes.reset_index(drop=True,inplace=True)
	extremes = extremes.sort_values(by=['height'], ascending=False, kind='mergesort') #note mergesort required to get same output as matlab, this has stable ordered duplicate values
	
	#compute range of pots
	pot_vals = np.percentile(good_data.height.values, gpd_pot_threshold)
	
	##### DECLUSTER EXTREME EVENTS
	# initialize declustered extremes with first extreme
	decl_extremes = extremes.iloc[[0]] #double [] makes pandas dataframe instead of series
	
	#check iteratively if new extreme is too close (in time) to declustered extremes
	for i in range(1,len(extremes)):	
		next_extreme = extremes.iloc[[i]] #pick new extreme
		next_extreme.reset_index(drop=True,inplace=True) #reset indices
		
		timediffs = decl_extremes['mytime'] - next_extreme.mytime[0] #calculate time difference between declustered extremes and next extreme
	
		# if none of the current extremes are less than 'cluster_lim' hours away from next extreme
		if all( abs(timediffs / np.timedelta64(3600, 's')) >= cluster_lim): 
			decl_extremes = pd.concat([decl_extremes,next_extreme]) #add next extreme
	
	decl_extremes.reset_index(drop=True,inplace=True)
	
	# calculate return periods in years	   
	decl_extremes['return_period'] = [(1 + len(good_data)/365.25/24)/ i for i in range(1,len(decl_extremes)+1)]	  
		
	# Return the information for this station
	return(station_name, station_lat, station_lon, pot_vals, nYears, good_data, decl_extremes)
	


def preprocess_gesla(gesladir, minDays, minYears, match_limit, center_year, pctPot, gpd_pot_threshold, cluster_lim, site_lats, site_lons, site_ids, pipeline_id):
	
	# Extract the gesla station information
	(station_names, station_lats, station_lons, station_filenames) = extract_gesla_locations(gesladir)
	
	# Match the target lats/lons with their nearest stations
	min_idx = [mindist(x,y,station_lats,station_lons, match_limit) for x,y in zip(site_lats, site_lons)]
	
	# Generate a list of input files and matched IDSs for the matched tide gauge data
	matched_filenames = []
	matched_ids = []
	for i in np.arange(len(min_idx)):
		if min_idx[i] is not None:
			matched_filenames.append([station_filenames[x] for x in min_idx[i]])
			matched_ids.append(site_ids[i])
	
	# If no matches are found, quit with error
	if not matched_filenames:
		raise Exception("No matches found within {} degrees for provided lat/lon list".format(match_limit))	
	
	# Define the output directory
	outdir = os.path.dirname(__file__)
	
	# Initialize station data dictionary to hold matched location data
	station_data = {}
	
	# Initialize variables to track the files that have been tested
	pass_files = {}
	fail_files = []
	
	# Loop over the matched files
	for i in np.arange(len(matched_filenames)):
		
		# This ID
		this_id = matched_ids[i]
		this_id_passed = False
		
		# Loop over the files within the match radius for this location
		for this_file in matched_filenames[i]:
			
			# This file was tested and passed
			if np.isin(this_file, list(pass_files.keys())):
				print("{0} already PASSED a previous check on data constraints. Mapping site ID {1} to site ID {2}.".format(this_file, pass_files[this_file], this_id))
				station_data[this_id] = station_data[pass_files[this_file]]
				this_id_passed = True
			
			# If this ID already has data, skip to the next ID
			if this_id_passed:
				continue
			
			# Move on to the next file if this file was already tested and failed
			if np.isin(this_file, fail_files):
				print("{0} already FAILED previous check on data constraints. Continuing with the next file".format(this_file))
				continue
			
			# This file has not been tested yet, go ahead and try to load it in.
			try:
				# Load this tide gauge data
				(match_name, match_lat, match_lon, pot_vals, nyears, obs, decl_extremes) = load_tg(os.path.join(gesladir, this_file), minDays, minYears, center_year, pctPot, gpd_pot_threshold, cluster_lim)
				
				# Put this entry into the station data dictionary
				station_data[this_id] = {'station_name':match_name, 'station_id':this_id, 'lon':match_lon, 'lat':match_lat, 'pot_vals':pot_vals, 'nyears':nyears, 'obs':obs, 'decl_extremes': decl_extremes}
				
				# This file passed, add it to the pass file dictionary
				pass_files[this_file] = this_id
				
				# This ID has data, set the flag to move onto the next ID
				this_id_passed = True
			
			except:
				# This file failed, add it to the fail list and continue with the next file
				print("{} did not pass the data constraints. Moving on to the next file.".format(this_file))
				fail_files.append(this_file)
				continue
		
		# Let the user know we didn't find a file that passes the data constraints
		if not this_id_passed:
			print("No locations within {0} degrees pass the data constraints for ID {1}.".format(match_limit, this_id))
	
	# Exit with error if we didn't find any gesla location suitable
	if not station_data:
		raise Exception("No data found for any requested locations")
	
	# Write station data to an output pickle
	outfile = open(os.path.join(outdir, "{}_station_data.pkl".format(pipeline_id)), 'wb')
	pickle.dump(station_data, outfile)
	outfile.close()
	
	return(list(station_data.keys()))




if __name__ == "__main__":
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the extreme sea-level workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--minDays', help="Minimum number of days in a valid year [default=250]", type=int, default=250)
	parser.add_argument('--minYears', help="Minimum number of years available [default=20]", type=int, default=20)
	parser.add_argument('--match_lim', help="Radius around requested locations to find a matching tide gauge in GESLA database", type=float, default=0.1)
	parser.add_argument('--center_year', help="This year +/- 9 years for centering data must be available", type=int, default=2000)
	parser.add_argument('--pctPot', help="Percentile for Point Over Threshold analysis [default=95]", type=int, default=95)
	parser.add_argument('--gpd_pot_threshold', help="Percentile for GPD analysis [default=99.7]", type=float, default=99.7)
	parser.add_argument('--cluster_lim', help="Maximum number of hours that define a cluster for extreme events [default=72]", type=int, default=72)
	parser.add_argument('--total_localsl_file', help="Total localized sea-level projection file. Site lats/lons are taken from this file and mapped to the GESLA database", default="total-workflow_localsl.nc")
	parser.add_argument('--target_years', help="Comma-delimited list of projection years of interest (e.g. 2050,2100)", default="2100")
	parser.add_argument('--gesla_dir', help="Directory containing GESLA database", default=os.path.join(os.path.dirname(__file__), "gesla_data"))
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Convert target_years argument from comma-delimited string to integer list
	target_years = [int(x) for x in args.target_years.split(",")]
	
	# Load the sea-level projection file
	nc = Dataset(args.total_localsl_file, 'r')
	
	#Extract data from NetCDF file
	site_lats = nc.variables['lat'][:].data
	site_lons = nc.variables['lon'][:].data
	site_ids = nc.variables['locations'][:].data
	proj_yrs = nc.variables['years'][:].data
	#proj_qnts = nc.variables['quantiles'][:].data
	#proj_slc_qnts = nc.variables['localSL_quantiles'][::,::,::].data
	proj_slc = nc.variables['sea_level_change'][::,::,::].data
	nc.close()
	
	# Make sure lats and lons are equal in length
	if not len(site_lats) == len(site_lons):
		raise Exception("Number of latitudes and longitudes not equal ({0} != {1})".format(len(site_lats), len(site_lons))) 
	
	# Test to make sure the list of site ids matches the length of the lats/lons
	if not len(site_lats) == len(site_ids):
		raise Exception("Number of site IDs not equal to number of locations provided in lat/lon list ({0} != {1})".format(len(site_ids), len(site_lats))) 
	
	# Find the target years that overlap the projection years
	slproj_targ_year_idx = np.flatnonzero(np.isin(proj_yrs, target_years))
	if len(slproj_targ_year_idx) == 0:
		raise Exception("Cannot find target years in projection years")
	
	# Define the gesla data directory
	gesladir = args.gesla_dir

	# Analysis options
	minDays = args.minDays
	minYears = args.minYears
	center_year = args.center_year
	pctPot = args.pctPot
	gpd_pot_threshold = args.gpd_pot_threshold
	cluster_lim = args.cluster_lim
	
	# Maximum distance from match allowed
	match_limit = args.match_lim
	
	# Run the GESLA preprocessing and return the matching site IDs
	matched_ids = preprocess_gesla(gesladir, minDays, minYears, match_limit, center_year, pctPot, gpd_pot_threshold, cluster_lim, site_lats, site_lons, site_ids, args.pipeline_id)
	
	# Find the site IDs that got a match in the GESLA database
	slproj_matched_ids_idx = np.flatnonzero(np.isin(site_ids, matched_ids))
	
	# Build the sea-level projection output dictionary
	slproj_output = {}
	for this_id_idx in slproj_matched_ids_idx:
		
		# Extract just the projections that match the target years and matched site ids
		# Note: Typecast the target year "idx" variable as a list in order to avoid singleton problems later.
		proj_slc_subset = proj_slc[::,[slproj_targ_year_idx],this_id_idx]
	
		# Generate the sea-level projection output dictionary
		slproj_output[site_ids[this_id_idx]] = {'site_lat': site_lats[this_id_idx], 'site_lon': site_lons[this_id_idx],\
			'site_id': site_ids[this_id_idx], 'proj_years': proj_yrs[slproj_targ_year_idx], \
			'proj_slc': proj_slc_subset}
	
	# Generate the sea-level projection output file
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_slproj_data.pkl".format(args.pipeline_id)), 'wb')
	pickle.dump(slproj_output, outfile)
	outfile.close()
	
	# Write the configuration output file
	config_output = {'minDays': minDays, 'minYears': minYears, 'match_lim': match_limit, \
		'center_year': center_year, 'pctPot': pctPot, 'cluster_lim': cluster_lim, \
		'total_localsl_file': args.total_localsl_file, 'target_years': target_years, \
		'gesla_dir': gesladir, 'pipeline_id': args.pipeline_id, 'gpd_pot_threshold': gpd_pot_threshold}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_config_data.pkl".format(args.pipeline_id)), 'wb')
	pickle.dump(config_output, outfile)
	outfile.close()
	
	exit()
