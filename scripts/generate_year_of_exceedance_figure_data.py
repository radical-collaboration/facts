import numpy as np
import os
import sys
import re
import fnmatch
import argparse
from netCDF4 import Dataset
from sample_from_quantiles import *


def YearsOfExceedance(sample, years, heights):

	# Get the sorted height indices for this sample
	#hsort_idx = np.argsort(sample)

	# Find the running maximum of this sample
	sample_max = np.maximum.accumulate(sample)

	# Find the years at which these heights are exceeded
	#exyears = np.interp(heights, sample[hsort_idx], years[hsort_idx])
	exyears = np.interp(heights, sample_max, years)

	# Return the years
	return(np.around(exyears))


def ExtractData(infile, outfile, heights):

	# Quantiles to extract from the data
	q = np.array([0.05, 0.17, 0.5, 0.83, 0.95])

	# Open the netCDF file
	nc = Dataset(infile, 'r')

	# Get the year and sample data
	ncyears = nc.variables['years'][...]
	#nclocalsl = nc.variables['localSL_quantiles'][...]
	nclocalsl = nc.variables['sea_level_change'][...]
	ncquantiles = nc.variables['quantiles'][...]
	#ncsiteids = nc.variables['id'][...]
	ncsiteids = nc.variables['locations'][...]
	nsamps = 5000

	# Get the attribute information
	nc_description = nc.__dict__['description']
	nc_history = nc.__dict__['history']
	nc_source = nc.__dict__['source']

	# Close the input netcdf file
	nc.close()

	# Initialize the milestone years output data structure
	exyears_q = []

	# Create the samples from the local distributions
	for j in np.arange(len(ncsiteids)):

		# Initialize samps
		samps = []

		# Loop over the years
		for i in np.arange(len(ncyears)):

			# Set the seed to ensure consistent time series
			np.random.seed(8071)

			# Generate the samples
			#samps.append(sample_from_quantiles(nclocalsl[:,j,i], ncquantiles, nsamps))
			samps.append(sample_from_quantiles(nclocalsl[:,i,j], ncquantiles, nsamps))

		# Convert samps to a numpy array
		samps = np.array(samps)

		# Find the exceedance years for this set of samples and heights
		this_exyears = np.apply_along_axis(YearsOfExceedance, 0, samps, years=ncyears, heights=heights)

		# Calculate the quantiles
		this_exyears_q = np.quantile(this_exyears, q, axis=1)

		# Append these quantiles to the output data structure
		exyears_q.append(this_exyears_q)

	# Convert the output data structure to a numpy array
	exyears_q = np.array(exyears_q)

	# Write the combined projections to a netcdf file
	rootgrp = Dataset(outfile, "w", format="NETCDF4")

	# Define Dimensions
	height_dim = rootgrp.createDimension("heights", len(heights))
	q_dim = rootgrp.createDimension("quantiles", len(q))
	site_dim = rootgrp.createDimension("sites", len(ncsiteids))

	# Populate dimension variables
	height_var = rootgrp.createVariable("heights", "i2", ("heights",))
	height_var.units = "mm"
	q_var = rootgrp.createVariable("quantiles", "f4", ("quantiles",))
	site_var = rootgrp.createVariable("sites", "i4", ("sites",))

	# Create a data variable
	exyearsq = rootgrp.createVariable("exceedance_years", "i2", ("sites", "quantiles", "heights"), zlib=True, complevel=4)

	# Assign attributes
	rootgrp.description = nc_description
	rootgrp.history = nc_history
	rootgrp.source = nc_source
	exyearsq.units = "-"

	# Put the data into the netcdf variables
	height_var[:] = heights
	q_var[:] = q
	site_var[:] = ncsiteids
	exyearsq[:,:] = exyears_q

	# Close the netcdf
	rootgrp.close()



def main(indir, outdir):

	# Get the list of workflows
	workflows = os.listdir(indir)
	workflows.sort()
	workflows = fnmatch.filter(workflows, "[!.]*")

	# Heights of interest
	heights = np.arange(100, 5001, 50)

	# Loop over the workflows
	for workflow in workflows:

		# Update user on which workflow is being processed
		print("Processing workflow: {}".format(workflow))

		# Get the scenarios for this workflow
		scenarios = os.listdir(os.path.join(indir, workflow))
		scenarios.sort()
		scenarios = fnmatch.filter(scenarios, "[!.]*[!.][!t][!x][!t]")

		# Loop through the scenarios for this workflow
		for scenario in scenarios:

			# Define the path to this workflow-scenario total file
			total_file = os.path.join(indir, workflow, scenario, "total-workflow.nc")

			# Generate the figure data for this file
			outfile = os.path.join(outdir, workflow, scenario, "{}_{}_milestone_figuredata.nc".format(workflow, scenario))
			ExtractData(total_file, outfile, heights)

	# Done
	return(None)


if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Produce the figure data files for the milestone plot for local projection workflows")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--indir', help="Input directory with global projection workflows")
	parser.add_argument('--outdir', help="Output directory for the figure data files")

	# Parse the arguments
	args = parser.parse_args()

	# Run the code with these options
	main(args.indir, args.outdir)



	sys.exit()
