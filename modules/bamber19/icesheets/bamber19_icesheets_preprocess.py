import argparse
import os
import sys
import pickle as p
import numpy as np
import scipy.io
from prelib import *

''' bamber19_preprocess_icesheets.py

This runs the preprocessing stage for the Bamber et al. 2019 ice sheet component of the IPCC AR6
workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code

'''

def bamber19_preprocess_icesheets(pyear_start, pyear_end, pyear_step, baseyear, scenario, pipeline_id, climate_data_file = ''):
	# Load the preprocess library from factslib/prelib.py
	prelib = PreProcess()

	# Define the target years
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)

	# Load the data
	filename = os.path.join(os.path.dirname(__file__), "SLRProjections190726core_SEJ_full.mat")
	mat = scipy.io.loadmat(filename)
	
	if len(climate_data_file) > 0:
		wais_sampsH, eais_sampsH, gis_sampsH = prelib.ExtractSamples(mat, 'corefileH', targyears, baseyear)
		wais_sampsL, eais_sampsL, gis_sampsL = prelib.ExtractSamples(mat, 'corefileL', targyears, baseyear)
        
		prelib.OutputDataAll(pipeline_id, eais_sampsH, wais_sampsH, gis_sampsH,  eais_sampsL, wais_sampsL, gis_sampsL, scenario, targyears, baseyear)
		
	else:
		scenario_map = {"rcp85": 'corefileH', 
				        "rcp26": 'corefileL',
						"tlim2.0win0.25": 'corefileL', 
						"tlim5.0win0.25": 'corefileH'
						}
		
		this_corefile = scenario_map[scenario]
		wais_samps, eais_samps, gis_samps = prelib.ExtractSamples(mat, this_corefile, targyears, baseyear)
		prelib.OutputDataScen(pipeline_id, eais_samps, wais_samps, gis_samps, scenario, targyears, baseyear)


if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the IPCC AR6 Bamber et al. 2019 ice sheet workflow",\
	epilog="Note: This is meant to be run as part of the IPCC AR6 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--pyear_start', help="Projection year start [default=2020]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Projection year end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Projection year step [default=10]", default=10, type=int)
	parser.add_argument('--scenario', help="Emissions scenario of interest [default=rcp85]", type=str, default="rcp85")
	parser.add_argument('--baseyear', help="Year to which projections are referenced [default = 2000]", default=2000, type=int)
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str, default="")
	
	# Parse the arguments
	args = parser.parse_args()
	
	if len(args.climate_data_file) == 0:
		bamber19_preprocess_icesheets(args.pyear_start, args.pyear_end, args.pyear_step, args.baseyear, args.scenario, args.pipeline_id)
	else:
		bamber19_preprocess_icesheets(args.pyear_start, args.pyear_end, args.pyear_step, args.baseyear, args.scenario, args.pipeline_id, args.climate_data_file)

	# Done
	sys.exit()