import os
import sys
import argparse
import pickle
import numpy as np
import fair
import pandas as pd


# Function that prepares the RCMIP emissions data
def prep_rcmip_emissions_conc(scenario, rcmip_file):

	from fair.constants import molwt
	from fair.constants.general import ppm_gtc
	from fair.constants.general import EARTH_RADIUS, SECONDS_PER_YEAR

	NTOA_ZJ = 4 * np.pi * EARTH_RADIUS**2 * SECONDS_PER_YEAR * 1e-21

	emis_all = pd.read_csv(rcmip_file)
	expt = scenario
	emis_subset = emis_all[(emis_all['Scenario']==expt)&(emis_all['Region']=='World')]

	species=['CO2|MAGICC Fossil and Industrial','CO2|MAGICC AFOLU','CH4','N2O','Sulfur','CO','VOC','NOx','BC','|OC','NH3',
			 'CF4','C2F6','C6F14','HFC23','HFC32','HFC4310mee','HFC125','HFC134a','HFC143a',
			 'HFC227ea','HFC245fa','SF6','CFC11','CFC12','CFC113','CFC114','CFC115','CCl4','CH3CCl3','HCFC22',
			 'HCFC141b','HCFC142b','Halon1211','Halon1202','Halon1301','Halon2402','CH3Br','|CH3Cl']
	emis = np.zeros((751,40))
	emis[:,0] = np.arange(1750,2501)
	for ie, specie in enumerate(species):
		try:
			tmp = emis_subset[emis_subset.Variable.str.endswith(specie)].loc[:,'1750':'2500'].values.squeeze()
			emis[:,ie+1] = pd.Series(tmp).interpolate().values
		except:
			emis[:,ie+1] = 0
	tmp = emis_subset[emis_subset['Variable']=='Emissions|NOx|MAGICC Fossil and Industrial|Aircraft'].loc[:,'1750':'2500'].values.squeeze()
	aviNOx = pd.Series(tmp).interpolate().values
	aviNOx_frac = aviNOx/emis[:,8]
	#
	unit_convert = np.ones(40)
	unit_convert[1]=0.001 * molwt.C/molwt.CO2
	unit_convert[2]=0.001 * molwt.C/molwt.CO2
	unit_convert[4]=0.001 * molwt.N2/molwt.N2O
	unit_convert[5]=molwt.S/molwt.SO2
	unit_convert[8]=molwt.N/molwt.NO2

	emis = emis * unit_convert

	## The section below will handle rcmip concentration files, which are not used
	## in this module.

	#conc_all = pd.read_csv("./rcmip/rcmip-concentrations-annual-means-v5-1-0.csv")
	#conc_subset = conc_all[(conc_all['Scenario']==expt)&(conc_all['Region']=='World')]
	#gases=['CO2','CH4','N2O','CF4','C2F6','C6F14','HFC23','HFC32','HFC4310mee','HFC125','HFC134a','HFC143a',
	#       'HFC227ea','HFC245fa','SF6','CFC11','CFC12','CFC113','CFC114','CFC115','CCl4','CH3CCl3','HCFC22',
	#       'HCFC141b','HCFC142b','Halon1211','Halon1202','Halon1301','Halon2402','CH3Br','CH3Cl']
	#conc = np.zeros((751,31))
	#for ig, gas in enumerate(gases):
	#    try:
	#        tmp = conc_subset[conc_subset.Variable.str.endswith(gas)].loc[:,'1750':'2500'].values.squeeze()
	#        conc[:,ig] = pd.Series(tmp).interpolate().values
	#    except:
	#        conc[:,ig] = 0

	#return(emis, conc)
	return (emis)


def fair_preprocess_temperature(scenario, rcmip_file, pipeline_id):

	# Definitions
	REFERENCE_YEAR = 1750

	# Load the emissions data for this scenario
	emis = prep_rcmip_emissions_conc(scenario, rcmip_file)

	# Save the preprocessed data to a pickle
	output = {"emis": emis, "REFERENCE_YEAR": REFERENCE_YEAR, "scenario": scenario, "rcmip_file": rcmip_file}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_preprocess.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile, protocol=-1)
	outfile.close()


	return(None)


if __name__ == "__main__":

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the preprocessing stage for the FAIR AR6 temperature module",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--rcmip_file', help="Full path to RCMIP emissions file", default="./rcmip/rcmip-emissions-annual-means-v5-1-0.csv")
	parser.add_argument('--scenario', help="SSP Emissions scenario",  default="ssp585")

	# Parse the arguments
	args = parser.parse_args()

	# Run the code
	fair_preprocess_temperature(scenario=args.scenario, rcmip_file=args.rcmip_file, pipeline_id=args.pipeline_id)


	sys.exit()
