import argparse
import os
import sys
import pickle as p
import numpy as np
import scipy.io

class PreProcess:

    def __init__(self) -> None:
        pass

    def ExtractSamples(self,mat, this_corefile, targyears, baseyear):
	    # Get the years available from the matlab core file
        mat_years = np.squeeze(mat[this_corefile][0,0][27])
	
        # Determine which matlab year indices match our target years
        mat_years_idx = np.isin(mat_years, targyears)
	
	    # Get the samples from the matlab core file
        samps = mat[this_corefile][0,0][21]
	
	    # Extract the samples for each ice sheet
        eais_samps = samps[:,20,:]
        wais_samps = samps[:,19,:]
        gis_samps = samps[:,18,:]
	
	    # Get the values for the baseyear of interest
        eais_refs = np.apply_along_axis(self.FindRefVals, axis=1, arr=eais_samps, years=mat_years, baseyear=baseyear)
        wais_refs = np.apply_along_axis(self.FindRefVals, axis=1, arr=wais_samps, years=mat_years, baseyear=baseyear)
        gis_refs = np.apply_along_axis(self.FindRefVals, axis=1, arr=gis_samps, years=mat_years, baseyear=baseyear)
	
	    # Center the projections to the reference period
        eais_samps -= eais_refs[:,np.newaxis]
        wais_samps -= wais_refs[:,np.newaxis]
        gis_samps -= gis_refs[:,np.newaxis]
	
	    # Subset for the target years
        eais_samps = eais_samps[:,mat_years_idx]
        wais_samps = wais_samps[:,mat_years_idx]
        gis_samps = gis_samps[:,mat_years_idx]

        return wais_samps, eais_samps, gis_samps
    
    def OutputDataAll(self,pipeline_id, eais_sampsH, wais_sampsH, gis_sampsH,  eais_sampsL, wais_sampsL, gis_sampsL, scenario, targyears, baseyear):

	    # Sum up the components to get total AIS samples
        ais_sampsH = eais_sampsH + wais_sampsH
        ais_sampsL = eais_sampsL + wais_sampsL

	    # Populate the output dictionary
        outdata = {'eais_sampsH': eais_sampsH, 
                   'wais_sampsH': wais_sampsH, 
                   'ais_sampsH': ais_sampsH,
                   'gis_sampsH': gis_sampsH, 
                   'eais_sampsL': eais_sampsL, 
                   'wais_sampsL': wais_sampsL,
                   'ais_sampsL': ais_sampsL, 
                   'gis_sampsL': gis_sampsL, 
                   'targyears': targyears, 
                   'baseyear': baseyear,
                   'scenario': scenario
                   }
	
	    # Define the data directory
        outdir = os.path.dirname(__file__)

	    # Write the rates data to a pickle file
        outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
        p.dump(outdata, outfile)
        outfile.close()

    def FindRefVals(self,timeseries, years, baseyear, append_yr=True):
	
	    # Append a zero to the beginning of the timeseries at year 2000
        if append_yr:
            timeseries = np.append(np.array([0.0]), timeseries)
            years = np.append(np.array([2000]), years)
	
	    # Interpolate to the appropriate base year
        ref_val = np.interp(baseyear, years, timeseries, left=0.0)
	
	    # Return the value
        return(ref_val)

