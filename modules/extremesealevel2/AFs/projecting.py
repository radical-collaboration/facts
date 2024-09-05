'''
@author: Tim Hermans
t(dot)h(dot)j(dot)hermans@uu(dot)nl
'''
import numpy as np

def compute_AFs(f,z_hist,slr,refFreq):
    '''computes amplification factor based on historical return curve(s) "z_hist" = F("f"), sea-level projections "slr" at slr.years and reference frequency "refFreq"'''    
    i_ref = np.argmin(np.abs(f-refFreq)) #find index of f closest to refFreq
    
    if (z_hist.ndim == 1) & (len(z_hist)==len(f)): #expand z_hist if only dimension is f (no samples)
        if np.isscalar(slr)==False:
            z_hist = np.repeat(z_hist[:,None],len(slr.samples),axis=1)
    
    i_z_min = np.nanargmin(z_hist,axis=0)[0] #find index of lowest defined z (highest defined frequency); this determines upper limit of AF. Same for all samples because no uncertainty in lower bound is considered.
    
    #z_fut = z_hist+slr
    z_hist_ = np.repeat(z_hist[...,None],len(slr.years),axis=-1) #expand z_hist matrix 
    z_fut = z_hist_ + slr.expand_dims({'f':f}) #compute future return curves by adding SLR for each frequency
    
    refZ_hist = z_hist_[i_ref] #historical reference height samples corresponding to refFreq
    
    if np.isnan(refZ_hist).all(): # when refFreq is higher or lower than where z_hist is defined (for all samples or for none because lower bound has no uncertainty), return nans   
        return [np.nan],np.nan,z_fut 
        
    iRefZ_hist_in_z_fut = np.nanargmin(np.abs(z_fut - refZ_hist),axis=0) #find future frequency corresponding to those heights
    
    AF = f[iRefZ_hist_in_z_fut]/refFreq #divide future frequency correspondong to refZ_hist by historical reference frequency (will evaluate to nan if refFreq=nan)
    
    max_AF = f[i_z_min]/f[i_ref] #keep track of what AF could maximally be given the return curves
    
    return AF,max_AF,z_fut


def compute_AF_timing(f,z_hist,slr,refFreq,AF):
    if AF<=1:
        print("Warning: requested AF<=1; returning projected timing=nan")
        return np.zeros(len(slr.samples)) * np.nan
    
    if np.isnan(refFreq): #refFreq undefined, for instance if location not sufficiently near protection information
        return np.zeros(len(slr.samples)) * np.nan
    
    slr_years=slr.years.values
    
    i_ref = np.argmin(np.abs(f-refFreq)) #find index of f closest to refFreq
    i_ref_AF = np.argmin(np.abs(f-refFreq*AF)) #find index of f closest to AF * refFreq
   
    if (z_hist.ndim == 1) & (len(z_hist)==len(f)): #expand z_hist if only dimension is f (no samples)
        z_hist = np.repeat(z_hist[:,None],len(slr.samples),axis=1)
        
    req_slr = z_hist[i_ref] - z_hist[i_ref_AF] #sea-level rise required to go from refFreq to AF*refFreq (#if one of these is undefined, that's the case for all samples?)
    
    if np.isnan(req_slr).all(): #reference return height and/or target return height not supported by z_hist
        return np.zeros(len(slr.samples)) * np.nan
    
    #find the first timestep at which SLR > required SLR:
    slr_minus_required = slr - np.repeat(req_slr[:,np.newaxis],len(slr_years),axis=1) #projected minus required sea-level rise at each timestep
    slr_minus_required = np.where(slr_minus_required<0,999,slr_minus_required) #minimize only consider where SLR > required SLR (so set negative to 999)
    
    imin = np.nanargmin(slr_minus_required,axis=-1) #find earliest timestep where SLR > required SLR

    #interpolate between earliest timestep and timestep before that (but use earliest timestep if it is the first timestep of the projections):
    slr_years_ = np.hstack((slr_years[0],slr_years)) #copy SLR info and preprend first timestep
    slr_ = np.hstack((np.tile(slr[:, [0]], 1),slr))
    
    #calculate dT/dZ slope
    timing = slr_years[imin]
    timing_tmin1 = slr_years_[imin]
    
    dTime = timing - timing_tmin1
    dZ = slr_[np.arange(len(slr.samples)),imin+1] - slr_[np.arange(len(slr.samples)),imin]
    
    #interpolate by subtracting dT/dZ * (SLR-required_SLR) from earliest timestep with SLR>required SLR
    timing = timing - dTime/(dZ+1e-8) * slr_minus_required[np.arange(len(slr.samples)),imin]
    
    timing = np.where(slr.isel(years=-1)<req_slr,9999,timing) #where not before 2150, set to 9999 (so that it counts as > last timestep in computing quantiles)
    
    return timing

