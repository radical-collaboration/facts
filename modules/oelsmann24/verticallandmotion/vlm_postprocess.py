import numpy as np 
import xarray as xr 
import pandas as pd 
import pickle as p
import os
import sys
import re
import argparse
from read_locationfile import ReadLocationFile
from scipy.spatial import cKDTree
from typing import Optional, Tuple
from scipy.interpolate import interp1d
from scipy.stats import norm
import time



''' vlm_preprocess.py

This runs the projection stage for the global vertical land motion component.

Parameters:
nsamps = Number of samples to draw
rng_seed = Seed value for the random number generator
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code

'''




def match_coord_vectors(site_lats, site_lons,lat_flat, lon_flat):
    """
    Match site coordinates to nearest grid points on a sphere.

    Parameters
    ----------
    site_lats, site_lons : array-like, shape (M,)
        Site latitudes/longitudes (deg).
    lat_flat, lon_flat   : array-like, shape (N,)
        Flattened grid lat/lon (deg), typically from 2D meshgrid.ravel().

    Returns
    -------
    idx : np.ndarray, shape (M,)
        Indices into (lat_flat, lon_flat) for the nearest grid point to each site.
    dist_km : np.ndarray, shape (M,)
        Great-circle distance (km) from each site to its matched grid point.
    """
    def sph2cart(lat_deg, lon_deg):
        lat = np.deg2rad(lat_deg)
        lon = np.deg2rad(lon_deg)
        coslat = np.cos(lat)
        x = coslat * np.cos(lon)
        y = coslat * np.sin(lon)
        z = np.sin(lat)
        return np.column_stack((x, y, z))

    # build tree from flattened grid points
    pts = sph2cart(lat_flat, lon_flat)          # shape (N,3)
    tree = cKDTree(pts)

    # convert sites and query
    sites = sph2cart(site_lats, site_lons)     # shape (M,3)
    d_euclid, idx = tree.query(sites, k=1)     # d_euclid = chord distances

    # convert to great-circle distance (using dot-product is numerically stable)
    # get the matched points' unit vectors
    matched = pts[idx]                          # shape (M,3)
    dot = np.einsum('ij,ij->i', matched, sites) # dot = cos(angle)
    dot = np.clip(dot, -1.0, 1.0)
    angle_rad = np.arccos(dot)                  # central angle (radians)

    EARTH_KM = 6371.0
    dist_km = EARTH_KM * angle_rad

    # results:
    # idx      -> indices into lat_flat/lon_flat nearest to each site
    # dist_km  -> great-circle distance in km for each site
    return idx,dist_km


def map_datasets_on_locationfile(
    dis_sel_2d,
    GIA_VLM_int,
    dis_sel_2d_grw_use,
    decay_vec_vlm,
    site_lats,
    site_lons,
    targyears,
    baseyear
):
    """
    Map gridded VLM/GIA fields to site locations and time grid.

    Parameters
    ----------
    dis_sel_2d : xr.Dataset
        Contains VLM trend/uncertainty on a (lat, lon) grid.
    GIA_VLM_int : xr.Dataset
        GIA fields on the same (lat, lon) grid. Must contain 'rsl_mean' and 'rsl_sterr'.
    dis_sel_2d_grw_use : xr.Dataset
        GRW ensemble std over (lat, lon, time) used to mimic present-day variability.
    decay_vec_vlm : np.ndarray
        Spatial-temporal weights (lat, lon, time) for blending; will be reshaped.
    site_lats, site_lons : array-like
        Site coordinates (deg).
    targyears : array-like
        Target projection years (absolute years).
    baseyear : int
        Base/reference year.

    Returns
    -------
    tuple of arrays:
        VLM_REC_trend_mapped                : (M,)   VLM→RSL trend at sites (mm/yr → mm after time mult)
        VLM_REC_trend_un_mapped             : (M,)   Uncertainty of VLM trend (mm/yr)
        VLM_GRW_trend_un_mapped             : (M,T5) Original GRW (5-year steps) mapped to sites (mm/yr)
        GIA_VLM_int_rsl_mean_mapped         : (M,)   GIA mean RSL trend (mm/yr)
        GIA_VLM_int_rsl_sterr_mapped        : (M,)   GIA sterr (mm/yr)
        VLM_GRW_trend_un_mapped_interp      : (M,Y)  GRW interpolated to targyears (mm/yr)
        WEIGHTS                              : (M,Y)  Time-dependent weights for each site

    Notes
    -----
    - VLM trend to RSL: uses `VLM_trend_coefficient_mean` (mm/yr) * (-1) + `GIA_VLM_int['asl']`.
    - GRW time axis is originally in 5-year steps; interpolated to `targyears`.
    - `decay_vec_vlm` is reshaped to (Npoints, Ny) then indexed by nearest-neighbor mapping.
    """

    # flatten datasets
    # VLM (GPS and reconstruction) trends + uncertainties
    VLM_REC_trend = (dis_sel_2d['VLM_trend_coefficient_mean'].fillna(0).values*1000*-1+GIA_VLM_int['asl'].values).flatten() # impact on RSL
    VLM_REC_trend_un = (dis_sel_2d['VLM_trend_coefficient_uncertainty'].fillna(0).values*1000).flatten()


    #dis_sel_2d_grw_use_matrix = dis_sel_2d_grw_use[:,:,:31]*decay_vec_vlm*1000.
    # VLM (GPS and reconstruction) GRW
    dim0, dim1, dim2 = dis_sel_2d_grw_use.dims
    VLM_GRW_trend_un= dis_sel_2d_grw_use.stack(points=(dim0, dim1)).T.values[:,1:] *1000.

    # VLM (GIA) trends + uncertainties
    GIA_VLM_int_rsl_mean = GIA_VLM_int['rsl_mean'].values.flatten()
    GIA_VLM_int_rsl_sterr = GIA_VLM_int['rsl_sterr'].values.flatten()

    # Interpolate
    lat1d = dis_sel_2d['lat'].values        # shape (720,)
    lon1d = dis_sel_2d['lon'].values        # shape (1440,)

    lon2d, lat2d = np.meshgrid(lon1d, lat1d)   # lon varies fastest in each row

    lat_flat = lat2d.ravel()   # length = 720*1440
    lon_flat = lon2d.ravel()

    idx,dist_km = match_coord_vectors(site_lats, site_lons,lat_flat, lon_flat)


    VLM_REC_trend_mapped = VLM_REC_trend[idx]
    VLM_REC_trend_un_mapped = VLM_REC_trend_un[idx]
    VLM_GRW_trend_un_mapped = VLM_GRW_trend_un[idx,:] # currently contains 60 timesteps, with 5 year stepsize
    GIA_VLM_int_rsl_mean_mapped = GIA_VLM_int_rsl_mean[idx]
    GIA_VLM_int_rsl_sterr_mapped = GIA_VLM_int_rsl_sterr[idx]
    n0, n1, n2 = decay_vec_vlm.shape
    WEIGHTS = decay_vec_vlm.reshape(n0 * n1, n2)[idx,:]
    # interpolate GRW timesteps on the new time vector

    arr = np.asarray(VLM_GRW_trend_un_mapped)   # shape (nspace, nt)
    nspace, nt = arr.shape

    orig_years = baseyear + np.arange(nt) * 5   # shape (nt,)

    # create interpolator along time axis (axis=1), allow extrapolation
    f = interp1d(orig_years, arr, axis=1, kind='linear',
                bounds_error=False, fill_value='extrapolate', assume_sorted=True)

    VLM_GRW_trend_un_mapped_interp = f(targyears)   # shape (nspace, ny)

    return VLM_REC_trend_mapped,VLM_REC_trend_un_mapped,VLM_GRW_trend_un_mapped,GIA_VLM_int_rsl_mean_mapped,GIA_VLM_int_rsl_sterr_mapped,VLM_GRW_trend_un_mapped_interp,WEIGHTS

def temp_decay(years,start_year,period1=2050,period2=2150):
    """
    Time-decay function blending VLM observations into GIA VLM trends.

    Parameters
    ----------
    years : array-like
        Absolute years (same grid as target years).
    start_year : int
        (Unused in current implementation; kept for signature stability)
    period1, period2 : int
        Blend window; cosine taper from 1 → 0 between [period1, period2].

    Returns
    -------
    decay_period_y : np.ndarray
        Cosine taper values during the blend window.
    decay_period_x : np.ndarray
        Years offset within the blend window (0..period2-period1).
    decay_vec : np.ndarray, shape (len(years),)
        Full time series of weights (1 before, cosine taper inside, 0 after).

    """
    decay_vec = np.ones(len(years))
    decay_period_x = years[(years>=period1) & (years<period2) ]
    decay_period_x=decay_period_x-decay_period_x[0]
    decay_period_y = ((np.cos(decay_period_x*2*np.pi/(2*(period2-period1)))+1)/2)
    decay_vec[(years>=period1) & (years<period2) ] = decay_period_y
    decay_vec[years>=period2]=0
    return decay_period_y,decay_period_x,decay_vec


def vlm_postprocess(
    nsamps,
    rng_seed,
    locationfile,
    baseyear,
    pyear_start,
    pyear_end,
    pyear_step,
    chunksize,
    pipeline_id
):
    """
    Run the VLM post-processing to produce local sea-level change samples.

    Steps
    -----
    1) Load `<pipeline_id>_data.pkl` (produced by preprocess stage)
    2) Read location file (ids, lats, lons) using `ReadLocationFile`
    3) Build target year vector; compute time-decay weights
    4) Map VLM & GIA fields to sites; interpolate GRW to target years
    5) Draw samples and construct (samples, years, locations) array
    6) Write NetCDF: `<pipeline_id>_localsl.nc`

    Outputs
    -------
    NetCDF file with variables:
      - sea_level_change[samples, years, locations] (mm)
      - lat[locations], lon[locations]
      - coords: years, locations, samples

    Notes
    -----
    - Uses two random permutations of inverse-normal quantiles for stratified sampling.
    - `ReadLocationFile` must return a tuple (_, site_ids, site_lats, site_lons).
    """
    # Read in the data from the preprocessing stage
    datafile = "{}_data.pkl".format(pipeline_id)
    try:
        f = open(datafile, 'rb')
    except:
        print("Cannot open datafile\n")
        sys.exit(1)

    # Extract the data from the file
    my_data = p.load(f)
    GIA_VLM_int = my_data['GIA_VLM_int']
    VLM_REC_MERGED = my_data['VLM_REC_MERGED']
    VLM_REC_MERGED_GRW = my_data['VLM_REC_MERGED_GRW']
    weights = my_data['weights']
    f.close()


    (_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)


    targyears = np.linspace(pyear_start,pyear_end,int((pyear_end-pyear_start)/pyear_step)+1)
    # targyears excludes baseyear here

    decay_period_y,decay_period_x,decay_vec = temp_decay(targyears,pyear_start,period1=2050,period2=2150)
    time_vector = targyears-baseyear

    # contains spatial and temporal weights to merge GIA and the other VLM observations
    decay_vec_vlm =(decay_vec * weights.values[:,:,None])
    #decay_vec_gia = (((decay_vec * w_d.values[:,:,None])*-1)+1)


    VLM_REC_trend_mapped,VLM_REC_trend_un_mapped,VLM_GRW_trend_un_mapped,GIA_VLM_int_rsl_mean_mapped,GIA_VLM_int_rsl_sterr_mapped,VLM_GRW_trend_un_mapped_interp,WEIGHTS = map_datasets_on_locationfile(VLM_REC_MERGED,GIA_VLM_int,VLM_REC_MERGED_GRW,decay_vec_vlm,site_lats, site_lons,targyears,baseyear)


    # Generate Samples
    rng = np.random.default_rng(rng_seed)
    rng2 = np.random.default_rng(rng_seed+1)
    x = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
    norm_inv = norm.ppf(x)
    norm_inv_perm = rng.permutation(norm_inv).astype(np.float32)
    norm_inv_perm2 = rng2.permutation(norm_inv).astype(np.float32)
    nc_missing_value = np.nan #np.iinfo(np.int16).min

    # Project
    years_mult = targyears - baseyear
    VLM_REC_trend_mapped_time = np.multiply.outer(years_mult,VLM_REC_trend_mapped)
    VLM_REC_trend_un_mapped_time = np.multiply.outer(years_mult,VLM_REC_trend_un_mapped)
    GIA_VLM_int_rsl_mean_mapped_time = np.multiply.outer(years_mult,GIA_VLM_int_rsl_mean_mapped)
    GIA_VLM_int_rsl_sterr_mapped_time = np.multiply.outer(years_mult,GIA_VLM_int_rsl_sterr_mapped) 


    sample_shape = VLM_REC_trend_mapped_time.shape   # (H, W)
    SAMPLES = np.empty((nsamps,) + sample_shape, dtype=np.int16)#*np.nan

    # Precompute / cast to float32 to reduce memory & speed up arithmetic
    VLM_REC_trend_mapped_time_ = VLM_REC_trend_mapped_time.astype(np.float32)
    VLM_REC_trend_un_mapped_time_ = VLM_REC_trend_un_mapped_time.astype(np.float32)
    gia_mean = GIA_VLM_int_rsl_mean_mapped_time.astype(np.float32)
    gia_sterr = GIA_VLM_int_rsl_sterr_mapped_time.astype(np.float32)
    grw = VLM_GRW_trend_un_mapped_interp.T.astype(np.float32)
    WEIGHTS_ = WEIGHTS.T.astype(np.float32)

    one_minus_weights = (WEIGHTS_ * -1.0) + 1.0   # shape (H, W)

    # Simple batch loop
    for start in range(0, nsamps, chunksize):
        end = min(start + chunksize, nsamps)
        b = end - start

        # shape (b, 1, 1) so broadcasting works with (H,W) arrays
        r1 = norm_inv_perm[start:end][:, None, None]
        r2 = norm_inv_perm2[start:end][:, None, None]

        # Broadcasted arithmetic -> shapes (b, H, W)
        VLM_obs_batch = WEIGHTS_[None, :, :] * ( VLM_REC_trend_mapped_time_[None, :, :] +
                                                VLM_REC_trend_un_mapped_time_[None, :, :] * r1 +
                                                grw[None, :, :] * r2 )

        # GIA_sample = (1 - WEIGHTS) * (gia_mean + gia_sterr*r1)
        GIA_batch = one_minus_weights[None, :, :] * ( gia_mean[None, :, :] +
                                                    gia_sterr[None, :, :] * r1 )

        # final samples: round then cast once
        VLM_batch = np.around(VLM_obs_batch + GIA_batch).astype(np.int16)

        SAMPLES[start:end, :, :] = VLM_batch
        print(start)


    ncvar_attributes = {"description": "Global VLM contributions to RSL change according to Global vlm workflow",
            "history": "Created " + time.ctime(time.time()),
            "source": "SLR Framework: Global vlm workflow",
            "scenario": 'NA',
            "baseyear": baseyear}
        

    vlm_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), SAMPLES, {"units":"mm", "missing_value":nc_missing_value}),
                            "lat": (("locations"), site_lats),
                            "lon": (("locations"), site_lons)},
        coords={"years": targyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)

    # Write the netcdf output file
    vlm_out.to_netcdf("{0}_localsl.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})
    return(None)



if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the oelsmann24 vertical land motion workflow",\
	epilog="Note: This is meant to be run as part of the oelsmann24 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator", default=1234, type=int)
	parser.add_argument('--baseyear', help="Base or reference year for projetions [default=2000]", default=2005, type=int)
	parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2100]", default=2150, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--chunksize', help="Number of locations to process at a time [default=50]", type=int, default=100)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Make sure the base year and target years are within data limits for this module
	if(args.baseyear < 2000):
		raise Exception("Base year cannot be less than year 2000: baseyear = {}".format(args.baseyear))
	if(args.baseyear > 2300):
		raise Exception("Base year cannot be greater than year 2300: baseyear = {}".format(args.baseyear))
	if(args.pyear_start < 2000):
		raise Exception("Projection year cannot be less than year 2000: pyear_start = {}".format(args.pyear_start))
	if(args.pyear_end > 2300):
		raise Exception("Projection year cannot be greater than year 2300: pyear_end = {}".format(args.pyear_end))

	# Make sure the target year stepping is positive
	if(args.pyear_step < 1):
		raise Exception("Projection year step must be greater than 0: pyear_step = {}".format(args.pyear_step))

	# Run the postprocessing stage
	vlm_postprocess(args.nsamps, args.seed, args.locationfile, args.baseyear, args.pyear_start, args.pyear_end, args.pyear_step, args.chunksize, args.pipeline_id)

	# Done
	exit()
	

	