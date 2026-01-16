import numpy as np 
import xarray as xr 
import pandas as pd 
import pickle as p
import os
import sys
import re
import argparse
from scipy.spatial import cKDTree
from typing import Optional, Tuple
import copy


''' vlm_preprocess.py

This runs the post-processing stage for the global vertical land motion component.

Parameters:
nsamps = Number of samples to draw
rng_seed = Seed value for the random number generator
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code

'''








EARTH_RADIUS_KM = 6371.0088

# ── Helper functions: spherical geometry & kernels ─────────────────────────────
def _to_xyz(lon_deg, lat_deg):
    """Convert lon/lat (deg) to unit-sphere 3D coordinates for fast KDTree queries.

    Parameters
    ----------
    lon_deg, lat_deg : array-like
        Longitudes and latitudes in degrees.

    Returns
    -------
    np.ndarray of shape (N, 3)
        Unit vectors (x, y, z) on the sphere.
    """
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    cl = np.cos(lat)
    x = cl * np.cos(lon)
    y = cl * np.sin(lon)
    z = np.sin(lat)
    return np.column_stack((x, y, z))



def _arc_km_from_chord(chord):
    """Convert chord length on the unit sphere to great-circle arc length (km)."""
    arc_rad = 2.0 * np.arcsin(np.clip(chord * 0.5, 0.0, 1.0))
    return EARTH_RADIUS_KM * arc_rad

def _kernel_weights(dist_km, L_km, kernel="cos", exponent=1.0):
    """
    Compute radial distance weights on [0, L_km]; outside the radius -> 0.

    Parameters
    ----------
    dist_km : array-like
        Distances in km from the evaluation point.
    L_km : float
        Kernel radius/scale in km.
    kernel : {'cos', 'tri', 'gaussian'}
        - 'cos'     : raised cosine (Hann) on [0, L]
        - 'tri'     : linear taper to 0 at L
        - 'gaussian': exp(-(d/L)^2), truncated at 3*L for speed
    exponent : float
        Optional exponent applied to weights (w ** exponent).

    Returns
    -------
    np.ndarray
        Weights aligned with `dist_km`.
    """
    d = np.asarray(dist_km)
    w = np.zeros_like(d, dtype=float)

    if kernel == "cos":  # raised cosine (Hann) on [0,L]
        inside = d <= L_km
        w[inside] = 0.5 * (1.0 + np.cos(np.pi * d[inside] / L_km))
    elif kernel == "tri":  # linear taper
        inside = d <= L_km
        w[inside] = 1.0 - d[inside] / L_km
    elif kernel == "gaussian":
        w = np.exp(-(d / L_km) ** 2)
        inside = d <= (3 * L_km)           # cut tail for speed
        w[~inside] = 0.0
    else:
        raise ValueError("kernel must be one of {'cos','tri','gaussian'}")

    if exponent != 1.0:
        w = w ** float(exponent)
    return w

def radial_smoothing(
    ds: xr.Dataset,
    *,
    var_name: str = "trend",
    uncert_name: str = "trend_un",     # if missing, uncertainty weighting is skipped
    dim: str = "x",
    lon_name: str = "lon",
    lat_name: str = "lat",
    L_km: float = 150.0,               # radius/scale (km)
    kernel: str = "cos",
    exponent: float = 1.0,
    use_uncert: bool = True,
    min_neighbors: int = 1
) -> xr.Dataset:
    """
    Radially smooth `ds[var_name]` using distance-decay weights and optional 1/σ² weights.

    Output variables added to a copy of `ds`:
      - `<var_name>_s`    : uncertainty- and distance-weighted mean
      - `<var_name>_s_un` : weighted standard deviation among neighbors

    Notes
    -----
    - Neighbor search uses a KDTree on unit-sphere 3D vectors with a chord-length cutoff
      equivalent to the specified arc radius L_km.
    - If `uncert_name` is absent or `use_uncert=False`, only distance weights are used.
    """

    if var_name not in ds:
        raise KeyError(f"{var_name!r} not found in dataset")

    # Extract arrays
    lon = np.asarray(ds[lon_name].values)
    lat = np.asarray(ds[lat_name].values)
    y   = np.asarray(ds[var_name].values).astype(float)

    # Optional uncertainty weights (1/variance)
    if use_uncert and (uncert_name in ds):
        sigma = np.asarray(ds[uncert_name].values).astype(float)
        sigma = np.where((sigma > 0) & np.isfinite(sigma), sigma, np.nan)
        w_unc = 1.0 / (sigma ** 2)
        w_unc = np.where(np.isfinite(w_unc), w_unc, 0.0)
        used_uncert = True
    else:
        w_unc = None
        used_uncert = False

    # KDTree on the sphere
    xyz = _to_xyz(lon, lat)
    tree = cKDTree(xyz)

    # Radius in chord space for neighbor search
    arc = L_km / EARTH_RADIUS_KM
    chord_radius = 2.0 * np.sin(arc / 2.0)

    y_s   = np.full_like(y, np.nan, dtype=float)
    y_s_un = np.full_like(y, np.nan, dtype=float)

    neighbors_list = tree.query_ball_point(xyz, r=chord_radius)

    for i, idxs in enumerate(neighbors_list):
        if len(idxs) < min_neighbors:
            continue

        # Distances (km) from point i to its neighbors
        euclid = np.linalg.norm(xyz[idxs] - xyz[i], axis=1)
        d_km = _arc_km_from_chord(euclid)

        # Distance-decay weights
        w = _kernel_weights(d_km, L_km, kernel=kernel, exponent=exponent)

        # Combine with 1/sigma^2 if provided
        if w_unc is not None:
            w = w * w_unc[idxs]

        # Exclude NaNs and zero weights
        valid = np.isfinite(y[idxs]) & (w > 0)
        if not np.any(valid):
            continue

        wv = w[valid]
        yv = y[idxs][valid]

        # Weighted mean
        wsum = np.sum(wv)
        ybar = np.sum(wv * yv) / wsum
        y_s[i] = ybar

        # Weighted standard deviation of neighbors
        # sqrt( sum(w*(x-mu)^2) / sum(w) )
        y_s_un[i] = np.sqrt(np.sum(wv * (yv - ybar) ** 2) / wsum)

    # Package outputs
    out = ds.copy()
    base_coords = {dim: ds[dim], lon_name: ds[lon_name], lat_name: ds[lat_name]}

    out[f"{var_name}_s"] = xr.DataArray(
        y_s, dims=(dim,), coords=base_coords,
        attrs={
            "description": f"Radially smoothed {var_name} (weighted mean)",
            "L_km": L_km, "kernel": kernel, "exponent": exponent,
            "uncert_weighting": used_uncert,
        },
    )
    out[f"{var_name}_s_un"] = xr.DataArray(
        y_s_un, dims=(dim,), coords=base_coords,
        attrs={
            "description": f"Weighted standard deviation of neighbors for {var_name}_s",
            "L_km": L_km, "kernel": kernel, "exponent": exponent,
            "uncert_weighting": used_uncert,
            "formula": "sqrt(sum(w*(x-mu)^2)/sum(w)) with w = decay * (1/sigma^2)",
        },
    )
    return out



def couple_short(
    a: xr.DataArray,
    b: xr.DataArray,
    *,
    dim: str = "x",
    lon_name: str = "lon",
    lat_name: str = "lat",
    id_coord_a: Optional[str] = None,   # ID for a (if coord present; else row index)
    id_coord_b: Optional[str] = None,   # ID for b
    limit_km: float = 100.0,
    return_distance: bool = True,
) -> Tuple[xr.DataArray, xr.DataArray, xr.DataArray, bool]:
    """
    Match the *shorter* dataset to the *other* by nearest neighbor on a sphere.

    Returns arrays of length = len(shorter set):
      - ids_short : IDs of the shorter dataset (for reference)
      - ids_other : IDs from the other dataset (NaN if nearest > limit_km)
      - dist_km   : great-circle distance (km) for each row (if return_distance=True)
      - short_is_a: True if `a` was the shorter set; False if `b` was shorter
    """

    # Extract lon/lat as numpy arrays
    lon_a = np.asarray(a[lon_name].values)
    lat_a = np.asarray(a[lat_name].values)
    lon_b = np.asarray(b[lon_name].values)
    lat_b = np.asarray(b[lat_name].values)

    # Helper to pick ID coord or fallback to row index
    def get_ids(da: xr.DataArray, id_coord: Optional[str]) -> np.ndarray:
        if id_coord is not None and id_coord in da.coords:
            return np.asarray(da[id_coord].values)
        return np.arange(da.sizes[dim], dtype=float)  # float so we can set NaN

    ids_a = get_ids(a, id_coord_a)
    ids_b = get_ids(b, id_coord_b)

    # Decide which is shorter → this will be the query set
    if a.sizes[dim] <= b.sizes[dim]:
        lon_short, lat_short, ids_short_np = lon_a, lat_a, ids_a
        lon_other, lat_other, ids_other_all = lon_b, lat_b, ids_b
        short_is_a = True
    else:
        lon_short, lat_short, ids_short_np = lon_b, lat_b, ids_b
        lon_other, lat_other, ids_other_all = lon_a, lat_a, ids_a
        short_is_a = False

    n = len(ids_short_np)

    # KDTree on the *other* set in 3D (unit sphere); query with the *short* set
    tree = cKDTree(_to_xyz(lon_other, lat_other))
    euclid, nn_idx = tree.query(_to_xyz(lon_short, lat_short), k=1)

    # Euclidean chord length -> great-circle arc length -> distance in km
    arc_rad = 2.0 * np.arcsin(np.clip(euclid * 0.5, 0.0, 1.0))
    dist_km = EARTH_RADIUS_KM * arc_rad

    # IDs of the nearest "other" points for each short point; NaN if beyond limit
    ids_other_for_short = ids_other_all[nn_idx].astype(float)
    ids_other_for_short[dist_km > limit_km] = np.nan

    # Build xarray outputs aligned to the shorter set; give lon/lat coords explicit dims
    base_coords = {
        dim: np.arange(n),
        lon_name: (dim, lon_short),
        lat_name: (dim, lat_short),
    }

    ids_short_da = xr.DataArray(
        ids_short_np.astype(float), dims=(dim,), coords=base_coords, name="id_short"
    )
    ids_other_da = xr.DataArray(
        ids_other_for_short, dims=(dim,), coords=base_coords, name="id_other"
    )

    if return_distance:
        dist_da = xr.DataArray(
            dist_km, dims=(dim,), coords=base_coords, name="distance_km"
        )
        return ids_short_da, ids_other_da, dist_da, short_is_a

    return ids_short_da, ids_other_da, None, short_is_a





def make_grw_grid(VLM_REC_MERGED_GRW_large_GRW,dis_sel_2d,ids_2):
    """
    Build a time-by-location dataset for GRW ensemble std on the target (lat, lon, time) grid.

    Parameters
    ----------
    VLM_REC_MERGED_GRW_large_GRW : xr.Dataset
        Dataset containing 'GRW_ensemble_std' with a time dimension.
    dis_sel_2d : xr.Dataset
        Target grid providing 'lat' and 'lon' coordinates (2D).
    ids_2 : np.ndarray
        Indices mapping from flattened target grid to source vector.

    Returns
    -------
    xr.Dataset
        Dataset with variable 'GRW_ensemble_std(lat, lon, time)' on the target grid.
    """
    all_datasets = []
    
    lat_a=dis_sel_2d['lat'].values
    lon_a=dis_sel_2d['lon'].values
    
    for time_id in np.arange(len(VLM_REC_MERGED_GRW_large_GRW['GRW_ensemble_std'].time)):
        data_set = np.empty([len(dis_sel_2d['VLM_trend_coefficient_mean'].values.flatten())])*np.nan
        np.put(data_set, list(ids_2.astype(int)),VLM_REC_MERGED_GRW_large_GRW['GRW_ensemble_std'][:,time_id])
        all_datasets.append(data_set)
    all_datasets = np.asarray(all_datasets).T
    dis_sel_2d_grw= xr.Dataset(
        data_vars=dict(
            GRW_ensemble_std=(["lat", "lon","time"], all_datasets.reshape(len(lat_a),len(lon_a),len(VLM_REC_MERGED_GRW_large_GRW['GRW_ensemble_std'].time))),
        ),
        coords=dict(
            lat=("lat", lat_a),
            lon=("lon", lon_a),

        ),
    )
    return dis_sel_2d_grw


def regrid_from_flat(data_in,variables,lon_a,lat_a):
    """
    Re-grid 1D scattered values back to a structured (lat, lon) grid.

    Parameters
    ----------
    data_in : xr.Dataset
        Must contain coords 'lon' and 'lat' and variables listed in `variables`.
    variables : list[str]
        Variables in `data_in` to re-grid.
    lon_a, lat_a : array-like
        Target 1D lon/lat coordinates.

    Returns
    -------
    (xr.Dataset, np.ndarray)
        - Dataset with each variable on shape (lat, lon).
        - Index map from source points to target flattened grid.
    """
    xx, yy = np.meshgrid(lat_a,lon_a)
    xx = xx.T.flatten()
    yy = yy.T.flatten()

    data_set = np.zeros([len(yy)])*0
    ids_ = np.empty([len(data_in.lon.values)])*0
    i=0

    for lond,latd in zip(data_in.lon.values,data_in.lat.values):    
        ids_[i] = int(np.argwhere((xx==latd) & (yy==lond))   )
        i=i+1

    variable = variables[0]
    data_set = np.zeros([len(yy)])*np.nan
    np.put(data_set, list(ids_.astype(int)), data_in[variable])

    ds = xr.Dataset(
        data_vars=dict(
            data=(["lat", "lon"], data_set.reshape(len(lat_a),len(lon_a))),
        ),
        coords=dict(
            lat=("lat", lat_a),
            lon=("lon", lon_a),

        ),
    )
    ds = ds.rename({'data':variable})
    if len(variables)>1:
        for var in variables[1:]:
            data_set = np.zeros([len(yy)])*np.nan
            np.put(data_set, list(ids_.astype(int)), data_in[var].values)
            ds[var]= (('lat','lon'), data_set.reshape(len(lat_a),len(lon_a)))

    return ds,ids_


def haversine(coord1, coord2):
    """
    Great-circle distance(s) using the haversine formula.

    Parameters
    ----------
    coord1 : tuple(float, float)
        (lat1, lon1) in degrees.
    coord2 : tuple(float, float) or arrays
        (lat2, lon2) in degrees. Can be arrays to get vectorized distances.

    Returns
    -------
    float or np.ndarray
        Distance(s) in meters.
    """
    R = EARTH_RADIUS_KM*1000  # Earth radius in meters
    lat1, lon1 = coord1
    lat2, lon2 = coord2

    phi1, phi2 = np.radians(lat1), np.radians(lat2) 
    dphi       = np.radians(lat2 - lat1)
    dlambda    = np.radians(lon2 - lon1)
    
    a = np.sin(dphi/2)**2 + \
        np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2
    if str(type(a))=="<class 'numpy.float64'>":
        return 2*R*math.atan2(math.sqrt(a), math.sqrt(1 - a))
    else:
        return 2*R*np.arctan(np.sqrt(a), np.sqrt(1 - a))


def compute_dist_between_datasets(data1,data2):
    """
    For each (lon, lat) in `data1`, find nearest neighbor in `data2` and return
    index and distance.

    Returns
    -------
    (np.ndarray, np.ndarray)
        - ids_  : indices into `data2` for nearest neighbors
        - dist_ : distances (meters) to those neighbors
    """
    lon_a= data1.lon.values
    lat_a = data1.lat.values

    lon_data =data2.lon.values
    lat_data =data2.lat.values

    ids_ = np.empty([len(lat_a)])*0
    dist_ = np.empty([len(lat_a)])*0
    coord2= [lat_data,lon_data]
    i=0
    for lond,latd in zip(lon_a,lat_a):    
        coord = (latd,lond)
        dist = haversine(coord,coord2)
        arg_ = np.argmin(dist)
        ids_[i] = arg_
        dist_[i] = dist[arg_]    
        i=i+1
    return ids_,dist_

def flatten_lon_lat(datain,name='sla',time=True):
    """
    Flatten a 2D (lat, lon) field into a 1D coordinate 'x' with preserved coords.

    Parameters
    ----------
    datain : xr.DataArray
    name : str
        Output variable name.
    time : bool
        If True, keep 'time' dimension in the output.

    Returns
    -------
    xr.Dataset
    """
    lon,lat = np.meshgrid(datain.lon,datain.lat)   
    flat_x = datain.values.flatten()
    coords={'lon': (['x'], lon.flatten()),'lat': (['x'], lat.flatten())}
    if time:
        coords['time'] = datain.time.values
        var_dict = {name: (['time','x'],  flat_x)}
    else:
        var_dict = {name: (['x'],  flat_x)}        

    ds = xr.Dataset(var_dict,coords=coords)  
    return ds    

def merge_gps_and_vlm_reconstruction(VLM_REC,GPS_VLM):
    """
    Merge GPS trends with VLM reconstruction using uncertainty-aware smoothing.

    Returns
    -------
    (xr.Dataset, xr.Dataset)
        - VLM_REC_MERGED     : concatenated dataset combining GPS trends and the VLM reconstruction
        - VLM_REC_MERGED_GRW : corresponding GRW ensemble std on merged set

    Notes
    -----
    - 'kind' = 0 for VLM reconstruction, 1 for GPS-smoothed contribution.
    - GRW_ensemble_std represents present-day variability (Gaussian random walks).
    - all units are in m/year or m
    """

    # match GPS and VLM_REC locations
    ids_short, ids_other, dist_km, short_is_a = couple_short(
        GPS_VLM.x, VLM_REC.x,
        dim="x", lon_name="lon", lat_name="lat",
        id_coord_a="ID",          # or None to use row indices
        id_coord_b="ID",          # adapt to your coord names
        limit_km=100.0,
        return_distance=True
    )

    # select GPS at locations not covered by the reconstruction
    GPS_VLM_sub=GPS_VLM.isel({'x':np.isnan(ids_other)})[['trend','trend_un','trend_CMR_corrected']]
    GPS_select = [
        [-35.68, 10.44, -20.37, 42.82],   # Canary Islands region (off NW Africa in the Atlantic)
        [-179.1, 7.0, -138.6, 40.2],      # Hawaii and central Pacific Ocean
        [51.9, -56.6, 92.4, -2.2],        # Western/central Indian Ocean (Madagascar, Mauritius, Seychelles)
        [27.0, -56.6, 92.4, -35.7],       # Southwestern Indian Ocean (Mozambique Channel, around Madagascar)
        [115.49, 1.01, 131.67, 21.04],    # Maritime Southeast Asia (Borneo, Philippines, Sulawesi)
        [-177.4, -53.1, -101.6, 6.7],     # South Pacific Ocean (Fiji, Tonga, Samoa region)
    ]

    select_gps_sub = np.asarray([False]*len(GPS_VLM_sub.lon))
    for box in GPS_select:
        select_gps_sub = select_gps_sub | ((GPS_VLM_sub.lon <box[2]) & (GPS_VLM_sub.lat <box[3]) & (GPS_VLM_sub.lon >box[0]) & (GPS_VLM_sub.lat >box[1]))

    GPS_VLM_sub = GPS_VLM_sub.sel({'x':select_gps_sub})
    GPS_VLM_sub['trend'][:] = GPS_VLM_sub['trend_CMR_corrected'][:].values

    # compute spatial averages based on GPS trends
    smoothed = radial_smoothing(
        GPS_VLM_sub,
        var_name="trend",
        uncert_name="trend_un",   # skip or set use_uncert=False if not available
        dim="x", lon_name="lon", lat_name="lat",
        L_km=100.0, kernel="cos", exponent=1.0, use_uncert=True
    )

    # match resolution with the 0.25 degree resolution of the VLM reconstruction
    smoothed['lon'][:] = np.round(smoothed['lon']*4)/4
    smoothed['lat'][:] = np.round(smoothed['lat']*4)/4
    smoothed = smoothed.drop({'x'}).rename({'trend_s':'VLM_trend_coefficient_mean','trend_s_un':'VLM_trend_coefficient_uncertainty'})/1000.
    smoothed['kind'] = copy.deepcopy(smoothed['VLM_trend_coefficient_uncertainty'])
    smoothed['kind'][:] =[1]*len(smoothed['kind'])


    VLM_REC['kind'] = copy.deepcopy(VLM_REC['VLM_trend_coefficient_uncertainty'])
    VLM_REC['kind'][:] =[0]*len(VLM_REC['kind'])
    VLM_REC_MERGED = xr.concat([VLM_REC[['VLM_trend_coefficient_mean','VLM_trend_coefficient_uncertainty','kind']],smoothed],dim='x')
    smoothed_grw = smoothed[['VLM_trend_coefficient_uncertainty']].expand_dims(dim={"time": VLM_REC['GRW_ensemble_std'].time}, axis=0)['VLM_trend_coefficient_uncertainty'].rename('GRW_ensemble_std').T.to_dataset()*0
    
    # GRW_ensemble_std is the ensemble std based on the Gaussian random walks, used to mimic present day variability
    VLM_REC_MERGED_GRW = xr.concat([VLM_REC[['GRW_ensemble_std']],smoothed_grw],dim='x')

    lat_merged = VLM_REC_MERGED.lat.values
    lon_merged = VLM_REC_MERGED.lon.values
    select_merged_sub = np.asarray([False]*len(lon_merged))

    # select regions with very little GNSS station coverage where GIA trends will be applied
    unselect_data_boxes = [[-119,-39,51,91], # Covers parts of North America
                        [-180,-38,68,91],    # Covers eastern Siberia to Greenland/Arctic Ocean
                        [52,77,10,28],       # Covers Middle East / Arabian Peninsula
                        [-20,50,69.1,80]]    # Covers Northern Europe and parts of the Arctic
    for box in unselect_data_boxes:
        select_merged_sub = select_merged_sub | ((VLM_REC_MERGED.lon >box[0]) & (VLM_REC_MERGED.lat <box[3]) & (VLM_REC_MERGED.lon <box[1]) & (VLM_REC_MERGED.lat >box[2]))
    VLM_REC_MERGED = VLM_REC_MERGED.sel({'x':~select_merged_sub.values})
    VLM_REC_MERGED_GRW = VLM_REC_MERGED_GRW.sel({'x':~select_merged_sub.values})
    return VLM_REC_MERGED,VLM_REC_MERGED_GRW


def interpolate_and_compute_weights(VLM_REC_MERGED,VLM_REC_MERGED_GRW,GIA_VLM,dist2coast):
    """
    Interpolate to the target projection grid and compute weights/masks.

    Returns
    -------
    (xr.Dataset, xr.Dataset, xr.DataArray, xr.DataArray)
        - GIA_VLM_int : GIA fields interpolated to target grid (+ 'asl')
        - dis_sel_2d  : structured grid with VLM trend/uncertainty and ancillary fields
        - dis_sel_2d_grw_use : GRW ensemble std (selected years) on grid
        - w_d         : weight mask blending GPS+recon vs GIA by coastal distance
    """

    lon_grid = np.arange((360*4))/4-180
    lat_grid = np.arange((180*4))/4-90

    distance_max = 150
    lon_ext = dist2coast['lon'].values
    lon_ext[lon_ext>180]=lon_ext[lon_ext>180]-360
    dist2coast.assign_coords(lon=dist2coast.lon*0+lon_ext)

    # interpolate to the same grid used in the projections
    dist2coast_int =dist2coast['z'][::10,::11].interp(lat=lat_grid,lon=lon_grid)
    dist2coast_int_flat = flatten_lon_lat(dist2coast_int,name='z',time=False)
    dis_sel = dist2coast_int_flat['z'][abs(dist2coast_int_flat['z'])<distance_max]
    dis_sel = dis_sel.to_dataset()

    # compute distances to match VLM_REC_MERGED with dis_sel
    ids_,dist_ = compute_dist_between_datasets(dis_sel,VLM_REC_MERGED)

    VLM_REC_MERGED_large = VLM_REC_MERGED.isel({'x':ids_.astype(int)})
    VLM_REC_MERGED_GRW_large_GRW = VLM_REC_MERGED_GRW.isel({'x':ids_.astype(int)})

    variables = ['VLM_trend_coefficient_mean','VLM_trend_coefficient_uncertainty']
    for var in variables:
        dis_sel[var] = copy.deepcopy(dis_sel['z'])
        dis_sel[var][:] = VLM_REC_MERGED_large[var].values    
    var='distance'
    dis_sel[var] = copy.deepcopy(dis_sel['z'])
    dis_sel[var][:] = dist_
    variables =  ['VLM_trend_coefficient_mean','VLM_trend_coefficient_uncertainty'] +['distance','z']
    dis_sel_2d,ids_2 = regrid_from_flat(dis_sel,variables,lon_grid,lat_grid)
    dis_sel_2d_grw = make_grw_grid(VLM_REC_MERGED_GRW_large_GRW,dis_sel_2d,ids_2)

    # select from year 25 because dis_sel_2d_grw['GRW_ensemble_std'] is zero before year 25
    dis_sel_2d_grw_use = xr.concat([dis_sel_2d_grw['GRW_ensemble_std'][:,:,0:1]*0,dis_sel_2d_grw['GRW_ensemble_std'][:,:,25::5],dis_sel_2d_grw['GRW_ensemble_std'][:,:,25::5]*0],dim='time').fillna(0)
    dist_thresh = 100000 # m
    # w_d contains weights == 1, where GPS and the VLM reconstruction is used, and == 0, where GIA data is used
    w_d = (dis_sel_2d['distance']<dist_thresh)*1*((np.cos(dis_sel_2d['distance'].values*2*np.pi/(2*dist_thresh))+1)/2)#(dis_sel_2d['distance']*-1+dist_thresh)/dist_thresh
    w_d = w_d.fillna(0)


    lon_prt =GIA_VLM['lon'].values
    lon_prt[lon_prt>180] = lon_prt[lon_prt>180]-360
    GIA_VLM['lon'] = lon_prt


    GIA_VLM_int = GIA_VLM.interp(lat=lat_grid,lon=lon_grid)
    GIA_VLM_int = GIA_VLM_int.interpolate_na(dim="lat",method='linear', fill_value="extrapolate").interpolate_na(dim="lon",method='linear', fill_value="extrapolate")
    GIA_VLM_int['asl'] = GIA_VLM_int['rsl_mean']+GIA_VLM_int['rad_mean'] # The contribution of GIA to RSL without the radial deformation (that is not captured by VLM observations)

    return GIA_VLM_int,dis_sel_2d,dis_sel_2d_grw_use,w_d

def load_datasets():
    """
    Load required datasets

    Returns
    -------
    (VLM_REC, GPS_VLM, GIA_VLM, dist2coast, example) : tuple of xr.Dataset
    """

    VLM_REC      = xr.open_dataset(os.path.join(os.path.dirname(__file__), "data/VLM_reconstruction.nc"))
    GPS_VLM      = xr.open_dataset(os.path.join(os.path.dirname(__file__), "data/NGL14_CMR_corrected.nc"))
    GIA_VLM      = xr.open_dataset(os.path.join(os.path.dirname(__file__), "data/GIA_stats.nc"))
    dist2coast   = xr.open_dataset(os.path.join(os.path.dirname(__file__), "data/dist2coast_1deg_v2.grd")).rename({'x':'lon','y':'lat'})
    example      = xr.open_dataset(os.path.join(os.path.dirname(__file__), "data/ssp585.2150.fair2.emuGLA.emulandice2.glaciers_quantiles.nc"))
    return VLM_REC,GPS_VLM,GIA_VLM,dist2coast,example

def global_preprocess_verticallandmotion(pipeline_id):
    """ Load, interpolate and combine datasets:

    VLM_REC: VLM reconstruction (Oelsmann et al., 2024)

    GPS: VLM trends from NGL

    GIA_VLM: GIA VLM and effects on RSL

    dist2coast: distance to coast grid

    example: File with the same output format as IPPC AR6 SL projections

    Parameters
    ----------
    pipeline_id : str
        Unique identifier for the pipeline run.

    """

    # Load required datasets
    VLM_REC,GPS_VLM,GIA_VLM,dist2coast,example = load_datasets()
    # Merge GPS and the VLM reconstruction
    VLM_REC_MERGED,VLM_REC_MERGED_GRW = merge_gps_and_vlm_reconstruction(VLM_REC,GPS_VLM)
    # Interpolate 
    GIA_VLM_int,VLM_REC_MERGED,VLM_REC_MERGED_GRW,weights= interpolate_and_compute_weights(VLM_REC_MERGED,VLM_REC_MERGED_GRW,GIA_VLM,dist2coast)

    # Populate the output dictionary
    outdata = {'GIA_VLM_int': GIA_VLM_int, 'VLM_REC_MERGED': VLM_REC_MERGED, 'VLM_REC_MERGED_GRW': VLM_REC_MERGED_GRW, 'weights': weights, 'example':example}
    
    # Define the data directory
    outdir = os.path.dirname(__file__)

    # Write the rates data to a pickle file
    outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
    p.dump(outdata, outfile)
    outfile.close()    



if __name__ == '__main__':

    # Initialize the command-line argument parser
    parser = argparse.ArgumentParser(description="Run the pre-processing stage for the global vertical land motion workflow",\
        epilog="Note: This is meant to be run as part of the oelsmann24 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
    parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

    # Parse the arguments
    args = parser.parse_args()

    # Run the preprocessing stage
    global_preprocess_verticallandmotion(args.pipeline_id)

    # Done
    exit()