##############################################################################################
#
# calculate_positive_rsl_percentage.r    09 July 2021
#
# AR5 states:
# "By the end of the 21st century, it is very likely that over about 95% of the world ocean, 
# regional sea level rise will be positive, and most regions that will experience a sea level 
# fall are located near current and former glaciers and ice sheets. About 70% of the global 
# coastlines are projected to experience a relative sea level change within 20% of the global 
# mean sea level change."
#
# We're going to make the same calculations with the new SLR projections.
#
##############################################################################################

# Area weight function
areaweight <- function(lons, lats, r=6367470){
  lat.bnds <- cbind(lats-0.5, lats+0.5)
  lon.bnds <- cbind(lons-0.5, lons+0.5)
  
  lat.bnds <- lat.bnds * pi/180
  lon.bnds <- lon.bnds * pi/180
  
  sapply(1:length(lats), function(x) r^2*(lon.bnds[,2]-lon.bnds[,1])*(sin(lat.bnds[x,2]) - sin(lat.bnds[x,1])))
  
}

# Directories
maindir <- "/Users/ggg46/research/slr_framework/code/ipcc_ar6_full_sample_regional_projections"
scriptdir <- sprintf("%s/scripts", maindir)
pbdir <- sprintf("%s/regional/pboxes/pb_1e", maindir)

# Load libraries and sourced data
library(ncdf4)
library(fields)
library(akima)
library(maps)
library(GSIF)

# Scenario for analysis
scenarios <- c("ssp119", "ssp126", "ssp245", "ssp370", "ssp585",
               "tlim1.5win0.25", "tlim2.0win0.25", "tlim3.0win0.25", "tlim4.0win0.25", "tlim5.0win0.25")

# Year of interest
targyear <- 2100

# Quantile of interest
this.q <- 0.1

# Landmask threshold Percentage (0-100).
# Values less than the first value is considered ocean.
# Values greater than the second value is considered land.
# Values in between are considered coastline.
mask.threshold <- c(10,90)

# Load the landmask dataset 
data(landmask)

# Make the landmask grid
lon.interp <- sort(unique(landmask$x))
lat.interp <- sort(unique(landmask$y))
mask.grid <- interp(landmask$x, landmask$y, landmask$mask, xo=lon.interp, yo=lat.interp)


# Loop over scenarios
for(scenario in scenarios){

  # Open the total slr projection file and load the meta fields
  ncfile <- sprintf("%s/%s/total-workflow.nc", pbdir, scenario)
  nc <- nc_open(ncfile)
  qs <- ncvar_get(nc, "quantiles")
  years <- ncvar_get(nc, "years")
  ids <- ncvar_get(nc, "locations")
  lats <- ncvar_get(nc, "lat")
  lons <- ncvar_get(nc, "lon")

  # Load the slr projection data that's consistent with the year and quantile of interest
  year.idx <- which.min(abs(years - targyear))
  q.idx <- which.min(abs(qs - this.q))
  slr <- ncvar_get(nc, "sea_level_change", start=c(1, year.idx, q.idx), count=c(-1,1,1))
  slr[is.na(slr)] <- -99999

  # Close the netcdf file
  nc_close(nc)

  # Separate the tide gauge indices from the global grid
  tg.idx <- which(ids < 100000)
  gg.idx <- which(ids > 100000)

  # Use simple interpolation to put the slr projections and land mask on the same grid
  slr.grid <- interp(lons[gg.idx], lats[gg.idx], slr[gg.idx], xo=lon.interp, yo=lat.interp)

  # Which grid indices are considered ocean and which are coastline
  ocean.idx <- which(mask.grid$z < mask.threshold[1])
  land.idx <- which(mask.grid$z > mask.threshold[2])
  coast.idx <- which(mask.grid$z >= mask.threshold[1] & mask.grid$z <= mask.threshold[2])

  ########### Calculations ################

  # What percent of the ocean area experiences greater than 0 regional SLR?
  #pct_positive_ocean_rsl <- (sum(slr.grid$z[ocean.idx] > 0.0, na.rm=T)/length(ocean.idx))*100
  pos.ocean <- slr.grid$z[ocean.idx] > 0.0
  pos.ocean[is.na(pos.ocean)] <- 0
  temp.weights <- areaweight(slr.grid$x, slr.grid$y)
  ocean.weights <- temp.weights[ocean.idx]/sum(temp.weights[ocean.idx])
  pct_positive_ocean_rsl <- sum(pos.ocean * ocean.weights) *  100

  print(sprintf("%s: Ocean area with RSL > 0: %.2f%%", scenario, pct_positive_ocean_rsl))
}

