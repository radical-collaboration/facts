# Quickly load and test localized projections
library(ncdf4)
library(fields)
library(akima)
library(maps)

# Load the netcdf file
dir <- "/Users/ggg46/research/slr_framework/code/ipcc_ar6_regional_projections/final_regional_projections/pboxes/pb_1e"
scenario <- "tlim2.0win0.25"
filename <- "verticallandmotion-kopp14-verticallandmotion_localsl.nc"
ncfile <- sprintf("%s/%s/%s", dir, scenario, filename)
nc <- nc_open(ncfile)

# Extract variables
lats <- ncvar_get(nc, "lat")
lons <- ncvar_get(nc, "lon")
vals <- ncvar_get(nc, "localSL_quantiles")
qs <- ncvar_get(nc, "quantiles")
nc_close(nc)

# Which quantiles to test?
#q.idx <- sapply(c(0.05, 0.17, 0.5, 0.83, 0.95), function(x) which.min(abs(qs - x)))
q.idx <- sapply(c(0.5), function(x) which.min(abs(qs - x)))

# Layout for maps
#layout.mat <- matrix(c(1,1,2,2,0,3,3,0,4,4,5,5), ncol=4, byrow=TRUE)
#layout(layout.mat)
par(mar=c(2,2,0,0)+1)

# Loop over the quantiles of interest
for(this.q in q.idx){
  
  # Interpolate for map values
  vals.map <- interp(lons, lats, vals[9,,this.q], duplicate="mean", nx=360, ny=180)
  vals.map$z[vals.map$z > 700] <- NA
  vals.map$z[vals.map$z < -300] <- NA
  
  # Generate the plot
  image(vals.map, zlim=c(-300,700), col=hcl.colors(12, "Blue-Red 3"), useRaster=TRUE)
  map(add=T, wrap=c(-180,180,-90), fill=TRUE, col="gray75")
}