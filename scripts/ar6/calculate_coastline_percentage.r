##############################################################################################
#
# calculate_coastline_percentage.r    08 February 2021
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
plotdir <- sprintf("%s/plots", maindir)
gmsldir <- "/Users/ggg46/research/slr_framework/code/ipcc_ar6_global_projections/final_gmsl_projections/pboxes/pb_1e"

# Load libraries and sourced data
library(ncdf4)
library(fields)
library(akima)
library(maps)
library(GSIF)

# Plot to a pdf
plot2pdf <- TRUE

# Scenarios
#scenarios <- c("ssp119", "ssp126", "ssp245", "ssp370", "ssp585")
scenarios <- c("tlim1.5win0.25", "tlim2.0win0.25", "tlim3.0win0.25", "tlim4.0win0.25", "tlim5.0win0.25")
scenario.cols <- c("blue", "green", "yellow", "orange", "red")

# Year of interest
targyear <- 2100

# Quantiles of interest
target.qs <- c(0.05, 0.17, 0.5, 0.83, 0.95)
#target.qs <- 0.5

# Landmask threshold Percentage (0-100).
# Values less than the first value is considered ocean.
# Values greater than the second value is considered land.
# Values in between are considered coastline.
mask.threshold <- c(10,90)

# Load the landmask dataset 
data(landmask)
lon.interp <- sort(unique(landmask$x))
lat.interp <- sort(unique(landmask$y))
if(!exists("mask.grid")) {
  mask.grid <- interp(landmask$x, landmask$y, landmask$mask, xo=lon.interp, yo=lat.interp)
}

# Loop over the quantiles of interest
for(this.q in target.qs){
  
  print(sprintf("q = %0.2f", this.q))

  # Open the plotting device if necessary
  plotfile <- sprintf("%s/coastline_fraction_%.0f_%.0f.pdf", plotdir, targyear, this.q*100)
  if(plot2pdf) pdf(plotfile, width=3.5, height=2.7, pointsize=7)

  # Open the table file
  table.file <- sprintf("%s/coastline_fraction_%.0f_%.0f.txt", plotdir, targyear, this.q*100)
  if(plot2pdf) table.filecon <- file(table.file, open = "wt")

  # Setup the plot
  par(mar=c(4,4,3,1)+0.1)
  plot.new()
  plot.window(xlim=c(-2,2), ylim=c(0,100))
  axis(1)
  axis(2)
  mtext("Multiples of GMSL [-]", side=1, line=2.5, cex=1.1, font=2)
  mtext("Coastal area [%]", side=2, line=2.5, cex=1.1, font=2)
  mtext(sprintf("Year: %.0f, Percentile: %.0f", targyear, this.q*100), side=3, line=1, cex=1.25, font=2)
  polygon(c(0.8, 0.8, 1.2, 1.2), c(0, 100, 100, 0), border=NA, col="gray85")
  abline(v=0, lty=2, col="gray")
  abline(h=c(0,100), lty=2, col="gray")
  box()

  # Start printing the table to screen
  gmsl.to.print <- seq(-2, 2, by=0.2)
  if(plot2pdf) writeLines(sprintf("Scenario, %s", paste(sprintf("%.1f", gmsl.to.print), collapse=", ")), table.filecon)

  # Initialize vector to hold coastline percentages
  all.coast.area <- array(NA, dim=length(scenarios))

  # Loop over the scenarios
  for(i in 1:(length(scenarios))){
  
    # This scenario
    scenario <- scenarios[i]

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
    
    # Load the GMSL pbox file
    gmsl.file <- sprintf("%s/%s/total-workflow.nc", gmsldir, scenario)
    nc <- nc_open(gmsl.file)
    gmsl <- ncvar_get(nc, "sea_level_change", start=c(1, year.idx, q.idx), count=c(-1,1,1))
    gmsl <- as.vector(gmsl)
    nc_close(nc)

    # Load the GMSL workflow files
    #wf1e.filename <- sprintf("%s/wf_1e/%s/total-workflow_globalsl.nc", gmsldir, scenario)
    #wf2e.filename <- sprintf("%s/wf_2e/%s/total-workflow_globalsl.nc", gmsldir, scenario)
    #nc <- nc_open(wf1e.filename)
    #wf1e.years <- ncvar_get(nc, "years")
    #wf1e.samps <- ncvar_get(nc, "samps")
    #nc_close(nc)
    #nc <- nc_open(wf2e.filename)
    #wf2e.years <- ncvar_get(nc, "years")
    #wf2e.samps <- ncvar_get(nc, "samps")
    #nc_close(nc)

    # Calculate the pbox quantiles for the gmsl data
    #wf1e.qs <- quantile(wf1e.samps[,which(wf1e.years==targyear)], probs=this.q, na.rm=T)
    #wf2e.qs <- quantile(wf2e.samps[,which(wf2e.years==targyear)], probs=this.q, na.rm=T)
    #if(this.q > 0.5){
    #  gmsl <- max(c(wf1e.qs, wf2e.qs))
    #}else if(this.q < 0.5){
    #  gmsl <- min(c(wf1e.qs, wf2e.qs))
    #}else{
    #  gmsl <- mean(c(wf1e.qs, wf2e.qs))
    #}
    
    

    # Separate the tide gauge indices from the global grid
    #tg.idx <- which(ids < 100000)
    gg.idx <- which(ids > 100000)

    # Use simple interpolation to put the slr projections and land mask on the same grid
    slr.grid <- interp(lons[gg.idx], lats[gg.idx], slr[gg.idx], xo=lon.interp, yo=lat.interp)

    # Which grid indices are considered ocean and which are coastline
    #ocean.idx <- which(mask.grid$z < mask.threshold[1])
    #land.idx <- which(mask.grid$z > mask.threshold[2])
    coast.idx <- which(mask.grid$z >= mask.threshold[1] & mask.grid$z <= mask.threshold[2])
  
    # Apply filters to this data
    slr.grid$z[-coast.idx] <- NA  # Set all non-coast grid cells to NA
    slr.grid$z[,slr.grid$y < -60] <- NA  # Remove Antarctica (south polar region)
    #slr.grid$z[,slr.grid$y > 60] <- NA  # Remove Greenland-ish (north polar region)
    slr.grid$z[slr.grid$z > 1e10] <- NA  # Set any netcdf fill values to NA
    
    # Calculate the grid weights
    grid.weights <- areaweight(slr.grid$x, slr.grid$y)

    # Plot the coastal area percentage
    #my.fun <- ecdf(slr.grid$z/gmsl)
    #lines(ecdf(slr.grid$z/gmsl), col=scenario.cols[i])
    total.coastal.area <- sum(grid.weights[coast.idx])
    coastal.order <- order(slr.grid$z/gmsl)
    plot.percentages <- (grid.weights[coastal.order]/total.coastal.area)
    lines((slr.grid$z/gmsl)[coastal.order], cumsum(plot.percentages)*100, col=scenario.cols[i])
    my.fun <- approxfun((slr.grid$z/gmsl)[coastal.order], cumsum(plot.percentages), method="constant")
    
  
    # Print values to screen
    out.string <- scenario
    for(this.factor in gmsl.to.print){
      out.string <- paste(out.string, sprintf("%.2f", my.fun(this.factor)*100), sep=", ")
    }
    if(plot2pdf) writeLines(out.string, table.filecon)
  
    # Calculate coast area within 20% of gmsl
    #temp.coast <- slr.grid$z[coast.idx]
    #pass.idx <- as.numeric(temp.coast > 0.8*gmsl & temp.coast < 1.2*gmsl)
    #pass.idx[is.na(pass.idx)] <- 0
    #coast.weights <- grid.weights[coast.idx]/sum(grid.weights[coast.idx])
    #coast.area <- sum(pass.idx * coast.weights) * 100
    #print(sprintf("%s: Coastal area within 20%% of GMSL: %.2f%%", scenario, coast.area))
    #all.coast.area[i] <- coast.area
    
    coast.area <- (my.fun(1.2) - my.fun(0.8))*100
    print(sprintf("%s: Coastal area within 20%% of GMSL: %.2f%%", scenario, coast.area))
    all.coast.area[i] <- coast.area
  }
  
  print("===========================================")
  
  # Legend
  legend("topleft", legend=sprintf("%s = %.1f%%", toupper(scenarios), all.coast.area), lty=1, col=scenario.cols, bg="white")

  # Turn of plotting device if necessary
  if(plot2pdf) dev.off()
  if(plot2pdf) close(table.filecon)
}