#############################################################
#
# venice_projections_timeseries.r   09 April 2021
#
# Produce time series bar plots of the NYNJ projections.
#
#############################################################

# Directories
maindir <- "/Users/ggg46/research/slr_framework/code/ipcc_ar6_nynj_projections"
datadir <- sprintf("%s/figure_data", maindir)
scriptdir <- sprintf("%s/scripts", maindir)
plotdir <- sprintf("%s/plots", maindir)

# Libraries and sourced code
library(ncdf4)
source(sprintf("%s/time_series_bars.r", scriptdir))
source(sprintf("%s/super_legend.r", scriptdir))

# Define SSP colors
ssp.reds <- c(30, 29, 234, 242, 132)
ssp.greens <- c(50, 51, 221, 17, 11)
ssp.blues <- c(132, 84, 61, 17, 34)
scenario.cols <- rgb(ssp.reds, ssp.greens, ssp.blues, maxColorValue = 255)
names(scenario.cols) <- c("ssp119", "ssp126", "ssp245", "ssp370", "ssp585")
scenario.print.names <- c("SSP1-1.9", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
names(scenario.print.names) <- c("ssp119", "ssp126", "ssp245", "ssp370", "ssp585")

# Available locations to plot
site.names <- c("NewYork", "AtlanticCity")
n.sites <- length(site.names)

# Which pbox and site to plot
my.pboxes <- c("pb_1e", "pb_2e", "pb_1f", "pb_2f")
n.pboxes <- length(my.pboxes)

for(this.pbox.idx in 1:n.pboxes){
for(this.site.idx in 1:n.sites){

this.pbox <- my.pboxes[this.pbox.idx]
pbox.print.names <- list("pb_1f" = "PB 1f: Med-Confidence Processes",
                         "pb_2f" = "PB 2f: Low-Confidence Processes",
                         "pb_1e" = "PB 1e: Med-Confidence Processes",
                         "pb_2e" = "PB 2e: Low-Confidence Processes")

# Available scenarios
pb.scenarios = list("pb_1f"=c("ssp119", "ssp126", "ssp245", "ssp370", "ssp585"),
                    "pb_2f"=c("ssp126", "ssp245", "ssp585"),
                    "pb_1e"=c("ssp119", "ssp126", "ssp245", "ssp370", "ssp585"),
                    "pb_2e"=c("ssp126", "ssp245", "ssp585"))
scenarios <- unlist(pb.scenarios[this.pbox])
n.scenarios <- length(scenarios)

# Plot years for pbox
long.plot.years <- seq(2020,2300,by=40)
short.plot.years <- seq(2020,2100,by=20)
pbox.plot.years <- list("pb_1e"=short.plot.years, "pb_2e"=short.plot.years,
                       "pb_1f"=long.plot.years, "pb_2f"=long.plot.years)

# Y-limits
pbox.ylims <- list("pb_1e"=c(0,2.0), "pb_2e"=c(0,2.0),
                   "pb_1f"=c(0,5), "pb_2f"=c(0,5))

# Options ============================
plot2dev <- TRUE
plot.quantiles <- c(0.05,0.17,0.83,0.95)
n.quantiles <- length(plot.quantiles)

plot.years <- as.numeric(unlist(pbox.plot.years[this.pbox]))
n.years <- length(plot.years)

# Output file name
out.file <- sprintf("%s/%s_%s_timeseries.pdf", plotdir, site.names[this.site.idx], this.pbox)

# ====================================

# Workflows to plot
n.workflows <- n.scenarios
dist.cols <- scenario.cols[scenarios]
workflow.names <- scenarios
workflow.names.print <- scenario.print.names[scenarios]

# Initialize the data structures
dist.vals <- array(NA, dim=c(n.quantiles, n.years, n.workflows))
dist.medians <- array(NA, dim=c(n.years, n.workflows))
dist.nmodels <- array(NA, dim=c(n.years, n.workflows))

# Loop over the workflows
for(i in 1:n.workflows){
  
  # Load the netcdf file
  nc <- nc_open(sprintf("%s/%s/%s/total-workflow_figuredata.nc", datadir, this.pbox, workflow.names[i]))
    
  # Extract the data
  nc.data <- ncvar_get(nc, "localSL_quantiles", collapse_degen = FALSE)/1000
  nc.years <- ncvar_get(nc, "years")
  nc.quantiles <- ncvar_get(nc, "quantiles")
  nc.lat <- ncvar_get(nc, "lat", start=this.site.idx, count=1)
  nc.lon <- ncvar_get(nc, "lon", start=this.site.idx, count=1)
  nc_close(nc)
    
  # Find the overlap in years and quantiles
  nc.years.idx <- which(nc.years %in% plot.years)
  nc.quantiles.idx <- sapply(plot.quantiles, function(x) which.min(abs(nc.quantiles - x)))
  nc.median.idx <- which.min(abs(nc.quantiles-0.5))
    
  # Load the data into the appropriate data structure
  dist.vals[,,i] <- t(nc.data[nc.years.idx,this.site.idx,nc.quantiles.idx])
  dist.medians[,i] <- t(nc.data[nc.years.idx,this.site.idx,nc.median.idx])
  
}

# Generate the plot
if(plot2dev) { pdf(out.file, height=2.7, width=4.2, pointsize=7)}
par(mar=c(4.5,4.8,3.5,1)+0.1)
TSDistBars(plot.years, dist.vals, dist.medians, dist.cols, ylim=unlist(pbox.ylims[this.pbox]),
          ylab="Sea-level change since 2005 [m]",
          title=sprintf("%s [%0.3f degN, %0.3f degE]\n%s", site.names[this.site.idx], nc.lat, nc.lon,pbox.print.names[this.pbox]))

# Apply the legend
LEGEND("topleft", workflow.names.print, lty=1, pch=21, cex=1.0,
       line.col=dist.cols, pt.bg=dist.cols, pt.col="black", bg="white", y.intersp = 1.1)

if(plot2dev){dev.off()}
}
}
