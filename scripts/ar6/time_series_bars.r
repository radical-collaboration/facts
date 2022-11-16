#############################################################################
#
# time_series_bars.r    07 April 2021
#
# Function that generates a time-series bar plot
#
#############################################################################

TSDistBars <- function(years, dist.vals, dist.medians, dist.cols, 
                       ylab="", title="", ylim=NA, ...){
  
  if(length(dim(dist.vals))==2){
    dim(dist.vals) <- c(dim(dist.vals), 1)
    dim(dist.medians) <- c(length(dist.medians), 1)
    }
  n.years <- length(years)
  n.workflows <- dim(dist.vals)[3]
  temp.vals <- c(dist.vals, dist.medians)
  
  if(any(is.na(ylim))){
    ylim=c(min(temp.vals, na.rm=T), max(temp.vals, na.rm=T))
  }
  
  scale.factor <- 0.04 * (ylim[2] - ylim[1])
  ylim <- c(ylim[1]-scale.factor, ylim[2]+scale.factor)
  
  # Calculate the x coordinates for the
  # background polygons
  poly.x <- apply(rbind(years[-1], years[-n.years]),2,mean)
  poly.x <- c(years[1] - (poly.x[1]-years[1]), 
              poly.x, 
              years[n.years] + (years[n.years]-poly.x[length(poly.x)]))
  poly.width <- diff(poly.x)

  # Start the plot
  par(xaxs="i", ...)
  plot.new()
  plot.window(xlim=range(poly.x), ylim=ylim)
  axis(1, at=years)
  axis(2)
  mtext("Year", 1, line=2.7, cex=1.0, font=2)
  mtext(ylab, 2, line=2.5, cex=1.0, font=2)
  mtext(title, 3, line=1.0, cex=1.2, font=2)
  poly.cols <- c("gray90", "white")
  for(k in 1:n.years){
    polygon(x=c(poly.x[k],poly.x[k],poly.x[k+1],poly.x[k+1]),
            y=c(ylim[1]-2*scale.factor,
                ylim[2]+2*scale.factor,
                ylim[2]+2*scale.factor,
                ylim[1]-2*scale.factor), 
            col=poly.cols[(k%%2)+1], border=NA)
  }
  grid(nx=NA, ny=NULL, col="gray75")
  box()
  
  # Initialize return variable for the x-coordinates of the workflows
  xcoords <- array(NA, dim=c(n.years, n.workflows))
  
  # Loop over the workflows
  for(j in 1:n.workflows){
    
    # Determine the x coordinate for this distribution
    this.x <- poly.x[-length(poly.x)] + j*(poly.width/(n.workflows+1))
    xcoords[,j] <- this.x
  
    # Plot the line segments
    segments(x0=this.x, y0=dist.vals[2,,j], y1=dist.vals[3,,j],
             col=dist.cols[j], lwd=3.2, lend=2)
    segments(x0=this.x, y0=dist.vals[1,,j], y1=dist.vals[4,,j],
             col=dist.cols[j], lwd=1.2, lend=2)
  
    # Plot the median as a point
    points(this.x, dist.medians[,j], pch=21, bg=dist.cols[j])
  }
  
  # Return the xcoordinates
  return(list("xcoords"=xcoords))
  
}