# Loads and constructs a Bathymetry grid for the area of interest.  
# Returns a grid of depths for the area of interest.
# Dependency: Have the ncdf module installed.


bathy <- function(inputFile, inputFileType, startX=0, startY=0, XDist, YDist, seriesName, debug=FALSE) {
	switch(inputFileType, 
		"netcdf" = {
		    if(file.exists(as.character(inputFile)) && require(ncdf)){
				library(ncdf)
		
	            ## open the netCDF file
	            ncdfObj = open.ncdf(inputFile)
	            
	            ## grab a slice (in grid form)
	            bGrid = get.var.ncdf(ncdfObj, 'z', start=c(startX, startY), 
	                    count=c( XDist, YDist))
		    }
			else {
				bGrid = simulate(XDist, YDist)
			}
		},
	    "arcgis" = {
			if(file.exists(as.character(inputFile)) && require(sp) && require(rgdal) 
					&& require(raster)){
				library(sp)
				library(raster)
				library(rgdal)
		        #For an arc/grid inputFile is the folder!:
		        bGrid = raster(inputFile)
	    	}
			else {
				bGrid = simulate(XDist, YDist)
			}
		})
    return(bGrid)
}

simulate <- function(XDist, YDist) {
	## Create test bGrid to use if real data unavailable
	nx = XDist
	ny = YDist
	x <- seq(-2*pi,2*pi,length=nx)
	X <- matrix(rep(x,ny),ny,nx,byrow=TRUE)
	y <- seq(-2*pi,2*pi,length=ny)
	Y <- matrix(rep(y,nx),ny,nx,byrow=FALSE)
	bGrid <- sin(X)*sin(Y)*abs(X)*abs(Y)-pi
	bGrid[bGrid>0] <- 0
	return(bGrid)
}