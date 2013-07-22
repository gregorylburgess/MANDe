# Loads and constructs a Bathymetry grid for the area of interest.  
# Returns a grid of depths for the area of interest.
# Dependency: Have the ncdf module installed.

library(sp)
library(raster)
library(rgdal)

bathy <- function(inputFile, inputFileType, startX=0, startY=0, XDist, YDist, seriesName, debug=FALSE) {

    if(file.exists(inputFile) && require(ncdf)){
		library(ncdf)
        switch(inputFileType, 
            "netcdf" = {
                ## open the netCDF file
                ncdfObj = open.ncdf(inputFile)
                
                ## grab a slice (in grid form)
                bGrid = get.var.ncdf(ncdfObj, 'z', start=c(startX, startY), 
                        count=c( XDist, YDist))
            },
            "arcgis" = {
                #For an arc/grid (albem_s1 is the folder!):
                bGrid = raster(inputFile)
            }
        )
    }
    else {
        ## Create test bGrid to use if real data unavailable
		nx = XDist
		ny = YDist
        x <- seq(-2*pi,2*pi,length=nx)
        X <- matrix(rep(x,ny),ny,nx,byrow=TRUE)
        y <- seq(-2*pi,2*pi,length=ny)
        Y <- matrix(rep(y,nx),ny,nx,byrow=FALSE)
        bGrid <- sin(X)*sin(Y)*abs(X)*abs(Y)-pi
        bGrid[bGrid>0] <- 0
        ## Plot bGrid
        ##image(x,y,bGrid)
        ##contour(x,y,bGrid,xlab='x',ylab='y',add=TRUE,nlevels=5)
    }

    ## return grid
    return(bGrid)
}
