
#' @name getBathy
#' @title Generates a Bathymetric Grid (BGrid) for the program to use.  Able to ingest NetCDF, 
#' ArcGIS, and ASC file formats.
#'
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param inputFile The relative path to the file to open.
#' @param inputFileType The type of file to open.  Vailid options are "netcdf", "arcgis", and "asc".
#' @param startX Starting index of the BGrid to take from the bathy file.
#' @param startY Starting index of the BGrid to take from the bathy file.
#' @param XDist The width of your desired BGrid.
#' @param YDist the height of your desired BGrid.
#' @param seriesName If you specified netcdf or arcgis, this is the name of the data series to use.
#' @param debug If enabled, turns on debug printing (console only).
#' @return A BGrid based on the parameters given.  If an error occurs, a default grid is provided.
getBathy <- function(inputFile, inputFileType, startX=0, startY=0, XDist, YDist, seriesName, debug=FALSE) {
	if (file.exists(as.character(inputFile))) {
		
		if(inputFileType == "netcdf" && require(ncdf)){
					library(ncdf)
		            ## open the netCDF file
		            ncdfObj = open.ncdf(inputFile)
		            
		            ## grab a slice (in grid form)
		            bGrid = get.var.ncdf(ncdfObj, 'z', start=c(startX, startY), 
		                    count=c( XDist, YDist))
	    }
		else if(inputFileType == "arcgis" && require(sp) && require(rgdal) && require(raster)){
					library(sp)
					library(raster)
					library(rgdal)
			        #For an arc/grid inputFile is the folder!:
			        bGrid = raster(inputFile)
    	}
		else if(inputFileType == "asc") {
					load(inputFile)
					bGrid = bath[startY:(startY+YDist),startX:(startX+XDist)]
		}
		else {
			bGrid = simulateBGrid(XDist,YDist)
		}
	}
	else {
		print("File not found.")
		bGrid = simulateBGrid(XDist,YDist)
	}

        ## Quick fix to get rid of NA in bGrid, should probably be interpolated (or something)
        if(any(is.na(bGrid))){
            bGrid[is.na(bGrid)] <- 0
        }

	#print(bGrid)
    return(bGrid)
}

#' @name simulateBGrid
#' @title Creates a default BGrid to use if real data could not be loaded.
#'
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param XDist The width of your desired BGrid.
#' @param YDist the height of your desired BGrid.
#' @return A default BGrid of dimensions XDist by YDist.
simulateBGrid <- function(XDist, YDist) {
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
