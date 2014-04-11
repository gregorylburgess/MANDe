#' Contains functions that read Bathymetry files and generate the Topography Grid


#' @name getBathy
#' @title Generate a topographyGrid.
#' @description Generates a Bathymetric Grid (topographyGrid) for the program to use.  Able to ingest NetCDF, 
#' ArcGIS, and ASC file formats.
#'
#' @param inputFile The relative path to the file to open.
#' @param inputFileType The type of file to open.  Vailid options are "netcdf", "arcgis", and "asc". The type "RData" can also be used and must refer to a file saved with the "save" command and containing a single variable, which must be a matrix containing bathymetry values.
#' @param startX Starting index of the topographyGrid to take from the bathy file.
#' @param startY Starting index of the topographyGrid to take from the bathy file.
#' @param XDist The width of your desired topographyGrid.
#' @param YDist the height of your desired topographyGrid.
#' @param seriesName If you specified netcdf or arcgis, this is the name of the data series to use.
#' @param debug If enabled, turns on debug printing (console only).
#' @param timestamp A timestamp of when the job was started (for logging purposes).
#' @return A topographyGrid based on the parameters given.  If an error occurs, a default grid is provided.
getBathy <- function(inputFile, inputFileType, startX=0, startY=0, XDist, YDist, seriesName, timestamp, debug=FALSE) {
	if (file.exists(as.character(inputFile))) {
			if(startX < 1 || startY < 1) {
				printError("topographyGrid x and y coordinates must be integers greater than 0.", timestamp)
			}
            if(inputFileType == "netcdf" && require(ncdf)){
                library(ncdf)
                ## open the netCDF file
                ncdfObj = open.ncdf(inputFile)
				#print(summary(ncdfObj))
                ## grab a slice (in grid form)
                topographyGrid = get.var.ncdf(ncdfObj, 'z', start=c(startX, startY), count=c( XDist, YDist))
	    }
   	    else if(inputFileType == "arcgis" && require(sp) && require(rgdal) && require(raster)){
                library(sp)
                library(raster)
                library(rgdal)
                ## For an arc/grid inputFile is the folder!:
                topographyGrid = raster(inputFile)
            }
	    else if(inputFileType == "asc") {
                bath = loadASC(inputFile)
                topographyGrid = bath[startY:(startY-1+YDist),startX:(startX-1+XDist)]
            }
            else if(inputFileType == "RData") {
                bathname = load(inputFile)
                topographyGrid = get(bathname)
                topographyGrid = topographyGrid[startY:(startY-1+YDist),startX:(startX-1+XDist)]
            } else {
                topographyGrid = simulatetopographyGrid(XDist,YDist)
            }
	} else {
            print("Bathymetry file not found.")
            topographyGrid = simulatetopographyGrid(XDist,YDist)
	}

        ## Check if all values in topographyGrid are NA
        if(all(is.na(topographyGrid))){
            printError("all values in topographyGrid are NA!", stop=TRUE)
        } else {
            ## Check if all values in topographyGrid are positive and if so make them negative
            if(all(topographyGrid >= 0)){
                print("Warning: No negative values found in topography grid and therefore nowhere to place sensors. Multiplying all values by -1. This may be inappropriate!!")
                topographyGrid <- -topographyGrid
            }
            ## Quick fix to get rid of NA in topographyGrid, should probably be interpolated (or something)
            if(any(is.na(topographyGrid))){
                print("Warning: NAs found in topographyGrid! setting to zero. This may be inappropriate so you may want to manually remove them.")
                topographyGrid[is.na(topographyGrid)] <- 0
        }
    }
	#print(topographyGrid)
    return(topographyGrid)
}

#' @name simulatetopographyGrid
#' @title Creates a default topographyGrid to use if real data could not be loaded.
#'
#' @param XDist The width of your desired topographyGrid.
#' @param YDist the height of your desired topographyGrid.
#' @return A default topographyGrid of dimensions XDist by YDist.
simulatetopographyGrid <- function(XDist, YDist) {
	nx = XDist
	ny = YDist
	x <- seq(-2*pi,2*pi,length=nx)
	X <- matrix(rep(x,ny),ny,nx,byrow=TRUE)
	y <- seq(-2*pi,2*pi,length=ny)
	Y <- matrix(rep(y,nx),ny,nx,byrow=FALSE)
	topographyGrid <- sin(X)*sin(Y)*abs(X)*abs(Y)-pi
	topographyGrid[topographyGrid>0] <- 0
	return(topographyGrid)
}

#' @name loadASC
#' @title Loads a topography grid from an ASC file.
#' @description The first five lines of the file must have this format:\cr
#' ncols         4026\cr
#' nrows         1039\cr
#' xllcorner     812625\cr
#' yllcorner     648140\cr
#' cellsize      5\cr
#' NODATA_value  -9999\cr
#'
#' Then come the data separated by " " (blank spaces).
#' @param inputFile Path to the asc file.
#' @return The loaded topograhy.
loadASC <- function(inputFile) {
    ncols <- as.numeric(scan(file=inputFile,what='character',nmax=2,sep='')[2])
    nrows <- as.numeric(scan(file=inputFile,what='character',skip=1,nmax=2,sep='')[2])
    xll <- as.numeric(scan(file=inputFile,what='character',skip=2,nmax=2,sep='')[2])
    yll <- as.numeric(scan(file=inputFile,what='character',skip=3,nmax=2,sep='')[2])
    dx <- as.numeric(scan(file=inputFile,what='character',skip=4,nmax=2,sep='')[2])
    nodata <- as.numeric(scan(file=inputFile,what='character',skip=5,nmax=2,sep='')[2])

    dat <- scan(file=inputFile,skip=6)

    dat[dat==nodata] <- NA
    bath <- matrix(dat,nrows,ncols,byrow=TRUE)
    bath <- bath[nrows:1,]
    bath <- t(bath)
    return(bath)
}

