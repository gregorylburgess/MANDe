#rm(list=ls()) 
source('src/Bathy.R')
#' @include src/Bathy.R
source('src/FishModel.R')
#' @include src/FishModel.R
source('src/Utility.R')
#' @include src/Utility.R

#' @name acousticRun
#' @title Design network with the provided parameters.
#' @description This is the main routine that calls several sub-routines when designing the network.
#' First the input parameters are checked for validity and default values are assigned if input 
#' values are missing. Then bathymetry is loaded and the fish distribution generated. Then the
#' goodness grid (sumGrid) is calculated depending on the design parameters (this is the heavy part).
#' When the sumGrid is finished sensors can be placed optimally and stats and figures generated.
#'
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param showPlots If TRUE plots are shown on the screen, if FALSE plots are stored in the img folder.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt If TRUE use vectorized R commands (faster).
#' @param save.inter If TRUE intermediary calculations are output as key inter.
#' @return A dictionary of return objects, see RETURN_DESCRIPTIONS.html for more info.
#' @export
acousticRun <- function(params, showPlots=FALSE, debug=FALSE, opt=FALSE, save.inter=FALSE){
    startTime = Sys.time()
    if(debug) {
        cat("\n[acousticRun]\n")
    }
    params = checkParams(params)

    ## Create/Load the Bathy grid for the area of interest
    bGrid = getBathy(params$inputFile, params$inputFileType, params$startX, params$startY, 
            params$XDist, params$YDist, params$seriesName, debug)
    bGrid = list("bGrid"=bGrid, "cellRatio"=params$cellSize)
    ## Convert parameter values from meters to number of grid cell 
    params = convertMetersToGrid(params,bGrid)
    ## Specify a standard scale of x and y axes if previously undefined
    if(!("x" %in% names(bGrid))) {
		bGrid$x = (1:dim(bGrid$bGrid)[1])*params$cellSize 
	}
    if(!("y" %in% names(bGrid))) {
		bGrid$y = (1:dim(bGrid$bGrid)[2])*params$cellSize
	}
    ## Calculate fish grid
    fGrid = fish(params, bGrid)

    ## Find good sensor placements
    sensors <- sensorFun(params$numSensors, bGrid, fGrid, params$range, params$bias, params, debug, opt, save.inter=save.inter)

    ## Stat analysis of proposed setup.
    statDict = getStats(params, bGrid, fGrid, sensors, debug, opt)
    ## Return Fish grid, Bathy grid, and Sensor Placements as a Dictionary.
    results = list("bGrid" = bGrid, "fGrid" = fGrid, "sumGrid"=sensors$sumGrid, "sensors" = sensors$sensorList, 
            "stats" = statDict, "params"=params)
    if(save.inter) results$inter = sensors$inter
    
    ## Graph results and make data file.
    results$filenames = graph(results,params,showPlots)
	
    endTime = Sys.time()
    results$runTime = endTime - startTime
    ## Email results
    if("userEmail" %in% names(params)) {
        from = "acousticwebapp@gmail.com"
        to = params$userEmail
        subject <- "Acoustic webapp results"
        body <- results                    
        mailControl=list(smtpServer="smtp.gmail.com")
	##sendmail(from=from,to=to,subject=subject,msg=body,control=mailControl)
    }
    return(results)
}

#' @name acousticTest
#' @title Executes a test run of the program, using default parameters.  No additional 
#' parameters are necessary.
#'
#' @param bias Choose between bias 1 (fish only), 2 (shadowing only) or 3 (fish and shadowing).
#' @param showPlots If TRUE plots are shown on the screen, if FALSE plots are stored in the img folder.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt If TRUE use vectorized R commands (faster).
#' @return A dictionary of return objects, see RETURN_DESCRIPTOINS.html for more info.
#' @export
acousticTest <- function(bias=1, showPlots=TRUE, debug=FALSE, opt=TRUE) {
	#### TEST RUN
	params = list()
	#notification option
	params$userEmail = "epy00n@hotmail.com"
	
	## Sensor variables
	params$numSensors = 4
	params$bias = bias
	params$sensorElevation <- 1
        params$shapeFcn <- 'shape.gauss'
	params$peak=.98 
        params$detectionRange <- 25
	
	# BGrid Variables
	params$inputFile = "src/palmyrabath.RData"
	params$inputFileType = "custom"
	params$seriesName = 'z'
	params$cellSize = 5 
	params$startX = 380
	params$XDist = 10
	params$startY = 1250
	params$YDist = 10
	
	## Suppression Variables
	params$suppressionRangeFactor = 2
	params$suppressionFcn = "detection.function"
	## This is only relevant with suppression.scale
	params$maxsuppressionValue = 1
	## This is only relevant with suppression.scale
	params$minsuppressionValue = .5 
	## Choose random walk type movement model
	params$fishmodel <- 'rw'
	if(params$fishmodel == 'ou'){
            ## OU parameter: center of home range
            params$mux <- 0.7 ## Proportion of x scale
            params$muy <- 0.5 ## Proportion of y scale
		
            ## ----- OU: Home range shape and size parameters: -----
            ## SD of home range in x direction, sdx > 0 in meters
            params$ousdx <- 25
            ## SD of home range in y direction, sdy > 0 in meters
            params$ousdy <- 25
            ## Correlation between directions, -1 < cor < 1
            params$oucor <- 0.7
	}
	## Apply vertical habitat range?
	vHabitatRange = TRUE
	if(vHabitatRange){
	    ## Minimum depth (shallowest depth)
	    params$mindepth <- -5
	    ## Maximum depth (deepest depth)
	    params$maxdepth <- -10
	}
	## Apply depth preference?
	depthPref = FALSE
	if(depthPref){
	    ## Depth preference of fish relative to bottom (in meters off the bottom)
	    params$dp <- 2
	    ## Strength of depth preference as a standard deviation, 95% of the time is spent within plus minus two dpsd
	    params$dpsd <- 2
	}
	return(acousticRun(params, showPlots=FALSE, debug, opt))
}

#acousticTest( bias=1)
