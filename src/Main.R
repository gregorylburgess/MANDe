# Main Method
# Supported operations:
#		1. User provides area of interest and sensor locations, asks for stats.
#		2. User provides area of interest and number of sensors, asks for optimal 
#				placement and stats.
#		3. User provides area of interest, fish behaviors, and number of sensors, 
#				asks for optimal placement and stats.
## Clear all variables
rm(list=ls()) 
source('src/Bathy.R')
source('src/FishModel.R')
source('src/Utility.R')


#' @name run
#' @title  the simulation with the provided parameters.
#'
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return A dictionary of return objects, see RETURN_DESCRIPTIONS.html for more info.
run <- function(params, debug=FALSE, opt=FALSE){
	startTime = Sys.time()
    if(debug) {
        cat("\n[run]\n")
    }
    params = checkParams(params)

    ## Create/Load the Bathy grid for the area of interest
    bGrid = getBathy(params$inputFile, params$inputFileType, params$startX, params$startY, 
            params$XDist, params$YDist, params$seriesName, debug)
    bGrid = list("bGrid"=bGrid, "cellRatio"=params$cellSize)
    ## Convert parameter values from meters to number of grid cell 
    params = convertMetersToGrid(params,bGrid)
    ## Specify a standard scale of x and y axes if previously undefined
    if(!('x' %in% names(bGrid))) {
		bGrid$x = (1:dim(bGrid$bGrid)[1])*params$cellSize 
	}
    if(!('y' %in% names(bGrid))) {
		bGrid$y = (1:dim(bGrid$bGrid)[2])*params$cellSize
	}
    ## Calculate fish grid
    fGrid = fish(params, bGrid)

    ## Find good sensor placements
    sensors <- sensorFun(params$numSensors, bGrid, fGrid, params$range, params$bias, params, debug, opt)

    ## Stat analysis of proposed setup.
    statDict = stats(params, bGrid, fGrid, sensors, debug, opt)
    ## Return Fish grid, Bathy grid, and Sensor Placements as a Dictionary.
    results = list("bGrid" = bGrid, "fGrid" = fGrid, "sumGrid"=sensors$sumGrid, "sensors" = sensors$sensorList, 
            "stats" = statDict, "params"=params)
	## Graph results
	results$filenames = graph(results,params)
	endTime = Sys.time()
	results$runTime = endTime - startTime
	## Email results
	if("userEmail" %in% names(params)) {
		from = "acousticwebapp@gmail.com"
		to = params$userEmail
		subject <- "Acoustic webapp results"
		body <- results                    
		mailControl=list(smtpServer="smtp.gmail.com")
		#sendmail(from=from,to=to,subject=subject,msg=body,control=mailControl)
	}
    return(results)
}

#' @name test
#' @title Executes a test run of the program, using default parameters.  No additional 
#' parameters are necessary.
#'
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return A dictionary of return objects, see RETURN_DESCRIPTOINS.html for more info.
test <- function(debug=FALSE, opt=FALSE) {
	#### TEST RUN
	params = list()
	#notification option
	params$userEmail = "epy00n@hotmail.com"
	
	## Sensor variables
	params$numSensors = 10
	params$bias = 1
	params$sensorElevation <- 1
        params$shapeFcn <- 'shape.gauss'
	params$peak=.98 
        params$detectionRange <- 100
	
	# BGrid Variables
	params$inputFile = "src/palmyrabath.RData"
	params$inputFileType = "asc"
	params$seriesName = 'z'
	params$cellSize = 5 
	params$startX = 350
	params$XDist = 100
	params$startY = 1200
	params$YDist = 100
	
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
	
	return(run(params, debug, opt))
}

##Rprof(tmp <- tempfile())
##asd <- test(opt=TRUE)
##Rprof()
##summaryRprof(tmp)

##system.time(result <- test(opt=TRUE))
#system.time(result <- test(debug=FALSE,opt=TRUE))

if(FALSE){
  print(result$stats$absRecoveryRate)
  print(result$stats$uniqRecoveryRate)
  ns <- length(result$sensors)
  sens <- matrix(unlist(result$sensors),ns,2,byrow=TRUE)
  xlab <- 'x dir'
  ylab <- 'y dir'
  plot.bathy <- TRUE
  graphics.off()
  plotGrid(result,type='fGrid',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
  #dev.new()
  #plotGrid(result,type='acousticCoverage',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
  dev.new()
  plotGrid(result,type='sumGrid',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
  dev.new()
  plotUniqueRR(result)
}
