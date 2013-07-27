# Main Method
# Supported operations:
#		1. User provides area of interest and sensor locations, asks for stats.
#		2. User provides area of interest and number of sensors, asks for optimal 
#				placement and stats.
#		3. User provides area of interest, fish behaviors, and number of sensors, 
#				asks for optimal placement and stats.

rm(list=ls()) ## Clear all variables
source('src/Bathy.R')
source('src/FishModel.R')
source('src/Utility.R')
library(sendmailR)

run <- function(params, debug=FALSE, opt=FALSE){
	startTime = Sys.time()
    if(debug) {
        cat("\n[run]\n")
    }
    params = checkParams(params)

    ## Create/Load the Bathy grid for the area of interest
    bGrid <- bathy(params$inputFile,
			params$inputFileType,
            params$startX, params$startY, 
            params$XDist, params$YDist,
            params$seriesName,
            debug)

    bGrid = list("bGrid"=bGrid, "cellRatio"=params$cellRatio)
    ## Specify a standard scale of x and y axes if previously undefined
    if(!('x' %in% names(bGrid))) bGrid$x <- seq(0,1,length=dim(bGrid$bGrid)[2])
    if(!('y' %in% names(bGrid))) bGrid$y <- seq(0,1,length=dim(bGrid$bGrid)[1])
    
    fGrid = fish(params, bGrid)
    ## Find good sensor placements
    sensors = sensors(params$numSensors, bGrid, fGrid, params$range, params$bias, params, debug, opt)
    
    ## Stat analysis of proposed setup.
    statDict = stats(params, bGrid, fGrid, sensors$sensorList)
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

test <- function(debug=FALSE, opt=FALSE) {
	#### TEST RUN
	params = list()
	#notification option
	params$userEmail = "epy00n@hotmail.com"
	
	## Array variables
	params$numSensors = 10 
	params$cellRatio = 1
	params$bias = 3
	
	## Receiver variables
	params$sd=10
	params$peak=.98 
	params$shapeFcn= "shape.gauss"
	params$range = 3*params$sd
	
	# BGrid Variables
	params$inputFile = "Pal_IKONOS\\pal_dball.asc"
	params$inputFileType = "custom"
	params$startX = 370
	params$startY = 1200
	params$XDist = 80
	params$YDist = 80
	params$seriesName = 'z'
	
	## Supression variables
	params$supressionFcn = "supression.scale"
	params$supressionRange = 20
	params$maxSupressionValue = 1
	params$minSupressionValue = .5
	## Mean squared displacement of fish (a proxy for movement capacity)
	params$msd <- 0.1
	## Sampling time step
	params$dt <- 1
	## Choose random walk type movement model
	params$fishmodel <- 'rw'
	## Set to TRUE if vertical habitat range is applied
	if(FALSE){
	    ## Minimum depth (shallowest depth)
	    params$mindepth <- -1
	    ## Maximum depth (deepest depth)
	    params$maxdepth <- -10
	}
	## Set to TRUE if depth preference should be applied
	if(FALSE){
	    ## Depth preference of fish relative to bottom (in meters off the bottom)
	    params$dp <- 2
	    ## Strength of depth preference as a standard deviation, 95% of the time is spent within plus minus two dpsd
	    params$dpsd <- 2
	}
	## Set to TRUE of Ornstein-Uhlenbeck (OU) movement should be applied
	if(FALSE){
	    ## Choose Ornstein-Uhlenbeck type movement model
	    params$fishmodel <- 'ou'
	    ## OU parameter: center of home range
	    params$mux <- 0.4
	    params$muy <- 0.2
	    ## OU: Attraction parameter, determines strength of attraction toward home range center
	    params$Bx <- 0.1
	    params$By <- 0.1
	    params$Bxy <- 0
	}
	
	## Print time stamp (to be able to check run time)
	result = run(params, debug, opt)
	return(result)
}

##Rprof(tmp <- tempfile())
##asd <- test(opt=TRUE)
##Rprof()
##summaryRprof(tmp)

#Rprof(tmp <- tempfile())
#asd <- test()
#Rprof()
#summaryRprof(tmp)
test()

