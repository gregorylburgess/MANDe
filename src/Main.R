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
##library(sendmailR)

run <- function(params, debug=FALSE, opt=FALSE){
	startTime = Sys.time()
    if(debug) {
        cat("\n[run]\n")
    }
    params = checkParams(params)

    ## Create/Load the Bathy grid for the area of interest
    bGrid <- bathy(params$inputFile, params$inputFileType, params$startX, params$startY, 
            params$XDist, params$YDist, params$seriesName, debug)
    bGrid = list("bGrid"=bGrid, "cellRatio"=params$cellRatio)
    ## Specify a standard scale of x and y axes if previously undefined
    if(!('x' %in% names(bGrid))) bGrid$x <- seq(0,1,length=dim(bGrid$bGrid)[1])
    if(!('y' %in% names(bGrid))) bGrid$y <- seq(0,1,length=dim(bGrid$bGrid)[2])
    
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
	params$numSensors = 4
	params$cellRatio = 1
	params$bias = 1
	
	## Receiver variables
	#params$sd=10 ## Palmyra
	params$sd = 5
	params$peak=.98 
	params$shapeFcn= "shape.gauss"
	params$range = 3*params$sd
	
	# BGrid Variables
	params$inputFile = "Pal_IKONOS\\pal_dball.asc"
	params$inputFileType = "arcgis"
	#params$startX = 370
	#params$startY = 1200
	#params$XDist = 80
	#params$YDist = 80
	params$startX = 1
	params$startY = 1
	params$XDist = 100
	params$YDist = 100
	params$seriesName = 'z'
	
	## suppression variables
        params$sparsity <- 0.5 ## This is a lower bound for sparsity
	params$suppressionFcn = "suppression.scale"
	##params$suppressionFcn = "detection.function"
        ## suppression range
        dists <- 1:max(c(params$XDist,params$YDist))
        dfvals <- do.call(params$shapeFcn, list(dists, params))
        params$detectionRange <- dists[min(which(dfvals<0.05))] ##This is different from range above as range is mostly used to cut out areas of grid, whereas detectionRange is closer to what we understand as the actual physical detection range, which is used in sparsity calculations
	params$suppressionRange = params$sparsity*2*params$detectionRange ## Using equation 8 in Pedersen & Weng 2013
	params$maxsuppressionValue = 1
	params$minsuppressionValue = .5
	## Mean squared displacement of fish (a proxy for movement capacity)
	params$msd <- 0.1
	## Sampling time step
	params$dt <- 1
	## Choose random walk type movement model
	params$fishmodel <- 'rw'
	## Set to TRUE if vertical habitat range is applied
	if(TRUE){
	    ## Minimum depth (shallowest depth)
	    params$mindepth <- -2
	    ## Maximum depth (deepest depth)
	    params$maxdepth <- -8
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
	    params$muy <- 0.4
	    ## OU: Attraction parameter, determines strength of attraction toward home range center
	    params$Bx <- 0.02
	    params$By <- 0.02
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

system.time(result <- test(opt=TRUE))

if(FALSE){
  print(result$stats$absRecoveryRate)
  print(result$stats$uniqRecoveryRate)
  ns <- length(result$sensors)
  graphics.off()
  ##image(result$bGrid$x,result$bGrid$y,result$sumGrid,main='sumGrid')
  image(result$bGrid$x,result$bGrid$y,result$fGrid,main='fGrid')
  contour(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,xlab='x',ylab='y',add=TRUE,nlevels=5)
  ##image(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,main='bGrid')
  for(i in 1:ns){
    points(result$bGrid$x[result$sensors[[i]]$r],result$bGrid$y[result$sensors[[i]]$c],pch=21,bg='blue',cex=4)
    text(result$bGrid$x[result$sensors[[i]]$r],result$bGrid$y[result$sensors[[i]]$c],i,col='white')
  }
  dev.new()
  ##image(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,main='bGrid')
  image(result$bGrid$x,result$bGrid$y,result$stats$acousticCoverage,main='Acoustic coverage')
  ##image(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid>=0,main='Land')
  contour(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,xlab='x',ylab='y',add=TRUE,nlevels=5)
  for(i in 1:ns){
    points(result$bGrid$x[result$sensors[[i]]$r],result$bGrid$y[result$sensors[[i]]$c],pch=21,bg='blue',cex=3)
    text(result$bGrid$x[result$sensors[[i]]$r],result$bGrid$y[result$sensors[[i]]$c],i,col='white')
    y <- result$sensors[[i]]$r - 0.5
    x <- result$sensors[[i]]$c - 0.5
    r <- 2*result$params$sd
    a <- seq(0,2*pi,length.out=100)
    X <- r*cos(a)+x
    Y <- r*sin(a)+y
    lines(Y/result$params$YDist,X/result$params$XDist,lty=3)
  }
}
