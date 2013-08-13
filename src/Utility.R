source("src/ShapeFunctions.R")
#' @title u
#' @name util
#' Finds a "good" set of sensor placements for a given setup [bGrid, fGrid, params].
#' Returns a list of locations as grid coordinates.
#' Bias cases:
#' 1. Fish density only.
#' 2. Visibility due to Bathy only.
#' 3. Detectable fish accounting for bathy.
#'
#' @param numSensors The number of sensors the program should place.
#' @param bGrid A valid BGrid.
#' @param fGrid A valid FGrid.
#' @param range The range of the sensor in bathymetric cells.
#' @param bias The goodness algorithm to use, choose 1, 2, or 3.  See above for descriptions.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return A dictionary of return objects, see RETURN_DESCRIPTIONS.html for more info.
sensorFun <- function(numSensors, bGrid, fGrid, range, bias, params, debug=FALSE, opt=FALSE) {
    if (debug) {
        cat("\n[sensorFun]\n")
        print("bGrid")
        print(bGrid)
        print("fGrid")
        print(fGrid)
        print(sprintf("bias=%g",bias))
        print("params")
        print(params)
    }
    
    sensorList = {}
    dims = dim(fGrid)
    rows = dims[1]
    cols = dims[2]
    grids = list("bGrid" = bGrid, "fGrid"=fGrid)
    
    # calculate the sumGrid
    grids = sumGridFun(grids, range, bias, params, debug, opt)
    sumGrid = grids$sumGrid
    # for each sensor, find a good placement
    print(paste('Suppression fcn:',params$suppressionFcn))
    for (i in 1:numSensors) {
        # find the max location 
        maxLoc = which.max(grids$sumGrid)
        # Switch the row/col vals since R references Grid coords as (y,x) instead of (x,y)
        c=ceiling(maxLoc/rows)
        r=(maxLoc %% rows)
        if (r==0) {
            r=rows
        }
        maxLoc = list(c=c,r=r)
        # append maxLoc to the sensor list.
        sensorList = c(sensorList, list(maxLoc))
        # down-weigh all near-by cells to discourage them from being chosen by the program
        if(params$suppressionFcn != 'detection.function.exact'){
          ##print('NOT using detection.function.exact')
          if(opt){
            grids$sumGrid = suppress.opt(grids$sumGrid, dim(fGrid), maxLoc, params, bGrid$bGrid, debug)
          }else{
            grids$sumGrid = supress(grids$sumGrid, dim(fGrid), maxLoc, params$suppressionFcn, 
                                    params$suppressionRange, params$minsuppressionValue, 
                                    params$maxsuppressionValue, params, debug)
          }
      }else{
          grids <- updateFGrid(maxLoc,grids,params,debug,opt)
          grids <- sumGridFun(grids, range, bias, params, debug, opt)
      }
    }
    return(list(sensorList=sensorList, sumGrid=sumGrid, sumGridSupp=grids$sumGrid))
}


#' Updates the FGrid after each sensor is placed to reflect which areas that are already covered by sensors.
#'
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return Returns the grids parameter, with an updated FGrid.
updateFGrid <- function(loc,grids,params,debug=FALSE,opt=FALSE){
  grid <- grids$fGrid
  bG <- grids$bGrid$bGrid
  dims <- dim(grid)
  rows <- dim(grid)[1]
  cols <- dim(grid)[2]
  ## Use range as in sumGrid calculation
  vals = getArea(loc, dims, params$range) 

  rind <- vals$rs:vals$re
  cind <- vals$cs:vals$ce
  nrows <- length(rind)
  ncols <- length(cind)
  Rind <- matrix(rep(rind,ncols),nrows,ncols)
  Cind <- matrix(rep(cind,nrows),nrows,ncols,byrow=TRUE)
  dist <- sqrt( (loc$c-Cind)^2 + (loc$r-Rind)^2 )
  ## Detection fun supp
  dgrid2 <- do.call(params$shapeFcn, list(dist, params)) 

  land <- bG >= 0
  sensorDepth <- bG + params$sensorElevation
  ng <- rows*cols
  nr <- rows
  ## If false then proportion of water column is calculated, if true depth preference is used
  dpflag <- FALSE 
  pctviz <- calc.percent.viz(loc$r,loc$c,rind,cind,ng,nr,bG,land,sensorDepth,dpflag,params)
  testmap <- matrix(0,rows,cols)
  testmap[pctviz$inds] <- pctviz$percentVisibility
  testmap[loc$r,loc$c] <- 1
  ## Line of sight supp
  dgrid1 <- testmap[rind,cind] 
  dgrid <- 1 - (dgrid1 * dgrid2)
  ## Downweigh observed region
  grid[rind,cind] <- grid[rind,cind] * dgrid 
  grids$fGrid <- grid
  return(grids)
}


#' Calculates the composite "goodness" grid for a particular bias.
#'
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param range The range of the sensor in bathymetric cells.
#' @param bias The goodness algorithm to use, choose 1, 2, or 3.  See package manual for more details.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return Returns the grids parameter, with an updated sumGrid.
sumGridFun <- function (grids, range, bias, params, debug=FALSE, opt=FALSE) {
    if (debug) {
        cat("\n[sumGrid]\n")
        print("bGrid")
        print(grids$bGrid)
        print("fGrid")
        print(grids$fGrid)
        print(sprintf("bias=%g", bias))
        print("params")
        print(params)
    }
    #Fish
    if (bias == 1) {
        if(opt){
            return(sumGrid.sumSimple.opt(grids, "fGrid", range, debug))
        }else{
            return(sumGrid.sumSimple(grids, "fGrid", range, debug))
        }
    }
    #Bathy
    else if (bias == 2) {
        if(opt){
            return(sumGrid.sumBathy.opt(grids, params, debug, opt))
        }else{
            return(sumGrid.sumBathy(grids, range, params$shapeFcn, params, debug))
        }
    }
    #Combo
    else if (bias == 3) {
        if(opt){
            ## Note sumGrid.sumBathy.opt also handles bias 3
            return(sumGrid.sumBathy.opt(grids, params, debug,opt))
        }else{
            return(sumGrid.sumProduct(grids, range, params$shapeFcn, params, debug))
        }
    }
    else {
        write("ERROR: Invalid Bias", stderr())
    }   
}


#' Simply sums the values within range of a cell, for each cell in the given grid.
#'
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param key A key to the dictionary provided in the 'grids' parameter specifying which grid should be summed.
#' @param range The range of the sensor in bathymetric cells.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns the grids parameter, with an updated sumGrid.
sumGrid.sumSimple <- function (grids, key, range, debug=FALSE) {
    tempGrid = get(key, grids)
    tempCopy = tempGrid
    rows = dim(tempGrid)[1]
    cols = dim(tempGrid)[2]
    
    for (i in 1:rows) {
        for(j in 1:cols) {
            vals = getArea(list(r=i,c=j), dim(tempGrid), range)
            tempGrid[i,j] = sum(tempCopy[vals$rs:vals$re, 
                            vals$cs:vals$ce], na.rm=TRUE)
        }
    }
    
    grids$sumGrid = tempGrid
    if(debug){
        cat("\n[sumGrid.sumSimple]\n")
        print("grids")
        print(grids)
    }
    return(grids)
}


#' Simply sums the values within range of a cell, for each cell in the given grid.  
#' [optimized, but gives different results than non-opt version, why?, see TestUtility.R for speed comparison]
#'
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param key A key to the dictionary provided in the 'grids' parameter specifying which grid should be summed.
#' @param range The range of the sensor in bathymetric cells.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns the grids parameter, with an updated sumGrid.
sumGrid.sumSimple.opt <- function (grids, key, range, debug=FALSE) {
	## Assume that range is an integer
    kernel <- rep(1,2*range+1) 
    tempGrid = get(key, grids)
    grids$sumGrid <- conv.2D(tempGrid,kernel,kernel)

    if(debug){
        cat("\n[sumGrid.sumSimple.opt]\n")
        print("grids")
        print(grids)
    }
    return(grids)
}


#' Sums the result of calling the detect() function on each cell within range of 
#' a target cell for each cell in the given grid.
#' 
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param range The range of the sensor in bathymetric cells.
#' @param shapeFcn The shape function that describes the attenuation of a sensor.  See ShapeFunctions.R for a list of supported functions.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns the grids parameter, with an updated sumGrid.
sumGrid.sumBathy <- function (grids, range, shapeFcn="shape.t", 
        params, debug=FALSE) {
    
    sumGrid0 = grids$bGrid$bGrid
    tempCpy = grids$bGrid$bGrid
    rows = dim(sumGrid0)[1]
    cols = dim(sumGrid0)[2]
    
    for (i in 1:rows) {
        for(j in 1:cols) {
            vals = getArea(list(r=i,c=j), dim(sumGrid0), range)
            rs = vals$rs
            re = vals$re
            cs = vals$cs
            ce = vals$ce
            visibilities = {}
            for (r in rs:re) {
                for (c in cs:ce) {
                    visibilities = c(
                            visibilities,
                            detect(tempCpy, sensorPos=list(r=i,c=j), tagPos=list(r=r,c=c), shapeFcn=shapeFcn,
                                    params, debug))
                }
            }
            sumGrid0[i,j] = sum(visibilities, na.rm=TRUE)
        }
    }
    
    grids$sumGrid = sumGrid0
    if(debug){
        cat("\n[sumGrid.sumBathy]\n")
        print("visibilities")
        print(visibilities)
        print("grids")
        print(grids)
    }
    return(grids)
}


#' # Sums the result of calling the detect() function on each cell within range of 
#' a target cell for each cell in the given grid.
#' [optimized version, which also supports bias 3]
#' 
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return Returns the grids parameter, with an updated sumGrid.
sumGrid.sumBathy.opt <- function (grids, params, debug=FALSE,opt=FALSE) {

    nr <- dim(grids$bGrid$bGrid)[1]
    nc <- dim(grids$bGrid$bGrid)[2]
    ng <- nr*nc
    bG <- grids$bGrid$bGrid
	## Allocate memory
    sumGrid <- matrix(0,nr,nc)
	## Round to integer range
    rng <- round(params$range)
    sensorDepth <- bG + params$sensorElevation
    belowSurf <- sensorDepth < 0
    land <- bG >= 0
    dpflag <- "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params
    usefGrid <- params$bias==3
    for(c in 1:nc){
        cind <- max(c(1,c-rng)):min(c(nc,c+rng))
        for(r in 1:nr){
			## Only calculate if sensor is below surface
            if(belowSurf[r,c]){ 
                rind <- max(c(1,r-rng)):min(c(nr,r+rng))
                pV <- calc.percent.viz(r,c,rind,cind,ng,nr,bG,land,sensorDepth,dpflag,params)
                probOfRangeDetection = do.call(params$shapeFcn, list(pV$dists, params))
                if(usefGrid) probOfRangeDetection <- probOfRangeDetection * grids$fGrid[pV$inds]
                sumGrid[r,c] = sum(probOfRangeDetection * pV$percentVisibility)
            }
        }
    }
    
    grids$sumGrid = sumGrid
    if(debug){
        cat("\n[sumGrid.sumBathy.opt]\n")
        print("grids")
        print(grids)
    }
    return(grids)
}

#' {{Martin}}
#' Calculates a matrix ()of proportion of water column visibile of the cells surrounding the current cell
#' The current cell looks at the surrounding cells within the detection range
#' and assigns a value to each of those cells, which is the proportion of the visible water column in that cell
#' Feel free to change this poor explanation 
#' 
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return Returns the grids parameter, with an updated sumGrid.
calc.percent.viz <- function(r,c,rind,cind,ng,nr,bG,land,sensorDepth,dpflag,params){
    ## Find rows and columns for relevant cells
    rvec <- rep(rind,length(cind))
    cvec <- sort(rep(cind,length(rind)))
	## Remove self cell
    tmp <- which(!(rvec==r & cvec==c))
    rvec <- rvec[tmp]
    cvec <- cvec[tmp]
	## Translate to single index
    inds <- sub2ind(rvec,cvec,nr) 
    ninds <- length(inds)
    ## Get depths, dists and slopes
    disttmp <- sort(sqrt((r-rvec)^2 + (c-cvec)^2),decreasing=TRUE,index=TRUE) ## This sorts after dist so longest dists are calculated first, then shorter ones might not be needed since they are already calculated for a long dist
    dists <- disttmp$x
    depths <- bG[inds[disttmp$ix]]
    ##print(c(r,c))
    ##print(dim(sensorDepth))
    slopes <- (depths-sensorDepth[r,c])/dists
    ibig2ismall <- rep(0,ng)
    ibig2ismall[inds[disttmp$ix]] <- 1:ninds
	## Assign small negative number to avoid problem with being exactly at the surface in pnorm
	vizDepths <- rep(-1e-4,ninds)
	remaining <- 1:ninds

	## Calculate visible depths
    while(length(remaining)>0){
        ii <- remaining[1]
        losinds <- getCells.new(list(r=r,c=c),list(r=rvec[disttmp$ix[ii]],c=cvec[disttmp$ix[ii]]), debug=FALSE, nr)
		## Get indices in small vectors (not whole grid)
		is <- ibig2ismall[losinds]
        d2 <- sort(dists[is],index=TRUE)
		## If LOS is blocked by land don't calculate for cells behind
        blocks <- land[losinds[d2$ix]] 
        if(any(blocks)){
            if(!all(blocks)){
                indsNoBlock <- 1:(min(which(blocks))-1)
            }else{
                indsNoBlock <- NULL
            }
        }else{
            indsNoBlock <- 1:length(losinds)
        }
            
        ## Here cummax ensures that the steepest slopes is used for calculating visible depth, it is important that the slopes are sorted in order of increasing distance from current cell to target cells, this is handled by d2$ix
        vizDepths[is[d2$ix[indsNoBlock]]] <- cummax(slopes[is[d2$ix[indsNoBlock]]])*d2$x[indsNoBlock] + sensorDepth[r,c]
        remaining <- setdiff(remaining,is)
        ##plot(c(c,cvec[disttmp$ix[ii]]),c(r,rvec[disttmp$ix[ii]]),typ='l',xlim=range(cvec),ylim=range(rvec),main=cc)
        ##points(cvec[disttmp$ix[remaining]],rvec[disttmp$ix[remaining]])
        ##Sys.sleep(0.3)
    }
	
	## Visible depths above water not valid (assign a number a little smaller than zero [just below surface])      
    vizDepths[vizDepths>0] <- -1e-4 
    indsNotLand <- depths<0
    ## if we have normal distribution data (depth preference), use it
    if(dpflag) {
        ## compute % fish visible from sensor to target cell
        mean = depths[indsNotLand] + params$depth_off_bottom
        mean[mean>0] <- 0
        sd = params$depth_off_bottom_sd
            
        ## Get cum probability for "available" water (between depth pref and surf)
        psurf <- pnorm(0,mean=mean,sd=sd)
        areaToCorrectFor <- psurf - pnorm(depths[indsNotLand],mean=mean,sd=sd)
        percentVisibility <- psurf - pnorm(vizDepths[indsNotLand],mean=mean,sd=sd)
        percentVisibility <- percentVisibility/areaToCorrectFor
    }else{
        ## if we don't have normal distribution data, assume equal distribution
        ## compute % visibility (of water column height) from sensor to target cell
        percentVisibility = vizDepths[indsNotLand] / depths[indsNotLand]
    }
    return(list(percentVisibility=percentVisibility,dists=dists[indsNotLand],inds=inds[disttmp$ix[indsNotLand]]))
}


#' Determines the likelihood of a tag at a given position is detectable by a sensor at a 
#' given position, using a specific shapeFunction.  This function considers Bathymetry and
#' sensor range.
#' Returns the percent chance of detection as a double between 0 [no chance of detection] 
#' and 1 [guaranteed detection].
#' [usable only in optimized version].
#' 
#' @param bGrid A valid BGrid.
#' @param sensorPos A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor's location on the BGrid.
#' @param tagPos A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen tag's location on the BGrid.
#' @param shapeFcn The shape function that describes the attenuation of a sensor.  See ShapeFunctions.R for a list of supported functions.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return Returns a value between 0 and 1 representing the probability of detecting a tag due bathymetry.
detect.new <- function(bGrid, sensorPos, tagPos, shapeFcn, params, debug=FALSE,opt=FALSE) {
    rvec = tagPos$rs:tagPos$re
    nr <- length(rvec)
    cvec = tagPos$cs:tagPos$ce
    nc <- length(cvec)
    rmat <- matrix(rep(rvec,nc),nr,nc)
    cmat <- matrix(rep(cvec,nr),nr,nc,byrow=TRUE)

    dist = sqrt((sensorPos$c - cmat)^2 + (sensorPos$r - rmat)^2)
    probOfRangeDetection = do.call(shapeFcn, list(dist, params))
	## Allocate memory
    probOfLOSDetection = matrix(0,nr,nc)
    for (r in 1:nr) {
        for (c in 1:nc) {
            probOfLOSDetection[r,c] = checkLOS(bGrid, sensorPos, list(r=rvec[r],c=cvec[c]), params, debug)
        }
    }
    
    probOfDetection = probOfRangeDetection * probOfLOSDetection
    if(debug) {
        cat("\n[detect]\n")
        print(sprintf("probOfLOSDetection=%g",probOfLOSDetection))
        print(sprintf("probOfRangeDetection=%g",probOfRangeDetection))
        print(sprintf("TotalProbOfDetection=%g",probOfDetection))
    }
    return(probOfDetection)
}


#' For each cell in a given grid, the function sums (the number of fish
#' times the probability of detection) for all cells within range.
#' 
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param range The range of the sensor in bathymetric cells.
#' @param shapeFcn The shape function that describes the attenuation of a sensor.  See ShapeFunctions.R for a list of supported functions.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns the grids parameter, with an updated sumGrid.
sumGrid.sumProduct <- function (grids, range, shapeFcn="shape.t", 
        params, debug=FALSE) {
    
    sumGrid0 = grids$bGrid$bGrid
    fGrid = grids$fGrid
    tempCpy = grids$bGrid$bGrid
    rows = dim(sumGrid0)[1]
    cols = dim(sumGrid0)[2]
    
    for (i in 1:rows) {
        for(j in 1:cols) {
            vals = getArea(list(r=i,c=j), dim(sumGrid0), range)
            rs = vals$rs
            re = vals$re
            cs = vals$cs
            ce = vals$ce
            visibilities = {}
            for (r in rs:re) {
                for (c in cs:ce) {
                    visibilities = c(
                            visibilities,
                            detect(tempCpy, sensorPos=list(r=i,c=j), tagPos=list(r=r,c=c),
									shapeFcn=shapeFcn, params, debug) * fGrid[r,c])
                }
            }
            sumGrid0[i,j] = sum(visibilities, na.rm=TRUE)
        }
    }
    
    grids$sumGrid = sumGrid0
    if(debug){
        cat("\n[sumGrid.sumProduct]\n")
        print("visibilities")
        print(visibilities)
        print("grids")
        print(grids)
    }
    return(grids)
}


#' Supresses the values of cells around a sensor using a specified suppressionFunction.
#' 
#' @param sumGrid A valid SumGrid.
#' @param dims The dimensions of the BGrid.  Just call dim() on the BGrid for this.
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param suppressionFcn The shape function that describes the attenuation of a sensor.  See ShapeFunctions.R for a list of supported functions.
#' @param suppressionRange How far out to apply suppression penalties, in bathymetric cells.
#' @param minsuppressionValue: The minimum allowable value to return.
#' @param maxsuppressionValue: The maximum allowable value to return (also the return value for suppression.static())
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns a suppressed sumGrid.
supress <- function(sumGrid, dims, loc, suppressionFcn, suppressionRange,
                    minsuppressionValue, maxsuppressionValue, params, debug=FALSE) {
    if(debug) {
        cat("\n[supress]\n")
        print(sprintf("suppressionFcn: %s", suppressionFcn))
        print(sprintf("loc: (%g,%g)",loc$c,loc$r))
        print("sumGrid")
        print(sumGrid)
    }
    vals = getArea(loc, dims, suppressionRange)
    mini = vals$rs
    maxi = vals$re
    minj = vals$cs
    maxj = vals$ce
    for (i in mini:maxi) {
        for (j in minj:maxj) {
                    dist = sqrt((loc$c - j)^2 + (loc$r - i)^2)
					sumGrid[i,j] = sumGrid[i,j] * do.call(suppressionFcn, list(dist, suppressionRange, 
                                        minsuppressionValue, maxsuppressionValue, params, debug))
            }
    }
    return(sumGrid)
}


#' Supresses the values of cells around a sensor using a specified suppressionFunction.
#' [Optimized using vectorization].
#' 
#' @param A valid SumGrid.
#' @param dims The dimensions of the BGrid.  Just call dim() on the BGrid for this.
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param A valid BGrid.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns a suppressed sumGrid.
suppress.opt <- function(sumGrid, dims, loc, params, bGrid, debug=FALSE) {
    if(debug) {
        cat("\n[suppress.opt]\n")
        print(sprintf("suppressionFcn: %s", suppressionFcn))
        print(sprintf("loc: (%g,%g)",loc$c,loc$r))
        print("sumGrid")
        print(sumGrid)
		print("bGrid")
		print(bGrid)
    }
    suppressionFcn <- params$suppressionFcn
    minsuppressionValue <- params$minsuppressionValue
    maxsuppressionValue <- params$maxsuppressionValue
    ## dfflag indicates wheter detection function should be used for suppression
    dfflag <- suppressionFcn=='detection.function' | suppressionFcn=='detection.function.shadow' | suppressionFcn=='detection.function.exact'
    rows <- dim(sumGrid)[1]
    cols <- dim(sumGrid)[2]
	
    if(dfflag){
      vals = getArea(loc, dims, params$range) ## Use range as in sumGrid calculation
    }else{
      vals = getArea(loc, dims, params$suppressionRange)
    }
    rind <- vals$rs:vals$re
    cind <- vals$cs:vals$ce
    nrows <- length(rind)
    ncols <- length(cind)
    Rind <- matrix(rep(rind,ncols),nrows,ncols)
    Cind <- matrix(rep(cind,nrows),nrows,ncols,byrow=TRUE)
    dist <- sqrt( (loc$c-Cind)^2 + (loc$r-Rind)^2 )

    if(suppressionFcn=='suppression.static'){
      supgrid <- matrix(maxsuppressionValue,nrows,ncols)
    }
    if(suppressionFcn=='suppression.scale'){
      sRange = minsuppressionValue - maxsuppressionValue
      supgrid = 1 - (sRange * (dist/params$suppressionRange) + maxsuppressionValue)
      supgrid[supgrid<0] <- 0
      supgrid[supgrid>1] <- 1
    }
    if(dfflag){
      supgrid2 <- do.call(params$shapeFcn, list(dist, params)) ## Detection fun supp
      if(suppressionFcn=='detection.function.shadow' | suppressionFcn=='detection.function.exact'){
        land <- bGrid >= 0
        sensorDepth <- bGrid + params$sensorElevation
        ng <- rows*cols
        nr <- rows
		## If false then proportion of water column is calculated, if true depth preference is used
		## {{Martin}} should this be a parameter?
        dpflag <- FALSE 
        pctviz <- calc.percent.viz(loc$r, loc$c, rind, cind, ng, nr, bGrid, land,
								   sensorDepth, dpflag, params)
        testmap <- matrix(0,rows,cols)
        testmap[pctviz$inds] <- pctviz$percentVisibility
        testmap[loc$r,loc$c] <- 1
		## Line of sight supp
        supgrid1 <- testmap[rind,cind]
        supgrid <- 1 - (supgrid1 * supgrid2)
      }else{
        supgrid <- 1 - supgrid2
      }
    }
	## Do suppression
    sumGrid[rind,cind] <- sumGrid[rind,cind] * supgrid

    return(sumGrid)
}

#' Returns a static value, effectively setting all cells within suppressionRange of a sensor to that number.
#' 
#' @param dist The distance between two cells on a grid.
#' @param suppressionRange How far out to apply suppression penalties, in bathymetric cells.
#' @param minsuppressionValue: The minimum allowable value to return.
#' @param maxsuppressionValue: The maximum allowable value to return (also the return value for suppression.static())
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns The value given in maxSuppression Value.
suppression.static <- function (dist, suppressionRange, minsuppressionValue, 
                               maxsuppressionValue, params, debug=FALSE) {
    return (maxsuppressionValue)
}


#' Returns a dynamic value based on distance from a chosen sensor. The returned value should be multiplied by the value to be scaled.
#'
#' @param dist The distance between two cells on a grid.
#' @param suppressionRange How far out to apply suppression penalties, in bathymetric cells.
#' @param minsuppressionValue: The minimum allowable value to return.
#' @param maxsuppressionValue: The maximum allowable value to return (also the return value for suppression.static())
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns The value given in maxSuppression Value.
suppression.scale <- function (dist, suppressionRange, minsuppressionValue, 
        maxsuppressionValue, params, debug=FALSE) {
    
    sRange = minsuppressionValue - maxsuppressionValue
    value = 1 - (sRange * (dist/suppressionRange) + maxsuppressionValue)
    value = max(0, value)
    value = min(1, value)
    if (debug) {
        print(sprintf("dist=%g", dist))
        print(sprintf("suppressionRange=%g", suppressionRange))
        print(sprintf("minsuppressionValue=%g", minsuppressionValue))
        print(sprintf("maxsuppressionValue=%g", maxsuppressionValue))
        print(sprintf("value=%g", value))
    }
    return (value)
}


#' Defines the "shape" of a sensor's range.
#' 
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param dims The dimensions of the BGrid.  Just call dim() on the BGrid for this.
#' @param range The range of the sensor in bathymetric cells.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns a dictionary of start/end indexes for rows and columns respectively named : {rs,re,cs,ce}.
getArea<-function(loc, dims, range, debug=FALSE) {
	# the row index for our central point
    r = loc$r
	# the col index for our central point
    c = loc$c
	# the max number of rows in the grid
    rows = dims[1]
	# the max number of cols in the grid
    cols = dims[2]
    
    # defines a square
    rs0 = max(1, r - range) 
    re0 = min(rows, r + range)
    cs0 = max(1 ,c - range)
    ce0 = min(cols, c + range)
    toRet = list(rs=rs0, re=re0, cs=cs0, ce=ce0)
    return(toRet)
}


#'Determines the likelihood of a tag at a given position is detectable by a sensor at a 
#' given position, using a specific shapeFunction.  This function considers Bathymetry and
#' sensor range.
#' 
#' @param bGrid A valid BGrid.
#' @param sensorPos A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor's location on the BGrid.
#' @param tagPos A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen tag's location on the BGrid.
#' @param shapeFcn The shape function that describes the attenuation of a sensor.  See ShapeFunctions.R for a list of supported functions.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return The percent chance of detection as a double between 0 [no chance of detection] and 1 [guaranteed detection].
detect <- function(bGrid, sensorPos, tagPos, shapeFcn, params, debug=FALSE) {

    dist = sqrt((sensorPos$c - tagPos$c)^2 + (sensorPos$r - tagPos$r)^2)
    probOfRangeDetection = do.call(shapeFcn, list(dist, params))
    
    probOfLOSDetection = checkLOS(bGrid, sensorPos, tagPos, params, debug)
    probOfDetection = probOfRangeDetection * probOfLOSDetection
    if(debug) {
        cat("\n[detect]\n")
        print(sprintf("probOfLOSDetection=%g",probOfLOSDetection))
        print(sprintf("probOfRangeDetection=%g",probOfRangeDetection))
        print(sprintf("TotalProbOfDetection=%g",probOfDetection))
    }
    return(probOfRangeDetection * probOfLOSDetection)
}



#' Returns the percent of the water column visible at a target cell from a starting cell.  
#' If 'depth_off_bottom' and 'depth_off_bottom_sd' are keys in 'params', then the algorithm 
#' assumes a normal distribution of fish within the specified zone and will return
#' the integral of the visible zone.  Otherwise, it assumes an equal distribution of fish 
#' throughout the watercolumn.
#' 
#' @param bGrid A valid BGrid.
#' @param startingCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor's location on the BGrid.
#' @param targetCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen tag's location on the BGrid.
#' @param shapeFcn The shape function that describes the attenuation of a sensor.  See ShapeFunctions.R for a list of supported functions.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return The percent of the watercolumn that is visible (as a double between 0 and 1).
checkLOS<- function(bGrid, startingCell, targetCell, params, debug=FALSE) {
    sensorElevation = params["sensorElevation"]
    dist = sqrt((startingCell$c - targetCell$c)^2 + (startingCell$r - targetCell$r)^2)
    if (dist ==0) {
        return(1)
    }
    # our sensor's z value
    sensorDepth = bGrid[startingCell$r, startingCell$c] + sensorElevation
    # retrieve list of intervening cells
    table = getCells(startingCell, targetCell, debug)

	# annotate each cell's z value from the bGrid
    table$z <-apply(table, 1, function(rows){ table$z = bGrid[rows[2],rows[1]]})
    # annotate each cell's percieved slope form our sensor to the cell
    table$m <-apply(table, 1, function(row) { 
                table$m = (row[3] - sensorDepth) / sqrt((row[1] - startingCell$c)^2 + (row[2] - startingCell$r)^2)
            })
    # take the max of all slopes as the limit on our LoS
    m = max(table$m)
    b = sensorDepth
    # y = mx + b
    targetCellsVisibleDepth = m*dist + b
    percentVisibility = 0
    
    # if we have normal distribution data, use it
    if( "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params) {
        # compute % fish visible from sensor to target cell
        mean = bGrid[targetCell$r,targetCell$c] + params$depth_off_bottom
        sd = depth_off_bottom_sd
		
        # pnorm gives the percent below the given point, so subtract from 1
        # to get the percent above the given point
		areaToCorrectFor = 1 - (pnorm(bGrid[targetCell$r,targetCell$c]))
        percentVisibility = 1 - (pnorm(targetCellsVisibleDepth,mean=mean,sd=sd))
		percentVisibility = percentVisibility/areaToCorrectFor
    }
    # if we don't have normal distribution data, assume equal distribution
    else {
        # compute % visibility (of water column height) from sensor to target cell
        percentVisibility = targetCellsVisibleDepth / bGrid[targetCell$r,targetCell$c]
    }
    
    percentVisibility = min(1, percentVisibility)
    percentVisibility = max(0, percentVisibility)
    if (debug) {
        cat("\n[checkLOS]\n")
        print(sprintf("sensorDepth=%g",sensorDepth))
        print(sprintf("dist/z: y=%gx+%g", m,b))
        print("Table:")
        print(table)
        print(sprintf("targetCellsVisibleDepth=%g", targetCellsVisibleDepth))
        print(sprintf("percentVisibility=%g",percentVisibility))
    }
    return(percentVisibility)
}


#' Returns the cells crossed by a beam from the starting cell to the target cell.
#' 
#' @param startingCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor's location on the BGrid.
#' @param targetCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen tag's location on the BGrid.
#' @param debug If enabled, turns on debug printing (console only).
#' @return A data.frame containing x and y indicies of cells crossed by a beam from the starting cell to the target cell.
getCells<-function(startingCell, targetCell, debug=FALSE) {
    sC=offset(startingCell)
    tC=offset(targetCell)
    m = (sC$r-tC$r)/(sC$c-tC$c)
    if(abs(m)== Inf) {
        if(m>0) {
            m = 999999
        }
        else {
            m = -999999
        }
    }
    # assume the sensor is in the middle of the cell
    b = sC$r-m*sC$c
    lowerX = min(startingCell$c, targetCell$c)
    upperX = max(startingCell$c, targetCell$c)
    lowerY = min(startingCell$r, targetCell$r)
    upperY = max(startingCell$r, targetCell$r)
    tx = {}
    ty= {}
    
    #STEEP SLOPES
    if(abs(m)>1) {
        startY = lowerY
        endY = upperY
        if(m<0){
            temp = startY
            startY = endY
            endY = temp
        }
        for( y in startY:endY) {
            x = (y-b)/m
            x1= ceiling(x)
            tx=c(tx,x1)
            ty =c(ty,y)
            if(y+1<=upperY) {
                tx=c(tx,x1)
                ty=c(ty,y+1)
            }
        }
    } else {
        #SLOW SLOPES
        startX = lowerX
        endX = upperX
        if(m<0){
            temp = lowerX
            startX = upperX
            endX = temp
        }
        for( x in startX:endX) {
            y = m * x + b
            y1 = ceiling(y)
            if(y == y1) {
                tx = c(tx,x)
                if (m<0) { 
                    ty= c(ty,y1+1)
                }
                else {
                    ty= c(ty,y1)
                }
            } else {
                tx=c(tx,x)
                ty =c(ty,y1)
                if(x+1<=upperX) {
                    tx=c(tx,x+1)
                    ty=c(ty,y1)
                }
            }
        }
    }
    start = list('x'=startingCell$c,'y'=startingCell$r)
    end = list('x'=targetCell$c,'y'=targetCell$r)
    grid = data.frame("x"=tx,"y"=ty)
    # return only unique values
    grid = unique(grid)
    # remove the starting cell
    grid = grid[!(grid$x == startingCell$c & grid$y == startingCell$r),]
    #grid = grid[!(grid$x == targetCell$c & grid$y == targetCell$r),]
    if(debug) {
        cat("\n[getCells]\n")
        print(sprintf("x/y: y = %gx + %g",m,b))
        print(sprintf("Starting Cell:(%g,%g)",startingCell$c,startingCell$r))
        print(sprintf("Target Cell: (%g,%g)",targetCell$c,targetCell$r))
        print("Table:")
        print(grid)
    }
    return(grid)
}


#' Returns the cells crossed by a beam from the starting cell to
#' the target cell. [Optimized version, see TestUtility.R to evaluate speed gain].
#' 
#' @param startingCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor's location on the BGrid.
#' @param targetCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen tag's location on the BGrid.
#' @param debug If enabled, turns on debug printing (console only).
#' @return A data.frame containing x and y indicies of cells crossed by a beam from the starting cell to the target cell.
getCells.opt <- function(startingCell, targetCell, debug=FALSE) {
  if(!(startingCell$r==targetCell$r & startingCell$c==targetCell$c)){
        sC <- offset(startingCell)
        tC <- offset(targetCell)
        e <- 1e-6
        a <- (tC$r-sC$r)/(tC$c-sC$c)
        if(abs(a)<=1){
            b <- sC$r - a*sC$c
            cPoints <- sort((startingCell$c):targetCell$c)-1
            cols <- ceiling(cPoints+e)
            rows <- ceiling(a*(cPoints+e) + b)
            if(a!=1){
                inds <- which(abs(diff(rows))>0)
                cols <- c(cols,cols[inds])
                rows <- c(rows,rows[inds+1])
            }
        }else{
            a <- 1/a
            b <- sC$c - a*sC$r
            rPoints <- sort((startingCell$r):targetCell$r)-1
            rows <- ceiling(rPoints+e)
            cols <- ceiling(a*(rPoints+e) + b)
            inds <- which(abs(diff(cols))>0)
            rows <- c(rows,rows[inds])
            cols <- c(cols,cols[inds+1])
        }
        ##crds <- data.frame("x"=cols,"y"=rows)
        useinds <- !(cols == startingCell$c & rows == startingCell$r)
        cellmat <- cbind(cols[useinds],rows[useinds])
        colnames(cellmat) <- c('x','y')
        return(as.data.frame(cellmat))
        ##return( crds[!(crds$x == startingCell$c & crds$y == startingCell$r),] )
        ##return(ret)
    }else{
        return(NULL)
    }
}


#{{Martin}} why are there two copies of this function?  Is there a significant difference?
#' Returns the cells crossed by a beam from the starting cell to
#' the target cell.
#' [Does not produce same output as getCells, same reported cells but order is different].
#' 
#' @param startingCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor's location on the BGrid.
#' @param targetCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen tag's location on the BGrid.
#' @param debug If enabled, turns on debug printing (console only).
#' @return A data.frame containing x and y indicies of cells crossed by a beam from the starting cell to the target cell.
getCells.new <- function(startingCell, targetCell, debug=FALSE, nr=NULL) {
    if(!(startingCell$r==targetCell$r & startingCell$c==targetCell$c)){
        sC <- offset(startingCell)
        tC <- offset(targetCell)
        e <- 1e-6
        a <- (tC$r-sC$r)/(tC$c-sC$c)
        if(abs(a)<=1){
            b <- sC$r - a*sC$c
            cPoints <- sort((startingCell$c):targetCell$c)-1
            cols <- ceiling(cPoints+e)
            rows <- ceiling(a*(cPoints+e) + b)
            if(a!=1){
                inds <- which(abs(diff(rows))>0)
                cols <- c(cols,cols[inds])
                rows <- c(rows,rows[inds+1])
            }
            if(!is.null(nr)) biginds <- sub2ind(rows,cols,nr)
        }else{
            a <- 1/a
            b <- sC$c - a*sC$r
            rPoints <- sort((startingCell$r):targetCell$r)-1
            rows <- ceiling(rPoints+e)
            cols <- ceiling(a*(rPoints+e) + b)
            inds <- which(abs(diff(cols))>0)
            rows <- c(rows,rows[inds])
            cols <- c(cols,cols[inds+1])
            if(!is.null(nr)) biginds <- sub2ind(rows,cols,nr)
        }
        useinds <- !(cols == startingCell$c & rows == startingCell$r)
        if(is.null(nr)){
          ##crds <- data.frame("x"=cols,"y"=rows)
          return( cbind(cols[useinds],rows[useinds]) )
          ##return( crds[!(crds$x == startingCell$c & crds$y == startingCell$r),] )
        }else{
          return( biginds[useinds] )
        }
    }else{
        return(NULL)
    }
}


#' Translates a cartesian point to the center of the gridcell it represents.  This is necessary
#' for drawing Line of Sight since we assume sensors are in the center of cells.
#' ex: the cartesian point (3,2) would be converted to (2.5, 1.5), which puts it in the
#' cell located at the third column, second row (aka the cell at (3,2) on a 1-based grid
#' system).
#' 
#' @param point A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the point to translate on the BGrid.
#' @return A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the translated point.
offset<- function(point){
    r= point$r
    c=point$c
    if(r>0) {
        r=r-.5
    } else {
        r=r+.5
    }
    if(c>0){
        c=c-.5
    } else {
        c=c+.5
    }
    return(list("r"=r,"c"=c))
}


#' Generates .png files for the visualizations of various grids.  You must have write access to the 
#' R working directory that the program is executed from.  Additionally, ensure that an 'img' folder exists there.
#' 
#' @param results A dictionary of return objects, the result of a successfull call to run() or sensorFun().
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param plot.bathy Specifies whether contour lines for bathymetry should be overlayed in the graphs.
#' @return A dictionary containing the filenames of the generated images.
graph <- function(result, params, plot.bathy=TRUE) {
	## Plotting
	graphics.off()
	filenames = {}
	time = as.numeric(Sys.time()) %% 1
        xlab <- 'x dir'
        ylab <- 'y dir'

	#{{Martin}} can we plot these dynamically instead of having a method for each one (perhaps with the exception of the RR graphs?}}
    
	## BGrid
	filenames$bGrid = sprintf("img/bGrid-%g.png", time)
	png(filenames$bGrid)
        plotBGrid(result,xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	dev.off()
	
	## FGrid
	filenames$fGrid = sprintf("img/fGrid-%g.png", time)
	png(filenames$fGrid)
        plotFGrid(result,xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	dev.off()
	
	## SumGrid
	filenames$sumGrid = sprintf("img/sumGrid-%g.png", time)
	png(filenames$sumGrid)
        plotSumGrid(result,xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	dev.off()
	
	## Acoustic Coverage
	filenames$acousticCoverage = sprintf("img/acousticCoverage-%g.png", time)
	png(filenames$acousticCoverage)
        plotAcousticCoverage(result,xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	dev.off()

    ## Unique Recovery Rate
	filenames$recoveryRates = sprintf("img/recoveryRates-%g.png", time)
	png(filenames$recoveryRates)
        plotUniqueRR(result)
	dev.off()

	return(filenames)
}


plotAcousticCoverage <- function(result,xlab='',ylab='',plot.bathy=TRUE){
    image(result$bGrid$x,result$bGrid$y,result$stats$acousticCoverage,main='Acoustic Coverage',xlab=xlab,ylab=ylab)
    if(plot.bathy) contour(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,add=TRUE,nlevels=5)
    plotSensors(result)
}


## Plot sumGrid
plotSumGrid <- function(result,xlab='',ylab='',plot.bathy=TRUE){
    image(result$bGrid$x,result$bGrid$y,result$sumGrid,main='Goodness grid',xlab=xlab,ylab=ylab)
    if(plot.bathy) contour(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,add=TRUE,nlevels=5)
    plotSensors(result)
}


## Plot fGrid
plotFGrid <- function(result,xlab='',ylab='',plot.bathy=TRUE){
    image(result$bGrid$x,result$bGrid$y,result$fGrid,main='fGrid',xlab=xlab,ylab=ylab)
    if(plot.bathy) contour(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,add=TRUE,nlevels=5)
    plotSensors(result)
}


## Plot bGrid
plotBGrid <- function(result,xlab='',ylab='',plot.bathy=TRUE){
    image(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,main='bGrid',xlab=xlab,ylab=ylab)
    if(plot.bathy) contour(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,add=TRUE,nlevels=5)
    plotSensors(result)
}


## Plot unique recovery rate as a function of number of sensors
plotUniqueRR <- function(result){
    ns <- length(result$sensors)
    nsmax <- length(result$stats$uniqRRs)
    par(mfrow=c(2,1),las=1)
    plot(0:ns,c(0,result$stats$uniqRRs[1:ns]),typ='l',xlab='Number of sensors',ylab='Unique recovery rate',ylim=c(0,1.02),xlim=c(0,nsmax))
    points(0:ns,c(0,result$stats$uniqRRs[1:ns]),pch=46,cex=3)
    lines(ns:length(result$stats$uniqRRs),result$stats$uniqRRs[ns:nsmax],lty=2)

    plot.intersect(ns,result$stats$uniqRecoveryRate,col='orange',lty=1)
    grid()
    text(0.05*length(result$stats$uniqRRs),result$stats$uniqRecoveryRate,round(result$stats$uniqRecoveryRate,digits=4),pos=3)
    legend('bottomright',c('Calculated','Projected','Requested'),lty=c(1,2,1),col=c(1,1,'orange'),bg='white')
    duRR <- diff(c(0,result$stats$uniqRRs))
    plot(1:ns,duRR[1:ns],typ='l',xlab='Number of sensors',ylab='Increase in unique RR',ylim=c(0,max(duRR)),xlim=c(0,nsmax))
    points(1:ns,duRR[1:ns],pch=46,cex=3)
    lines(ns:nsmax,duRR[ns:nsmax],lty=2)
    grid()
}


## Adds sensors to current plot (an existing plot is required)
plotSensors <- function(result,circles=TRUE,circlty=3){
  ns <- length(result$sensors)
  r <- result$params$detectionRange ## Radius of circle
  a <- seq(0,2*pi,length.out=100)
  sensx <- result$bGrid$x[result$stats$sensorMat[1:ns,2]] ## Cols
  sensy <- result$bGrid$x[result$stats$sensorMat[1:ns,1]] ## Rows
  points(sensx,sensy,pch=21,bg='blue',cex=3)
  text(sensx,sensy,1:ns,col='white')
  for(i in 1:ns){
    ##text(sensx[i],sensy[i],i,col='white')
    if(circles){
      y <- result$sensors[[i]]$r - 0.5
      x <- result$sensors[[i]]$c - 0.5
      X <- r*cos(a)+x
      Y <- r*sin(a)+y
      lines(Y/result$params$YDist,X/result$params$XDist,lty=circlty)
    }
  }
}


#' Provides Statistical data on detection, given a particular bGrid, fGrid, and sensor arrangement.
#' Returns a dictionary of staistical values.
#' 
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param bGrid A valid BGrid.
#' @param fGrid A valid FGrid.
#' @param sensors The result of a successful call to sensorFun().
#' @return A dictionary of statistical values.
stats <- function(params, bGrid, fGrid, sensors, debug=FALSE, opt=FALSE) {
    statDict <- list()
    numSensors <- length(sensors$sensorList)
    rows <- dim(fGrid)[1]
    cols <- dim(fGrid)[2]
    
    ## Calculate the value of each sensor (as the increase in unique recoveryrate)
    numProj <- 2*numSensors
	## sumGrid suppressed by numSensors
    sumGridSupp <- sensors$sumGridSupp
    sensorList <- sensors$sensorList
    ## Calculate locations of projected sensors
    for (i in 1:(numProj-numSensors)) { 
      ## find the max location 
      maxLoc = which.max(sumGridSupp)
      ## Switch the row/col vals since R references Grid coords differently
      c=ceiling(maxLoc/rows)
      r=(maxLoc %% rows)
      if (r==0) {
        r=rows
      }
      maxLoc = list(c=c,r=r)
      ## append maxLoc to the sensor list.
      sensorList = c(sensorList, list(maxLoc))
      ## down-weigh all near-by cells to discourage them from being chosen by the program
      if(opt){
        sumGridSupp = suppress.opt(sumGridSupp, dim(fGrid), maxLoc, params, bGrid$bGrid, debug)
      }else{
        sumGridSupp = supress(sumGridSupp, dim(fGrid), maxLoc, params$suppressionFcn, 
          params$suppressionRange, params$minsuppressionValue, 
          params$maxsuppressionValue, params, debug)
      }
    }

    xSens <- rep(0,numProj)
    ySens <- rep(0,numProj)
    for(i in 1:numProj){
        xSens[i] <- sensorList[[i]]$c
        ySens[i] <- sensorList[[i]]$r
    }
    ## Calculate distance matrix needed to calculate sparsity
    distMat <- matrix(0,numSensors,numSensors)
    for(i in 1:numSensors){
        distMat[i,] <- sqrt((xSens[i]-xSens[1:numSensors])^2 + (ySens[i]-ySens[1:numSensors])^2)
    }
    
    ## a is the median of the distances between the receivers
    print(a <- median(distMat[upper.tri(distMat)]))
    ## delta is a sparsity measure (see Pedersen & Weng 2013)
    statDict$delta <- a/(2*params$detectionRange) 
    ## phi is a dimensionless indicator of movement capacity relative to detection range, it can also be viewed as a signal to noise ratio
    ##statDict$phi <- params$msd/params$detectionRange
    ## Distance maps (the distance from any grid cell to a receiver)
    rows <- dim(fGrid)[1]
    cols <- dim(fGrid)[2]
    X <- matrix(rep(1:cols,rows),rows,cols,byrow=TRUE)
    Y <- matrix(rep(1:rows,cols),rows,cols,byrow=FALSE)
    dimap <- array(0,dim=c(rows,cols,numProj))
    for(i in 1:numProj) dimap[,,i] <- sqrt( (X-xSens[i])^2 + (Y-ySens[i])^2 ) ## Distance to receiver
    
    ## Horizontal detection maps using detection function
    demap <- array(0,dim=c(rows,cols,numProj))
    for(i in 1:numProj) demap[,,i] <- do.call(params$shapeFcn, list(dimap[,,i], params))
    ##filled.contour(gy,gx,demap[,,2],color.palette=rainbow,main=c('Detection map'),xlab='x',ylab='y',plot.axes = {axis(1); axis(2); points(r$x,r$y)})

    ## Incorporate vertical detection probability using line of sight
    ng <- rows*cols
    bG <- bGrid$bGrid
    nr <- dim(bG)[1]
    land <- bG >= 0
    sensorDepth <- bG + params$sensorElevation
    rng <- params$range
	#{{Martin}} should this be a parameter?
	## If false then proportion of water column is calculated, if true depth preference is used
    dpflag <- FALSE 
    for(i in 1:numProj){
      r <- ySens[i]
      c <- xSens[i]
      cind <- max(c(1,c-rng)):min(c(cols,c+rng))
      rind <- max(c(1,r-rng)):min(c(rows,r+rng))
      pctviz <- calc.percent.viz(ySens[i],xSens[i],rind,cind,ng,nr,bG,land,sensorDepth,dpflag,params)
      testmap <- matrix(0,rows,cols)
      testmap[pctviz$inds] <- pctviz$percentVisibility
      testmap[r,c] <- 1
      demap[,,i] <- demap[,,i] * testmap
    }

    ## Coverage map
    cover <- matrix(1,rows,cols)
    covertmp <- array(0,dim=c(rows,cols,numProj))
    uniqRRs <- rep(0,numProj)
    for(i in 1:numProj){
		## Probability of no detection
        cover <- cover * (1 - demap[,,i]) 
        covertmp[,,i] <- 1-cover
        uniqRRs[i] <- sum(covertmp[,,i] * fGrid)
        ##if(i==numSensors)  statDict$acousticCoverage <- covertmp ## Save coverage map for numSensors
    }
    duniqRRs <- diff(c(0,uniqRRs))
    print(duniqRRs)
    srt <- sort(duniqRRs,index=TRUE,decreasing=TRUE) ## Sort list so best sensors come first
    sensorMat <- matrix(unlist(sensorList),numProj,2,byrow=TRUE)
    
    statDict$sensorMat <- sensorMat[srt$ix,]
    statDict$uniqRRs <- cumsum(srt$x)
    statDict$acousticCoverage <- 1-apply(1-demap[,,srt$ix[1:numSensors]],c(1,2),prod) ## Acoustic coverage for best sensors
    
    ## Absolute recovery rate (here we don't care about getting the same ping multiple times)
    demapmat <- apply(demap[,,srt$ix[1:numSensors]],c(1,2),sum)
    statDict$absRecoveryRate <- sum(demapmat*fGrid)
	
    ## Calculate recovery rate of unique detections    
    statDict$uniqRecoveryRate <- sum(statDict$acousticCoverage * fGrid)
    
    return(statDict)
}


#' Provides default parameter values if none are provided.
#' 
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @return The 'params' parameter, populated with default values where necessary.
checkParams <- function(params) {

    names = names(params)
	## Cast all possible strings to numbers (JSON makes everything strings)
	for (name in names) {
		if(!is.na((as.numeric(params[name])))) {
			params[name] = as.numeric(params[name])
		}
	}
    if(!('numSensors' %in% names)) {
        write("Error: 'numSensors' is required", stderr())
    }
    if(!('bias' %in% names)) {
        write("Error: 'bias' value is required.")
    }
	if('dp' %in% names  && !('inputFile' %in% names)) {
		write("Error: Using dp option without a known input file may be bad!.
				For example, if the generated habitat grid contains no cells near
				the depth specified, no fish will be generated.", stderr())
	}

    # suppression Function Defaults
    if(!('suppressionFcn' %in% names)) {
        params$suppressionFcn = "suppression.static"
        params$suppressionRange = 2
        params$maxsuppressionValue = 0
        params$minsuppressionValue = 0
    }	else {
		params$suppressionFcn = as.character(params$suppressionFcn)
	}
    
    # Shape Function Defaults
    if(!('shapeFcn' %in% names)) {
        params$shapeFcn= "shape.gauss"
        params$sd=.3334
        params$peak=.75 
    }	else {
		params$shapeFcn = as.character(params$shapeFcn)
	}
	if(!('range' %in% names)) {
		params$range = 3*params$sd
	}
    if(!('sensorElevation' %in% names)){
        params$sensorElevation <- 1
    }   else {
                params$sensorElevation = as.numeric(params$sensorElevation)
        }
    
    # Bathymetry defaults
	if(('inputfile' %in% names)) {
		params$inputfile = as.character(params$inputfile)
	}
    if(!('cellRatio' %in% names)) {
        params$cellRatio = 1
    }
    if(!('startX' %in% names)) {
        params$startX = 9000
    }
    if(!('startY' %in% names)) {
        params$startY = 8000
    }
    if(!('XDist' %in% names)) {
        params$startY = 10
    }
    if(!('YDist' %in% names)) {
        params$startY = 10
    }
    if(!('seriesName' %in% names)) {
        params$seriesName = 'z'
    }	else {
		params$seriesName = as.character(params$seriesName)
	}
    
    # Fish Modeling
    if(!('fishmodel' %in% names)) {
        params$fishmodel <- 'rw'
    }	else {
		if(params$fishmodel == "True") {
			params$fishmodel = 'ou'
		}
		if(params$fishmodel == "False") {
			params$fishmodel = 'rw'
		}
		params$fishmodel = as.character(params$fishmodel)
	}
    
    return(params)
}


#' Performs the convolution operation in 1D.
#' 
#' @param fun {{Martin}} please fill these in. 
#' @param kern something.
#' @return something.
conv.1D <- function(fun,kern){
  lk <- length(kern)
  kern <- kern[lk:1]
  lk2 <- 0.5*(lk-1)
  lf <- length(fun)
  test <- convolve(fun,kern,type='o')
  test[(lk2+1):(lf+lk-lk2-1)]
}


#' Performs the convolution operation in 2D (using two 1D kernels though)
#' 
#' @param mat {{Martin}} please fill these in.
#' @param kx something.
#' @param ky something.
#' @return something.
conv.2D <- function(mat,kx,ky){
  dimmat <- dim(mat)
  matout <- matrix(0,dimmat[1],dimmat[2])
  for(i in 1:dimmat[1]) matout[i,] <- conv.1D(mat[i,],kx)
  for(i in 1:dimmat[2]) matout[,i] <- conv.1D(matout[,i],ky)
  matout
}


#' Converts a row, col index to a linear index within a matrix.
#' 
#' @param row The row index.
#' @param col The col index.
#' @param dims The dimensions of the BGrid.  Just call dim() on the parent matrix for this.
#' @return The translated linear index.
sub2ind <- function(row,col,dims){
    (col-1)*dims[1] + row
}


#' Plots a circle with radius r and center (x,y).
#' 
#' @param x The x coordinate of the center of the circle.
#' @param y The y coordinate of the center of the circle.
#' @param r The radius of the circle.
#' @param lty Line type option.  See R documentation on the lines() function for more details.
#' @param col Line color option.  See R documentation on the lines() function for more details.
#' @return The drawn circle.
plot.circle <- function(x,y,r,lty=1,col=1){
  a <- seq(0,2*pi,length.out=100)
  X <- r*cos(a)+x; Y <- r*sin(a)+y
  lines(X,Y,lty=lty,col=col)
}


#' Plots line intersects (useful for plotting recovery rate as a function of sensors).
#' 
#' @param n An x value.
#' @param rate A y value.
#' @param col Line color option.  See R documentation on the lines() function for more details.
#' @param lty Line type option.  See R documentation on the lines() function for more details.
#' @return The drawn lines of intersection.
plot.intersect <- function(n,rate,col=1,lty=1){
  lines(rep(n,2),c(0,rate),col=col,lty=lty)
  lines(c(0,n),rep(rate,2),col=col,lty=lty)
}
