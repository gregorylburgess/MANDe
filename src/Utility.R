#' @include src/ShapeFunctions.R
source('src/ShapeFunctions.R')




#' @name sensorFun
#' @title Calls functions to generate a 'goodness' grid and choose sensor locations.
#' @details Finds a "good" set of sensor placements for a given setup [topographyGrid, behaviorGrid, params].
#' Returns a list of locations as grid coordinates, and the calculated goodnessGrid (goodness grid).
#' Bias value controls the 'goodness' algorithm that gets called.
#' Bias cases:
#' 1. Fish density only.
#' 2. Visibility due to Bathy only.
#' 3. Detectable fish accounting for bathy.
#'
#' @param numSensors The number of sensors the program should place.
#' @param topographyGrid A valid topographyGrid.
#' @param behaviorGrid A valid behaviorGrid.
#' @param range The range of the sensor in bathymetric cells.
#' @param bias The goodness algorithm to use, choose 1, 2, or 3.  See above for descriptions.
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param save.inter If TRUE intermediary calculations are output as key inter.
#' @param silent If set to TRUE, disables status printing.
#' @param multi If set to TRUE, uses multicore to parallelize calculations.
#' @return A dictionary of return objects, see RETURN_DESCRIPTIONS.html for more info.
sensorFun = function(numSensors, topographyGrid, behaviorGrid, range, bias, params, debug=FALSE, silent = FALSE, save.inter=FALSE, multi=FALSE) {
    if (debug) {
        cat("\n[sensorFun]\n")
        print("topographyGrid")
        print(topographyGrid)
        print("behaviorGrid")
        print(behaviorGrid)
        print(sprintf("bias=%g",bias))
        print("params")
        print(params)
    }

    ## This is used to collect intermediary grids (warning: may use up lots of memory)
    if(save.inter) inter = list()
    
    sensorList = {}
    dims = dim(behaviorGrid)
    rows = dims[1]
    cols = dims[2]
    grids = list("topographyGrid" = topographyGrid, "behaviorGrid"=behaviorGrid)
    
    # calculate the goodnessGrid
    print("Calculating initial goodness grid")
    grids = goodnessGridFun(grids, range, bias, params, debug=debug, silent=silent, multi=multi)
    goodnessGrid = grids$goodnessGrid
    if(save.inter) inter[[1]] = grids
    ## place user-defined sensors, and down weigh them
    if("sensorList" %in% names(params) && length(params$sensorList) > 0) {
        i = length(params$sensorList)
        while(i > 0) {
            print(paste("Placing predefined sensor number",length(params$sensorList)-i+1),"of",length(params$sensorList))
            loc = params$sensorList[[i]]
            ## invert the incoming r/c values
            placement = list(c=loc$r, r=loc$c)
            grids = sensorFun.suppressHelper(placement, grids, range, bias, params, debug=debug, multi=multi)
            sensorList = c(sensorList, list(placement))
            i = i - 1
        }
    }
    
    ## for each sensor, find a good placement
    i = numSensors
    while(i > 0) {
        print(paste("Placing sensor number",numSensors-i+1,"of",numSensors))
        ## find the max location 
        maxLoc = which.max(grids$goodnessGrid)
        ## Switch the row/col vals since R references Grid coords as (y,x) instead of (x,y)
        c = ceiling(maxLoc/rows)
        r = (maxLoc %% rows)
        if (r==0) {
            r=rows
        }
        maxLoc = list(c=c,r=r)
        ##print(paste('Placed sensor',i,'at: ',maxLoc$c,maxLoc$c))
        
        ## append maxLoc to the sensor list.
        sensorList = c(sensorList, list(maxLoc))
        ## down-weigh all near-by cells to discourage them from being chosen by the program
        grids = sensorFun.suppressHelper(maxLoc, grids, range, bias, params, debug=debug, multi=multi)
        if(save.inter){
            ## Save intermediary grids
            inter[[numSensors-i+1+1]] = grids
        }
        i = i - 1
    }
    print("Done placing sensors")
    if(save.inter){
        return(list(sensorList=sensorList, goodnessGrid=goodnessGrid, goodnessGridSupp=grids$goodnessGrid, inter=inter))
    }else{
        return(list(sensorList=sensorList, goodnessGrid=goodnessGrid, goodnessGridSupp=grids$goodnessGrid))
    }
}

#' @title Performs suppression on the goodnessGrid/behaviorGrid depending upon the suppressionFcn specified.
#' @description The suppressionFcn 'detection.function.exact' indicates that the behaviorGrid must be suppressed, then the entire 
#' goodnessGrid recalculated based upon the new behaviorGrid.  Any other suppressionFcn simply calls suppress.opt, which just suppresses
#' the goodnessGrid, avoiding the need to recalculate the whole goodnessGrid (a very expensive operation).  The latter approach can be
#' considered a 'fast approximation', while the former an 'exact count'.
#'
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param grids A dictionary containing the keys 'topographyGrid', 'behaviorGrid', and 'goodnessGrid', which hold a valid topographyGrid, behaviorGrid and goodnessGrid.
#' @param range The range of the sensor in bathymetric cells.
#' @param bias The goodness algorithm to use, choose 1, 2, or 3.  See above for descriptions.
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param multi If set to TRUE, uses multicore to parallelize calculations.
#' @return Returns the grids parameter, with an updated behaviorGrid.
sensorFun.suppressHelper = function(loc, grids, range, bias, params, debug=FALSE, multi=FALSE) {
	if(params$suppressionFcn != 'detection.function.exact'){
		grids$goodnessGrid = suppress.opt(grids$goodnessGrid, dim(grids$behaviorGrid), loc, params, grids$topographyGrid$topographyGrid, debug)
	}else{
            print("Updating goodness grid")
            ## Downweigh observed region
            grids = updatebehaviorGrid(loc,grids,params,debug)
            grids = goodnessGridFun(grids, range, bias, params, debug=debug, multi=multi)
	}
	return(grids)
}

#' @title Updates the behaviorGrid after each sensor is placed to reflect which areas that are already covered by sensors.
#' @description When a sensor is placed the behaviorGrid must be updated to reflect where unique signals are emitted.
#' This is done by calculating the coverage of the sensor given by loc and then downweighing the behaviorGrid
#' using this coverage such that locations that are well covered by the sensor at loc is downweighed more
#' that poorly covered locations.
#'
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param grids A dictionary containing the keys 'topographyGrid', 'behaviorGrid', and 'goodnessGrid', which hold a valid topographyGrid, behaviorGrid and goodnessGrid.
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns the grids parameter, with an updated behaviorGrid.
updatebehaviorGrid = function(loc, grids, params, debug=FALSE){
  if(debug){
      cat("\n[updatebehaviorGrid]\n")
      print("loc")
      print(loc)
  }
  grid = grids$behaviorGrid
  bG = grids$topographyGrid$topographyGrid
  dims = dim(grid)
  rows = dim(grid)[1]
  cols = dim(grid)[2]
  ## get rows and column indices for relevant area
  vals = getArea(loc, dims, params$range) 

  rind = vals$rs:vals$re
  cind = vals$cs:vals$ce
  nrows = length(rind)
  ncols = length(cind)
  ## Make matrices with row and column indices to allow vectorized calculations
  Rind = matrix(rep(rind,ncols),nrows,ncols)
  Cind = matrix(rep(cind,nrows),nrows,ncols,byrow=TRUE)
  ## Calculate a matrix containing the distances from loc to all cells
  dist = sqrt( (loc$c-Cind)^2 + (loc$r-Rind)^2 )
  ## Calculate the detection function value at all grid points
  dgrid2 = do.call(params$shapeFcn, list(dist, params)) 
  ## Create a matrix where land cells have value TRUE
  land = bG >= 0
  ## Create the depth value of a sensor placed at loc
  sensorDepth = bG[loc$r,loc$c] + params$sensorElevation
  ## If dpflag is false then proportion of water column is calculated, if true depth preference is used
  dpflag = "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params
  ## Calculate the proportion of signals in each of the surrounding cell that can be detected by a sensor at loc
  pctviz = calc.percent.viz(loc$r, loc$c, rind, cind, bG, land, sensorDepth, dpflag, params, debug)
  ## testmap is a matrix with size as the full grid containing the percentage visibility of each cell
  ## Initialize
  testmap = matrix(0,rows,cols)
  ## Insert values at correct indices
  testmap[pctviz$linearIndex] = pctviz$percentVisibility
  ## 100% detected in self cell
  testmap[loc$r,loc$c] = 1
  ## Copy relevant area to dgrid1
  dgrid1 = testmap[rind,cind]
  ## dgrid contains downweighing values
  dgrid = 1 - (dgrid1 * dgrid2)
  ## Downweigh observed region
  grid[rind,cind] = grid[rind,cind] * dgrid
  grids$behaviorGrid = grid
  return(grids)
}


#' @name goodnessGridFun
#' @title Calculates the composite "goodness" grid for a particular bias.
#' @description Calls a particular goodnessGrid function based on the bias and opt values.  Actual work
#' 			is done by the called function.
#' @param bias The goodness algorithm to use, choose 1, 2, or 3.  See above for descriptions.
#' @param grids A dictionary containing the keys 'topographyGrid', 'behaviorGrid', and 'goodnessGrid', which hold a valid topographyGrid, behaviorGrid and goodnessGrid.
#' @param range The range of the sensor in bathymetric cells.
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param silent If set to TRUE, disables status printing.
#' @param multi If set to TRUE, uses multicore to parallelize calculations.
#' @return Returns the grids parameter, with an updated goodnessGrid.
goodnessGridFun = function (grids, range, bias, params, debug=FALSE, silent=FALSE, multi=FALSE) {
    if (debug) {
        cat("\n[goodnessGrid]\n")
        print("topographyGrid")
        print(grids$topographyGrid)
        print("behaviorGrid")
        print(grids$behaviorGrid)
        print(sprintf("bias=%g", bias))
        print("params")
        print(params)
    }
    status [toString(params$timestamp)] <<- 0
    topographyGrid = grids$topographyGrid$topographyGrid
    ## Remove all NAs from the Grids
    topographyGrid[is.na(topographyGrid)] = 0
    grids$topographyGrid$topographyGrid = topographyGrid
    grids$behaviorGrid[is.na(grids$behaviorGrid)] = 0
	
    ##Fish
    if (bias == 1) {
		if(debug) {
			print("bias=1; Calling goodnessGrid.sumSimple")
		}
        return(goodnessGrid.sumSimple.opt(grids, "behaviorGrid", params, debug, silent))
    }
    #Bathy
    else if (bias == 2) {
		if(debug) {
			print("bias=2; Calling goodnessGrid.sumBathy.opt")
		}
                if(multi) {
                    require(multicore)
                    return(goodnessGrid.sumBathy.multi(grids, params, debug, silent))
                }else{
                    return(goodnessGrid.sumBathy.opt(grids, params, debug, silent))
                }
    }
    #Combo
	else if (bias == 3) {
        ## Note goodnessGrid.sumBathy.opt also handles bias 3
		if(debug) {
			print("bias=3; Calling goodnessGrid.sumBathy.opt")
		}
                if(multi) {
                    require(multicore)
                    return(goodnessGrid.sumBathy.multi(grids, params, debug, silent))
                }else{
                    return(goodnessGrid.sumBathy.opt(grids, params, debug, silent))
                }
    }
    else {
        stop("ERROR: Invalid Bias")
    }   
}

#' @title Simply sums the values within range of a cell, for each cell in the given grid.
#' @description This is a speed optimized version, which uses the convolution operation (which
#' mainly gains it speed from using FFT) as an alternative to running through tedious
#' R for loops.
#'
#' @param grids A dictionary containing the keys 'topographyGrid', 'behaviorGrid', and 'goodnessGrid', which hold a valid topographyGrid, behaviorGrid and goodnessGrid.
#' @param key A key to the dictionary provided in the 'grids' parameter specifying which grid should be summed.
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param silent If set to TRUE, disables status printing.
#' @return Returns the grids parameter, with an updated goodnessGrid.
goodnessGrid.sumSimple.opt = function (grids, key, params, debug=FALSE, silent=FALSE) {
	
    ## Create a vector of distances to the cells that can be sensed by a sensor in the current cell
    subdists = 1:params$range
    ## Concatenate the "other" side and add self cell
    dists = c(rev(subdists),0,subdists)
    ## Calculate the detection function value at all distances and that is our kernel
    kernel = do.call(params$shapeFcn, list(dists, params))
    ## Check that the length of the kernel is as it should be
    if(length(kernel) != 2*params$range+1) {
		printError(paste('[goodnessGrid.sumSimple.opt]: length of kernel was:',
					length(kernel),'expected:',2*params$range+1))
	}
    ## Extract relevant grid as given by key
    tempGrid = get(key, grids)
    ## Do convolution. This operation is identical to the for loop in goodnessGrid.sumSimple
    ## For more general information about how the convolution operation is defined google it! wikipedia has a decent explanation.
	if(debug) {
		print("Calling conv.2D")
		print("tempGrid")
		print(tempGrid)
		print("Kernel")
		print(kernel)
	}
    grids$goodnessGrid = conv.2D(tempGrid,kernel,kernel, params$timestamp, silent)
    ## Calculate a matrix containing the depth of hypothetical sensors placed in each cell as an offset from the bottom
    sensorDepth = grids$topographyGrid$topographyGrid + params$sensorElevation
    ## Calculate a matrix where TRUE values indicate that a grid cell could contain a sensor below the surface
    belowSurf = sensorDepth < 0
    ## Set the goodness to zero in cells where a sensor would not be below the surface (it would stick out of the water)
    grids$goodnessGrid[!belowSurf] = 0
    
    if(debug){
        cat("\n[goodnessGrid.sumSimple.opt]\n")
        print("grids")
        print(grids)
    }
    return(grids)
}


#' @title Calculates the goodnessGrid when a line of sight bias is chosen (bias 2 or 3).
#' @description Loops through all cells where sensor placement is valid (where sensor would be below surface)
#' and calculates goodness. If bias is 2 only bathymetry (line of sight) is used to calculate goodness, whereas if
#' bias is 3 both bathymetry and fish distribution (behaviorGrid) are used. This function uses vectorized calculations.
#' 
#' @param grids A dictionary containing the keys 'topographyGrid', 'behaviorGrid', and 'goodnessGrid', which hold a valid topographyGrid, behaviorGrid and goodnessGrid.
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param silent If set to TRUE, turns off status printing.
#' @return Returns the grids parameter, with an updated goodnessGrid.
goodnessGrid.sumBathy.opt = function (grids, params, debug=FALSE, silent=FALSE) {
    nr = dim(grids$topographyGrid$topographyGrid)[1]
    nc = dim(grids$topographyGrid$topographyGrid)[2]
    ## Calculate the number of cells in the topographyGrid
    ng = nr*nc
    ## Make a copy of the topographyGrid to make code look nicer (could get rid of to save memory)
    bG = grids$topographyGrid$topographyGrid
    ## Initialize the goodnessGrid matrix (allocate memory)
    goodnessGrid = matrix(0,nr,nc)
    ## Make sure we use an integer range
    rng = round(params$range)
    ## Calculate a matrix containing the depth of hypothetical sensors placed in each cell as an offset from the bottom
    sensorDepth = bG + params$sensorElevation
    ## Calculate a matrix where TRUE values indicate that a grid cell could contain a sensor below the surface
    belowSurf = sensorDepth < 0
    ## Create a matrix where land cells have value TRUE
    land = bG >= 0
    ## If dpflag is false then proportion of water column is calculated, if true depth preference is used
    dpflag = "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params
    usebehaviorGrid = params$bias==3
    for(c in 1:nc){
        comp = c/nc
        if(!silent) {
            print(sprintf("LOS progress: %g", comp))
        }
        status[toString(params$timestamp)] <<- comp
        goodnessGrid[,c] <- goodness.multi.helper(c, nc, nr, rng, belowSurf, bG, land, sensorDepth, dpflag, params, usebehaviorGrid, grids, debug=debug)
        ## Column indices
        #cind = max(c(1,c-rng)):min(c(nc,c+rng))
        #for(r in 1:nr){
            ## Only calculate if sensor is below surface
            ## {{Patch}}
        #    cell = belowSurf[r,c]
        #    if(!is.na(cell) && cell){
                ## Row indices
        #        rind = max(c(1,r-rng)):min(c(nr,r+rng))
                ## Calculate the proportion of signals in each of the surrounding cell that can be detected by a sensor at (r,c)
				#pV has the following keys: percentVisibility, dists, linearIndex
        #        pV = calc.percent.viz(r, c, rind, cind, bG, land, sensorDepth[r,c], dpflag, params, debug)
                ## Calculate the detection function value at all grid points
        #        probOfRangeDetection = do.call(params$shapeFcn, list(pV$dists, params))
                ## If bias == 3 include the behaviorGrid in the calculations, if not just use bathymetry and detection function
	#			if(usebehaviorGrid) {
	#				probOfRangeDetection = probOfRangeDetection * grids$behaviorGrid[pV$linearIndex]
	#			}
                ## Calculate goodness of cell (r,c) by summing detection probabilities of all visible cells
        #        goodnessGrid[r,c] = sum(probOfRangeDetection * pV$percentVisibility)
        #    }
        #}
    }
	status[toString(params$timestamp)] <<- 1
    grids$goodnessGrid = goodnessGrid
    if(debug){
        cat("\n[goodnessGrid.sumBathy.opt]\n")
        print("grids")
        print(grids)
    }
    return(grids)
}


#' @title Calculates the goodnessGrid when a line of sight bias is chosen (bias 2 or 3) using multicore!
#' @description Loops through all cells where sensor placement is valid (where sensor would be below surface)
#' and calculates goodness. If bias is 2 only bathymetry (line of sight) is used to calculate goodness, whereas if
#' bias is 3 both bathymetry and fish distribution (behaviorGrid) are used. This function uses vectorized and parallel calculations using the multicore library.
#' 
#' @param grids A dictionary containing the keys 'topographyGrid', 'behaviorGrid', and 'goodnessGrid', which hold a valid topographyGrid, behaviorGrid and goodnessGrid.
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param silent If set to TRUE, turns off status printing.
#' @return Returns the grids parameter, with an updated goodnessGrid.
goodnessGrid.sumBathy.multi = function (grids, params, debug=FALSE, silent=FALSE) {
    ## Open fifo progress file
    progfile <- fifo(tempfile(), open="w+b", blocking=T)
    if (inherits(fork(), "masterProcess")) {
        ## Child
        progress <- 0.0
        while (progress < 1 && !isIncomplete(progfile)) {
            msg <- readBin(progfile, "double")
            progress <- progress + as.numeric(msg)
            cat(sprintf("Multi LOS progress: %.2f%%\r", progress * 100))
        } 
        exit()
    }
    nr = dim(grids$topographyGrid$topographyGrid)[1]
    nc = dim(grids$topographyGrid$topographyGrid)[2]
    ## Calculate the number of cells in the topographyGrid
    ng = nr*nc
    ## Make a copy of the topographyGrid to make code look nicer (could get rid of to save memory)
    bG = grids$topographyGrid$topographyGrid
    ## Initialize the goodnessGrid matrix (allocate memory)
    goodnessGrid = matrix(0,nr,nc)
    ## Make sure we use an integer range
    rng = round(params$range)
    ## Calculate a matrix containing the depth of hypothetical sensors placed in each cell as an offset from the bottom
    sensorDepth = bG + params$sensorElevation
    ## Calculate a matrix where TRUE values indicate that a grid cell could contain a sensor below the surface
    belowSurf = sensorDepth < 0
    ## Create a matrix where land cells have value TRUE
    land = bG >= 0
    ## If dpflag is false then proportion of water column is calculated, if true depth preference is used
    dpflag = "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params
    usebehaviorGrid = params$bias==3

    CS <- 1:nc
    ## mclapply is a function of the multicore package parallelizing the LOS calculation
    ## Need to implement a progress indicator, see: http://stackoverflow.com/questions/10984556/is-there-way-to-track-progress-on-a-mclapply
    goodnesslist <- mclapply(CS,goodness.multi.helper, nc, nr, rng, belowSurf, bG, land, sensorDepth, dpflag, params, usebehaviorGrid, grids, progfile, debug)
    ## goodnesslist is a list of vectors, which must be combined as a matrix
    goodnessGrid <- matrix(unlist(goodnesslist),nr,nc)
    status[toString(params$timestamp)] <<- 1
    grids$goodnessGrid = goodnessGrid
    if(debug){
        cat("\n[goodnessGrid.sumBathy.multi]\n")
        print("grids")
        print(grids)
    }
    ## Close progress file
    close(progfile)
    return(grids)
}



#' @title Helper function for multicore LOS calculation.
#' @description This is the function that goes into mclapply
#' 
#' @param r Row of the current cell in the topographyGrid.
#' @param c Column of the current cell in the topographyGrid.
#' @param rind Row indices of the topographyGrid to calculate visibility percentage.
#' @param cind Column indices of the topographyGrid to calculate visibility percentage.
#' @param topographyGrid A valid topographyGrid.
#' @param land Matrix containing logicals (TRUE = land cell) indicating wheter a cell 
#' in the topographyGrid is a land cell.
#' @param debug If enabled, turns on debug printing (console only).
#' @param sensorDepth Depth of sensor in current cell.
#' @param dpflag If TRUE depth preference is used meaning that the percentage of visible 
#' fish is calculated, if FALSE visible water column is calculated.
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @return Returns a dictionary with three keys (all vectors): percentVisibility contains
#' the percentage fish/signals visible in the surrounding cells, linearIndex contains the linear
#' indices in the topographyGrid to which the visibilities pertain, dists contains the distance
#' from the current cell to each of the returned cells as given by linearIndex.
goodness.multi.helper = function(c, nc, nr, rng, belowSurf, bG, land, sensorDepth, dpflag, params, usebehaviorGrid, grids, progfile=NA, debug=FALSE){
    goodnessVec = rep(0,nr)

    ## Column indices
    cind = max(c(1,c-rng)):min(c(nc,c+rng))
    for(r in 1:nr){
        ## Only calculate if sensor is below surface
        ## {{Patch}}
        cell = belowSurf[r,c]
        if(!is.na(cell) && cell){
            ## Row indices
            rind = max(c(1,r-rng)):min(c(nr,r+rng))
            ## Calculate the proportion of signals in each of the surrounding cell that can be detected by a sensor at (r,c)
            ##pV has the following keys: percentVisibility, dists, linearIndex
            pV = calc.percent.viz(r, c, rind, cind, bG, land, sensorDepth[r,c], dpflag, params, debug)
            ## Calculate the detection function value at all grid points
            probOfRangeDetection = do.call(params$shapeFcn, list(pV$dists, params))
            ## If bias == 3 include the behaviorGrid in the calculations, if not just use bathymetry and detection function
            if(usebehaviorGrid) {
                probOfRangeDetection = probOfRangeDetection * grids$behaviorGrid[pV$linearIndex]
            }
            ## Calculate goodness of cell (r,c) by summing detection probabilities of all visible cells
            goodnessVec[r] = sum(probOfRangeDetection * pV$percentVisibility)
        }
    }
    ## Send progress update
    if(class(progfile)[1]=='fifo') writeBin(1/nc, progfile)
    return(goodnessVec)
}


#' @name calc.percent.viz
#' @title Calculates the percentage of the water column in the surrounding cells that is visible to a sensor placed in the current cell.
#' @description Calculates a matrix centered around the current cell (r,c) and containing the percentage
#' of the water column in the surrounding cells that is visible to a sensor placed in the
#' current cell.
#'
#' 
#' @param r Row of the current cell in the topographyGrid.
#' @param c Column of the current cell in the topographyGrid.
#' @param rind Row indices of the topographyGrid to calculate visibility percentage.
#' @param cind Column indices of the topographyGrid to calculate visibility percentage.
#' @param topographyGrid A valid topographyGrid.
#' @param land Matrix containing logicals (TRUE = land cell) indicating wheter a cell 
#' in the topographyGrid is a land cell.
#' @param debug If enabled, turns on debug printing (console only).
#' @param sensorDepth Depth of sensor in current cell.
#' @param dpflag If TRUE depth preference is used meaning that the percentage of visible 
#' fish is calculated, if FALSE visible water column is calculated.
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @return Returns a dictionary with three keys (all vectors): percentVisibility contains
#' the percentage fish/signals visible in the surrounding cells, linearIndex contains the linear
#' indices in the topographyGrid to which the visibilities pertain, dists contains the distance
#' from the current cell to each of the returned cells as given by linearIndex.
calc.percent.viz = function(r, c, rind, cind, topographyGrid, land, sensorDepth, dpflag, params, debug=FALSE){
	if(debug){
		print("[calc.percent.viz]")
		print(sprintf("point: (%g,%g)",r,c))
	}
    rows = dim(topographyGrid)[1]
    cols = dim(topographyGrid)[2]
    nr = rows

    ## Find rows and columns for relevant cells
    rvec = rep(rind,length(cind))
    cvec = sort(rep(cind,length(rind)))
    ## Remove self cell
    tmp = which(!(rvec==r & cvec==c))
	#print(tmp)
    rvec = rvec[tmp]
    cvec = cvec[tmp]
	#print(rvec)
	#print(cvec)
    ## Translate from row col index to to single index
    linearIndex = sub2ind(rvec,cvec,nr) 
    linearIndexLength = length(linearIndex)
    ## Calculate distances from the current cell to the surrounding cells within range
    ## This sorts after dist so longest dists are calculated first, then shorter ones might not be needed since they are already calculated for a long dist
	## sort() returns two items, x and ix;  ix is the index of x in the sorted array.
    disttmp = sort(sqrt((r-rvec)^2 + (c-cvec)^2), decreasing=TRUE, index=TRUE)
	#print("Disttmp")
	#print(disttmp)
    ## Save actual distances in the dists vector
    dists = disttmp$x
    ## Get depths at the sorted cells by using disttmp$ix, which contains the sorted indices
    depths = topographyGrid[linearIndex[disttmp$ix]]
    ## Calculate the line of sight slopes to each of the sorted cells
    slopes = (depths-sensorDepth)/dists
    ## Create a vector, which can be used to easily map an index in the sorted vector to an index in the rind by cind matrix
    ibig2ismall = rep(0, rows*cols)
    ibig2ismall[linearIndex[disttmp$ix]] = 1:linearIndexLength

    ## Initialize vizDepths vector, this will be filled with visible depths below
    ## Assign small negative number to avoid problem with being exactly at the surface in pnorm
    vizDepths = rep(-1e-4,linearIndexLength)
    ## Initialize the remaining vector, which contains the indices of the cells for which the vizDepth has not yet been calculated
    remaining = 1:linearIndexLength

    ## Calculate visible depths as long as uncalculated cells remain
    while(length(remaining)>0){
        ## Since cells are sorted by decreasing distance taking the first index of remaining always gives the farthest uncalculated cell
        targetCell = remaining[1]
        ## Find the indices of the cells within the line of sight from r,c to targetCell, loslinearIndex are indices of big grid
        cells = getCells.opt(list(r=r,c=c),list(r=rvec[disttmp$ix[targetCell]],c=cvec[disttmp$ix[targetCell]]), debug, nr)
		if (all(is.na(cells))) {
			remaining = remaining[-1]
		}
        ## Get indices in small sorted vector (not whole grid)
		sortedIndexes = ibig2ismall[cells]
        ## Sort the distances within LOS
        d2 = sort(dists[sortedIndexes],index=TRUE)
        ## Find indices of obstacles (land areas) in LOS. If LOS is blocked by land don't calculate for cells behind.
        blocks = land[cells[d2$ix]]
        
        if(any(blocks, na.rm=TRUE)){
            if(!all(blocks)){
                ## Find indices that are not blocked
                linearIndexNoBlock = 1:(min(which(blocks), na.rm = TRUE)-1)
            }else{
                linearIndexNoBlock = NULL
            }
        }else{
            ## If theres no obstacles all cells in LOS are actually visible
            linearIndexNoBlock = 1:length(cells)
        }

        ## Vectorized calculation of the visible depths of the unblocked cells with LOS using the simple formula for a line y = ax + b
        ## a is cummax(slopes[is[d2$ix[linearIndexNoBlock]]])
        ## Here cummax ensures that the steepest slopes between the current cell (r,c) cell along the LOS is used for calculating visible depth, it is important that the slopes are sorted in order of increasing distance from current cell to target cells, this is handled by d2$ix
        ## x is d2$x[linearIndexNoBlock], the sorted distances from (r,c)
        ## b is sensorDepth, which is the intercept of the line
        vizDepths[sortedIndexes[d2$ix[linearIndexNoBlock]]] = cummax(slopes[sortedIndexes[d2$ix[linearIndexNoBlock]]])*d2$x[linearIndexNoBlock] + sensorDepth
        ## Remove the cells from remaining for which calculations are done
        remaining = setdiff(remaining, sortedIndexes)
    }
	
    ## Visible depths above water not valid (assign a number a little smaller than zero [just below surface])      
    vizDepths[vizDepths>0] = -1e-4
    ## Find indices that are not land cells
    linearIndexNotLand = depths<0
    ## if we have normal distribution data (depth preference), use it
    if(dpflag) {
        ## compute % fish visible from sensor to target cell
        ## calculate the mean fish depth in all cells
        mean = depths[linearIndexNotLand] + params$depth_off_bottom
        ## Values above water are set to be at the surface
        mean[mean>0] = 0
        ## Save SD in sd for prettier code
        sd = params$depth_off_bottom_sd
            
        ## Get cum probability for "available" water (between depth pref and surf)
        ## Cumulative probability below surface
        psurf = pnorm(0,mean=mean,sd=sd)
        ## Calculate the percentage of visible fish (the area under the visible part of the normal curve)
        percentVisibility = psurf - pnorm(vizDepths[linearIndexNotLand],mean=mean,sd=sd)
        ## Probability within available water (from bottom to surface), if this value is different from 1 it means that the vertical fish distribution extends below the bottom and/or above the surface, in which case we need to normalise the probability between bottom and surface so it sums to one.
        areaToCorrectFor = psurf - pnorm(depths[linearIndexNotLand],mean=mean,sd=sd)
        ## Correct for the distribution extending out of bounds
        percentVisibility = percentVisibility/areaToCorrectFor
    }else{
        ## if we don't have normal distribution data, assume equal distribution
        ## compute % visibility (of water column height) from sensor to target cell
        percentVisibility = vizDepths[linearIndexNotLand] / depths[linearIndexNotLand]
    }
    return(list(percentVisibility=percentVisibility, dists=dists[linearIndexNotLand], linearIndex=linearIndex[disttmp$ix[linearIndexNotLand]]))
}

#' @title Suppresses the values of cells around a sensor using a specified suppressionFunction.
#' @description This is an optimized version, which uses vectorization.
#' 
#' @param goodnessGrid A valid goodnessGrid.
#' @param dims The dimensions of the topographyGrid.  Just call dim() on the topographyGrid for this.
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param topographyGrid valid topographyGrid.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns a suppressed goodnessGrid.
suppress.opt = function(goodnessGrid, dims, loc, params, topographyGrid, debug=FALSE) {
    if(debug) {
        cat("\n[suppress.opt]\n")
        print(sprintf("suppressionFcn: %s", params$suppressionFcn))
        print(sprintf("loc: (%g,%g)",loc$c,loc$r))
        print("goodnessGrid")
        print(goodnessGrid)
        print("topographyGrid")
        print(topographyGrid)
    }
    suppressionFcn = params$suppressionFcn
    minsuppressionValue = params$minsuppressionValue
    maxsuppressionValue = params$maxsuppressionValue
    ## dfflag indicates whether detection function variant should be used for suppression
    dfflag = suppressionFcn=='detection.function' ||
			 suppressionFcn=='detection.function.shadow' ||
			 suppressionFcn=='detection.function.exact'
    rows = dim(goodnessGrid)[1]
    cols = dim(goodnessGrid)[2]
	
    ## Find the values that are affected by suppression (are within suppression range)
    vals = getArea(loc, dims, params$suppressionRange)

    ## Construct vectors containing row and col indices
    rind = vals$rs:vals$re
    cind = vals$cs:vals$ce
    nrows = length(rind)
    ncols = length(cind)
    ## Create matrices containing row and col indices and use for vectorized calculation of distances from loc to surrounding cells
    Rind = matrix(rep(rind,ncols),nrows,ncols)
    Cind = matrix(rep(cind,nrows),nrows,ncols,byrow=TRUE)
    dist = sqrt( (loc$c-Cind)^2 + (loc$r-Rind)^2 )

    if(suppressionFcn=='suppression.static'){
      supgrid = 1-matrix(maxsuppressionValue,nrows,ncols)
    }
    if(suppressionFcn=='suppression.scale'){
      sRange = minsuppressionValue - maxsuppressionValue
      supgrid = 1 - (sRange * (dist/params$suppressionRange) + maxsuppressionValue)
      supgrid[supgrid<0] = 0
      supgrid[supgrid>1] = 1
    }
    ## Use detection function variant
    if(dfflag){
      ## Save a temporary copy of params
      partmp = params
      ## Alter the SD value such that the shape function uses the suppression SD and not the detection SD since the shape function looks for the sd key
      partmp$sd = partmp$suppsd 
      ## Detection fun suppression
      supgrid2 = do.call(params$shapeFcn, list(dist, partmp))
      ## If shadowing should be accounted for in the suppression
      if(suppressionFcn=='detection.function.shadow' | suppressionFcn=='detection.function.exact'){
        ## Create a matrix where land cells have value TRUE
        land = topographyGrid >= 0
        ## Create the depth value of a sensor placed at loc
        sensorDepth = topographyGrid[loc$r,loc$c] + params$sensorElevation
        ## If dpflag is false then proportion of water column is calculated, if true depth preference is used
        dpflag = "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params
        ## Calculate the proportion of signals in each of the surrounding cell that can be detected by a sensor at loc
        pctviz = calc.percent.viz(loc$r, loc$c, rind, cind, topographyGrid, land, sensorDepth, dpflag, params, debug)
        ## testmap is a matrix with size as the full grid containing the percentage visibility of each cell
        ## Initialize
        testmap = matrix(0,rows,cols)
        ## Insert values at correct indices
        testmap[pctviz$linearIndex] = pctviz$percentVisibility
        ## 100% detected in self cell
        testmap[loc$r,loc$c] = 1
        ## Copy relevant area to dgrid1
        supgrid1 = testmap[rind,cind]
        ## supgrid contains suppression values
        supgrid = 1 - (supgrid1 * supgrid2)
      }else{
        ## supgrid contains suppression values
        supgrid = 1 - supgrid2
      }
    }
    ## Do suppression at the relevant indices given by rind and cind
    goodnessGrid[rind,cind] = goodnessGrid[rind,cind] * supgrid

    return(goodnessGrid)
}


#' @title Constant suppression.
#' @description Returns a static value, effectively setting all cells within suppressionRange of a sensor to that number.
#' 
#' @param dist The distance between two cells on a grid.
#' @param suppressionRange How far out to apply suppression penalties, in bathymetric cells.
#' @param minsuppressionValue The minimum allowable value to return.
#' @param maxsuppressionValue The maximum allowable value to return (also the return value for suppression.static()).
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns The value given in maxSuppression Value.
suppression.static = function (dist, suppressionRange, minsuppressionValue, 
                               maxsuppressionValue, params, debug=FALSE) {
    return (maxsuppressionValue)
}


#' @title Linear suppression.
#' @description Returns a dynamic value based on distance from a chosen sensor. The returned value should be multiplied by the value to be scaled.
#'
#' @param dist The distance between two cells on a grid.
#' @param suppressionRange How far out to apply suppression penalties, in bathymetric cells.
#' @param minsuppressionValue The minimum allowable value to return.
#' @param maxsuppressionValue The maximum allowable value to return (also the return value for suppression.static()).
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns The value given in maxSuppression Value.
suppression.scale = function (dist, suppressionRange, minsuppressionValue, 
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


#' @title Get rows and cols of area.
#' @description Defines the "shape" of a sensor's range.
#' 
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param dims The dimensions of the topographyGrid.  Just call dim() on the topographyGrid for this.
#' @param range The range of the sensor in bathymetric cells.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns a dictionary of start/end indexes for rows and columns respectively named : {rs,re,cs,ce}.
getArea=function(loc, dims, range, debug=FALSE) {
    ## the row index for our central point
    r = loc$r
    ## the col index for our central point
    c = loc$c
    ## the max number of rows in the grid
    rows = dims[1]
    ## the max number of cols in the grid
    cols = dims[2]
    
    ## defines a square
    rs0 = max(1, r - range) 
    re0 = min(rows, r + range)
    cs0 = max(1 ,c - range)
    ce0 = min(cols, c + range)
    toRet = list(rs=rs0, re=re0, cs=cs0, ce=ce0)
    return(toRet)
}

#' @title Returns the cells crossed by a beam from the starting cell to the target cell.
#' @description This is an optimized using vectorization. Returns the cells crossed by a beam from the starting cell to the target cell. Note that the starting cell is omitted from the result set.
#' 
#' @param startingCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor's location on the topographyGrid.
#' @param targetCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen tag's location on the topographyGrid.
#' @param debug If enabled, turns on debug printing (console only).
#' @param nr Number of rows in topographyGrid.
#' @return If nr != NULL a vector is returned containing the linear indices
#' (not row and col) of cells in the topographyGrid crossed by a beam from the starting
#' cell to the target cell. If nr == NULL a matrix with a row column and a col
#' column is returned.
getCells.opt = function(startingCell, targetCell, debug=FALSE, nr=NULL) {
	if(debug) {
		print("[getCells.opt]")
		print("Starting Cell:")
		print(startingCell)
		print("Target Cell:")
		print(targetCell)
		print(is.na(startingCell))
		print(is.na(targetCell))
		print(length(targetCell$c) == 0)
		print(length(targetCell$r) == 0)
	}
	if(is.na(startingCell) || is.na(targetCell) || 
			length(targetCell$c) == 0 || length(targetCell$r) == 0 || 
			length(startingCell$c) == 0 || length(startingCell$r) == 0) {
		print("something is null")
		return(NA)
	}
    if (!(startingCell$r == targetCell$r && startingCell$c == targetCell$c)) {
        ## Offset starting and target cells
        sC = offset(startingCell)
        tC = offset(targetCell)
        ## Define a value slightly larger than -1
        e = 1e-6 - 1
        ## Slope of a beam (line) from start to target cell
        ## Note: this is a slope in the horizontal not vertial plane
        ## rows act as y vals, cols act as x vals
        a = (tC$r-sC$r)/(tC$c-sC$c)
		if(abs(a)== Inf) {
			if(a>0) {
				a = 999999
			}
			else {
				a = -999999
			}
		}
		if(debug) print(sprintf("slope=%g",a))
        ## If absolute slope is less than 1
        if(abs(a)<=1){
            ## Intercept of beam
            b = sC$r - a*sC$c
            ## Column indices touched by beam (shift by e and then use ceiling to be sure to get the right column, this may unnecessary)
            cPoints = sort(startingCell$c:targetCell$c) + e
            cols = ceiling(cPoints)
            ## Find rows touched by beam
            rows = ceiling(a*cPoints + b)
            ## More than one row can be touched per column if -1 < a < 1
            if(a!=1){
                ## inds contains the indices where within one column two rows were touched 
                inds = which( abs(diff(rows))>0 )
                ## Find the rows and add them to the rows vector
                rows = c(rows,rows[inds+1])
                ## Add the corresponding columns also
                cols = c(cols,cols[inds])
            }
            ## Convert to linear (not row col) indices in the topographyGrid
            if(!is.null(nr)) {
				biginds = sub2ind(rows,cols,nr)
			}
        }else{
            ## If absolute slope is larger than 1 make rows act as x vals and cols as y vals
            ## Change slope accordingly
            a = 1/a
            ## Find intercept
            b = sC$c - a*sC$r
            ## Row indices touched by beam (shift by e and then use ceiling, this may unnecessary)
            rPoints = sort((startingCell$r):targetCell$r) + e
            rows = ceiling(rPoints)
            ## Find columns touched by beam
            cols = ceiling(a*rPoints + b)
            ## If within a row multiple columns were touched add them and their rows
            inds = which(abs(diff(cols))>0)
            rows = c(rows,rows[inds])
            cols = c(cols,cols[inds+1])
            ## Convert to linear (not row col) indices in the topographyGrid
            if(!is.null(nr)) biginds = sub2ind(rows,cols,nr)
        }
        ## Indices to use are all of them except the starting cell
        useinds = !(cols == startingCell$c & rows == startingCell$r)
        if(is.null(nr)){
          ## If nr is not input return a matrix with rows and cols
          return( cbind(cols[useinds],rows[useinds]) )
        }else{
          ## If nr was input return the linear indices in the topographyGrid
          return( biginds[useinds] )
        }
    }else{
        ## Return NULL if target cell is the starting cell
        return(NULL)
    }
}


#' @title Offset rows and cols to cell midpoint.
#' @description Translates a cartesian point to the center of the gridcell it represents.  This is necessary
#' for drawing Line of Sight since we assume sensors are in the center of cells.
#' ex: the cartesian point (3,2) would be converted to (2.5, 1.5), which puts it in the
#' cell located at the third column, second row (aka the cell at (3,2) on a 1-based grid
#' system).
#' 
#' @param point A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the point to translate on the topographyGrid.
#' @return A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the translated point.
offset= function(point){
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


#' @title Generate acoustic plots.
#' @description Generates .png files for the visualizations of various grids.  You must have write access to the 
#' R working directory that the program is executed from.  Additionally, ensure that an 'img' folder exists there.
#' 
#' @param result A dictionary of return objects, the result of a successfull call to run() or sensorFun().
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param showPlots If TRUE plots are shown on the screen, if FALSE plots are stored in the img folder.
#' @param plot.bathy Specifies whether contour lines for bathymetry should be overlayed in the graphs.
#' @return A dictionary containing the filenames of the generated images.
graph = function(result, params, showPlots, plot.bathy=TRUE) {
	time = "1"
	path = ""
        if(!showPlots) {
			if('timestamp' %in%  names(params)) {
				# Prevent R from using scientific notation (messes up filenames on windows)
				time = as.numeric(params$timestamp)
			}
			if(!file.exists(paste(path,"img", sep=""))) {
				dir.create(paste(path,"img", sep=""))
			}
			if(!file.exists(paste(path,"txt", sep=""))) {
				dir.create(paste(path,"txt", sep=""))
			}
			if(!file.exists(paste(path,"zip", sep=""))) {
				dir.create(paste(path,"zip", sep=""))
			}
		}
	## Plotting
	graphics.off()
	filenames = {}
        xlab = 'x dir'
        ylab = 'y dir'
	
	## topographyGrid
	if(!showPlots){
          filenames$topographyGrid = paste(path, "img/", time, "-TopographyGrid.png", sep="")
          png(filenames$topographyGrid)
        }else{
            dev.new()
        }
        plotGrid(result,type='topographyGrid',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	if(!showPlots) dev.off()
	
	## behaviorGrid
	if(!showPlots){
          filenames$behaviorGrid = paste(path, "img/", time, "-BehaviorGrid.png", sep="")
          png(filenames$behaviorGrid)
        }else{
            dev.new()
        }
        plotGrid(result,type='behaviorGrid',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	if(!showPlots) dev.off()
	
	## goodnessGrid
	if(!showPlots){
          filenames$goodnessGrid = paste(path, "img/", time, "-GoodnessGrid.png", sep="")
          png(filenames$goodnessGrid)
        }else{
            dev.new()
        }
        plotGrid(result,type='goodnessGrid',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	if(!showPlots) dev.off()
	
	## Acoustic Coverage
	if(!showPlots){
            filenames$coverageGrid = paste(path, "img/", time, "-CoverageGrid.png", sep="")
            png(filenames$coverageGrid)
        }else{
            dev.new()
        }
        plotGrid(result,type='coverageGrid',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	if(!showPlots) dev.off()

        ## Unique Recovery Rate
	if(!showPlots){
            filenames$recoveryRates = paste(path, "img/", time, "-RecoveryRates.png", sep="")
            png(filenames$recoveryRates)
        }else{
            dev.new()
        }
	print("[graph]")
	##print(result$stats)
        plotUniqueRR(result)
	if(!showPlots) dev.off()

	
	filenames = writeFiles(filenames, result, path, time, zip=TRUE)
	##print(filenames)
	return(filenames)
}


#' @title Adds txt and zip file paths to the filenames attribute, then writes the files.
#' @description Writes the txt dump of raw data used in the job, png files of the
#' graphs generated by the job, and optionally a zip file containing all of the above.  All
#' filenames include a timestamp value so that associated files are recognizable.
#' @param result A dictionary of return objects, the result of a successfull call to run() or sensorFun().
#' @param filenames A dictionary containing file paths; keys are: topographyGrid, behaviorGrid, goodnessGrid, coverageGrid, and
#' recoveryRates.
#' @param path A path to prepend to all the filenames.
#' @param time The timestamp label to include in all filenames.
#' @param zip If true, writes a zip file containing the txt dump, and images of a run.
#' @return A dictionary containing the filenames of the generated images.
writeFiles = function(filenames, result, path, time, zip=TRUE) {
	if (grepl("Rcheck",getwd())) {
		setwd("..")
	}
	print(getwd())	
	## Write results to a text file
	filename = paste(path, "txt/", time, "-Results.txt", sep="")
	jsonFile = paste(path, "txt/", time, "-Results.json", sep="")
        rdataFile = paste(path, "txt/", time, "-Results.RData", sep="")
        shortresFile = paste(path, "txt/", time, "-shortResults.txt", sep="")
	file.create(filename)
	capture.output(print(result), file=filename)
	filenames$txt = filename
	filenames$jsonFile = jsonFile
	filenames$rdataFile = rdataFile
	filenames$shortresFile = shortresFile
	# If true, write a zipped copy of files
	if (zip) {
		# Zip the text results and image files
		filename = paste(path, "zip/", time, ".zip", sep="")
		zip(zipfile=filename, files=filenames, flags="-r9X", extras="", zip=Sys.getenv("R_ZIPCMD", "zip"))
		filenames$zip = filename
	}
        ## Save in compressed R format (this should probably be save in a different folder, but will do for now)

        save('result',file=rdataFile)
        ## Save formatted text files with statistics

        cat(paste('- Acoustic network design results, generated:',Sys.time(),'\n\n'),file=shortresFile)
        if(is.null(result$errors)){
            cat(paste('Number of sensors placed:',params$numSensors,'\n'),file=shortresFile,append=TRUE)
            cat(paste('Detection range:',params$detectionRange,'meters\n'),file=shortresFile,append=TRUE)
            cat(paste('Using line of sight?:',!params$bias==1,'\n'),file=shortresFile,append=TRUE)
            cat(paste('Using fish behavior?:',!params$bias==2,'\n'),file=shortresFile,append=TRUE)
            cat(paste('Number of grid cells (row,col): ',params$XDist*params$YDist,' (',params$YDist,',',params$XDist,')\n',sep=''),file=shortresFile,append=TRUE)
            cat(paste('Suppression function:',params$suppressionFcn,'\n'),file=shortresFile,append=TRUE)
            cat(paste('Run time:',round(as.numeric(result$runTime, units='mins'),2),'mins\n'),file=shortresFile,append=TRUE)
            sensors <- result$stats$sensorMat[,c(2,1)]
            sensors <- cbind(sensors,sensors+matrix(c(rep(params$startX,params$numSensors),rep(params$startY,params$numSensors)),params$numSensors,2),round(result$stats$uniqRRs,3),round(c(result$stats$uniqRRs[1],diff(result$stats$uniqRRs)),3))
            colnames(sensors) <- c('loc_row','loc_col','glob_row','glob_col','Uniq_RR','d_Uniq_RR')
            rownames(sensors) <- paste('Sensor_',1:params$numSensors,sep='')
            cat(paste('\nOptimal sensor indices:\n'),file=shortresFile,append=TRUE)
            capture.output(sensors,file=shortresFile,append=TRUE)
            cat(paste('\nNetwork sparsity (delta):',round(result$stats$delta,3),'\n'),file=shortresFile,append=TRUE)
            cat(paste('Absolute recovery rate:',round(result$stats$absRecoveryRate,3),'\n'),file=shortresFile,append=TRUE)
            cat(paste('Unique recovery rate:',round(result$stats$uniqRecoveryRate,3),'\n'),file=shortresFile,append=TRUE)
            ##write.table(sensors,file=shortresFile,append=TRUE,sep='\t',row.names=FALSE)
        }else{
            cat(paste('Errors:',result$errors,'\n'),file=shortresFile,append=TRUE)
        }

	# Append the new txt and zip file locations, and clear all the old data from the result set.
	result$filenames = filenames
	result$topographyGrid = NULL
	result$behaviorGrid = NULL
	result$goodnessGrid = NULL
	result$coverageGrid = NULL
	return(filenames)
}


#' @title Plots the grid specified by the input type.
#' @description In addition to the grid itself sensor locations are also plotted along
#' with numbers indicating the order in which sensors were placed. Furthermore, bathymetry
#' contours can be overlayed using the plot.bathy flag.
#' 
#' @param result A dictionary of return objects, the result of a successfull call to run() or sensorFun().
#' @param type Character specifying grid type. Available grids: 'topographyGrid', 'behaviorGrid', 'goodnessGrid', or 'coverageGrid'.
#' @param main Set title of plot.
#' @param xlab Set label of x axis.
#' @param ylab Set label of y axis.
#' @param plot.bathy Specifies whether bathymetry contour lines should be added to plots.
#' @param plot.sensors Specifies if sensors should be added to plot.
#' @param ... Additional parameters to image, see ?image.
#' @return Nothing.
plotGrid = function(result,type='topographyGrid',main=type,xlab='',ylab='',plot.bathy=TRUE,plot.sensors=TRUE,...){
    ## n is number of colors in palette
    n = 24
    col = heat.colors(n)
    if(type=='topographyGrid'){
        ##col = colorRampPalette(c("darkviolet","navy","white"))(n)
        col = colorRampPalette(c("navy","white"))(n)
        grid = result$topographyGrid$topographyGrid
    }
    if(type=='behaviorGrid'){
        col = colorRampPalette(c("white","red", "yellow", "forestgreen"))(n)
        grid = result$behaviorGrid
    }
    if(type=='goodnessGrid'){
        col = colorRampPalette(c("white","red", "yellow", "forestgreen"))(n)
        grid = result$goodnessGrid
    }
    if(type=='coverageGrid'){
        col = colorRampPalette(c("white","red", "yellow", "forestgreen"))(n)
        grid = result$stats$coverageGrid
    }
    ## Plot the actual grid as an image
    image(result$topographyGrid$x,result$topographyGrid$y,grid,main=main,xlab=xlab,ylab=ylab,col=col,...)
    if(plot.bathy) {
        ## Add bathymetry contour
        contour(result$topographyGrid$x,result$topographyGrid$y,result$topographyGrid$topographyGrid,add=TRUE,nlevels=5)
    }
    if(plot.sensors){
        ## Add sensors and their numbers
        plotSensors(result)
    }
}


#' @title Plot unique recovery rate as a function of placed sensors.
#' @description Two-way plot showing cumulative unique recovery rate (top) and gained
#' unique recovery rate (bottom) per sensor as a function of number of sensors (in the
#' order they were placed). This is useful when assessing the value of adding more sensors
#' and when assessing the number of sensors required to obtain a certain recovery rate.
#' 
#' @param result A dictionary of return objects, the result of a successfull call to run() or sensorFun().
#' @return Nothing.
plotUniqueRR = function(result){
    ## Find number of placed sensors
    ns = length(result$sensors)
	##print(result$sensors)
	##print(ns)
    ## Find number of sensors for which unique recovery rate was calculated
    nsmax = length(result$stats$uniqRRs)
    ## Calculate the max value to use for the y-axis in TOP PLOT
    ## It looks good to show 1 as the max y val, but only if we are relatively ce (above 0.7)
    ymax = ifelse(max(result$stats$uniqRRs)>0.7,1.02,max(result$stats$uniqRRs)+0.2)
    ## Make two way plot
    par(mfrow=c(2,1),las=1)
    ## TOP PLOT
    plot(0:ns,c(0,result$stats$uniqRRs[1:ns]),typ='l',xlab='Number of sensors',ylab='Unique recovery rate',ylim=c(0,ymax),xlim=c(0,nsmax))
    points(0:ns,c(0,result$stats$uniqRRs[1:ns]),pch=46,cex=3)
    lines(ns:length(result$stats$uniqRRs),result$stats$uniqRRs[ns:nsmax],lty=2)
    plotIntersect(ns,result$stats$uniqRecoveryRate,col='orange',lty=1)
    grid()
    text(0.05*length(result$stats$uniqRRs),result$stats$uniqRecoveryRate,round(result$stats$uniqRecoveryRate,digits=4),pos=3)
    ##legend('bottomright',c('Calculated','Requested','Projected'),lty=c(1,1,2),col=c(1,'orange',1),bg='white')
    duRR = diff(c(0,result$stats$uniqRRs))
    ## BOTTOM PLOT
    plot(1:ns,duRR[1:ns],typ='l',xlab='Number of sensors',ylab='Increase in unique RR',ylim=c(0,max(duRR)),xlim=c(0,nsmax))
    points(1:ns,duRR[1:ns],pch=46,cex=3)
    lines(ns:nsmax,duRR[ns:nsmax],lty=2)
    grid()
}


#' @title Adds sensors to current plotGrid plot.
#' @description This function assumes a plot created by plotGrid exists.
#' 
#' @param result A dictionary of return objects, the result of a successfull call to run() or sensorFun().
#' @param circles If TRUE circles with radius equal to the detection range are drawn around sensors.
#' @param circlty Line type for circles.
#' @return Nothing.
plotSensors = function(result,circles=TRUE,circlty=3){
  ns = dim(result$stats$sensorMat)[1]
  ##print("[plotSensors]")
  ##print("result$stats")
  ##print(result$stats)
  #if (is.null(ns)) {
#	  ns = matrix(result$stats$sensorMat,nrow=length(result$stats$sensorMat)/2,ncol=2)
 # }
  ##print(length(result$stats$sensorMat))
  ##print(result$stats$sensorMat)
  ##print(dim(result$sensorMat))
  ##print("ns")
  ##print(ns)
  ## Radius of circle
  r = result$params$detectionRange
  ## Radian values for a full circle
  a = seq(0, 2 * pi, length.out=100)
  ## Cols
  sensx = result$topographyGrid$x[result$stats$sensorMat[1:ns, 2]]
  ## Rows
  sensy = result$topographyGrid$y[result$stats$sensorMat[1:ns, 1]]
  ## Plot sensor range as circles
  if(circles){
	i = ns
    while(i > 0){
      X = r*cos(a) + sensx[ns-i+1]
      Y = r*sin(a) + sensy[ns-i+1]
      lines(X,Y,lty=circlty)
	  i = i - 1
    }
  }
  ## Plot sensor locations
  points(sensx,sensy,pch=21,bg='blue',cex=3)
  text(sensx,sensy,1:ns,col='white')
}


#' @title Provides statistical data on detection, given a particular topographyGrid, behaviorGrid, and sensor arrangement.
#' @description This function calculates placements of projected sensors and the value
#' of each sensor (requested and projected) as given in increase in unique recovery rate.
#' Additionally, the acoustic coverage map, unique recovery rate, absolute recovery rate, and sparsity
#' are calculated and returned after placing the requested number of sensors.
#' 
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param topographyGrid A valid topographyGrid.
#' @param behaviorGrid A valid behaviorGrid.
#' @param sensors The result of a successful call to sensorFun().
#' @param debug If enabled, turns on debug printing (console only).
#' @return A dictionary (list) of statistical values containing the keys: "delta", "sensorMat"         
#' "uniqRRs", "coverageGrid", "absRecoveryRate", "uniqRecoveryRate".
getStats = function(params, topographyGrid, behaviorGrid, sensors, debug=FALSE) {
    if(debug) {
        print("[getstats]")
    }
    statDict = list()
    ## Number of requested sensors
    numSensors = length(sensors$sensorList)
    rows = dim(behaviorGrid)[1]
    cols = dim(behaviorGrid)[2]
    rng = params$range
    
    ## Calculate the value of each sensor (as the increase in unique recoveryrate)
    numProj = numSensors + params$projectedSensors
    ## goodnessGrid suppressed by numSensors
    goodnessGridSupp = sensors$goodnessGridSupp
    sensorList = sensors$sensorList
    ## Calculate locations of projected sensors
    ns = numProj-numSensors
    i = ns
    while (i > 0) { 
        ## find the max location 
        maxLoc = which.max(goodnessGridSupp)
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
        goodnessGridSupp = suppress.opt(goodnessGridSupp, dim(behaviorGrid), maxLoc, params, topographyGrid$topographyGrid, debug)
        i = i - 1
    }

    xSens = rep(0,numProj)
    ySens = rep(0,numProj)
    i = numProj
    while(i > 0){
        k = numProj-i+1
        xSens[k] = sensorList[[k]]$c
        ySens[k] = sensorList[[k]]$r
        i = i - 1
    }
    
    ## Distance maps (the distance from any grid cell to a receiver)
    rows = dim(behaviorGrid)[1]
    cols = dim(behaviorGrid)[2]
    X = matrix(rep(1:cols,rows),rows,cols,byrow=TRUE)
    Y = matrix(rep(1:rows,cols),rows,cols,byrow=FALSE)
    dimap = array(0,dim=c(rows,cols,numProj))
    i = numProj
    while(i > 0) {
        k = numProj-i+1
        ## Distance to receiver
        dimap[,,k] = sqrt( (X-xSens[k])^2 + (Y-ySens[k])^2 )
        i = i - 1
    }
    
    ## Horizontal detection maps using detection function
    demap = array(0,dim=c(rows,cols,numProj))
    i = numProj
    while(i > 0) {
        k = numProj-i+1
        demap[,,k] = do.call(params$shapeFcn, list(dimap[,,k], params))
        i = i - 1
    } 

    ## Incorporate vertical detection probability using line of sight
    if(params$bias!=1){
        bG = topographyGrid$topographyGrid
        ## Create a matrix where land cells have value TRUE
        land = bG >= 0
        ## Calculate a matrix containing the depth values of a sensor placed in each grid cell
        sensorDepth = bG + params$sensorElevation
        ## If dpflag is false then proportion of water column is calculated, if true depth preference is used
        dpflag = "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params
        i = numProj
        while(i > 0) {
	    k = numProj-i+1
            r = ySens[k]
            c = xSens[k]
            cind = max(c(1,c-rng)):min(c(cols,c+rng))
            rind = max(c(1,r-rng)):min(c(rows,r+rng))
            ## Calculate the proportion of signals in each of the surrounding cell that can be detected by a sensor at loc
            pctviz = calc.percent.viz(ySens[k], xSens[k], rind, cind, bG, land, sensorDepth[ySens[k], xSens[k]], dpflag, params, debug)
            ## testmap is a matrix with size as the full grid containing the percentage visibility of each cell
            ## Initialize
            testmap = matrix(0,rows,cols)
            ## Insert values at correct indices
            testmap[pctviz$linearIndex] = pctviz$percentVisibility
            ## 100% detected in self cell
            testmap[r,c] = 1
            ## Update demap (detection map)
            demap[,,k] = demap[,,k] * testmap
            i = i - 1
        }
    }

    ## Coverage map
    cover = matrix(1,rows,cols)
    uniqRRs = rep(0,numProj)
	i = numProj
	while(i > 0) {
		k = numProj-i+1
        r = ySens[k]
        c = xSens[k]
        cind = max(c(1,c-rng)):min(c(cols,c+rng))
        rind = max(c(1,r-rng)):min(c(rows,r+rng))
        ## Probability of no detection
        cover[rind,cind] = cover[rind,cind] * (1 - demap[rind,cind,k]) 
        covertmp = 1-cover
        uniqRRs[k] = sum(covertmp * behaviorGrid)
		i = i - 1
    }
    duniqRRs = diff(c(0,uniqRRs))
    ##print(uniqRRs)
    ##print(duniqRRs)
    ## Sort list so best sensors come first
    srt = sort(duniqRRs,index=TRUE,decreasing=TRUE) 
    sensorMat = matrix(unlist(sensorList),numProj,2,byrow=TRUE)
    ## Store the unsorted sensor list (this order may not have a decreasing value per sensor)
    ## That is, the sensor placed as say 10 may be better than sensor 9 (this can happen when
    ## not using detection.function.exact as suppression)
    statDict$sensorMatNotSorted = sensorMat
    ## Order sensor such that best are placed first
    statDict$sensorMat = sensorMat[srt$ix,]
    ## Re-calculate unique recovery rate (OBS: this should really be re-calculated by iteratively
    ## placing sensors and calculating unique RR as above)
    statDict$uniqRRs = cumsum(srt$x)
    ## Acoustic coverage for best sensors
    statDict$coverageGrid = 1-apply(1-demap[,,srt$ix[1:numSensors]],c(1,2),prod)

    ## Calculate distance matrix needed to calculate sparsity
    distVec = rep(0,numSensors)
    i = numProj
    while(i > 0) {
        k = numProj-i+1
        dists = sqrt((topographyGrid$x[xSens[srt$ix[k]]]-topographyGrid$x[xSens[srt$ix[1:numSensors]]])^2 + (topographyGrid$y[ySens[srt$ix[k]]]-topographyGrid$y[ySens[srt$ix[1:numSensors]]])^2)
        distVec[k] = min(dists[dists>0])
        i = i - 1
    }
    ## a is the median of the distances between the receivers
    if(debug){
        print("distVec")
        print(distVec)
    }
    a = median(distVec)
    ## delta is a sparsity measure (see Pedersen & Weng 2013)
    statDict$delta = a/(2*params$detectionRange)
    
    ## Absolute recovery rate (here we don't care about getting the same ping multiple times)
    demapmat = apply(demap[,,srt$ix[1:numSensors]],c(1,2),sum)
    statDict$absRecoveryRate = sum(demapmat*behaviorGrid)
	
    ## Calculate recovery rate of unique detections    
    statDict$uniqRecoveryRate = sum(statDict$coverageGrid * behaviorGrid)
    
    return(statDict)
}


#' @title Check model parameters.
#' @description Provides default parameter values if none are provided.
#' 
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param stop If TRUE, stops the program when an error occurs.
#' @return The 'params' parameter, populated with default values where necessary.
checkParams = function(params, stop=TRUE) {
    names = names(params)
	## Cast all possible strings to numbers (JSON makes everything strings)
	## Additionally, check for NA, NaN, +/-INF values.
	for (name in names) {
		if(!is.na(suppressWarnings(as.numeric(params[name])))) {
			if(!is.finite(as.numeric(params[name]))) {
				printError(sprintf("Found NA, NaN, or +/-INF values in %s.", name), stop=stop)
			}
			params[name] = as.numeric(params[name])
		}
	}
	# timestamp
	if(!('timestamp' %in% names)) {
		params$timestamp = -1
	}
	
	# numSensors
    if(!('numSensors' %in% names)) {
        params$numSensors = 0
    }
	else {
		checkForMin("numSensors", as.numeric(params$numSensors), 0, stop)
	}
	
	# projectedSensors
	if(!('projectedSensors' %in% names)) {
		params$projectedSensors = 0
	}
	else {
		checkForMin("projectedSensors", as.numeric(params$projectedSensors), 0, stop)
	}
	
	#userSensorList
	if('userSensorList' %in% names) {
		# if params$userSensorList exists, clean it and parse it into params$sensorList
		cleaned = gsub("\\s", "", params$userSensorList)
		rawPointList = strsplit(cleaned, ",")[[1]]
		points = {}
		len = floor(length(rawPointList)/2)
		i=len
		while(i > 0) {
			r = as.numeric(rawPointList[2])
			c = as.numeric(rawPointList[1])
			
			if (!(is.finite(r) && is.finite(c))) {
				printError("A user-defined point containing NaN, NA, or Inf was found.", stop)
			}
			if (r <= 0 || c <= 0) {
				printError("A user-defined point is out of bounds.", stop)
			}
			point = list(r=r,c=c )
			points = c(points, list(point))
			rawPointList = rawPointList[-2]
			rawPointList = rawPointList[-1]
			i = i - 1
		}
		params$sensorList = points
		
		# if not enough sensors were specified, throw an error!
		if((params$numSensors + length(points) + params$projectedSensors  <= 1)) {
			printError("Please specify/allow the program to place/project a total of at least two sensors.", stop)
		}
	}
	# if not enough sensors were specified, throw an error!
	else if((params$numSensors + params$projectedSensors  <= 1)) {
		printError("Please specify/allow the program to place/project a total of at least two sensors.", stop)
	}
	
	# bias
    if(!('bias' %in% names)) {
        params$bias = 1
    }
	if(!(params$bias %in% c(1,2,3))) {
		printError("Error: Bias value must be 1, 2, or 3.", stop)
	}
	
	# Warn users if they use depth pref without knowing what the map looks like
	if('depth_off_bottom' %in% names  && !('inputFile' %in% names)) {
		printError("Error: Using dp option without a known input file may be bad!.
				For example, if the generated habitat grid contains no cells near
				the depth specified, no fish will be generated.", stop)
	}
	
    # inputFile
    if(('inputFile' %in% names)) {
        params$inputFile = as.character(params$inputFile)
    }
	else {
		# default to the 1km grid
		params$inputFile = "src/himbsyn.bathytopo.1km.v19.grd/himbsyn.bathytopo.1km.v19.grd"
	}
	
	# inputFileType
    if(!('inputFileType' %in% names)) {
		params$inputFileType = 'ncdf'
	}
	else {
		supportedFileTypes = c("netcdf", "arcgis", "asc")
		if(!(params$inputFileType %in% supportedFileTypes)) {
			printError("Invalid 'inputFileType' value.", stop)
		}
        params$inputFileType = as.character(params$inputFileType)
    }
	
	# cellSize
    if(!('cellSize' %in% names)) {
        params$cellSize = 1000
    }
	else {
		checkForMin("cellSize", as.numeric(params$cellSize), 1, stop)
	}
	
	# detectionRange
	if(!('detectionRange' %in% names)) {
		params$detectionRange = 2000
	}
	else {
		checkForMin("detectionRange", as.numeric(params$detectionRange), 1, stop)
	}
	if(params$detectionRange <= params$cellSize) {
		printError("Detection Range of a sensor must greater than a cell's width.", stop)
	}
	
	# startX
    if(!('startX' %in% names)) {
        params$startX = 308
    }
	else {
		checkForMin("startX", as.numeric(params$startX), 1, stop)
	}
	
	# startY
    if(!('startY' %in% names)) {
        params$startY = 452
    }
	else {
		checkForMin("startY", as.numeric(params$startY), 1, stop)
	}
	
	#XDist
    if(!('XDist' %in% names)) {
        params$XDist = 50
    }
	else {
		checkForMin("XDist", as.numeric(params$XDist), 1, stop)
	}
	
	#YDist
    if(!('YDist' %in% names)) {
        params$YDist = 50
    }
	else {
		checkForMin("YDist", as.numeric(params$YDist), 1, stop)
	}
	
	#seriesName
    if(!('seriesName' %in% names)) {
        params$seriesName = 'z'
    }	else {
		params$seriesName = as.character(params$seriesName)
	}
	

    # shapeFcn
    if(!('shapeFcn' %in% names)) {
        params$shapeFcn = "shape.gauss"
    }	
	else {
		# Currently, only gauss is defined.
		if (params$shapeFcn != "shape.gauss") {
			printError("Currently, the only valid value for shapeFcn is 'shape.gauss'.", stop)
		}
	}
	if(!('peak' %in% names)) {
		peak = .98
	}
	else {
		checkForMin("peak", as.numeric(params$peak), 0.00001, stop)
		checkForMax("peak", as.numeric(params$peak), 1, stop)
	}
	
	#sensorElevation
    if(!('sensorElevation' %in% names)){
        params$sensorElevation = 1
    }   
	else {
		checkForMin("sensorElevation", as.numeric(params$sensorElevation), 0, stop)
    }
    
    ## Calculate an approximate upper bound for suppression range
    ## First calculate the area of the study region
    studyArea = params$XDist * params$YDist * params$cellSize^2
    ## Then find the available area per sensor
    areaPerSensor = studyArea/params$numSensors
    ## Find the radius of a circle covering this area
    ## (then if sensors had this range the total area covered is equal to the
    ## area of the study region.) 
    requiredDetectionRange = sqrt(areaPerSensor/pi)
    ## However, since a suppressionRangeFactor of 2 is the default value use this
    ## as the minimum max value. Also multiply requiredDetectionRange by 2 because
    ## suppression resulting in non-overlapping detection zones should have a
    ## suppressionRangeFactor, which is twice the detectionRange. Using 2 may, however,
    ## result in too large suppression so a trade-off is 1.5. (we could think more about this)
    maxSuppressionRangeFactor = max(c(2,1.5*requiredDetectionRange/params$detectionRange))
    
    # SuppressionFcn
    if(!('suppressionFcn' %in% names)) {
        params$suppressionFcn = "detection.function"
    }	
	else {
		validFcns = c("suppression.static", "suppression.scale", "detection.function",
					  "detection.function.shadow", "detection.function.exact")
		if (!(as.character(params$suppressionFcn) %in% validFcns)) {
			printError("Invalid 'suppressionFcn' value.", stop)
		}
		params$suppressionFcn = as.character(params$suppressionFcn)
	}
	# 0 < max/min suppressionValue < 1
	if(!('maxsuppressionValue' %in% names)) {
		params$maxsuppressionValue = .5
	}
	else {
		checkForMin("maxsuppressionValue", params$maxsuppressionValue, 0, stop)
		checkForMax("maxsuppressionValue", params$maxsuppressionValue, 1, stop)
	}
	if(!('minsuppressionValue' %in% names)) {
		params$minsuppressionValue = .1
	}
	else {
		checkForMin("minsuppressionValue", params$minsuppressionValue, 0, stop)
		checkForMax("minsuppressionValue", params$minsuppressionValue, 1, stop)
	}
	# minsuppressionValue < maxsuppressionValue
	if(!(params$minsuppressionValue < params$maxsuppressionValue)) {
		printError("'minsuppressionValue' must be greater than 'maxsuppressionValue'.")
	}
	
	# SuppressionRange Factor
    if(!('suppressionRangeFactor' %in% names)) {
        params$suppressionRangeFactor = 2
    }	
	else {
		checkForMin("suppressionRangeFactor", as.numeric(params$suppressionRangeFactor), 0, stop)
		params$suppressionRangeFactor = as.numeric(params$suppressionRangeFactor)
	}
    if(params$suppressionRangeFactor > maxSuppressionRangeFactor){
        params$suppressionRangeFactor = maxSuppressionRangeFactor
    }
    
    # FishModel
    if(!('fishmodel' %in% names)) {
        params$fishmodel = 'rw'
    }	
	else {
		params$fishmodel = as.character(params$fishmodel)
		if(params$fishmodel == "True" | params$fishmodel == "ou") {
			params$fishmodel = 'ou'
			#OU vals
			required = c('mux','muy','ousdx','ousdy','oucor')
			if(!all(required %in% names)) {
				printError("Missing OU model variables.  Required values are 'mux', 'muy', 'ousdx', 'ousdy', and 'oucor'")
			}
			#validation for ou vars
			else {
				#mux, muy must be porportions of the x/y axis, between 0 and 1)
				checkForMin('mux', params$mux, 0, stop)
				checkForMin('muy', params$muy, 0, stop)
				checkForMax('mux', params$mux, 1, stop)
				checkForMax('muy', params$muy, 1, stop)
				
				#ousdx and ousdy need only be non negative
				checkForMin('ousdx', params$ousdx, 0, stop)
				checkForMin('ousdy', params$ousdy, 0, stop)
				
				# -1 < oucor < 1
				checkForMin('oucor', params$oucor, -1, stop)
				checkForMax('oucor', params$oucor, 1, stop)
			}
		}
		else if(params$fishmodel == "False" | params$fishmodel == "rw") {
			params$fishmodel = 'rw'
		}
		else {
			printError("Invalid 'fishmodel' value.", stop)
		}
	}
	
	#Vertical Habitat Range
	if('mindepth' %in% names && 'maxdepth' %in% names) {
		#mindepth must be non-positive
		checkForMax('mindepth', params$mindepth, 0, stop)
	}
	
	#Depth Preference
	if('depth_off_bottom' %in% names && 'depth_off_bottom_sd' %in% names) {
		#must be non-negative
		checkForMin('depth_off_bottom', params$depth_off_bottom, 0, stop)
		checkForMin('depth_off_bottom_sd', params$depth_off_bottom_sd, 0, stop)
	}
	
    return(params)
}

#' @name checkForMin
#' @title Checks that value is greater than the provided minVal, throwing an error if it isn't.
#' @description Checks taht value is greater than minVal, printing an error with the provided name
#' if it isn't.
#' 
#' @param name A text name for the value variable (this is what gets printed in the error message.
#' @param value The value to test.
#' @param minVal The minimum value for value.
#' @param stop If TRUE, stops the program when an error occurs.
#' @return Nothing.
checkForMin = function(name, value, minVal, stop=TRUE){
	if(value < minVal ){
		printError(sprintf("'%s' value must be at least %g, recieved %g.", name, minVal, value), stop)
	}
}

#' @name checkForMax
#' @title Checks that value is less than the provided maxVal, throwing an error if it isn't.
#' @description Checks that value is less than maxVal, printing an error with the provided name
#' if it isn't.
#' 
#' @param name A text name for the value variable (this is what gets printed in the error message.
#' @param value The value to test.
#' @param maxVal The max value for value.
#' @param stop If TRUE, stops the program when an error occurs.
#' @return Nothing.
checkForMax = function(name, value, maxVal, stop=TRUE){
	if(value > maxVal ){
		printError(sprintf("'%s' value must be less than %g, recieved %g.", name, maxVal, value), stop)
	}
}

#' @title Converts input parameters from meters to grid cells.
#' @description These are used for internal calculations and are invisible to the user.
#' 
#' @param params A dictionary of parameters, see ?acousticRun for more info.
#' @param topographyGrid A valid topographyGrid.
#' @return The 'params' parameter, populated with necessary internal variables ("sd",
#' "suppsd", "range", "suppressionRange") with unit 'grid cells'.
convertMetersToGrid = function(params, topographyGrid=NA){
  ## Cell size in meters
  cellSize = params$cellSize
  ## Detection range with SD=1, dx=1
  detectrng1 = abs(qnorm(0.05 / 2 /params$peak))
  ## SD in grid cells
  params$sd = params$detectionRange / detectrng1 / cellSize
  ## SD used for suppression in grid cells
  params$suppsd = params$suppressionRangeFactor * params$sd
  ## Range in grid cells, used to cut out sub grids from large grid
  params$range = round(3 * params$sd)
  ## Using equation 8 in Pedersen & Weng 2013
  params$suppressionRange = round(params$suppressionRangeFactor * params$detectionRange / cellSize) 
  return(params)
}

#' @title Performs the convolution operation in 1D.
#' 
#' @param fun Vector containing values to perform convolution operation on.
#' @param kern One-dimensional convolution kernel given as a vector preferably with odd-numbered length.
#' @return A vector containing the result of the convolution operation with same length as fun.
conv.1D = function(fun, kern){
  ## Length of kernel
  lk = length(kern)
  ## Prepare kern for input to R's built-in convolve
  kern = kern[lk:1]
  ## Do convolution, because of the "open" type the returned vector is longer than fun
  test = convolve(fun,kern,type='open')
  ## Calculate variable needed when extracting the relevant information from test
  lk2 = 0.5 * (lk - 1)
  lf = length(fun)
  ## Extract and return relevant information from test
  return(test[(lk2 + 1):(lf + lk - lk2 - 1)])
}


#' @title Performs the convolution operation in 2D.
#' @description This is a simple generalization of conv.1D to two dimensions in the sense that
#' it performs two 1D convolutions using two 1D kernels (which can be different). This
#' function is not able to do an actual 2D convolution using a convolution matrix (this functionality
#' is not needed here anyway).
#' 
#' @param mat Matrix containing values to perform convolution operation on
#' @param kx Convolution kernel in x direction
#' @param ky Convolution kernel in y direction
#' @param timestamp A timestamp reference for status updates.
#' @param silent If set to TRUE, disables status printing.
#' @return A matrix containing the result of the convolution operation with same dimensions as mat.
conv.2D = function(mat, kx, ky, timestamp=1, silent=FALSE){
	
  dimmat = dim(mat)
  ## Initialize return matrix
  matout = matrix(0, dimmat[1], dimmat[2])
  x = dimmat[1]
  y = dimmat[2]
  ## Convolve in x-direction
  for(i in 1:x) {
	  matout[i,] = conv.1D(mat[i,],kx)
	  pct = (i/(2*x))
	  status[toString(timestamp)] <<- pct
	  if(!silent) {
		  print(sprintf("LoS Progress:%g", pct))
	  }
  }
  ## Convolve in y-direction (important to use matout here)
  for(i in 1:y) {
	  matout[,i] = conv.1D(matout[,i],ky)
	  pct = (i/(2*x) + .5)
	  status[toString(timestamp)] <<- pct
	  if(!silent) {
		  print(sprintf("LoS Progress:%g", pct))
	  }
  }
  if(!silent) {
	  print(sprintf("LoS Progress:%g", 1))
  }
  status[toString(timestamp)] <<- 1
  return(matout)
}


#' @title Converts a row, col index to a linear index within a matrix.
#' 
#' @param row The row index.
#' @param col The col index.
#' @param dims The dimensions of the topographyGrid.  Just call dim() on the parent matrix for this.
#' @return The translated linear index.
sub2ind = function(row, col, dims){
    (col-1) * dims[1] + row
}


#' @title Plots line intersects on an existing plot, which must exist.
#' @description A horizontal line from zero to rate and a vertial line
#' from zero to n are drawn. This function is useful for highlighting an
#' important value when plotting recovery rate as a function of sensors.
#' 
#' @param n An x value.
#' @param rate A y value.
#' @param col Line color option.  See R documentation on the lines() function for more details.
#' @param lty Line type option.  See R documentation on the lines() function for more details.
#' @return No returned value, lines are drawn on existing plot.
plotIntersect = function(n, rate, col=1, lty=1){
  lines(rep(n, 2), c(0, rate), col=col, lty=lty)
  lines(c(0, n), rep(rate, 2), col=col, lty=lty)
}

#' @title Prints errors.
#'
#' @param msg The message to print
#' @param stop If TRUE, stops the program when an error occurs.
#' @return The passed message.
printError = function(msg, stop=TRUE) {
	print(msg)
	traceback()
	if(stop) {
		stop(msg)
	}
}

#' @title Returns the percent completion of the calculation of the goodnessGrid for a given job.
#'
#' @param id The id of the job to query.  This is always the associated timestamp for the job.
#' @return The status of the job as a percentage between zero and one.
checkStatus = function(id) {
	return(status[toString(id)])
}

#' @title Sets the percent completion of the calculation of the goodnessGrid for a given job.
#'
#' @param id The id of the job to query.  This is always the associated timestamp for the job.
#' @param value The percent completion of the job.
#' @return None.
setStatus = function(id, value) {
	status[toString(id)] <<- value
	if(require('multicore')) {
		data = raw(2)
		data[1] = id
		data[2] = value
		sendMaster(data)
	}
}
