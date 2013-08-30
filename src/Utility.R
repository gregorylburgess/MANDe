#' @include src/ShapeFunctions.R
source('src/ShapeFunctions.R')
library('rjson')

#' @name sensorFun
#' @title Calls functions to generate a 'goodness' grid and choose sensor locations.
#' @details Finds a "good" set of sensor placements for a given setup [bGrid, fGrid, params].
#' Returns a list of locations as grid coordinates, and the calculated sumGrid (goodness grid).
#' Bias value controls the 'goodness' algorithm that gets called.
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
sensorFun = function(numSensors, bGrid, fGrid, range, bias, params, debug=FALSE, opt=FALSE) {
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
    for (i in 1:numSensors) {
        # find the max location 
        maxLoc = which.max(grids$sumGrid)
        # Switch the row/col vals since R references Grid coords as (y,x) instead of (x,y)
        c = ceiling(maxLoc/rows)
        r = (maxLoc %% rows)
        if (r==0) {
            r=rows
        }
        maxLoc = list(c=c,r=r)
        print(paste('Placed sensor',i))
        ##print(paste('Placed sensor',i,'at: ',maxLoc$c,maxLoc$c))
        # append maxLoc to the sensor list.
        sensorList = c(sensorList, list(maxLoc))
        # down-weigh all near-by cells to discourage them from being chosen by the program
        if(params$suppressionFcn != 'detection.function.exact'){
          ##print('NOT using detection.function.exact')
          if(opt){
            grids$sumGrid = suppress.opt(grids$sumGrid, dim(fGrid), maxLoc, params, bGrid$bGrid, debug)
          }else{
            grids$sumGrid = suppress(grids$sumGrid, dim(fGrid), maxLoc, params$suppressionFcn, 
                                    params$suppressionRange, params$minsuppressionValue, 
                                    params$maxsuppressionValue, params, debug)
          }
        }else{
            grids = updateFGrid(maxLoc,grids,params,debug,opt)
            grids = sumGridFun(grids, range, bias, params, debug, opt)
        }
    }
    return(list(sensorList=sensorList, sumGrid=sumGrid, sumGridSupp=grids$sumGrid))
}


#' @title Updates the FGrid after each sensor is placed to reflect which areas that are already covered by sensors.
#' @description When a sensor is placed the FGrid must be updated to reflect where unique signals are emitted.
#' This is done by calculating the coverage of the sensor given by loc and then downweighing the FGrid
#' using this coverage such that locations that are well covered by the sensor at loc is downweighed more
#' that poorly covered locations.
#'
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return Returns the grids parameter, with an updated FGrid.
updateFGrid = function(loc,grids,params,debug=FALSE,opt=FALSE){
  if(debug){
      cat("\n[updateFGrid]\n")
      print("loc")
      print(loc)
  }
  grid = grids$fGrid
  bG = grids$bGrid$bGrid
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
  pctviz = calc.percent.viz(loc$r,loc$c,rind,cind,bG,land,sensorDepth,dpflag,params)
  ## testmap is a matrix with size as the full grid containing the percentage visibility of each cell
  ## Initialize
  testmap = matrix(0,rows,cols)
  ## Insert values at correct indices
  testmap[pctviz$inds] = pctviz$percentVisibility
  ## 100% detected in self cell
  testmap[loc$r,loc$c] = 1
  ## Copy relevant area to dgrid1
  dgrid1 = testmap[rind,cind]
  ## dgrid contains downweighing values
  dgrid = 1 - (dgrid1 * dgrid2)
  ## Downweigh observed region
  grid[rind,cind] = grid[rind,cind] * dgrid
  grids$fGrid = grid
  return(grids)
}


#' @name sumGridFun
#' @title Calculates the composite "goodness" grid for a particular bias.
#' @description Calls a particular sumGrid function based on the bias and opt values.  Actual work
#' 			is done by the called function.
#'
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param range The range of the sensor in bathymetric cells.
#' @param bias The goodness algorithm to use, choose 1, 2, or 3.  See package manual for more details.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return Returns the grids parameter, with an updated sumGrid.
sumGridFun = function (grids, range, bias, params, debug=FALSE, opt=FALSE) {
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
            return(sumGrid.sumSimple.opt(grids, "fGrid", params, debug))
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


#' @title Calculates sumGrid for bias 1.
#' @description Simply sums the values within range of a cell, for each cell in the given grid.
#'
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param key A key to the dictionary provided in the 'grids' parameter specifying which grid should be summed.
#' @param range The range of the sensor in bathymetric cells.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns the grids parameter, with an updated sumGrid.
sumGrid.sumSimple = function (grids, key, range, debug=FALSE) {
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


#' @title Simply sums the values within range of a cell, for each cell in the given grid.
#' @description This is a speed optimized version, which uses the convolution operation (which
#' mainly gains it speed from using FFT) as an alternative to running through tedious
#' R for loops.
#'
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param key A key to the dictionary provided in the 'grids' parameter specifying which grid should be summed.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns the grids parameter, with an updated sumGrid.
sumGrid.sumSimple.opt = function (grids, key, params, debug=FALSE) {
    ## Create a vector of distances to the cells that can be sensed by a sensor in the current cell
    subdists = 1:params$range
    ## Concatenate the "other" side and add self cell
    dists = c(rev(subdists),0,subdists)
    ## Calculate the detection function value at all distances and that is our kernel
    kernel = do.call(params$shapeFcn, list(dists, params))
    ## Check that the length of the kernel is as it should be
    if(length(kernel) != 2*params$range+1) {
		print(paste('[sumGrid.sumSimple.opt]: length of kernel was:',
					length(kernel),'expected:',2*params$range+1))
	}
    ## Extract relevant grid as given by key
    tempGrid = get(key, grids)
    ## Do convolution. This operation is identical to the for loop in sumGrid.sumSimple
    ## For more general information about how the convolution operation is defined google it! wikipedia has a decent explanation.
    grids$sumGrid = conv.2D(tempGrid,kernel,kernel)

    if(debug){
        cat("\n[sumGrid.sumSimple.opt]\n")
        print("grids")
        print(grids)
    }
    return(grids)
}


#' @title Calculates sumGrid that accounts for bathymetric shadowing.
#' @description Sums the result of calling the detect() function on each cell within range of 
#' a target cell for each cell in the given grid.
#' 
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param range The range of the sensor in bathymetric cells.
#' @param shapeFcn The shape function that describes the attenuation of a sensor.  See ShapeFunctions.R for a list of supported functions.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns the grids parameter, with an updated sumGrid.
sumGrid.sumBathy = function (grids, range, shapeFcn="shape.t", 
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


#' @title Calculates the sumGrid when a line of sight bias is chosen (bias 2 or 3).
#' @description Loops through all cells where sensor placement is valid (where sensor would be below surface)
#' and calculates goodness. If bias is 2 only bathymetry (line of sight) is used to calculate goodness, whereas if
#' bias is 3 both bathymetry and fish distribution (fGrid) are used. This function uses vectorized calculations.
#' 
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return Returns the grids parameter, with an updated sumGrid.
sumGrid.sumBathy.opt = function (grids, params, debug=FALSE,opt=FALSE) {

    nr = dim(grids$bGrid$bGrid)[1]
    nc = dim(grids$bGrid$bGrid)[2]
    ## Calculate the number of cells in the bGrid
    ng = nr*nc
    ## Make a copy of the bGrid to make code look nicer (could get rid of to save memory)
    bG = grids$bGrid$bGrid
    ## Initialize the sumGrid matrix (allocate memory)
    sumGrid = matrix(0,nr,nc)
    ## Make sure we use an integer range
    rng = round(params$range)
    ## Calculate a matrix containing the depth of hypothetical sensors placed in each cell as an offset from the bottom
    sensorDepth = bG + params$sensorElevation
    ## Calculate a matrix where grid cells containing TRUE values would containg a sensor below the surface
    belowSurf = sensorDepth < 0
    ## Create a matrix where land cells have value TRUE
    land = bG >= 0
    ## If dpflag is false then proportion of water column is calculated, if true depth preference is used
    dpflag = "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params
    usefGrid = params$bias==3
    for(c in 1:nc){
        comp = c/nc
        print(sprintf("completed:%g", comp))
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
                pV = calc.percent.viz(r,c,rind,cind,bG,land,sensorDepth[r,c],dpflag,params)
                ## Calculate the detection function value at all grid points
                probOfRangeDetection = do.call(params$shapeFcn, list(pV$dists, params))
                ## If bias == 3 include the fGrid in the calculations, if not just use bathymetry and detection function
                if(usefGrid) probOfRangeDetection = probOfRangeDetection * grids$fGrid[pV$inds]
                ## Calculate goodness of cell (r,c) by summing detection probabilities of all visible cells
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

#' @name calc.percent.viz
#' @title Calculates the percentage of the water column in the surrounding cells that is visible to a sensor placed in the current cell.
#' @description Calculates a matrix centered around the current cell (r,c) and containing the percentage
#' of the water column in the surrounding cells that is visible to a sensor placed in the
#' current cell.
#'
#' 
#' @param r Row of the current cell in the bGrid.
#' @param c Column of the current cell in the bGrid.
#' @param rind Row indices of the bGrid to calculate visibility percentage.
#' @param cind Column indices of the bGrid to calculate visibility percentage.
#' @param bGrid A valid bGrid.
#' @param land Matrix containing logicals (TRUE = land cell) indicating wheter a cell 
#' in the bGrid is a land cell.
#' @param sensorDepth Depth of sensor in current cell.
#' @param dpflag If TRUE depth preference is used meaning that the percentage of visible 
#' fish is calculated, if FALSE visible water column is calculated.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @return Returns a dictionary with three keys (all vectors): percentVisibility contains
#' the percentage fish/signals visible in the surrounding cells, inds contains the linear
#' indices in the bGrid to which the visibilities pertain, dists contains the distance
#' from the current cell to each of the returned cells as given by inds.
calc.percent.viz = function(r,c,rind,cind,bGrid,land,sensorDepth,dpflag,params){
    rows = dim(bGrid)[1]
    cols = dim(bGrid)[2]
    ng = rows*cols
    nr = rows

    ## Find rows and columns for relevant cells
    rvec = rep(rind,length(cind))
    cvec = sort(rep(cind,length(rind)))
    ## Remove self cell
    tmp = which(!(rvec==r & cvec==c))
    rvec = rvec[tmp]
    cvec = cvec[tmp]
    ## Translate from row col index to to single index
    inds = sub2ind(rvec,cvec,nr) 
    ninds = length(inds)
    ## Calculate distances from the current cell to the surrounding cells within range
    ## This sorts after dist so longest dists are calculated first, then shorter ones might not be needed since they are already calculated for a long dist
    disttmp = sort(sqrt((r-rvec)^2 + (c-cvec)^2),decreasing=TRUE,index=TRUE)
    ## Save actual distances in the dists vector
    dists = disttmp$x
    ## Get depths at the sorted cells by using disttmp$ix, which contains the sorted indices
    depths = bGrid[inds[disttmp$ix]]
    ## Calculate the line of sight slopes to each of the sorted cells
    slopes = (depths-sensorDepth)/dists
    ## Create a vector, which can be used to easily map an index in the sorted vector to an index in the rind by cind matrix
    ibig2ismall = rep(0,ng)
    ibig2ismall[inds[disttmp$ix]] = 1:ninds

    ## Initialize vizDepths vector, this will be filled with visible depths below
    ## Assign small negative number to avoid problem with being exactly at the surface in pnorm
    vizDepths = rep(-1e-4,ninds)
    ## Initialize the remaining vector, which contains the indices of the cells for which the vizDepth has not yet been calculated
    remaining = 1:ninds

    ## Calculate visible depths as long as uncalculated cells remain
    while(length(remaining)>0){
        ## Since cells are sorted by decreasing distance taking the first index of remaining always gives the farthest uncalculated cell
        ii = remaining[1]
        ## Find the indices of the cells within the line of sight from r,c to ii, losinds are indices of big grid
        losinds = getCells.opt(list(r=r,c=c),list(r=rvec[disttmp$ix[ii]],c=cvec[disttmp$ix[ii]]), debug=FALSE, nr)
        ## Get indices in small sorted vector (not whole grid)
        is = ibig2ismall[losinds]
        ## Sort the distances within LOS
        d2 = sort(dists[is],index=TRUE)
        ## Find indices of obstacles (land areas) in LOS. If LOS is blocked by land don't calculate for cells behind.
        blocks = land[losinds[d2$ix]]
        
        if(any(blocks, na.rm=TRUE)){
            if(!all(blocks)){
                ## Find indices that are not blocked
                indsNoBlock = 1:(min(which(blocks))-1)
            }else{
                indsNoBlock = NULL
            }
        }else{
            ## If theres no obstacles all cells in LOS are actually visible
            indsNoBlock = 1:length(losinds)
        }

        ## Vectorized calculation of the visible depths of the unblocked cells with LOS using the simple formula for a line y = ax + b
        ## a is cummax(slopes[is[d2$ix[indsNoBlock]]])
        ## Here cummax ensures that the steepest slopes between the current cell (r,c) cell along the LOS is used for calculating visible depth, it is important that the slopes are sorted in order of increasing distance from current cell to target cells, this is handled by d2$ix
        ## x is d2$x[indsNoBlock], the sorted distances from (r,c)
        ## b is sensorDepth, which is the intercept of the line
        vizDepths[is[d2$ix[indsNoBlock]]] = cummax(slopes[is[d2$ix[indsNoBlock]]])*d2$x[indsNoBlock] + sensorDepth
        ## Remove the cells from remaining for which calculations are done
        remaining = setdiff(remaining,is)
    }
	
    ## Visible depths above water not valid (assign a number a little smaller than zero [just below surface])      
    vizDepths[vizDepths>0] = -1e-4
    ## Find indices that are not land cells
    indsNotLand = depths<0
    ## if we have normal distribution data (depth preference), use it
    if(dpflag) {
        ## compute % fish visible from sensor to target cell
        ## calculate the mean fish depth in all cells
        mean = depths[indsNotLand] + params$depth_off_bottom
        ## Values above water are set to be at the surface
        mean[mean>0] = 0
        ## Save SD in sd for prettier code
        sd = params$depth_off_bottom_sd
            
        ## Get cum probability for "available" water (between depth pref and surf)
        ## Cumulative probability below surface
        psurf = pnorm(0,mean=mean,sd=sd)
        ## Calculate the percentage of visible fish (the area under the visible part of the normal curve)
        percentVisibility = psurf - pnorm(vizDepths[indsNotLand],mean=mean,sd=sd)
        ## Probability within available water (from bottom to surface), if this value is different from 1 it means that the vertical fish distribution extends below the bottom and/or above the surface, in which case we need to normalise the probability between bottom and surface so it sums to one.
        areaToCorrectFor = psurf - pnorm(depths[indsNotLand],mean=mean,sd=sd)
        ## Correct for the distribution extending out of bounds
        percentVisibility = percentVisibility/areaToCorrectFor
    }else{
        ## if we don't have normal distribution data, assume equal distribution
        ## compute % visibility (of water column height) from sensor to target cell
        percentVisibility = vizDepths[indsNotLand] / depths[indsNotLand]
    }
    return(list(percentVisibility=percentVisibility,dists=dists[indsNotLand],inds=inds[disttmp$ix[indsNotLand]]))
}


#' @title Updates sumGrid with probability of fish detections.
#' @description For each cell in a given grid, the function sums the number of fish
#' times the probability of detection) for all cells within range.
#' 
#' @param grids A dictionary containing the keys 'bGrid', 'fGrid', and 'sumGrid', which hold a valid BGrid, FGrid and SumGrid.
#' @param range The range of the sensor in bathymetric cells.
#' @param shapeFcn The shape function that describes the attenuation of a sensor.  See ShapeFunctions.R for a list of supported functions.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns the grids parameter, with an updated sumGrid.
sumGrid.sumProduct = function (grids, range, shapeFcn="shape.t", 
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


#' @title Suppress when a sensor is placed.
#' @description Suppresses the values of cells around a sensor using a specified suppressionFunction.
#' 
#' @param sumGrid A valid SumGrid.
#' @param dims The dimensions of the BGrid.  Just call dim() on the BGrid for this.
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param suppressionFcn The shape function that describes the attenuation of a sensor.  See ShapeFunctions.R for a list of supported functions.
#' @param suppressionRange How far out to apply suppression penalties, in bathymetric cells.
#' @param minsuppressionValue The minimum allowable value to return.
#' @param maxsuppressionValue The maximum allowable value to return (also the return value for suppression.static()).
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns a suppressed sumGrid.
suppress = function(sumGrid, dims, loc, suppressionFcn, suppressionRange,
                    minsuppressionValue, maxsuppressionValue, params, debug=FALSE) {
    if(debug) {
        cat("\n[suppress]\n")
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


#' @title Suppresses the values of cells around a sensor using a specified suppressionFunction.
#' @description This is an optimized version, which uses vectorization.
#' 
#' @param sumGrid A valid SumGrid.
#' @param dims The dimensions of the BGrid.  Just call dim() on the BGrid for this.
#' @param loc A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor location.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param bGrid valid BGrid.
#' @param debug If enabled, turns on debug printing (console only).
#' @return Returns a suppressed sumGrid.
suppress.opt = function(sumGrid, dims, loc, params, bGrid, debug=FALSE) {
    if(debug) {
        cat("\n[suppress.opt]\n")
        print(sprintf("suppressionFcn: %s", params$suppressionFcn))
        print(sprintf("loc: (%g,%g)",loc$c,loc$r))
        print("sumGrid")
        print(sumGrid)
        print("bGrid")
        print(bGrid)
    }
    suppressionFcn = params$suppressionFcn
    minsuppressionValue = params$minsuppressionValue
    maxsuppressionValue = params$maxsuppressionValue
    ## dfflag indicates whether detection function variant should be used for suppression
    dfflag = suppressionFcn=='detection.function' ||
			 suppressionFcn=='detection.function.shadow' ||
			 suppressionFcn=='detection.function.exact'
    rows = dim(sumGrid)[1]
    cols = dim(sumGrid)[2]
	
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
        land = bGrid >= 0
        ## Create the depth value of a sensor placed at loc
        sensorDepth = bGrid[loc$r,loc$c] + params$sensorElevation
        ## If dpflag is false then proportion of water column is calculated, if true depth preference is used
        dpflag = "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params
        ## Calculate the proportion of signals in each of the surrounding cell that can be detected by a sensor at loc
        pctviz = calc.percent.viz(loc$r, loc$c, rind, cind, bGrid, land, sensorDepth, dpflag, params)
        ## testmap is a matrix with size as the full grid containing the percentage visibility of each cell
        ## Initialize
        testmap = matrix(0,rows,cols)
        ## Insert values at correct indices
        testmap[pctviz$inds] = pctviz$percentVisibility
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
    sumGrid[rind,cind] = sumGrid[rind,cind] * supgrid

    return(sumGrid)
}


#' @title Constant suppression.
#' @description Returns a static value, effectively setting all cells within suppressionRange of a sensor to that number.
#' 
#' @param dist The distance between two cells on a grid.
#' @param suppressionRange How far out to apply suppression penalties, in bathymetric cells.
#' @param minsuppressionValue The minimum allowable value to return.
#' @param maxsuppressionValue The maximum allowable value to return (also the return value for suppression.static()).
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
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
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
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
#' @param dims The dimensions of the BGrid.  Just call dim() on the BGrid for this.
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


#' @title Checks if a tag is detectable by a sensor.
#' @description Determines the likelihood of a tag at a given position is detectable by a sensor at a 
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
detect = function(bGrid, sensorPos, tagPos, shapeFcn, params, debug=FALSE) {

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



#' @title Returns the percent of the water column visible at a target cell from a starting cell.  
#' @description If 'depth_off_bottom' and 'depth_off_bottom_sd' are keys in 'params', then the algorithm 
#' assumes a normal distribution of fish within the specified zone and will return
#' the integral of the visible zone.  Otherwise, it assumes an equal distribution of fish 
#' throughout the watercolumn.
#' 
#' @param bGrid A valid BGrid.
#' @param startingCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor's location on the BGrid.
#' @param targetCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen tag's location on the BGrid.
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param debug If enabled, turns on debug printing (console only).
#' @return The percent of the watercolumn that is visible (as a double between 0 and 1).
checkLOS= function(bGrid, startingCell, targetCell, params, debug=FALSE) {
    sensorElevation = params$sensorElevation
    # find the linear distance between the cells
    dist = sqrt((startingCell$c - targetCell$c)^2 + (startingCell$r - targetCell$r)^2)
    if (dist ==0) {
        return(1)
    }
    # our sensor's z value
    sensorDepth = bGrid[startingCell$r, startingCell$c] + as.numeric(sensorElevation)
    # retrieve list of intervening cells
    table = getCells(startingCell, targetCell, debug)

	# annotate each cell's z value from the bGrid
    table$z =apply(table, 1, function(rows){ table$z = bGrid[rows[2],rows[1]]})
    # annotate each cell's percieved slope form our sensor to the cell
    table$m =apply(table, 1, function(row) { 
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
        sd = params$depth_off_bottom_sd
		
        # pnorm gives the percent below the given point, so subtract from 1
        # to get the percent above the given point
		# {{GREG}}
		# emulate below line format/params
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


#' @title Find cell crossed by a beam.
#' @description Returns the cells crossed by a beam from the starting cell to the target cell.
#' 
#' @param startingCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor's location on the BGrid.
#' @param targetCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen tag's location on the BGrid.
#' @param debug If enabled, turns on debug printing (console only).
#' @return A data.frame containing x and y indicies of cells crossed by a beam from the starting cell to the target cell.
getCells=function(startingCell, targetCell, debug=FALSE) {
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


#' @title Returns the cells crossed by a beam from the starting cell to
#' the target cell.
#' @description This is an optimized using vectorization. It does not
#' produce exactly the same output as getCells: the same cells are reported,
#' but the cell order is different so this function should only be used
#' together with other optimized code. 
#' 
#' @param startingCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen sensor's location on the BGrid.
#' @param targetCell A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the chosen tag's location on the BGrid.
#' @param debug If enabled, turns on debug printing (console only).
#' @param nr Number of rows in bGrid.
#' @return If nr != NULL a vector is returned containing the linear indices
#' (not row and col) of cells in the bGrid crossed by a beam from the starting
#' cell to the target cell. If nr == NULL a matrix with a row column and a col
#' column is returned.
getCells.opt = function(startingCell, targetCell, debug=FALSE, nr=NULL) {
    if(!(startingCell$ r== targetCell$r && startingCell$c == targetCell$c)){
        ## Offset starting and target cells
        sC = offset(startingCell)
        tC = offset(targetCell)
        ## Define a value slightly larger than -1
        e = 1e-6 - 1
        ## Slope of a beam (line) from start to target cell
        ## Note: this is a slope in the horizontal not vertial plane
        ## rows act as y vals, cols act as x vals
        a = (tC$r-sC$r)/(tC$c-sC$c)
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
            ## Convert to linear (not row col) indices in the bGrid
            if(!is.null(nr)) biginds = sub2ind(rows,cols,nr)
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
            ## Convert to linear (not row col) indices in the bGrid
            if(!is.null(nr)) biginds = sub2ind(rows,cols,nr)
        }
        ## Indices to use are all of them except the starting cell
        useinds = !(cols == startingCell$c & rows == startingCell$r)
        if(is.null(nr)){
          ## If nr is not input return a matrix with rows and cols
          return( cbind(cols[useinds],rows[useinds]) )
        }else{
          ## If nr was input return the linear indices in the bGrid
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
#' @param point A dictionary containing the keys 'r' and 'c', which hold the row and column indicies of the point to translate on the BGrid.
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
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param showPlots If TRUE plots are shown on the screen, if FALSE plots are stored in the img folder.
#' @param plot.bathy Specifies whether contour lines for bathymetry should be overlayed in the graphs.
#' @return A dictionary containing the filenames of the generated images.
graph = function(result, params, showPlots, plot.bathy=TRUE) {
	time = "1"
        if(!showPlots) {
			if('timestamp' %in%  names(params)) {
				print("Using timestamp")
				time = params$timestamp
			}
			if(!file.exists("img")) {
				dir.create("img")
			}
		}
	## Plotting
	graphics.off()
	filenames = {}
        xlab = 'x dir'
        ylab = 'y dir'
		## Write results to a text file
		filename = sprintf("txt/%s.txt", time)
		file.create(filename)
		jsonFile = sprintf("txt/%s.json", time)
		file.create(jsonFile)
		
		cat(toJSON(result), file=jsonFile, append=FALSE)
		capture.output(print(result), file=filename)
		filenames$txt = filename
	## BGrid
	filenames$bGrid = sprintf("img/bGrid-%s.png", time)
	if(!showPlots){
            png(filenames$bGrid)
        }else{
            dev.new()
        }
        plotGrid(result,type='bGrid',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	if(!showPlots) dev.off()
	
	## FGrid
	filenames$fGrid = sprintf("img/fGrid-%s.png", time)
	if(!showPlots){
            png(filenames$fGrid)
        }else{
            dev.new()
        }
        plotGrid(result,type='fGrid',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	if(!showPlots) dev.off()
	
	## SumGrid
	filenames$sumGrid = sprintf("img/sumGrid-%s.png", time)
	if(!showPlots){
            png(filenames$sumGrid)
        }else{
            dev.new()
        }
        plotGrid(result,type='sumGrid',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	if(!showPlots) dev.off()
	
	## Acoustic Coverage
	filenames$acousticCoverage = sprintf("img/acousticCoverage-%s.png", time)
	if(!showPlots){
            png(filenames$acousticCoverage)
        }else{
            dev.new()
        }
        plotGrid(result,type='acousticCoverage',xlab=xlab,ylab=ylab,plot.bathy=plot.bathy)
	if(!showPlots) dev.off()

        ## Unique Recovery Rate
	filenames$recoveryRates = sprintf("img/recoveryRates-%s.png", time)
	if(!showPlots){
            png(filenames$recoveryRates)
        }else{
            dev.new()
        }
        plotUniqueRR(result)
	if(!showPlots) dev.off()
	
	filename = paste("zip/", time, ".zip", sep="")
	zip(zipfile=filename, files=filenames, flags="-r9X", extras="", zip=Sys.getenv("R_ZIPCMD", "zip"))
	filenames$zip = filename
	return(filenames)
}

#' @title Plots the grid specified by the input type.
#' @description In addition to the grid itself sensor locations are also plotted along
#' with numbers indicating the order in which sensors were placed. Furthermore, bathymetry
#' contours can be overlayed using the plot.bathy flag.
#' 
#' @param result A dictionary of return objects, the result of a successfull call to run() or sensorFun().
#' @param type Character specifying grid type. Available grids: 'bGrid', 'fGrid', 'sumGrid', or 'acousticCoverage'.
#' @param main Set title of plot.
#' @param xlab Set label of x axis.
#' @param ylab Set label of y axis.
#' @param plot.bathy Specifies whether contour lines for bathymetry should be overlayed in the graphs.
#' @return Nothing.
plotGrid = function(result,type='bGrid',main=type,xlab='',ylab='',plot.bathy=TRUE){
    if(type=='bGrid'){
        grid = result$bGrid$bGrid
    }
    if(type=='fGrid'){
        grid = result$fGrid
    }
    if(type=='sumGrid'){
        grid = result$sumGrid
    }
    if(type=='acousticCoverage'){
        grid = result$stats$acousticCoverage
    }
    ## Plot the actual grid as an image
    image(result$bGrid$x,result$bGrid$y,grid,main=main,xlab=xlab,ylab=ylab)
    if(plot.bathy) {
        ## Add bathymetry contour
        contour(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,add=TRUE,nlevels=5)
    }
    ## Add sensors and their numbers
    plotSensors(result)
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
    ## Find number of sensors for which unique recovery rate was calculated
    nsmax = length(result$stats$uniqRRs)
    ## Calculate the max value to use for the y-axis in TOP PLOT
    ## It looks good to show 1 as the max y val, but only if we are relatively close (above 0.7)
    ymax = ifelse(max(result$stats$uniqRRs)>0.7,1.02,max(result$stats$uniqRRs))
    ## Make two way plot
    par(mfrow=c(2,1),las=1)
    ## TOP PLOT
    plot(0:ns,c(0,result$stats$uniqRRs[1:ns]),typ='l',xlab='Number of sensors',ylab='Unique recovery rate',ylim=c(0,ymax),xlim=c(0,nsmax))
    points(0:ns,c(0,result$stats$uniqRRs[1:ns]),pch=46,cex=3)
    lines(ns:length(result$stats$uniqRRs),result$stats$uniqRRs[ns:nsmax],lty=2)
    plotIntersect(ns,result$stats$uniqRecoveryRate,col='orange',lty=1)
    grid()
    text(0.05*length(result$stats$uniqRRs),result$stats$uniqRecoveryRate,round(result$stats$uniqRecoveryRate,digits=4),pos=3)
    legend('bottomright',c('Calculated','Requested','Projected'),lty=c(1,1,2),col=c(1,'orange',1),bg='white')
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
  ns = length(result$sensors)
  ## Radius of circle
  r = result$params$detectionRange
  ## Radian values for a full circle
  a = seq(0, 2 * pi, length.out=100)
  ## Cols
  sensx = result$bGrid$x[result$stats$sensorMat[1:ns, 2]]
  ## Rows
  sensy = result$bGrid$y[result$stats$sensorMat[1:ns, 1]] 
  points(sensx,sensy,pch=21,bg='blue',cex=3)
  text(sensx,sensy,1:ns,col='white')
  if(circles){
    for(i in 1:ns){
      X = r*cos(a) + sensx[i]
      Y = r*sin(a) + sensy[i]
      lines(X,Y,lty=circlty)
    }
  }
}


#' @title Provides statistical data on detection, given a particular bGrid, fGrid, and sensor arrangement.
#' @description This function calculates placements of projected sensors and the value
#' of each sensor (requested and projected) as given in increase in unique recovery rate.
#' Additionally, the acoustic coverage map, unique recovery rate, absolute recovery rate, and sparsity
#' are calculated and returned after placing the requested number of sensors.
#' 
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param bGrid A valid BGrid.
#' @param fGrid A valid FGrid.
#' @param sensors The result of a successful call to sensorFun().
#' @param debug If enabled, turns on debug printing (console only).
#' @param opt Tells the program to use vectorized R commands.
#' @return A dictionary of statistical values containing the keys: "delta", "sensorMat"         
#' "uniqRRs", "acousticCoverage", "absRecoveryRate", "uniqRecoveryRate".
getStats = function(params, bGrid, fGrid, sensors, debug=FALSE, opt=FALSE) {
	if(debug) {
		print("[getstats]")
	}
    statDict = list()
    ## Number of requested sensors
    numSensors = length(sensors$sensorList)
    rows = dim(fGrid)[1]
    cols = dim(fGrid)[2]
    rng = params$range
    
    ## Calculate the value of each sensor (as the increase in unique recoveryrate)
    numProj = 2*numSensors
    ## sumGrid suppressed by numSensors
    sumGridSupp = sensors$sumGridSupp
    sensorList = sensors$sensorList
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
        sumGridSupp = suppress(sumGridSupp, dim(fGrid), maxLoc, params$suppressionFcn, 
          params$suppressionRange, params$minsuppressionValue, 
          params$maxsuppressionValue, params, debug)
      }
    }

    xSens = rep(0,numProj)
    ySens = rep(0,numProj)
    for(i in 1:numProj){
        xSens[i] = sensorList[[i]]$c
        ySens[i] = sensorList[[i]]$r
    }
    
    ## Distance maps (the distance from any grid cell to a receiver)
    rows = dim(fGrid)[1]
    cols = dim(fGrid)[2]
    X = matrix(rep(1:cols,rows),rows,cols,byrow=TRUE)
    Y = matrix(rep(1:rows,cols),rows,cols,byrow=FALSE)
    dimap = array(0,dim=c(rows,cols,numProj))
    for(i in 1:numProj) {
		## Distance to receiver
		dimap[,,i] = sqrt( (X-xSens[i])^2 + (Y-ySens[i])^2 )
	}
    
    ## Horizontal detection maps using detection function
    demap = array(0,dim=c(rows,cols,numProj))
    for(i in 1:numProj) demap[,,i] = do.call(params$shapeFcn, list(dimap[,,i], params))

    ## Incorporate vertical detection probability using line of sight
    if(params$bias!=1){
      bG = bGrid$bGrid
      ## Create a matrix where land cells have value TRUE
      land = bG >= 0
      ## Calculate a matrix containing the depth values of a sensor placed in each grid cell
      sensorDepth = bG + params$sensorElevation
      ## If dpflag is false then proportion of water column is calculated, if true depth preference is used
      dpflag = "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params
      for(i in 1:numProj){
        r = ySens[i]
        c = xSens[i]
        cind = max(c(1,c-rng)):min(c(cols,c+rng))
        rind = max(c(1,r-rng)):min(c(rows,r+rng))
        ## Calculate the proportion of signals in each of the surrounding cell that can be detected by a sensor at loc
        pctviz = calc.percent.viz(ySens[i],xSens[i],rind,cind,bG,land,sensorDepth[ySens[i],xSens[i]],dpflag,params)
        ## testmap is a matrix with size as the full grid containing the percentage visibility of each cell
        ## Initialize
        testmap = matrix(0,rows,cols)
        ## Insert values at correct indices
        testmap[pctviz$inds] = pctviz$percentVisibility
        ## 100% detected in self cell
        testmap[r,c] = 1
        ## Update demap (detection map)
        demap[,,i] = demap[,,i] * testmap
      }
    }

    ## Coverage map
    cover = matrix(1,rows,cols)
    uniqRRs = rep(0,numProj)
    for(i in 1:numProj){
        r = ySens[i]
        c = xSens[i]
        cind = max(c(1,c-rng)):min(c(cols,c+rng))
        rind = max(c(1,r-rng)):min(c(rows,r+rng))
        ## Probability of no detection
        cover[rind,cind] = cover[rind,cind] * (1 - demap[rind,cind,i]) 
        covertmp = 1-cover
        uniqRRs[i] = sum(covertmp * fGrid)
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
    statDict$acousticCoverage = 1-apply(1-demap[,,srt$ix[1:numSensors]],c(1,2),prod)

    ## Calculate distance matrix needed to calculate sparsity
    distVec = rep(0,numSensors)
    for(i in 1:numSensors){
        dists = sqrt((bGrid$x[xSens[srt$ix[i]]]-bGrid$x[xSens[srt$ix[1:numSensors]]])^2 + (bGrid$y[ySens[srt$ix[i]]]-bGrid$y[ySens[srt$ix[1:numSensors]]])^2)
        distVec[i] = min(dists[dists>0])
    }
    ## a is the median of the distances between the receivers
    a = median(distVec)
    ## delta is a sparsity measure (see Pedersen & Weng 2013)
    statDict$delta = a/(2*params$detectionRange)
    
    ## Absolute recovery rate (here we don't care about getting the same ping multiple times)
    demapmat = apply(demap[,,srt$ix[1:numSensors]],c(1,2),sum)
    statDict$absRecoveryRate = sum(demapmat*fGrid)
	
    ## Calculate recovery rate of unique detections    
    statDict$uniqRecoveryRate = sum(statDict$acousticCoverage * fGrid)
    
    return(statDict)
}


#' @title Check model parameters.
#' @description Provides default parameter values if none are provided.
#' 
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @return The 'params' parameter, populated with default values where necessary.
checkParams = function(params) {

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

    # Bathymetry defaults
    if(('inputFile' %in% names)) {
        params$inputFile = as.character(params$inputFile)
    }
    if(('inputFileType' %in% names)) {
        params$inputFileType = as.character(params$inputFileType)
    }   else {
                params$inputFileType = 'custom'
        }
    if(!('cellSize' %in% names)) {
        params$cellSize = 1
    }
    if(!('startX' %in% names)) {
        params$startX = 1
    }
    if(!('startY' %in% names)) {
        params$startY = 1
    }
    if(!('XDist' %in% names)) {
        params$XDist = 21
    }
    if(!('YDist' %in% names)) {
        params$YDist = 21
    }
    if(!('seriesName' %in% names)) {
        params$seriesName = 'z'
    }	else {
		params$seriesName = as.character(params$seriesName)
	}

    # Shape Function Defaults
    if(!('shapeFcn' %in% names)) {
        params$shapeFcn = "shape.gauss"
        params$detectionRange = 3*params$cellSize
        params$peak = 0.98
    }	else {
		params$shapeFcn = "shape.gauss" ##as.character(params$shapeFcn) ## Always use Gauss for simplicity
	}
	##if(!('range' %in% names)) {
	##	params$range = 3*params$sd
	##}
    if(!('sensorElevation' %in% names)){
        params$sensorElevation = 1
    }   else {
                params$sensorElevation = as.numeric(params$sensorElevation)
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
    
    # Suppression Defaults
    if(!('suppressionFcn' %in% names)) {
        params$suppressionFcn = "detection.function"
    }	else {
		params$suppressionFcn = as.character(params$suppressionFcn)
	}
    if(!('suppressionRangeFactor' %in% names)) {
        params$suppressionRangeFactor = 2
    }	else {
		params$suppressionRangeFactor = as.numeric(params$suppressionRangeFactor)
                if(params$suppressionRangeFactor > maxSuppressionRangeFactor){
                    params$suppressionRangeFactor = maxSuppressionRangeFactor
                }
	}
    
    # Fish Modeling
    if(!('fishmodel' %in% names)) {
        params$fishmodel = 'rw'
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

#' @title Converts input parameters from meters to grid cells.
#' @description These are used for internal calculations and are invisible to the user.
#' 
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.
#' @param bGrid A valid bGrid.
#' @return The 'params' parameter, populated with necessary internal variables ("sd",
#' "suppsd", "range", "suppressionRange") with unit 'grid cells'.
convertMetersToGrid = function(params,bGrid){
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
  params$suppressionRange = round(params$suppressionRangeFactor *
				 params$detectionRange / cellSize) 
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
#' @return A matrix containing the result of the convolution operation with same dimensions as mat.
conv.2D = function(mat, kx, ky){
  dimmat = dim(mat)
  ## Initialize return matrix
  matout = matrix(0, dimmat[1], dimmat[2])
  ## Convolve in x-direction
  for(i in 1:dimmat[1]) matout[i,] = conv.1D(mat[i,],kx)
  ## Convolve in y-direction (important to use matout here)
  for(i in 1:dimmat[2]) matout[,i] = conv.1D(matout[,i],ky)
  return(matout)
}


#' @title Converts a row, col index to a linear index within a matrix.
#' 
#' @param row The row index.
#' @param col The col index.
#' @param dims The dimensions of the BGrid.  Just call dim() on the parent matrix for this.
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
