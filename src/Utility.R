# Provides utility functions.

source("src/ShapeFunctions.R")

# Finds a "good" set of sensor placements for a given setup [bGrid, fGrid, params].
# Returns a list of locations as grid coordinates.
# Bias cases:
# 1. fish only
# 2. bathy only
# 3. detectable fish due to bathy
sensors <- function(numSensors, bGrid, fGrid, range, bias, params, debug=FALSE, opt=FALSE) {
    if (debug) {
        cat("\n[sensors]\n")
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
    grids = sumGrid(grids, range, bias, params, debug, opt)
    sumGrid = grids$sumGrid
    # for each sensor, find a good placement
    for (i in 1:numSensors) {
        # find the max location 
        maxLoc = which.max(grids$sumGrid)
        # Switch the row/col vals since R references Grid coords differently
        c=ceiling(maxLoc/rows)
        r=(maxLoc %% rows)
        if (r==0) {
            r=rows
        }
        maxLoc = list(c=c,r=r)
        # append maxLoc to the sensor list.
        sensorList = c(sensorList, list(maxLoc))
        # down-weigh all near-by cells to discourage them from being chosen by the program
        grids$sumGrid = supress(grids$sumGrid, dim(fGrid), maxLoc, params$supressionFcn, 
                                params$supressionRange, params$minSupressionValue, 
                                params$maxSupressionValue, params, debug)
    }
    return(list(sensorList=sensorList, sumGrid=sumGrid))
}


# Calculates the composite "goodness" grid for a particular bias.
sumGrid<- function (grid, range, bias, params, debug=FALSE, opt=FALSE) {
    if (debug) {
        cat("\n[sumGrid]\n")
        print("bGrid")
        print(grid$bGrid)
        print("fGrid")
        print(grid$fGrid)
        print(sprintf("bias=%g", bias))
        print("params")
        print(params)
    }
    #Fish
    if (bias == 1) {
        if(opt){
            return(sumGrid.sumSimple.opt(grid, "fGrid", range, debug))
        }else{
            return(sumGrid.sumSimple(grid, "fGrid", range, debug))
        }
    }
    #Bathy
    else if (bias == 2) {
        if(opt){
            return(sumGrid.sumBathy.opt(grid, params, debug, opt))
        }else{
            return(sumGrid.sumBathy(grid, range, params$shapeFcn, params, debug))
        }
    }
    #Combo
    else if (bias == 3) {
        if(opt){
            ## Note sumGrid.sumBathy.opt also handles bias 3
            return(sumGrid.sumBathy.opt(grid, params, debug,opt))
        }else{
            return(sumGrid.sumProduct(grid, range, params$shapeFcn, params, debug))
        }
    }
    else {
        write("ERROR: Invalid Bias", stderr())
    }   
}

# Simply sums the values within range of a cell for each cell in the given grid.
sumGrid.sumSimple <- function (grid, key, range, debug=FALSE) {
    tempGrid = get(key, grid)
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
    
    grid$sumGrid = tempGrid
    if(debug){
        cat("\n[sumGrid.sumSimple]\n")
        print("grid")
        print(grid)
    }
    return(grid)
}

# Simply sums the values within range of a cell for each cell in the given grid.
# [optimized, but gives different results than non-opt version, why?, see
# TestUtility.R for speed comparison]
sumGrid.sumSimple.opt <- function (grid, key, range, debug=FALSE) {
    kernel <- rep(1,2*range+1) ## Assume that range is an integer
    tempGrid = get(key, grid)
    grid$sumGrid <- conv.2D(tempGrid,kernel,kernel)

    if(debug){
        cat("\n[sumGrid.sumSimple]\n")
        print("grid")
        print(grid)
    }
    return(grid)
}

# Sums the result of calling the detect() function on each cell within range of 
# a target cell for each cell in the given grid.
sumGrid.sumBathy <- function (grid, range, shapeFcn="shape.t", 
        params, debug=FALSE) {
    
    sumGrid0 = grid$bGrid$bGrid
    tempCpy = grid$bGrid$bGrid
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
    
    grid$sumGrid = sumGrid0
    if(debug){
        cat("\n[sumGrid.sumBathy]\n")
        print("visibilities")
        print(visibilities)
        print("grid")
        print(grid)
    }
    return(grid)
}

# Sums the result of calling the detect() function on each cell within range of 
# a target cell for each cell in the given grid.
# [optimized version, which also supports bias 3]
sumGrid.sumBathy.opt <- function (grid, params, debug=FALSE,opt=FALSE) {

    ##grid <- list(bGrid=bGrid)
    nr <- dim(grid$bGrid$bGrid)[1]
    nc <- dim(grid$bGrid$bGrid)[2]
    ng <- nr*nc
    bG <- grid$bGrid$bGrid
        
    sumGrid <- matrix(0,nr,nc) ## Allocate memory
    rng <- round(params$range) ## Round to integer range
    params$sensorElevation <- 1
    sensorDepth <- bG + params$sensorElevation
    belowSurf <- sensorDepth < 0
    land <- bG >= 0
    dpflag <- "depth_off_bottom" %in% params && "depth_off_bottom_sd" %in% params
    usefGrid <- params$bias==3
    for(c in 1:nc){
        cind <- max(c(1,c-rng)):min(c(nc,c+rng))
        for(r in 1:nr){
            if(belowSurf[r,c]){ ## Only calculate if sensor is below surface
                rind <- max(c(1,r-rng)):min(c(nr,r+rng))
                pV <- calc.percent.viz(r,c,rind,cind,ng,nr,bG,land,sensorDepth,dpflag,params)
                probOfRangeDetection = do.call(params$shapeFcn, list(pV$dists, params))
                if(usefGrid) probOfRangeDetection <- probOfRangeDetection * grid$fGrid[pV$inds]
                sumGrid[r,c] = sum(probOfRangeDetection * pV$percentVisibility)
            }
        }
    }
    
    grid$sumGrid = sumGrid
    if(debug){
        cat("\n[sumGrid.sumBathy]\n")
        print("grid")
        print(grid)
    }
    return(grid)
}


calc.percent.viz <- function(r,c,rind,cind,ng,nr,bG,land,sensorDepth,dpflag,params){
    ## Find rows and columns for relevant cells
    rvec <- rep(rind,length(cind))
    cvec <- sort(rep(cind,length(rind)))
    tmp <- which(!(rvec==r & cvec==c)) ## Remove self cell
    rvec <- rvec[tmp]
    cvec <- cvec[tmp]
    inds <- sub2ind(rvec,cvec,nr) ## Translate to single index
    ninds <- length(inds)
    ## Get depths, dists and slopes
    disttmp <- sort(sqrt((r-rvec)^2 + (c-cvec)^2),decreasing=TRUE,index=TRUE) ## This sorts after dist so longest dists are calculated first, then shorter ones might not be needed since they are already calculated for a long dist
    dists <- disttmp$x
    depths <- bG[inds[disttmp$ix]]
    slopes <- (depths-sensorDepth[r,c])/dists
    ibig2ismall <- rep(0,ng)
    ibig2ismall[inds[disttmp$ix]] <- 1:ninds
    vizDepths <- rep(-1e-4,ninds) ## Assign small negative number to avoid problem with being exactly at the surface in pnorm
    remaining <- 1:ninds

    while(length(remaining)>0){ ## Calculate visible depths
        ii <- remaining[1]
        losinds <- getCells.new(list(r=r,c=c),list(r=rvec[disttmp$ix[ii]],c=cvec[disttmp$ix[ii]]), debug=FALSE, nr)
        is <- ibig2ismall[losinds] ## Get indices in small vectors (not whole grid)
        d2 <- sort(dists[is],index=TRUE)
        blocks <- land[losinds[d2$ix]] ## If LOS is blocked by land don't calculate for cells behind
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
          
    vizDepths[vizDepths>0] <- -1e-4 ## Visible depths above water not valid (assign a number a little smaller than zero [just below surface])
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

# Determines the likelihood of a tag at a given position is detectable by a sensor at a 
# given position, using a specific shapeFunction.  This function considers Bathymetry and
# sensor range.
# Returns the percent chance of detection as a double between 0 [no chance of detection] 
# and 1 [guaranteed detection].
# [usable only in optimized version]
detect.new <- function(bGrid, sensorPos, tagPos, shapeFcn, params, debug=FALSE,opt=FALSE) {
    rvec = tagPos$rs:tagPos$re
    nr <- length(rvec)
    cvec = tagPos$cs:tagPos$ce
    nc <- length(cvec)
    rmat <- matrix(rep(rvec,nc),nr,nc)
    cmat <- matrix(rep(cvec,nr),nr,nc,byrow=TRUE)

    dist = sqrt((sensorPos$c - cmat)^2 + (sensorPos$r - rmat)^2)
    probOfRangeDetection = do.call(shapeFcn, list(dist, params))

    probOfLOSDetection = matrix(0,nr,nc) ## Allocate memory
    for (r in 1:nr) {
        for (c in 1:nc) {
            probOfLOSDetection[r,c] = checkLOS(bGrid, sensorPos, list(r=rvec[r],c=cvec[c]), params, debug)
        }
    }
    
    ##probOfLOSDetection = checkLOS(bGrid, sensorPos, tagPos, params, debug)
    probOfDetection = probOfRangeDetection * probOfLOSDetection
    if(debug) {
        cat("\n[detect]\n")
        print(sprintf("probOfLOSDetection=%g",probOfLOSDetection))
        print(sprintf("probOfRangeDetection=%g",probOfRangeDetection))
        print(sprintf("TotalProbOfDetection=%g",probOfDetection))
    }
    ##return(probOfRangeDetection * probOfLOSDetection)
    return(probOfDetection)
}


# For each cell in a given grid, the function sums (the number of fish
# times the probability of detection) for all cells within range
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
                            detect(tempCpy, sensorPos=list(r=i,c=j), tagPos=list(r=r,c=c), shapeFcn=shapeFcn,
                                    params, debug) * fGrid[r,c])
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

# Supresses the values of cells around a sensor using various supressionFunctions.
# minSupression: The minimum value to return
# maxSupression: The maximum value to return (also the return value for supression.static())
supress <- function(grid, dims, loc, supressionFcn, supressionRange,
                    minSupressionValue, maxSupressionValue, params, debug=FALSE) {
    if(debug) {
        cat("\n[supress]\n")
        print(sprintf("supressionFcn: %s", supressionFcn))
        print(sprintf("loc: (%g,%g)",loc$c,loc$r))
        print("grid")
        print(grid)
    }
    vals = getArea(loc, dims, supressionRange)
    mini = vals$rs
    maxi = vals$re
    minj = vals$cs
    maxj = vals$ce
    for (i in mini:maxi) {
        for (j in minj:maxj) {
                    dist = sqrt((loc$c - j)^2 + (loc$r - i)^2)
                    grid[i,j] = grid[i,j] * do.call(supressionFcn, list(dist, supressionRange, 
                                        minSupressionValue, maxSupressionValue, params, debug))
            }
    }
    return(grid)
}

# Returns a static value defined by 'maxSupressionValue'
supression.static <- function (dist, supressionRange, minSupressionValue, 
                               maxSupressionValue, params, debug=FALSE) {
    return (maxSupressionValue)
}

# Returns a dynamic value based on distance from a point.
# Returned value should be multiplied by the value to be scaled.
supression.scale <- function (dist, supressionRange, minSupressionValue, 
        maxSupressionValue, params, debug=FALSE) {
    
    sRange = minSupressionValue - maxSupressionValue
    value = 1 - (sRange * (dist/supressionRange) + maxSupressionValue)
    value = max(0, value)
    value = min(1, value)
    if (debug) {
        print(sprintf("dist=%g", dist))
        print(sprintf("supressionRange=%g", supressionRange))
        print(sprintf("minSupressionValue=%g", minSupressionValue))
        print(sprintf("maxSupressionValue=%g", maxSupressionValue))
        print(sprintf("value=%g", value))
    }
    return (value)
}

# Defines the "shape" of a sensor's range, returns an set of start/end indexes
# for rows and columns respectively named : {rs,re,cs,ce}.
getArea<-function(loc, dim, range, debug=FALSE) {
    r = loc$r # the row index for our central point
    c = loc$c # the col index for our central point
    rows = dim[1] # the max number of rows in the grid
    cols = dim[2] # the max number of cols in the grid
    
    # defines a square
    rs0 = max(1, r - range) 
    re0 = min(rows, r + range)
    cs0 = max(1 ,c - range)
    ce0 = min(cols, c + range)
    toRet = list(rs=rs0, re=re0, cs=cs0, ce=ce0)
    return(toRet)
}


# Determines the likelihood of a tag at a given position is detectable by a sensor at a 
# given position, using a specific shapeFunction.  This function considers Bathymetry and
# sensor range.
# Returns the percent chance of detection as a double between 0 [no chance of detection] 
# and 1 [guaranteed detection].
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


# Returns the percent of the water column visible at a target cell from a
# starting cell.
checkLOS<- function(bGrid, startingCell, targetCell, params, debug=FALSE) {
    sensorElevation = params["sensorElevation"]
    sensorElevation = 1
    dist = sqrt((startingCell$c - targetCell$c)^2 + (startingCell$r - targetCell$r)^2)
    if (dist ==0) {
        return(1)
    }
    # our sensor's z value
    sensorDepth = bGrid[startingCell$r, startingCell$c] + sensorElevation
    # retrieve list of intervening cells
    table = getCells(startingCell, targetCell, debug) ######getCells returns nothing because the cells are adjacent...

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

# Returns the cells crossed by a beam from the starting cell to
# the target cell.
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
    # uniques
    grid = unique(grid)
    # remove start and end cells
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

# Returns the cells crossed by a beam from the starting cell to
# the target cell.
# [Optimized version, see TestUtility.R to evaluate speed gain]
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


# Returns the cells crossed by a beam from the starting cell to
# the target cell.
# [Does not produce same output as getCells]
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

# Offsets a cartesian point towards the center of the gridcell it represents.
# ex: the cartesian point (3,2) would be converted to (2.5, 1.5), which puts it in the
# cell located at the third column, second row (aka the cell at (3,2) on a 1-based grid
# system).
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

graph <- function(result, params) {
	## Plotting
	graphics.off()
	filenames = {}
	time = as.numeric(Sys.time()) %% 1
	
	## BGrid
	filenames$bGrid = sprintf("img/bGrid-%g.png", time)
	png(filenames$bGrid)
	image(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,main='bGrid')
	for(i in 1:params$numSensors) points(result$bGrid$x[result$sensors[[i]]$r],result$bGrid$y[result$sensors[[i]]$c])
	contour(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,xlab='x',ylab='y',add=TRUE,nlevels=5)
	dev.off()
	
	## FGrid
	filenames$fGrid = sprintf("img/fGrid-%g.png", time)
	png(filenames$fGrid)
	image(result$bGrid$x,result$bGrid$y,result$fGrid,main='fGrid')
	for(i in 1:params$numSensors) points(result$bGrid$x[result$sensors[[i]]$r],result$bGrid$y[result$sensors[[i]]$c])
	dev.off()
	
	## SumGrid
	filenames$sumGrid = sprintf("img/sumGrid-%g.png", time)
	png(filenames$sumGrid)
	image(params$startY:(params$startY+params$YDist),params$startX:(params$startX+params$XDist),result$sumGrid,main='sumGrid')
	#image(1:params$XDist, 1:params$YDist ,result$sumGrid,main='sumGrid')
	for(i in 1:params$numSensors) points(result$sensors[[i]]$r,result$sensors[[i]]$c)
	dev.off()
	
	## acoustic coverage
	filenames$acousticCoverage = sprintf("img/acousticCoverage-%g.png", time)
	png(filenames$acousticCoverage)
	image(params$startY:(params$startY+params$YDist),params$startX:(params$startX+params$XDist),result$stats$acousticCoverage,main='acousticCoverage')
	
	#image(1:params$XDist, 1:params$YDist ,result$stats$acousticCoverage,main='acousticCoverage')
	for(i in 1:params$numSensors) points(result$sensors[[i]]$r,result$sensors[[i]]$c)
	dev.off()
	return(filenames)
}
# Provides Statistical data on detection, given a particular bGrid, fGrid, and sensor 
# arrangement.
# Returns a dictionary of staistical values.
stats <- function(params, bGrid, fGrid, sensors) {
    statDict <- list()
    numSensors <- length(sensors)
    xSens <- rep(0,numSensors)
    ySens <- rep(0,numSensors)
    for(i in 1:numSensors){
        xSens[i] <- bGrid$x[sensors[[i]]$c]
        ySens[i] <- bGrid$y[sensors[[i]]$r]
    }
    distMat <- matrix(0,numSensors,numSensors)
    for(i in 1:numSensors){
        distMat[i,] <- sqrt((xSens[i]-xSens)^2 + (ySens[i]-ySens)^2)
    }
    ## a is the median of the distances between the receivers
    a <- median(distMat[upper.tri(distMat)])
    ## delta is a sparsity measure (see Pedersen & Weng 2013)
    statDict$delta <- a/(2*params$range) 
    ## phi is a dimensionless indicator of movement capacity relative to detection range, it can also be viewed as a signal to noise ratio
    statDict$phi <- params$msd/params$range
    ## Distance maps (the distance from any grid cell to a receiver)
    rows <- dim(fGrid)[1]
    cols <- dim(fGrid)[2]
    X <- matrix(rep(bGrid$x,rows),rows,cols,byrow=TRUE)
    Y <- matrix(rep(bGrid$y,cols),rows,cols,byrow=FALSE)
    dimap <- array(0,dim=c(rows,cols,numSensors))
    for(i in 1:numSensors) dimap[,,i] <- sqrt( (X-xSens[i])^2 + (Y-ySens[i])^2 ) ## Distance to receiver
    ##filled.contour(gy,gx,dimap[,,2],color.palette=rainbow,main=c('Distance map'),xlab='x',ylab='y',plot.axes = {axis(1); axis(2); points(r$x,r$y)})
    
    ## Detection maps
    demap <- array(0,dim=c(rows,cols,numSensors))
    for(i in 1:numSensors) demap[,,i] <- do.call(params$shapeFcn, list(dimap[,,i], params))
    ##filled.contour(gy,gx,demap[,,2],color.palette=rainbow,main=c('Detection map'),xlab='x',ylab='y',plot.axes = {axis(1); axis(2); points(r$x,r$y)})
    
    ## Coverage map
    cover <- matrix(1,rows,cols)
    for(i in 1:numSensors){
        cover <- cover * (1 - demap[,,i]) ## Probability of no detection
    }
    cover <- 1 - cover ## Detection probability at location
    statDict$acousticCoverage <- cover
    
    return(statDict)
}

# Provides default parameter values if none are provided.
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

    # Supression Function Defaults
    if(!('supressionFcn' %in% names)) {
        params$supressionFcn = "supression.static"
        params$supressionRange = 2
        params$maxSupressionValue = 0
        params$minSupressionValue = 0
    }	else {
		params$supressionFcn = as.character(params$supressionFcn)
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
		params$fishmodel = as.character(params$fishmodel)
	}
    
    return(params)
}

## Performs the convolution operation in 1D
conv.1D <- function(fun,kern){
  lk <- length(kern)
  kern <- kern[lk:1]
  lk2 <- 0.5*(lk-1)
  lf <- length(fun)
  test <- convolve(fun,kern,type='o')
  test[(lk2+1):(lf+lk-lk2-1)]
}

## Performs the convolution operation in 2D (using two 1D kernels though)
conv.2D <- function(mat,kx,ky){
  dimmat <- dim(mat)
  matout <- matrix(0,dimmat[1],dimmat[2])
  for(i in 1:dimmat[1]) matout[i,] <- conv.1D(mat[i,],kx)
  for(i in 1:dimmat[2]) matout[,i] <- conv.1D(matout[,i],ky)
  matout
}

## Convert from row col to linear index
sub2ind <- function(row,col,dim){
    (col-1)*dim[1] + row
}


##checkParams({})
