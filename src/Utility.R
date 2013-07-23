# Provides utility functions.

source("src/ShapeFunctions.R")

# Finds a "good" set of sensor placements for a given setup [bGrid, fGrid, params].
# Returns a list of locations as grid coordinates.
# Bias cases:
# 1. fish only
# 2. bathy only
# 3. detectable fish due to bathy
sensors <- function(numSensors, bGrid, fGrid, range, bias, params, debug=FALSE) {
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
    grids = sumGrid(grids, range, bias, params, debug)
    sumGrid = grids$sumGrid
	print(fGrid)
	print(sumGrid)
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
sumGrid<- function (grid, range, bias, params, debug=FALSE) {
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
        return(sumGrid.sumSimple(grid, "fGrid", range, debug))
    }
    #Bathy
    else if (bias == 2) {
        return(sumGrid.sumBathy(grid, range, params$shapeFcn, params, debug))
    }
    #Combo
    else if (bias == 3) {
        return(sumGrid.sumProduct(grid, range, params$shapeFcn, params, debug))
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
    if(TRUE){
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
        percentVisibility = 1 - (pnorm(targetCellsVisibleDepth,mean=mean,sd=sd))
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
	filenames$bGrid = sprintf("img/bGrid %g.png", time)
	png(filenames$bGrid)
	image(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,main='bGrid')
	contour(result$bGrid$x,result$bGrid$y,result$bGrid$bGrid,xlab='x',ylab='y',add=TRUE,nlevels=5)
	dev.off()
	#dev.new()
	filenames$fGrid = sprintf("img/fGrid %g.png", time)
	png(filenames$fGrid)
	image(result$bGrid$x,result$bGrid$y,result$fGrid,main='fGrid')
	numSensors <- length(result$sensors)
	for(i in 1:numSensors) points(result$bGrid$x[result$sensors[[i]]$r],result$bGrid$y[result$sensors[[i]]$c])
	dev.off()
	#dev.new()
	filenames$sumGrid = sprintf("img/sumGrid %g.png", time)
	png(filenames$sumGrid)
	image(1:params$XDist, 1:params$YDist ,result$sumGrid,main='sumGrid')
	dev.off()
	result$filenames = filenames
	return(result)
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
	print(names)
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
    if(!('range' %in% names)) {
        params$range = 1
    }
    # Supression Function Defaults
    if(!('supressionFcn' %in% names)) {
        params$supressionFcn = "supression.static"
        params$supressionRange = 1
        params$maxSupressionValue = 0
        params$minSupressionValue = 0
    }
    
    # Shape Function Defaults
    if(!('shapeFcn' %in% names)) {
        params$shapeFcn= "shape.t"
        params$sd=1
        params$peak=.75 
    }
    # Bathymetry defaults
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
    }
    # Fish Modeling
    if(!('fishmodel' %in% names)) {
        params$fishmodel <- 'rw'
    }
    return(params)
}



dist = 0
supressionRange = 100
minSupressionValue = .25
maxSupressionValue = .75
params = {}
print(supression.scale(dist, supressionRange, minSupressionValue, 
        maxSupressionValue, params))

