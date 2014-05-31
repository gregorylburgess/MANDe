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
#' goodness grid (goodnessGrid) is calculated depending on the design parameters (this is the heavy part).
#' When the goodnessGrid is finished sensors can be placed optimally and stats and figures can be generated.
#' Source code for building the package and installing a web-based GUI can be downloaded from https://github.com/gregorylburgess/acoustic-deploy.
#' @section INPUT DETAILS params:
#' - Behaviour parameters
#' 
#' $fishmodel: defines the fish model algorithm. Possible values are random walk ('rw') or Ornstein-Uhlenbeck ('ou'). Default: RW. The choice here should reflect the prior knowledge about the long-term distribution of fish within the study region.
#'
#' 'rw' means purely diffusive movements with no locational or directional preference. This results in a uniform distribution of fish over the study region. A random walk is useful when no knowledge about fish behaviour is available or can be assumed.
#'
#' 'ou' results in movements within a limited region mimicking fish with a home range. This results in a normal distribution with two parameters: 1) the home range center (mean) specified by $mux and $muy as proportions of the X and Y axes. 2) the home range extents (variance) and shape (covariance) specified by the standard deviations in X and Y directions in meters, $ousdx and $ousdy respectively, and by the correlation between X and Y specified by $oucor. This is useful if the study species is known to prefer a specific geographic region. A rule of thumb says that for each direction in isolation approximately 95% of the time is spent within plus minus two standard deviations. The correlation is used if the elliptical home range shape is angled relative to the X,Y coordinate system. The correlation must be in the range [-1; 1]. Note: multiple home range centers can be specified by assigning vectors containing values for each home range, e.g. $mux = c(0.3, 0.6) etc..
#' 
#' Optional: Vertical habitat restriction. Useful if the fish is known to live in a certain vertical habitat say from -10 to -50 meters. If unspecified the species is assumed to be able to utilize the whole water column. The minimum and maximum depth must be specified in meters using $mindepth (shallow) and $maxdepth (deep), for example $mindepth = -5 and $maxdepth = -10. Only areas within the vertical habitat are considered in the network design.
#'
#' Optional: Preferred depth relative to bottom. This is specified by a normal distribution with a mean preferred height relative to the bottom ($depth_off_bottom) given in meters off the bottom, and a standard deviation ($depth_off_bottom_sd) given in meters.
#'
#' - Sensor parameters
#'
#' $numSensors: Specifies the number of sensors to be placed. Positive integer values are accepted.
#'
#' $projectedSensors: Specifies the number of sensors to project. This is can be used to see the development in sensor value if more sensors were placed.
#'
#' $userSensorList: Userdefined sensor locations specified in grid cell coordinates. For example, pre-placed sensors in row1=3 col1=4 and row2=7 col2=9 are specified as $userSensorList <- '3, 4, 7, 9'.
#' 
#' $bias: Specifies how the "goodness" of a cell is determined. Possible values are 1, 2, or 3. Default: 1. $bias = 1 indicates that a "good" cell has a high number of animals within detection range (ignoring line of sight). This is useful for sensors not restricted to line-of-sight detection. $bias = 2 indicates that a "good" cell has the best visibility (taking into account bathymetry and shadowing, but completely ignoring fish density). This is useful for networks restricted to line-of-sight detection and having no prior knowledge of animal habitat. $bias = 3 indicates that a "good" cell has a high number of visible fish (incorporating both bathymetry and animal density). This is useful for networks restricted to line-of-sight detection, and having some idea of animal habitat.
#'
#' $sensorElevation: Specifies the sensor vertical placement height in meters over the bottom. This value is important if accounting for line of sight ($bias = 2 or 3).
#' 
#' $shapeFcn: Sensor detection function (or shape function). Currently the only possible value is 'shape.gauss'. The detection function determines which functional shape to represent horizontal acoustic attenuation in the water. The detection function specifies how the chance of signal detection declines as a function of distance to sensor. Ranging experiments should preferably be carried out locally at the study site to determine appropriate values of the two detection function parameters:
#' 
#'  $peak: The probability of detecting a fish located right next to the sensor. Specifies a maximum value for the shape function. Values should be a decimal between 0.05 and 1.
#' 
#'  $detectionRange: The distance in meters from the sensor where the chance of detecting a signal is 0.05.
#'
#' - Bathymetry parameters
#' 
#' $inputFile: Path of a bathymetry file relative to the current working directory.
#' 
#' $inputFileType: Specifies the file type of the bathymetry file. Possible values are 'netcdf, 'arcgis', and 'asc'.
#' 
#' $seriesName: For $inputFileType = 'netcdf' or 'arcgis', this specifies the name of the data series to use.
#' 
#' $cellSize: Specifies the size of the grid cell in meters.
#' 
#' $startX: Specifies the starting X cell-index of the area of interest for the file chosen. Valid values are non-negative integers.
#' 
#' $startY: Specifies the starting Y cell-index of the area of interest for the file chosen. Valid values are non-negative integers.
#' 
#' $XDist: Specifies how far the area of interest extends along the X axis. Units are in number of bathymetric grid cells. Valid values are non-negative integers.
#' 
#' $YDist: Specifies how far the area of interest extends along the Y axis. Units are in number of bathymetric grid cells. Valid values are non-negative integers.
#'
#' - Suppression parameters
#' 
#' To prevent the program from placing sensors too near or on top of each other, a suppression function is used to reduce the goodness of cells around already placed sensors. Another way to think of this mechanism is that the program tries to enforce a minimum bound on the sparsity of the sensor network with the result that sensors are more spread out.
#'
#' $suppressionFcn: Specify the suppression function to use. Possible values are 'suppression.static', 'suppression.scale', 'detection.function', 'detection.function.shadow', and 'detection.function.exact'. Default is 'detection.function'.
#'
#' 'suppression.static' will replace all cells within range of a sensor with the value specified by $minsuppressionValue.
#'
#' 'suppression.scale' will multiply the values of cells within range of a sensor by a scaling factor according to the cell's distance from the sensor. Nearby cells receive a higher scaling factor, and more distant cells receive a lower scaling factor. The scaling factor is linearly related to the distance between the sensor and cell. Minimum and maximum values for scaling factors are specified by $minsuppressionValue and $maxsuppressionValue.
#'
#' 'detection.function' will use the inverse of the detection function (that is 'one minus the detection function') to down scale goodness of grid cells near the sensor, however it does not account for objects that block signals. Grid cells' goodness will increase as a function of distance to the sensor.
#'
#' 'detection.function.shadow' will use the inverse of the detection function to down scale goodness of grid cells near the sensor, and does take objects that block signal into account. This means that sensors on opposite sides of a blocking wall will not affect each other's goodness. Unblocked grid cells' goodness will increase as a function of distance to the sensor.
#'
#' 'detection.function.exact' the above suppression functions do not recalculate the goodness grid after placing a new sensor and therefore only provide an approximately optimal solution. This suppression function provides an exact solution by iteratively recalculating the goodness after each sensor placement, and is therefore slower by a factor equal to the number of sensors.
#'
#' $suppressionRangeFactor: Specifies the range of suppression in multiples of the detection range. Valid values are non-negative real numbers. Default is 2. A suppression factor of 2 enforces a suppression range of two times the detection range. Cells within the suppression range of a sensor will be subject to the specified suppression function. To deter sensor detection overlap, a good value to use is equal to or higher than the sensor detection range specified above.
#' @param params A dictionary (list) of parameters, see below for more details.
#' @param showPlots If TRUE plots are shown on the screen, if FALSE plots are stored in the img folder.
#' @param debug If enabled, turns on debug printing (console only).
#' @param save.inter If TRUE intermediary calculations are output as key inter.
#' @param multi If set to TRUE, uses multicore to parallelize calculations.
#' @param silent If set to TRUE, the program will not print status updates for LOS calculation (which may take a very long time).
#' @return A dictionary (list) containing the following return objects:
#'
#' $stats: A dictionary (list) containing summary statistics of the designed network.
#'
#' $stats$delta: The network sparsity is a measure of sensor closeness. If the sparsity is smaller than 1, the network mostly has detection functions that overlap, whereas a sparsity is larger than 1 implies a sparser network with mostly non-overlapping detection functions. Thus, for sparsity < 1 the spatial density of sensors is high, which will make the network suited for estimating detailed animal movements. Sparser arrays (sparsity > 1) will result in higher uncertainty of location estimates, but will in turn have larger spatial extents for a fixed number of sensors thus improving the ability of the network to detect extreme movements.
#'
#' $stats$absRecoveryRate: The absolute recovery rate is the expected total number of detections relative to emitted signals of the network. The absolute recovery rate can exceed a value of 1. This happens if a fish often is detected simultaneously on multiple sensors.
#'
#' $stats$uniqRecoveryRate: The unique recovery rate is an indicator of the overall performance of the network. In short, the unique recovery rate is calculated as the expected proportion of all emitted signals that would be detected by at least one sensor. The unique recovery rate is bounded between zero and one ranging from no coverage to perfect coverage.
#'
#' $stats$uniqRRs: A vector of length equal to the number of sensors indicating the unique recovery rate after each sensor is placed. This is useful to see the improvement in unique recovery rate when placing additional sensors.
#'
#' $stats$coverageGrid: A matrix indicating the acoustic coverage of the study area represented as the probability of at least one sensor detecting a fish located in a given cell.
#'
#' $stats$sensorMat: A matrix containing the indices found by the algorithm where sensors should be placed in the grid matrix.
#' 
#' $runTime: Clock time spent designing the network.
#'
#' $params: parameters used to design network. This also contains parameters that were not necessarily defined by the user and therefore were assigned default values.
#'
#' $errors: errors encountered during run time.
#'
#' $filenames: list of filenames where output was stored.
#'
#' $topographyGrid: A matrix containing the topography of the study area.
#'
#' $behaviorGrid: A matrix highlighting the concentration of animals. This is a result of the selected animal movement model.
#'
#' $goodnessGrid: A matrix containing in each cell the goodness of placing the first sensor in that cell. This is useful to get an impression of potential shadowing effects if params$bias = 2 or 3.
#' @examples 
#' result <- acousticTest() # check the code of acousticTest for details
#' @export
acousticRun <- function(params, showPlots=FALSE, debug=FALSE, save.inter=FALSE, silent=FALSE, multi=FALSE) {
    startTime = Sys.time()
    print(save.inter)
    if(debug) {
        cat("\n[acousticRun]\n")
    }
    if (!exists("acousticStatus", where = -1, mode = "any",inherits = TRUE)) {
        acousticStatus <<- {}
    }
    acousticErrors <<- {}
    topographyGrid = {}
    behaviorGrid = {}
    goodnessGrid = {}
    sensors = {}
    sensors$goodnessGrid = {}
    sensors$sensorList = {}
    sensors$inter = {}
    statDict = {}
    results = {}
    filenames = {}
    
    tryCatch({
        if(!silent) {
            cat("\n  Checking params\n")
        }
        params = checkParams(params=params)
        if(!silent) {
            cat("\n  Parameters checked!\n")
        }
	
        ## Create/Load the Bathy grid for the area of interest
        topographyGrid = getBathy(params$inputFile, params$inputFileType, params$startX, params$startY, 
            params$XDist, params$YDist, params$seriesName, debug)
        topographyGrid = list("topographyGrid"=topographyGrid, "cellRatio"=params$cellSize)
        ## Convert parameter values from meters to number of grid cell 
        params = convertMetersToGrid(params,topographyGrid)
        ## Specify a standard scale of x and y axes if previously undefined
        if(!("x" %in% names(topographyGrid))) {
            topographyGrid$x = (1:dim(topographyGrid$topographyGrid)[1])*params$cellSize 
        }
        if(!("y" %in% names(topographyGrid))) {
            topographyGrid$y = (1:dim(topographyGrid$topographyGrid)[2])*params$cellSize
        }
        ## Calculate fish grid
        print('Calculate fish grid')
        behaviorGrid = fish(params, topographyGrid)
	
        ## Find good sensor placements
        sensors <- sensorFun(params$numSensors, topographyGrid, behaviorGrid, params$range, params$bias, params, debug, save.inter=save.inter, silent=silent, multi=multi)
	
        ## Stat analysis of proposed setup.
        statDict = getStats(params, topographyGrid, behaviorGrid, sensors, debug)
        
        ## Return Fish grid, Bathy grid, and Sensor Placements as a Dictionary.
        results = list("topographyGrid" = topographyGrid, "behaviorGrid" = behaviorGrid, "goodnessGrid"=sensors$goodnessGrid, "sensors" = sensors$sensorList, 
            "stats" = statDict, "params"=params, "errors"=acousticErrors[toString(params$timestamp)])
        
        if(save.inter) {
            print("Saving intermediary calculations in key inter")
            results$inter = sensors$inter
        }
        endTime = Sys.time()
        results$runTime = endTime - startTime
        
        ## Graph results and make data file.
        results$filenames = graph(results,params,showPlots=showPlots, debug=debug)
        print(names(results))
        invisible(results)
        
    }, error = function(e) {
        print("Error")
        print(e)
        appendError(e, toString(params$timestamp))
    }, finally = function(e){}  )
    
    ## only params and errors should actually have values
    results = list("topographyGrid" = topographyGrid, "behaviorGrid" = behaviorGrid, "goodnessGrid"=sensors$goodnessGrid, "sensors" = sensors$sensorList, "stats" = statDict, "inter" = sensors$inter, "filenames"=filenames, "params"=params, "errors"=acousticErrors[toString(params$timestamp)])
    
    ## writeFiles returns json and txt file locations, but only if showPlots == FALSE
    if(!showPlots) results$filenames = writeFiles(filenames, results, path="", as.numeric(params$timestamp), showPlots=showPlots, zip=FALSE, debug)
    ##print(results$filenames)
    ## Return results invisibly (don't print to screen if unassigned because they are usually very long)
    invisible(results)
}

#' @name acousticTest
#' @title Executes a test run of the program, using default parameters.
#' @description Executes a test run of the program, using default parameters.  No additional 
#' parameters are necessary. The code for this function can be used as a template for new projects.
#' @param bias Determines whether to account for species behavior and/or detection shadows when designing network. Choose between bias 1 (behavior only), 2 (shadowing only) or 3 (behavior and shadowing).
#' @param real If TRUE real topographical data are DOWNLOADED (from ftp://ftp.soest.hawaii.edu/pibhmc/website/data/pria/bathymetry/Pal_IKONOS.zip) and used for analysis, if FALSE a made-up topography is used.
#' @param exact If TRUE use exact calculations (slow because goodness grid is updated with line of sight after each sensor placement), if FALSE use approximate calculations (faster).
#' @param multi If TRUE use multicore package to speed up calculations, if FALSE don't. This is useful if your CPU has multiple cores, however it will probably slow down all other processes on your system.
#' @param showPlots If TRUE plots are shown on the screen, if FALSE plots are stored in the img folder.
#' @param paper If TRUE runs analysis that is shown as example in the scientific paper describing this method. If exact=TRUE example 1 is run, otherwise example 2 is run. Note: to run these examples real topography must be downloaded. The user will be asked at runtime if this is OK.
#' @param debug If enabled, turns on debug printing (console only).
#' @param silent If set to TRUE, the program will not print status updates for LOS calculation (which may take a very long time).
#' @return A dictionary of return objects, see ?acousticRun for all details.
#' @export
acousticTest <- function(bias=3, real=FALSE, exact=FALSE, multi=FALSE, showPlots=TRUE, paper=FALSE, silent=FALSE, debug=FALSE) {
	acousticStatus <<- {}
	#### TEST RUN
	params = list()

        if(paper){
            params = list()
            params$timestamp = 00
            params$numSensors = 6
            params$bias = 3
            params$sensorElevation = 1
            params$shapeFcn = "shape.gauss"
            params$peak = 0.98
            params$detectionRange <- 120
            ##params$inputFile = "pal_dbmb.asc"
            ##params$inputFileType = "asc"
            params$inputFile = "pal_dbmb.RData"
            params$inputFileType = "RData"
            params$seriesName = 'bath'
            if(!file.exists(params$inputFile)){
                answ <- readline('Did not find required topography file! do you want to download it? [y/n] ')
                if(answ=='yes' | answ=='y' | answ=='Y'){
                    print("Downloading real topography...")
                    dest <- 'Pal_IKONOS.zip'
                    download.file(url='ftp://ftp.soest.hawaii.edu/pibhmc/website/data/pria/bathymetry/Pal_IKONOS.zip',destfile=dest)
                    unzip(zipfile=dest)
                    filename = "pal_dbmb.asc"
                    bath = loadASC(filename)
                    save(bath,file=params$inputFile)
                } else {
                    stop('Stopping, cannot run paper examples without downloading real topography data!')
                }
            }
            params$cellSize = 5 ## meters
            ##params$maxsuppressionValue = 1
            ##params$minsuppressionValue = 0.5
            params$fishmodel <- "ou"
            params$mux <- 0.5
            params$muy <- 0.5
            params$ousdx <- 300
            params$ousdy <- 300
            params$oucor <- 0
            ## Example 1+2
            params$startX = 290
            params$XDist = 230
            params$startY = 1125
            params$YDist = 240
            params$depth_off_bottom <- 0.5
            params$depth_off_bottom_sd <- 1.5

            if(exact){
                ## --- Scenario 1 ---
                cat("  Running paper scenario 1. (set exact=FALSE if you want to run scenario 2)\n")
                cat("  \n")
                params$suppressionFcn = "detection.function.exact"
            } else {
                ## --- Scenario 2 ---
                cat("  Running paper scenario 2. (set exact=TRUE if you want to run scenario 1)\n")
                params$suppressionRangeFactor = 3
                params$suppressionFcn = "detection.function"
            }
        } else {
            ## Sensor variables
            params$timestamp = 0
            params$numSensors = 4
            params$bias = bias
            params$sensorElevation = 1
            params$shapeFcn = 'shape.gauss'
            params$peak = .98
            params$detectionRange <- 20
            
            ## topographyGrid Variables
            if(real){
                params$inputFile = "pal_dbmb.RData"
                params$inputFileType = "RData"
                params$seriesName = 'bath'
                if(!file.exists(params$inputFile)){
                    print("Downloading real topography...")
                    dest <- 'Pal_IKONOS.zip'
                    download.file(url='ftp://ftp.soest.hawaii.edu/pibhmc/website/data/pria/bathymetry/Pal_IKONOS.zip',destfile=dest)
                    unzip(zipfile=dest)
                    filename = "pal_dbmb.asc"
                    bath = loadASC(filename)
                    save(bath,file=params$inputFile)
                }
            }
            params$cellSize = 5
            params$startX = 290
            params$XDist = 30
            params$startY = 1125
            params$YDist = 30
            
            ## Suppression Variables
            if(exact){
                params$suppressionFcn = "detection.function.exact"
            } else {
                params$suppressionRangeFactor = 1
                params$suppressionFcn = "detection.function"
            }

            ## Behaviour variables
            params$fishmodel <- 'ou'
            params$mux <- .3 
            params$muy <- .3 
            params$ousdx <- 45
            params$ousdy <- 45
            params$oucor <- 0

            ## Depth preference
            if(real){
                params$depth_off_bottom <- 3
                params$depth_off_bottom_sd <- 3
            }
        }
	return(acousticRun(params=params, showPlots=showPlots, silent=silent, debug=debug, multi=multi))
}

#' @name appendError
#' @title Appends execution errors to a global acousticErrors dictionary, with the 'time' variable as the key.
#' @description Creates an entry in the global 'acousticErrors' dictionary with the 'time' value as a key, and
#' msg as a value, then prints msg to the terminal.
#' @param msg The error message to print.
#' @param time A timestamp to use as an identifier for associating related outputs.
#' @return The error message that was passed in.
appendError = function(msg, time) {
	acousticErrors[toString(time)] <<- msg[1]
	print(acousticErrors[toString(time)])
}
##results = acousticTest(showPlots=FALSE, silent=TRUE)
