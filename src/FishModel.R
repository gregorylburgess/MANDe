#' Generates the Behaviour Grid

#' @name fish
#' @title Generates a Fish Location Grid (behaviorGrid) for the program to use.  
#' @description Values in the cells of this grid are expressed as a percentage of 
#' the total number of fish on the grid.
#' 
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.  
#' @param topographyGrid A valid topographyGrid.
#' @return An behaviorGrid of the same dimensions as the provided topographyGrid.
fish <- function(params, topographyGrid) {
    require(mvtnorm)
    rows <- dim(topographyGrid$topographyGrid)[1]
    cols <- dim(topographyGrid$topographyGrid)[2]
    land <- topographyGrid$topographyGrid>=0
    switch(params$fishmodel,
            ## Random walk case
            rw={ 
                print("Using RW model")
                behaviorGrid <- matrix(1,rows,cols)
            },
            ## Ornstein-Uhlenbeck case
            ou={
                print("Using OU model")
                Y <- matrix(rep(topographyGrid$y,rows),rows,cols,byrow=TRUE)
                X <- matrix(rep(topographyGrid$x,cols),rows,cols,byrow=FALSE)
                XY <- cbind(as.vector(X),as.vector(Y))
                nocenters <- length(params$mux)
                behaviorGrid <- matrix(0,rows,cols)
                for(i in 1:nocenters){
                    mux <- min(topographyGrid$x) + diff(range(topographyGrid$x))*params$mux[i]
                    muy <- min(topographyGrid$y) + diff(range(topographyGrid$y))*params$muy[i]
                    varx <- params$ousdx[i]^2
                    vary <- params$ousdy[i]^2
                    covxy <- params$oucor[i] * params$ousdx[i] * params$ousdy[i]
                    hrCov <- matrix(c(vary,covxy,covxy,varx),2,2)
                    hrVals <- dmvnorm(XY,c(mux,muy),hrCov)
                    hrVals <- hrVals/max(hrVals)
                    behaviorGridtmp <- matrix(hrVals,rows,cols,byrow=FALSE)
                    ## Add contribution from center i to behaviorGrid
                    behaviorGrid <- behaviorGrid + behaviorGridtmp
                }
            }
    )
	if('mindepth' %in% names(params) && 'maxdepth' %in% names(params)){
		behaviorGrid[!verticalHabitat(params$mindepth,params$maxdepth,topographyGrid$topographyGrid)] <- 0
	}
    ## Set land areas to zero
    behaviorGrid[land] <- 0
    ## Make sure behaviorGrid sums to one
    behaviorGrid <- behaviorGrid/sum(behaviorGrid) 
    return (behaviorGrid)
}

#' @name verticalHabitat
#' @title A vectorized helper function that indicates if a cell has a depth 
#' between mindepth and maxdepth.  
#' @description Returns boolean grid.  Cells between the given depths are set to TRUE,
#' and all other cells are set to FALSE.
#' @param mindepth The shallowest depth a fish will visit.
#' @param maxdepth The deepest depth a fish will visit.
#' @param topographyGrid A valid topographyGrid.
#' @return A boolean grid with values set to TRUE if the cell is between the depth restrictions,
#' and FALSE otherwise.
verticalHabitat <- function(mindepth,maxdepth,topographyGrid) {
	topographyGrid < mindepth & topographyGrid > maxdepth
} 
