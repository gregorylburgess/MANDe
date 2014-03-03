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
                mux <- min(topographyGrid$x) + diff(range(topographyGrid$x))*params$mux
                muy <- min(topographyGrid$y) + diff(range(topographyGrid$y))*params$muy
                varx <- params$ousdx^2
                vary <- params$ousdy^2
                covxy <- params$oucor * params$ousdx * params$ousdy
                hrCov <- matrix(c(vary,covxy,covxy,varx),2,2)
                Y <- matrix(rep(topographyGrid$y,rows),rows,cols,byrow=TRUE)
                X <- matrix(rep(topographyGrid$x,cols),rows,cols,byrow=FALSE)
                XY <- cbind(as.vector(X),as.vector(Y))
                hrVals <- dmvnorm(XY,c(mux,muy),hrCov)
                behaviorGrid <- matrix(hrVals,rows,cols,byrow=FALSE)
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
