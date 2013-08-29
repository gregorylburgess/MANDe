#' @name fish
#' @title Generates a Fish Location Grid (FGrid) for the program to use.  
#' @description Values in the cells of this grid are expressed as a percentage of 
#' the total number of fish on the grid.
#' 
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.  
#' @param bGrid A valid BGrid.
#' @return An fGrid of the same dimensions as the provided bGrid.
fish <- function(params, bGrid) {
    require(mvtnorm)
    rows <- dim(bGrid$bGrid)[1]
    cols <- dim(bGrid$bGrid)[2]
    land <- bGrid$bGrid>=0
    switch(params$fishmodel,
            ## Random walk case
            rw={ 
                print("Using RW model")
                fGrid <- matrix(1,rows,cols)
                if('mindepth' %in% names(params) & 'maxdepth' %in% names(params)){
                    print('RW: Using vertical habitat to calculate fGrid')
                    fGrid[!verticalHabitat(params$mindepth,params$maxdepth,bGrid$bGrid)] <- 0
                }
            },
            ## Ornstein-Uhlenbeck case
            ou={
                print("Using OU model")
                mux <- min(bGrid$x) + diff(range(bGrid$x))*params$mux
                muy <- min(bGrid$y) + diff(range(bGrid$y))*params$muy
                varx <- params$ousdx^2
                vary <- params$ousdy^2
                covxy <- params$oucor * params$ousdx * params$ousdy
                hrCov <- matrix(c(vary,covxy,covxy,varx),2,2)
                Y <- matrix(rep(bGrid$y,rows),rows,cols,byrow=TRUE)
                X <- matrix(rep(bGrid$x,cols),rows,cols,byrow=FALSE)
                XY <- cbind(as.vector(X),as.vector(Y))
                hrVals <- dmvnorm(XY,c(mux,muy),hrCov)
                fGrid <- matrix(hrVals,rows,cols,byrow=FALSE)
                if('mindepth' %in% names(params) & 'maxdepth' %in% names(params)){
                    print('OU: Using vertical habitat to calculate fGrid')
                    fGrid[!verticalHabitat(params$mindepth,params$maxdepth,bGrid$bGrid)] <- 0
                }
            }
    )
    ## Set land areas to zero
    fGrid[land] <- 0
    ## Make sure fGrid sums to one
    fGrid <- fGrid/sum(fGrid) 
    return (fGrid)
}

#' @name verticalHabitat
#' @title A vectorized helper function that indicates if a cell has a depth 
#' between mindepth and maxdepth.  
#' @description Returns boolean grid.  Cells between the given depths are set to TRUE,
#' and all other cells are set to FALSE.
#' @param mindepth The shallowest depth a fish will visit.
#' @param maxdepth The deepest depth a fish will visit.
#' @param bGrid A valid BGrid.
#' @return A boolean grid with values set to TRUE if the cell is between the depth restrictions,
#' and FALSE otherwise.
verticalHabitat <- function(mindepth,maxdepth,bGrid) {
	bGrid < mindepth & bGrid > maxdepth
}
