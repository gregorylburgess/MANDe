# Models Fish Behavior.  Outputs a Fish grid with the same dimensions as the Bathymetry grid,
# containing the percentage of transmissions sent per cell.
library(mvtnorm)

#' Generates a Fish Location Grid (FGrid) for the program to use.  Values in the cells of this grid 
#' are expressed as a percentage of the total number of fish on the grid.
#' 
#' @param params A dictionary of parameters, see PARAMETER_DESCRIPTIONS.html for more info.  
#' @param bGrid A valid BGrid.
#' @return An FGrid of the same dimensions as the provided BGrid.
#' @export
fish <- function(params, bGrid) {
    rows <- dim(bGrid$bGrid)[1]
    cols <- dim(bGrid$bGrid)[2]
    land <- bGrid$bGrid>=0
    switch(params$fishmodel,
            rw={ ## Random walk case
                fGrid <- matrix(1,rows,cols)
                if('mindepth' %in% names(params) & 'maxdepth' %in% names(params)){
                    print('RW: Using vertical habitat to calculate fGrid')
                    fGrid[!verticalHabitat(params$mindepth,params$maxdepth,bGrid$bGrid)] <- 0
                }
            },
            ou={ ## Ornstein-Uhlenbeck case
                mux <- min(bGrid$x) + diff(range(bGrid$x))*params$mux
                muy <- min(bGrid$y) + diff(range(bGrid$y))*params$muy
                varx <- params$ousdx^2
                vary <- params$ousdy^2
                covxy <- params$oucor * params$ousdx * params$ousdy
                hrCov <- matrix(c(vary,covxy,covxy,varx),2,2)
                ##sigma <- params$msd/sqrt(params$dt)
                ##B <- matrix(c(params$By,params$Bxy,params$Bxy,params$Bx),2,2)
                ##hrCov <- 0.5*sigma^2*solve(B)
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
    fGrid[land] <- 0 ## Set land areas to zero
    fGrid <- fGrid/sum(fGrid) ## Make sure fGrid sums to one
    return (fGrid)
}

#' A helper function that indicates if a cell has a depth between mindepth and maxdepth.  
#' 
#' @param mindepth The shallowest depth a fish will visit.
#' @param maxdepth The deepest depth a fish will visit.
#' @param bGrid A valid BGrid.
#' @return TRUE if the cell is between the depth restrictions, FALSE otherwise.
#' @export
verticalHabitat <- function(mindepth,maxdepth,bGrid) bGrid < mindepth & bGrid > maxdepth
