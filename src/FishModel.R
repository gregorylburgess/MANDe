# Models Fish Behavior.  Outputs a Fish grid with the same dimensions as the Bathymetry grid,
# containing the percentage of transmissions sent per cell.
library(mvtnorm)

fish <- function(params, bGrid) {
    rows <- dim(bGrid$bGrid)[1]
    cols <- dim(bGrid$bGrid)[2]
    land <- bGrid$bGrid>=0
    switch(params$fishmodel,
            rw={ ## Random walk case
                print('rw')
                fGrid <- matrix(1,rows,cols)
                if('mindepth' %in% names(params) & 'maxdepth' %in% names(params)){
                    print('RW: Using vertical habitat to calculate fGrid')
                    fGrid[!verticalHabitat(params$mindepth,params$maxdepth,bGrid$bGrid)] <- 0
                }
            },
            ou={ ## Ornstein-Uhlenbeck case
                print('ou')
                sigma <- params$msd/sqrt(params$dt)
                B <- matrix(c(params$By,params$Bxy,params$Bxy,params$Bx),2,2)
                hrCov <- 0.5*sigma^2*solve(B)
                X <- matrix(rep(bGrid$x,rows),rows,cols,byrow=TRUE)
                Y <- matrix(rep(bGrid$y,cols),rows,cols,byrow=FALSE)
                XY <- cbind(as.vector(X),as.vector(Y))
                hrVals <- dmvnorm(XY,c(params$muy,params$mux),hrCov)
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

verticalHabitat <- function(mindepth,maxdepth,bGrid) bGrid < mindepth & bGrid > maxdepth
