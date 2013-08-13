source("src/Utility.R")
# Use a non-square grid to ensure that columns and rows
# are being correctly referenced
r=5
c=2
bGrid = matrix(c(1,2,3,4,5,6,7,8,9,10), 
		nrow=r, 
		ncol=c) 
cellRatio = 1
fGrid = bGrid
bGrid = list(bGrid=bGrid, cellRatio=cellRatio)
grids = list(bGrid=bGrid, fGrid=fGrid)
range = 1
params = checkParams(list(numSensors=0, shapeFcn="shape.t", range=1, sd=1, peak=.75))



SpeedTests.getCells.opt <- function(){
	startingCell <- list(r=3,c=3)
	targetCell <- list(r=7,c=7)
	print('Old fun:')
	print(system.time(for(i in 1:10000) getCells(startingCell, targetCell)))
	print('Optimized fun:')
	print(system.time(for(i in 1:10000) getCells.opt(startingCell, targetCell)))
}

SpeedTests.sumGrid.sumSimple <- function() {
	print('Old fun:')
	ng <- 500
	grid <- list(fGrid=matrix(1:ng^2,ng,ng))
	at <- system.time(a <- sumGrid.sumSimple(grid, 'fGrid', params$range, debug=FALSE)$sumGrid)
	print(at)
	print('Optimized fun:')
	bt <- system.time(b <- sumGrid.sumSimple.opt(grid, 'fGrid', params$range, debug=FALSE)$sumGrid)
	print(bt)
}

SpeedTests.getCells.opt()
SpeedTests.sumGrid.sumSimple()