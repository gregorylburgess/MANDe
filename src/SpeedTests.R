source("src/Utility.R")
# Use a non-square grid to ensure that columns and rows
# are being correctly referenced
r=5
c=2
topographyGrid = matrix(c(1,2,3,4,5,6,7,8,9,10), 
		nrow=r, 
		ncol=c) 
cellRatio = 1
behaviorGrid = topographyGrid
topographyGrid = list(topographyGrid=topographyGrid, cellRatio=cellRatio)
grids = list(topographyGrid=topographyGrid, behaviorGrid=behaviorGrid)
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

SpeedTests.goodnessGrid.sumSimple <- function() {
	print('Old fun:')
	ng <- 500
	grid <- list(behaviorGrid=matrix(1:ng^2,ng,ng))
	at <- system.time(a <- goodnessGrid.sumSimple(grid, 'behaviorGrid', params$range, debug=FALSE)$goodnessGrid)
	print(at)
	print('Optimized fun:')
	bt <- system.time(b <- goodnessGrid.sumSimple.opt(grid, 'behaviorGrid', params$range, debug=FALSE)$goodnessGrid)
	print(bt)
}

SpeedTests.getCells.opt()
SpeedTests.goodnessGrid.sumSimple()