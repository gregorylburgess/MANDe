rm(list=ls()) 
source("src/Utility.R")

#todo
#check Params
#update fGrid

#tells the program that warnings are errors
options(warn=2)
# global size for Grids created by resetGrids()
r=7
c=5
# global error tolerance
tolerance = 1e-10

#' Resets the test parameters to a default state.
#' 
#' @return Returns a dictionary of reset parameters.
#' @export
resetParams <- function () {
	return(checkParams(list(numSensors=0, shapeFcn="shape.gauss", bias=1, range=1, sd=1, peak=.75, timestamp=1)))
	
}

#' Resets the test Grids to a default state.
#' 
#' @return Returns a dictinary of reset grids containing the keys 'bGrid', 'fGrid', and 'sumGrid'.
#' @export
resetGrids = function () {
	# Use a non-square grid to ensure that columns and rows
	# are being correctly referenced
	status <<- {}
	fGrid = matrix(c(1:35), 
			nrow=r, 
			ncol=c) 
	bGrid = -1 * fGrid
	sumGrid = {}
	bGrid = list(bGrid=bGrid, cellRatio=1)
	grids = list(bGrid=bGrid, fGrid=fGrid, sumGrid=sumGrid)
	return(grids)
}

#' Tests the sumGrid function and all of its bias options.
#' 
#' @return N/A
#' @export
TestUtility.sumGrid <- function (debug=FALSE) {
	sumGrid = {}
	
	na = matrix(rep(NA,35), nrow=r, ncol=c)
	zero = matrix(rep(0,35), nrow=r, ncol=c) 
	
	solGrid = list(
	    s1 = matrix(c(0, 9.28503479100948, 11.2849196505924, 13.2848045101754, 15.2846893697583, 17.2845742293412, 13.4510883353333, 
						16.7541167440698, 24.7943045429025, 27.5492272698917, 30.3041499968809, 33.05907272387, 35.8139954508592, 
						27.243350166755, 30.7533107611504, 44.0787636318267, 46.8336863588159, 49.5886090858051, 52.3435318127942, 
						55.0984545397834, 41.2425441838356, 44.752504778231, 63.3632227207509, 66.1181454477401, 68.8730681747293, 
						71.6279909017184, 74.3829136287076, 55.2417382009161, 38.8129620667501, 54.7112807156447, 56.7111655752276, 
						58.7110504348105, 60.7109352943935, 62.7108201539764, 46.4274257702694),
	            nrow=r, 
	            ncol=c),
		
	    s2 = matrix(c(0, 1.67509226284183, 1.70191680542725, 1.72337643949558, 1.74093432191513, 1.75556589059808, 1.03713887305137, 
						1.64061736071105, 2.51857211434967, 2.56211423021532, 2.61729597058962, 2.65566171177554, 2.68420475826093, 
						1.81994479280308, 1.80966619908725, 2.73927863823349, 2.75183759213253, 2.76262146462847, 2.77200465668237, 
						2.78025963086771, 1.84498103254947, 1.85641173840996, 2.80006540752516, 2.80542792693555, 2.81031642748013, 
						2.81479395567453, 2.81891247919612, 1.85970823240029, 1.14389589552535, 1.87644554097683, 1.87804824518218, 
						1.87952766444866, 1.88089749710281, 1.88216948456737, 1.18570557044753),
	            nrow=r, 
	            ncol=c),
		
	    s3 = matrix(c(0, 8.94867932186257, 10.8651924679732, 12.7817056140838, 14.6982187601943, 16.6147319063049, 10.753417380961, 
					15.7870307546273, 23.5500045060423, 26.3393134371621, 29.2625437398143, 32.1857740424666, 35.1090043451188, 
					23.8930998228565, 28.7927918249904, 43.8786952530755, 46.8019255557277, 49.7251558583799, 52.6483861610321, 
					55.5716164636844, 37.3086918456306, 42.2083838477645, 64.3413073716411, 67.2645376742933, 70.1877679769455, 
					73.1109982795977, 76.03422858225, 50.7242838684048, 29.03899356579, 49.4867517643161, 51.4032649104267, 
					53.3197780565373, 55.2362912026479, 57.1528043487585, 35.6532343603592),
	   			nrow=r, 
	            ncol=c)
	)
    # for each bias option,
    for( i in 1:length(solGrid)) {
		params = resetParams()
		params$bias = i
		grids = resetGrids()
		
		
		# Calculate the sumgrid of the Grids provided by resetGrids() and the given bias
        sumGrid[[i]] = sumGridFun(grids, params$range, i, params, debug, silent=TRUE)$sumGrid
        # compare the results to our solGrid
        if (isTRUE(all.equal(as.character(solGrid[[i]]), as.character(sumGrid[[i]], tolerance=tolerance)))) {
            print(sprintf("[SumGrid:bias %g]: Pass", i))
        }
		# print any inconsistencies
        else {
			print("Expected:")
			print(solGrid[[i]])
			print("Recieved:")
			print(sumGrid[[i]])
            stop(sprintf("[SumGrid bias=%g]: FAIL", i))
        }
		
		
		# Calculate the sumgrid of a grid full of NAs and a grid full of zeroes
		invalidGrids = list(na=na, zero=zero)
		gridNames = names(invalidGrids)
		k = 1
		for(grid in invalidGrids) {
			params = resetParams()
			grids = resetGrids()
			# set bGrid and fGrid to either a grid of NAs or a grid of zeroes
			grids$fGrid = grid
			grids$bGrid$bGrid = grid
			# call the function with obviously bad inputs
			result = sumGridFun(grids, params$range, i, params, debug=FALSE, silent=TRUE)$sumGrid
			# we should get back a grid of all zeroes for both cases.
			if (isTRUE(all.equal(zero, result, tolerance=tolerance))) {
				print(sprintf("[SumGrid:bias %g - %s ]: Pass", i, gridNames[[k]]))
			}
			else {
				print("Expected:")
				print(zero)
				print("Recieved:")
				print(result)
				stop(sprintf("[SumGrid bias=%g - %s]: FAIL", i, gridNames[[k]]))
			}
			k = k + 1
		}
		
		
    }

}


#' Tests the zeroOut function and the getArea function.
#' 
#' @return N/A
#' @export
TestUtility.suppress.scale <- function() {
    tests = list(
        list(d=0, min=.25, max=.75, ans=.25),
        list(d=100, min=.25, max=.75, ans=.75),
        list(d=25, min=.25, max=.75, ans=.375),
        list(d=0, min=0, max=1, ans=0),
        list(d=100, min=0, max=1, ans=1),
        list(d=25, min=0, max=1, ans=.25)
    )
    suppressionRange = 100
    params = {}
    error = FALSE
	
    for( test in tests) {
        val = suppression.scale(test$d, suppressionRange, test$min, test$max, params)
        if (val != test$ans) {
			stop(sprintf("Error: [suppress.scale] incorect value.  Expected %g, recieved %g",test$ans, val ))
			error = TRUE
			return()
		}
    }
	if (!error) {
		print("[suppress.scale]: Pass")
	}
}


#' Tests the getCells functions.
#' 
#' @param opt If TRUE, tests the vectorized version, else tests the unvectorized version.
#' @return N/A
#' @export
TestUtility.getCells <- function() {
	grids = resetGrids()
    params = resetParams()
	cells = {}
	# Test the conditions where slope is negative, positive, infinite, and the start and end cells are the same.
    startCells = list(
					list(r=1,c=1),
					list(r=1,c=1),
					list(r=5,c=2)
				)
				
    endCells = list(
					list(r=5,c=2),
					list(r=1,c=1),
					list(r=1,c=1)
				)
				
	points = list(  
					matrix(c(1,1,2,2,2,2,3,3,4,5),
							nrow=5, 
							ncol=2
				 	),
					NULL,
					matrix(c(1,1,1,2,2,1,2,3,3,4),
							nrow=5, 
							ncol=2
				)
			)
	#call getCells.opt()
	for (i in 1:length(startCells)) {
	    result = (getCells.opt(startCells[[i]], endCells[[i]], debug=FALSE))
			#validate result size
			if(length(points[[i]]) != length(result)) {
				print("Expected:")
				print(points[[i]])
				print("Recieved")
				print(result)
				stop(sprintf("[getCells %g]: FAIL: missing or duplicate cells! ", i))
			}
			# validate result content
			for(point in points[i]) {
				if(length(which(result[,1]==point[,1] & result[,2]==point[,2])) == 1) {
					print("Expected:")
					print(point)
					print("Recieved")
					print(result)
					stop(sprintf("[getCells %g]: FAIL: missing cells! ", i))
				}
			}
	}
	print("[getCells.opt]: Pass")
}

#' Tests the getArea function.
#' 
#' @param opt If TRUE, tests the vectorized version, else tests the unvectorized version.
#' @return N/A
#' @export
TestUtility.getArea = function () {
	grid = matrix(c(1:100),nrow=20,ncol=5)
	dims = dim(grid)
	
	vals = list(# top left corner, bottom Right corner, top right corner, bottom left corner, center
				list(loc=(list(r=1,c=1)), range=1), 
				list(loc=(list(r=20,c=5)), range=1),
				list(loc=(list(r=1,c=5)), range=1),
				list(loc=(list(r=20,c=1)), range=1),
				list(loc=(list(r=10,c=3)), range=1)
		
			)
	sols = list(
				list(rs=1, re=2, cs=1, ce=2),
				list(rs=19, re=20, cs=4, ce=5),
				list(rs=1, re=2, cs=4, ce=5),
				list(rs=19, re=20, cs=1, ce=2),
				list(rs=9, re=11, cs=2, ce=4)
			)
	i = 1
	for (test in vals) {
		result = getArea(test$loc, dims, test$range, debug=FALSE)
		for(key in names(sols[[i]])) {
			if(result[[key]] != sols[[i]][[key]]) {
				print("Expected:")
				print(sols[[i]])
				print("Recieved:")
				print(result)
				stop("[TestUtility.getArea: FAIL")
			}
		}
		i = i + 1
	}
	print("[getArea]: Pass")
}

#' Tests the getArea function.
#' 
#' @param opt If TRUE, tests the vectorized version, else tests the unvectorized version.
#' @return N/A
#' @export
TestUtility.offset = function() {
	vals = list(
			list(r=1,c=1), 
			list(r=-1,c=1),
			list(r=1,c=-1),
			list(r=-1,c=-1),
			list(r=0,c=0)
	)
	
	sols = list(
			list(r=.5,c=.5), 
			list(r=-.5,c=.5),
			list(r=.5,c=-.5),
			list(r=-.5,c=-.5),
			list(r=.5,c=.5)
	)
	i = 1
	for (test in vals) {
		result = offset(test)
		for(key in names(sols[[i]])) {
			if(result[[key]] != sols[[i]][[key]]) {
				print("Expected:")
				print(sols[[i]])
				print("Recieved:")
				print(result)
				stop(sprintf("[TestUtility.offset %i]: FAIL",i))
			}
		}
		i = i + 1
	}
	print("[offset]: Pass")
}

TestUtility.calc.percent.viz = function() {
	sensorElevation = 1
	# our test cell can see to the corner of each 5x5 test grid
	range = 3
	# our test cell is in the center of each 5x5 test grid
	params = resetParams()
	fGrid = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), nrow=5, ncol=5)
	loc = list(c=3,r=3)
	params$params$depth_off_bottom = 4
	params$depth_off_bottom_sd = 1
	
	# each grid is 5x5
	vals = 	list(
			# a cell surrounded by a wall
			c1 = matrix(c(5,5,5,5,5,5,0,0,0,5,5,0,5,0,5,5,0,0,0,5,5,5,5,5,5), nrow=5, ncol=5),
			
			# a cell on a step-cliff (wall on one side)
			c2 = matrix(c(5,5,5,5,5,4,4,4,4,4,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0), nrow=5, ncol=5),
			
			# a cell on a hill
			c3 = matrix(c(5,5,5,5,5,5,5,5,5,5,5,5,2,5,5,5,5,5,5,5,5,5,5,5,5), nrow=5, ncol=5),
			
			# a cell on a hill with one diagonal obstruction
			c4 = matrix(c(5,5,5,5,5,5,5,5,5,5,5,5,2,5,5,5,0,5,5,5,5,5,5,5,5), nrow=5, ncol=5),
	
			# a cell surrounded by a low wall
			c5 = matrix(c(5,5,5,5,5,5,4,4,4,5,5,4,5,4,5,5,4,4,4,5,5,5,5,5,5), nrow=5, ncol=5)
	)
	
	sols = list(
				c1 = list(
						 percentVisibility = c(2e-05, 2e-05, 2e-05, 2e-05, 2e-05, 2e-05, 2e-05,
								 2e-05, 2e-05, 2e-05, 2e-05, 2e-05, 2e-05, 2e-05, 2e-05, 2e-05),
					  	  dists = c(2.8284271247462, 2.8284271247462, 2.8284271247462, 2.8284271247462, 
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998, 
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998, 
								 2.0000000000000, 2.0000000000000, 2.0000000000000, 2.0000000000000),
					 	  linearIndex = c(1, 5, 21, 25, 2, 4, 6, 10, 16, 20, 22, 24, 3, 11, 
								 15, 23)
					 ),
				c2 = list(
						 percentVisibility = c(1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000,
								 1.0000000, 1.0000000, 1.0000000, 1.0000000, 0.85355339059327, 
								 0.85355339059327, 1.0000000, 1.0000000, 1.0000000), 
						  dists = c(2.8284271247462, 2.8284271247462, 2.2360679774998, 2.2360679774998,
								 2.2360679774998, 2.2360679774998, 2.0000000000000, 2.0000000000000,
								 2.0000000000000, 1.4142135623731, 1.4142135623731, 1.0000000000000,
								 1.0000000000000, 1.0000000000000),
						  linearIndex = c( 1, 5, 2, 4, 6, 10, 3, 11, 15, 7, 9, 8, 12, 14)
			  		 ),
				c3 = list(
						 percentVisibility = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
								 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
						  dists = c(2.8284271247462, 2.8284271247462, 2.8284271247462, 2.8284271247462,
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
								 2.0000000000000, 2.0000000000000, 2.0000000000000, 2.0000000000000, 
								 1.4142135623731, 1.4142135623731, 1.4142135623731, 1.4142135623731,
								 1.0000000000000, 1.0000000000000, 1.0000000000000, 1.0000000000000),
						  linearIndex = c(1, 5, 21, 25, 2, 4, 6, 10, 16, 20, 22, 24,
								 3, 11, 15, 23, 7, 9, 17, 19, 8, 12, 14, 18)
			   		 ),
			   c4 = list(
					    percentVisibility = c(1e+00, 1e+00, 2e-05, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00,
							     2e-05, 1e+00, 2e-05, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00,
								 1e+00, 1e+00, 1e+00, 1e+00, 1e+00),
						dists = c(2.8284271247462, 2.8284271247462, 2.8284271247462, 2.8284271247462, 
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
								 2.0000000000000, 2.0000000000000, 2.0000000000000, 2.0000000000000,
								 1.4142135623731, 1.4142135623731, 1.4142135623731, 1.0000000000000,
								 1.0000000000000, 1.0000000000000, 1.0000000000000),
						linearIndex = c(1, 5, 21, 25, 2, 4, 6, 10, 16, 20, 22, 24,
								 3, 11, 15, 23, 7, 9, 19, 8, 12, 14, 18)
					 ),
			  
			   c5 = list(
					    percentVisibility = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
								 0.8, 0.8, 0.8, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
					    dists = c(2.8284271247462, 2.8284271247462, 2.8284271247462, 2.8284271247462,
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
								 2.0000000000000, 2.0000000000000, 2.0000000000000, 2.0000000000000,
								 1.4142135623731, 1.4142135623731, 1.4142135623731, 1.4142135623731,
								 1.0000000000000, 1.0000000000000, 1.0000000000000, 1.0000000000000),
					    linearIndex = c(1, 5, 21, 25, 2, 4, 6, 10, 16, 20, 22, 24, 3, 11, 15, 23, 7,
								 9, 17, 19, 8, 12, 14, 18)
			   ),
			   c6 = list(
					    percentVisibility = c(1.5910485584147e-08, 1.5910485584147e-08,
						 	     1.5910485584147e-08, 1.5910485584147e-08, 1.5910485584147e-08,
							     1.5910485584147e-08, 1.5910485584147e-08, 1.5910485584147e-08,
							     1.5910485584147e-08, 1.5910485584147e-08, 1.5910485584147e-08,
							     1.5910485584147e-08, 1.5910485584147e-08, 1.5910485584147e-08,
							     1.5910485584147e-08, 1.5910485584147e-08),
					    dists = c(2.8284271247462, 2.8284271247462, 2.8284271247462, 2.8284271247462, 
							     2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998, 
							     2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998, 
							     2.0000000000000, 2.0000000000000, 2.0000000000000, 2.0000000000000),
					    linearIndex = c(1, 5, 21, 25, 2, 4, 6, 10, 16, 20, 22, 24, 3, 11, 
							     15, 23)
			   ),
			   c7 = list(
					    percentVisibility = c(1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000,
							     1.0000000, 1.0000000, 1.0000000, 1.0000000, 0.78487534586011, 
							 	 0.78487534586011, 1.0000000, 1.0000000, 1.0000000), 
					    dists = c(2.8284271247462, 2.8284271247462, 2.2360679774998, 2.2360679774998,
							 	 2.2360679774998, 2.2360679774998, 2.0000000000000, 2.0000000000000,
							 	 2.0000000000000, 1.4142135623731, 1.4142135623731, 1.0000000000000,
							 	 1.0000000000000, 1.0000000000000),
					    linearIndex = c( 1, 5, 2, 4, 6, 10, 3, 11, 15, 7, 9, 8, 12, 14)
			   ),
			   c8 = list(
					    percentVisibility = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
								 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
					    dists = c(2.8284271247462, 2.8284271247462, 2.8284271247462, 2.8284271247462,
							 	 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
							 	 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
							 	 2.0000000000000, 2.0000000000000, 2.0000000000000, 2.0000000000000, 
							 	 1.4142135623731, 1.4142135623731, 1.4142135623731, 1.4142135623731,
							 	 1.0000000000000, 1.0000000000000, 1.0000000000000, 1.0000000000000),
					    linearIndex = c(1, 5, 21, 25, 2, 4, 6, 10, 16, 20, 22, 24,
							 	 3, 11, 15, 23, 7, 9, 17, 19, 8, 12, 14, 18)
			   ),
			   c9 = list(
					    percentVisibility = c(1.0000000000000e+00, 1.0000000000000e+00, 
							     1.5910485584147e-08, 1.0000000000000e+00, 1.0000000000000e+00,
							     1.0000000000000e+00, 1.0000000000000e+00, 1.0000000000000e+00,
							     1.5910485584147e-08, 1.0000000000000e+00, 1.5910485584147e-08,
							     1.0000000000000e+00, 1.0000000000000e+00, 1.0000000000000e+00,
							     1.0000000000000e+00, 1.0000000000000e+00, 1.0000000000000e+00,
							     1.0000000000000e+00, 1.0000000000000e+00, 1.0000000000000e+00,
							     1.0000000000000e+00, 1.0000000000000e+00, 1.0000000000000e+00),
					    dists = c(2.8284271247462, 2.8284271247462, 2.8284271247462, 2.8284271247462, 
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
								 2.0000000000000, 2.0000000000000, 2.0000000000000, 2.0000000000000,
								 1.4142135623731, 1.4142135623731, 1.4142135623731, 1.0000000000000,
							 	 1.0000000000000, 1.0000000000000, 1.0000000000000),
					    linearIndex = c(1, 5, 21, 25, 2, 4, 6, 10, 16, 20, 22, 24,
							 	 3, 11, 15, 23, 7, 9, 19, 8, 12, 14, 18)
			   ),
			   c10 = list(
					    percentVisibility = c(0.59427143559031, 0.59427143559031, 0.59427143559031,
								 0.59427143559031, 0.59427143559031, 0.59427143559031, 0.59427143559031,
								 0.59427143559031, 0.59427143559031, 0.59427143559031, 0.59427143559031,
								 0.59427143559031, 0.59427143559031, 0.59427143559031, 0.59427143559031,
								 0.59427143559031, 1.00000000000000, 1.00000000000000, 1.00000000000000,
								 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000,
								 1.00000000000000),
					    dists = c(2.8284271247462, 2.8284271247462, 2.8284271247462, 2.8284271247462, 
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
								 2.2360679774998, 2.2360679774998, 2.2360679774998, 2.2360679774998,
								 2.0000000000000, 2.0000000000000, 2.0000000000000, 2.0000000000000,
								 1.4142135623731, 1.4142135623731, 1.4142135623731, 1.4142135623731,
								 1.0000000000000, 1.0000000000000, 1.0000000000000, 1.0000000000000),
					    linearIndex = c(1, 5, 21, 25, 2, 4, 6, 10, 16, 20, 22, 24, 3, 11, 15, 23, 7,
								 9, 17, 19, 8, 12, 14, 18)
			   )
				
	)

	#Result has the following keys: percentVisibility, dists, linearIndex

	j = 1
	for (dpflag in c(FALSE, TRUE)) {
		params$dpflag = dpflag
		for (bGrid in vals) {
			bGrid = bGrid * -1
			land = bGrid >= 0
			sensorDepth = bGrid + sensorElevation
			cells = getArea(loc, dim(bGrid), range)
			rind = cells$rs:cells$re
			cind = cells$cs:cells$ce
			result = calc.percent.viz(loc$r, loc$c, rind, cind, bGrid, land, sensorDepth[loc$r,loc$c], dpflag, params, debug=FALSE)
			for(key in c("percentVisibility", "dists", "linearIndex")) {
				#print(!all(result[[key]] == sols[[j]][[key]]))
				options(digits=14)
				if(!isTRUE(all.equal(result[[key]], sols[[j]][[key]], tolerance=tolerance))) {
					print(key)
					print(j)
					print("Expected:")
					print(sols[[j]][[key]])
					print("Recieved:")
					print(result[[key]])
					print("Diff:")
					print(result[[key]]-sols[[j]][[key]])
					stop(sprintf("[TestUtility.calc.percent.viz %i]: FAIL",j))
				}
			}
			j = j+1
		}
	}
	print("[calc.percent.viz]: Pass")
}




#' Runs a battery of tests.
#' 
#' @param opt If TRUE, tests the vectorized version, else tests the unvectorized version.
#' @return N/A
#' @export
runTests = function () {
	TestUtility.sumGrid()
	TestUtility.getCells()
	TestUtility.suppress.scale()
	TestUtility.getArea()
	TestUtility.offset()
	TestUtility.calc.percent.viz()
	print("SUCCESS! All Tests Pased.")
}


#Executes the tests, and stops if an error occurs
tryCatch({
	runTests()
	}, warning = function(e) {
	traceback()
	print(e)
	}, error = function(e) {
	traceback()
	print(e)
	}, finally = function(e){})
