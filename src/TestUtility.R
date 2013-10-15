rm(list=ls()) 
source("src/Utility.R")

r=7
c=5

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
    s1 = matrix(c(0, 9.28503479100948, 11.2849196505924, 13.2848045101754, 15.2846893697583, 17.2845742293412, 13.4510883353333, 
					16.7541167440698, 24.7943045429025, 27.5492272698917, 30.3041499968809, 33.05907272387, 35.8139954508592, 
					27.243350166755, 30.7533107611504, 44.0787636318267, 46.8336863588159, 49.5886090858051, 52.3435318127942, 
					55.0984545397834, 41.2425441838356, 44.752504778231, 63.3632227207509, 66.1181454477401, 68.8730681747293, 
					71.6279909017184, 74.3829136287076, 55.2417382009161, 38.8129620667501, 54.7112807156447, 56.7111655752276, 
					58.7110504348105, 60.7109352943935, 62.7108201539764, 46.4274257702694),
            nrow=r, 
            ncol=c)
	
    s2 = matrix(c(0, 1.67509226284183, 1.70191680542725, 1.72337643949558, 1.74093432191513, 1.75556589059808, 1.03713887305137, 
					1.64061736071105, 2.51857211434967, 2.56211423021532, 2.61729597058962, 2.65566171177554, 2.68420475826093, 
					1.81994479280308, 1.80966619908725, 2.73927863823349, 2.75183759213253, 2.76262146462847, 2.77200465668237, 
					2.78025963086771, 1.84498103254947, 1.85641173840996, 2.80006540752516, 2.80542792693555, 2.81031642748013, 
					2.81479395567453, 2.81891247919612, 1.85970823240029, 1.14389589552535, 1.87644554097683, 1.87804824518218, 
					1.87952766444866, 1.88089749710281, 1.88216948456737, 1.18570557044753),
            nrow=r, 
            ncol=c)
s3 = matrix(c(0, 8.94867932186257, 10.8651924679732, 12.7817056140838, 14.6982187601943, 16.6147319063049, 10.753417380961, 
				15.7870307546273, 23.5500045060423, 26.3393134371621, 29.2625437398143, 32.1857740424666, 35.1090043451188, 
				23.8930998228565, 28.7927918249904, 43.8786952530755, 46.8019255557277, 49.7251558583799, 52.6483861610321, 
				55.5716164636844, 37.3086918456306, 42.2083838477645, 64.3413073716411, 67.2645376742933, 70.1877679769455, 
				73.1109982795977, 76.03422858225, 50.7242838684048, 29.03899356579, 49.4867517643161, 51.4032649104267, 
				53.3197780565373, 55.2362912026479, 57.1528043487585, 35.6532343603592),
   			nrow=r, 
            ncol=c)
    solGrid = list(s1=s1, s2=s2, s3=s3)
    sumGrid = {}
    for( i in 1:length(solGrid)) {
		params = resetParams()
		params$bias = i
		grids = resetGrids()
        sumGrid[[i]] = sumGridFun(grids, params$range, i, params, debug)$sumGrid
        
        if (isTRUE(all.equal(as.character(solGrid[[i]]), as.character(sumGrid[[i]])))) {
            print(sprintf("[SumGrid:bias %g]: Pass", i))
        }
        else {
			print("Expected:")
			print(solGrid[[i]])
			print("Recieved:")
			print(sumGrid[[i]])
			#print("Diff")
			#print(sumGrid[[i]]-solGrid[[i]])
			#print(sum(sumGrid[[i]])-sum(solGrid[[i]]))
			#print(split(as.character(sumGrid[[i]]), rep(1:ncol(sumGrid[[i]]), each = nrow(sumGrid[[i]]))))
            stop(sprintf("[SumGrid bias=%g]: FAIL", i))
        }
    }

}


#' Tests the zeroOut function and the getArea function.
#' 
#' @return N/A
#' @export
TestUtility.zeroOut <- function () {
    params = resetParams()
    dims = dim(fGrid)
    loc = list(r=dims[1], c=dims[2])
    value = 0
    # zeroing out the grid around the bottom right corner
    # should result in the bottom two cells being '0'
    solGrid = matrix(c(1,2,3,0,0,6,7,8,0,0),
            nrow=r, 
            ncol=c)
    
    fGrid = zeroOut(fGrid, dims, loc, range, value)
    if (isTRUE(all.equal(fGrid, solGrid))) {
        print("[zeroOut: %s]: Pass")
    }
    else {
		stop(sprintf("[zeroOut: %s]: FAIL",i))
        print("solGrid:")
        print(solGrid)
        print("result:")
        print(fGrid)
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
    params = resetParams()
    startCell = list(r=1,c=1)
    endCell = list(r=5,c=2)
	optLabel = ""

    cells = (getCells.opt(startCell, endCell, debug=FALSE))
	colnames(cells) <- c("x","y")
	cells = as.data.frame(cells)
    points = list(  list(x=1,y=2),
                    list(x=1,y=3),
                    list(x=2,y=3),
                    list(x=2,y=4),
                    list(x=2,y=5))
    if (dim(cells)[1] > 5) {
		stop("[getCells]: FAIL: too many results, duplicates exist!")
    }
	print(cells)
    for (point in points) {
        if(dim(subset(cells, x==point$x & y==point$y))[1] != 1) {
			stop(sprintf("[getCells]: FAIL: (%g,%g) missing or duplicate! ",point$x,point$y))
        }
    }
    print(sprintf("[getCells.opt]: Pass"))
}

runTests = function () {
	TestUtility.sumGrid()
	TestUtility.getCells()
TestUtility.suppress.scale()
}

tryCatch({
	runTests()
}, error = function(e) {
	traceback()
	print("Error")
	print(e)
	}, finally = function(e){})
