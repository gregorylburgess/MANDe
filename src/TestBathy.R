source('src/Bathy.R')
ncdf = list("ncdf")
arcgis = list("raster")
##arcgis = list("rgdal", "raster", "sp")
prereqs = list(ncdf=ncdf, arcgis=arcgis)

## Check that required packages installed
TestBathy.checkPackages = function() {
	for (prereq in names(prereqs)) {
		for (requirement in prereqs[[prereq]]) {
			if(!(require(requirement, character.only=TRUE))) {
				print(sprintf("Package %s not installed!  Unable to parse %s files!", requirement, prereq))
			}
			else {
				print(sprintf("Package %s was successfully loaded.",requirement))
			}
		}
	}
}

TestBathy.getBathy = function() {
	inputFile = "src/palmyra_40m.grd"
	linkLocation = "ftp://ftp.soest.hawaii.edu/pibhmc/website/data/pria/bathymetry/Palmyra_40m.grd.zip"
	inputFileType = "netcdf"
	startX = 100
	startY = 100
	XDist = 5
	YDist = 5
	seriesName = 'z'
	timestamp = 0
	debug = FALSE
	
	sol = matrix(c(-1915.0788574219, -1896.7038574219, -1885.7969970703, -1878.0434570312,
				   -1874.9871826172, -1892.5827636719, -1881.6988525391, -1869.2602539062,
				   -1861.2454833984, -1852.3382568359, -1873.7113037109, -1866.0177001953,
				   -1851.7595214844, -1845.0474853516, -1840.6657714844, -1859.3164062500,
				   -1851.4396972656, -1839.8825683594, -1827.2858886719, -1817.6904296875,
				   -1842.9913330078, -1841.9722900391, -1832.9112548828, -1820.4696044922,
				   -1806.1882324219),
					 nrow=5,
				 	 ncol=5)

	topographyGrid = getBathy(inputFile, inputFileType, startX, startY, XDist, YDist, seriesName, timestamp, debug)
	
	#The default simulated grid
	fakeGrid = simulatetopographyGrid(XDist, YDist)
	
	#Ensure that the system is indeed pulling the results from our file and not simply
	#simulating them.
	if(isTRUE(all.equal(topographyGrid,fakeGrid))) {
		printError(sprintf("Missing Test File: %s.  Please download it here: %s, unzip it, and place it in the src/ directory.", inputFile, linkLocation))
	}
	
	#check that we got the right data
	if(!isTRUE(all.equal(topographyGrid,sol))){
		print("Expected:")
		print(sol)
		print("Recieved:")
		print(topographyGrid)
		stop("[TestBathy.getBathy %i]: FAIL")
	}
	print("[getBathy]: Pass")
}

TestBathy.simulatetopographyGrid = function() {
	XDist = 5
	YDist = 5
	fakeGrid = simulatetopographyGrid(XDist, YDist)
	sol = c(XDist,YDist)
	if(!(all.equal(dim(fakeGrid), sol))) {
		printError("[simulatetopographyGrid]: FAIL!")
	}
	print("[simulatetopographyGrid]: Pass")
}

#' Runs a battery of tests.
#' 
#' @return N/A
#' @export
runBathyTests = function () {
	print("--- Testing Bathy Functions ---")
	TestBathy.checkPackages()
	TestBathy.simulatetopographyGrid()
	TestBathy.getBathy()
	print("Success! All Bathy Tests Passed.")
}


#Executes the tests, and stops if an error occurs
tryCatch({
			runBathyTests()
		}, warning = function(e) {
			print(e[1])
		}, error = function(e) {
			print(e[1])
		}, finally = function(e){})
