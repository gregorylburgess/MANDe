source('src/Main.R')
test = function (timestamp, opt, lonLatUserSensorList) {
	#### TEST RUN
	params = list()
	## Sensor variables
	params$timestamp = timestamp

	#turn n/s/e/w boundaries into start points & distance for area of interest
	# 'Palmyra_40m.grd' specific data
	#maximum boundaries of the netcdf file
	n= 6.040206
	s= 5.772677
	w= -162.357356
	e= -161.916751
	
	
	#the user's area of interest
	myn = 5.9
	mys = 5.8560
	myw = -162.1755
	mye = -161.9961
		
	#degrees per cell for Row and Column ... don't touch!
	dpcR =  0.0003635
	dpcC =  0.000362
	xoffset = 0*dpcC
	yoffset= 0*dpcR
	params$startX = floor((myw-w)/dpcC)  		#Western most x index (in grid cells) for the area of interest
	params$XDist = floor((mye-myw)/dpcC)		#the extent (in grid cells) for the area of interest in the x direction
	params$startY = floor((mys-s)/dpcR)    		#Southern most y index (in grid cells) for the area of interest
	params$YDist = floor((myn-mys)/dpcR)  		#the extent (in grid cells) for the area of interest in the y direction
	
	#strip whitespace
	cleaned = gsub("\\s", "", lonLatUserSensorList)
	#split on commas
	data = strsplit(cleaned, ",")
	userSensorList = ""
	len=floor(length(data[[1]])/2)
	data = data[[1]]
	if(opt) {
		params$numSensors = len
	}
	else {
		params$numSensors=0
		while(len > 0) {
			y=as.numeric(data[2])
			x=as.numeric(data[1])
			r=floor((y + yoffset - mys)/dpcR) + 1	#the y coordinate of the test point
			c=floor((x + xoffset -myw)/dpcC) + 1 # the x coordinate of the test point
			#delete the points from data
			data = data[-2]
			data = data[-1]
			userSensorList = paste(userSensorList,c, ",", r, sep="")
			if(len != 1) {
				userSensorList = paste(userSensorList, ",", sep="")
			}
			# decriment length
			len = len - 1
		}
		print("userList")
		print(userSensorList)
		params$userSensorList = userSensorList
	}
	
	params$bias = 3
	params$sensorElevation = 10
	params$shapeFcn = 'shape.gauss'
	params$peak = .98 
	params$detectionRange = 220
	
	# topographyGrid Variables
	params$inputFile = "src/palmyra_40m.grd"
	params$inputFileType = "netcdf"
	params$seriesName = 'z'
	params$cellSize = 40 
	
	## Suppression Variables
	params$suppressionRangeFactor = 2
	params$suppressionFcn = "suppression.scale"
	## This is only relevant with suppression.scale
	params$maxsuppressionValue = 1
	## This is only relevant with suppression.scale
	params$minsuppressionValue = .5 
	## Choose random walk type movement model
	params$fishmodel = 'rw'
	## Apply vertical habitat range?
	vHabitatRange = FALSE
	if(vHabitatRange){
		## Minimum depth (shallowest depth)
		params$mindepth = -5
		## Maximum depth (deepest depth)
		params$maxdepth = -10
	}
	## Apply depth preference?
	depthPref = FALSE
	if(depthPref){
		## Depth preference of fish relative to bottom (in meters off the bottom)
		params$dp = 2
		## Strength of depth preference as a standard deviation, 95% of the time is spent within plus minus two dpsd
			params$dpsd = 2
	}
	
	print(params)
	return(acousticRun(params=params, showPlots=FALSE))
}



inputFile = "src/PW_VR2LocationsForWebAppSorryGregIOweYouABeer.csv"
con  = file(inputFile, open = "r")
fileContent =  readLines(con, n = -1)
data = (strsplit(fileContent, "\n"))
for (line in data) {
	nameEnd = regexpr(',', line, TRUE)
	stringLen = nchar(line)
	name = substr(line,1,nameEnd-1)
	sensorList = substr(line,nameEnd+1,stringLen)
	print(name)
	print(paste(name, sep=""))
	tryCatch({
		test(timestamp=name, FALSE, sensorList)
		test(timestamp=paste("-",name, sep=""), TRUE, sensorList)
	},
	error = function(e) {
			print(e)
	})
}
