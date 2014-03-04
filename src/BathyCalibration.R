
input50m = "src/himbsyn.bathy.v19.grd/himbsyn.bathy.v19.grd"
input1km ="src/himbsyn.bathytopo.1km.v19.grd/himbsyn.bathytopo.1km.v19.grd"
inputFileType = "netcdf"
XDist= 1
YDist= 1
seriesName="z"
timestamp=-1
debug=FALSE

sim50m = function() {
	nrow=41
	ncol=3
	probs = matrix(c(3441.900, 7477.267, 338,
					3459.533, 7481.867, 380,
					3417.600, 7435.233, 376,
					9348.800, 7614.467, 330,
					9555.567, 8061.900, 325,
					7817.300, 8841.700, 200,
					8018.067, 9217.500, 182,
					8005.700, 9305.867, 118,
					7910.267, 9174.500, 253,
					7716.967, 8529.967, 290,
					8674.700, 8818.433, 338,
					8474.933, 7946.400, 343,
					9095.567, 8545.533, 320,
					8699.267, 8298.767, 320.6,
					8529.600, 8218.500, 327.8,
					8836.700, 8417.533, 302.5,
					9013.633, 11481.567, 313,
					8006.100, 9306.000, 116,
					7910.267, 9173.467, 256,
					7978.967, 9289.800, 279,
					7817.033, 8841.800, 220,
					8005.500, 9306.533, 116,
					7716.867, 8529.867, 285,
					7978.967, 9289.800, 279,
					8524.433, 8811.367, 100,
					8851.067, 8610.400, 115,
					8511.967, 8862.467, 377,
					8878.000, 8612.733, 354,
					9287.500, 7637.367, 330,
					9588.900, 8061.900, 325,
					8524.900, 8811.367, 100,
					8997.567, 8506.833, 313,
					8474.733, 7946.400, 343,
					9095.567, 7454.467, 320,
					8529.600, 7780.500, 330,
					8699.267, 7701.233, 322,
					8836.700, 7582.467, 316,
					8468.667, 11626.933, 348,
					7910.267, 10825.500, 261,
					7916.967, 11470.033, 285,
					7978.967, 10710.200, 279), 
			 		nrow=ncol,ncol=nrow)
	
	rslt = matrix(1,nrow=nrow, ncol=1)
	for(i in 1:dim(probs)[2]) {
		startY = ceiling(probs[1,i])
		startX = ceiling(probs[2,i])
		cell = getBathy(input50m, inputFileType, startX, startY, XDist, YDist, seriesName, timestamp, debug)
		rslt[i] = abs(cell[1] + probs[3,i])
	}
	print(rslt)
}

sim1km = function() {
	nrow=41
	ncol=3
	probs = matrix(c(172.095, 373.863, 338,
					172.977, 374.093, 380,
					170.880, 371.762, 376,
					467.440, 380.723, 330,
					477.778, 403.095, 325,
					390.865, 442.085, 200,
					400.903, 460.875, 182,
					400.285, 465.293, 118,
					395.513, 458.725, 253,
					385.848, 426.498, 290,
					433.735, 440.922, 338,
					423.747, 397.320, 343,
					454.778, 427.277, 320,
					434.963, 414.938, 320.6,
					426.480, 410.925, 327.8,
					441.835, 420.877, 302.5,
					450.682, 574.078, 313,
					400.305, 465.300, 116,
					395.513, 458.673, 256,
					398.948, 464.490, 279,
					390.852, 442.090, 220,
					400.275, 465.327, 116,
					385.843, 426.493, 285,
					398.948, 464.490, 279,
					426.222, 440.568, 100,
					442.553, 430.520, 115,
					425.598, 443.123, 377,
					443.900, 430.637, 354,
					464.375, 381.868, 330,
					479.445, 403.095, 325,
					426.245, 440.568, 100,
					449.878, 425.342, 313,
					423.737, 397.320, 343,
					454.778, 372.723, 320,
					426.480, 389.025, 330,
					434.963, 385.062, 322,
					441.835, 379.123, 316,
					423.433, 581.347, 348,
					395.513, 541.275, 261,
					395.848, 573.502, 285,
					398.948, 535.510, 279),
					nrow=ncol,ncol=nrow)
					
					rslt = matrix(1,nrow=nrow, ncol=1)
					for(i in 1:dim(probs)[2]) {
						startY = ceiling(probs[1,i])
						startX = ceiling(probs[2,i])
						print(probs[1,i])
						print(probs[2,i])
						print("====")
						cell = getBathy(input1km, inputFileType, startX, startY, XDist, YDist, seriesName, timestamp, debug)
						rslt[i] = abs(cell[1] + probs[3,i])
					}
					print(rslt)
}

#' @name testplot
#' @title Generates a topography Graph.
#' @description Creates a test plot to help improve targeting on a netcdf dataset.
#' @param x The longitude of a test point.
#' @param y The latitude of a test point.
#' @param myn The test Northern limit of the area to graph .
#' @param mys The test Southern limit of the area to graph.
#' @param mye The test Eastern limit of the area to graph.
#' @param myw The test Western limit of the area to graph.
#' @param timestamp A unique string/integer to prevent graph file overwrites.
testplot = function(x, y, myn, mys, mye, myw, timestamp, ...) {
	    #Data specific to the 1km grid file.
		inputFile = "src/himbsyn.bathytopo.1km.v19.grd/himbsyn.bathytopo.1km.v19.grd" 	#path to the 1km v19 bathy data grid file
		inputFileType = "netcdf" 		#tells getBathy that this is a netcdf file 
		seriesName = 'z' 		#netcdf variable name
		n = 25; 			#Northern limit of the netcdf file
		s = 17; 			# Southern limit of hte netcdf file
		e = -153; 		#Eastern limit of the netcdf file
		w = -162;		#Western limit of the netcdf file
		dpcC = .00999;		#Degrees per cell in the Y direction, as given by the bathy dataset
		dpcR = .00999;		#Degrees per cell in the X direction, as given by the bathy dataset
		
		startX = floor((myw-w)/dpcC)  		#Western most x index (in grid cells) for the area of interest
		XDist = floor((mye-myw)/dpcC)		#the extent (in grid cells) for the area of interest in the x direction
		startY = floor((mys-s)/dpcR)    		#Southern most y index (in grid cells) for the area of interest
		YDist = floor((myn-mys)/dpcR)  		#the extent (in grid cells) for the area of interest in the y direction
		xoffset = 2*dpcC
		yoffset= 2*dpcR
		r=(y + yoffset - mys)/dpcR 		#the y coordinate of the test point
		c=(x + xoffset -myw)/dpcC	# the x coordinate of the test point
		
		#Print out the calulation of the test point's x/y coordinates
		print(paste("r=(",y," - ",mys,")/",dpcR,"=",r))
		print(paste("c=(",x," - ",myw,")/",dpcC,"=",c))
		
		print(paste("startx=",startX))
		print(paste("starty=",startY))
	    # Create/Load the Bathy grid for the area of interest
		library(ncdf)
		ncdfObj = open.ncdf(inputFile)		# open the netCDF file
		topographyGrid = get.var.ncdf(ncdfObj, seriesName, start=c(startX, startY), count=c( XDist, YDist))		# grab a slice (in grid form)
		
	    # Specify a standard scale of x and y axes if previously undefined
		xaxis = (1:dim(topographyGrid)[1])		#Create a range of numbers from 1 to the number of columns in the area of interest
		yaxis = (1:dim(topographyGrid)[2])		#Create a range of numbers from 1 to the number of columns in the area of interest
		
		#graphing stuff
		col = colorRampPalette(c("navy","white"))(24)
		png(paste("topographyGrid-",timestamp, ".png"))
		image(x=xaxis, y=yaxis, z=topographyGrid, zlim=c(-6000,0),xlab='x',ylab='y',col=col,...)
		
		contour(xaxis,yaxis,topographyGrid,add=TRUE,nlevels=5)		# Add bathymetry contour
		sensx = list(c)		# Graphing library expects a vector containng x values
		sensy = list(r)		# Graphing library expects a vector containng y values
		points(sensx,sensy,pch=20,bg='red',cex=3)		# Plot the test point
		dev.off()
}
#plot a test point around Big Island
x= -155.003
y=19.333
myn=20
mys=19
mye=-154.66
myw=-156.0
timestamp="BI"
testplot(x, y, myn, mys, mye, myw, timestamp)

#plot a test point around the western most tip of Oahu (Kaena point)
x= -158.28
y=21.575
myn=21.76
mys=21.13
mye=-157.7
myw=-158.39
timestamp="Oahu"
testplot(x, y, myn, mys, mye, myw, timestamp)

#plot a test point around the south western most tip of Niihau (pueo point)
x=-160.072
y=21.896
myn=22.02
mys=21.767
myw=-160.3
mye=-160.02
timestamp="NI"
testplot(x, y, myn, mys, mye, myw, timestamp)

#plot a test point around a south western sea mount
x=-161.335
y=17.718
myn=17.78
mys=17.7
myw=-161.35
mye=-161.257
timestamp="SW"
testplot(x, y, myn, mys, mye, myw, timestamp)

#plot a test point around a south eastern sea mount
x=-154.083
y=17.147
myn=17.267
mys=17.032
myw=-154.192
mye=-153.942
timestamp="SE"
testplot(x, y, myn, mys, mye, myw, timestamp)

#plot a test point around a north eastern sea mount
x=-153.212
y=24.885
myn=24.937
mys=24.825
myw=-153.281
mye=-153.151
timestamp="NE"
testplot(x, y, myn, mys, mye, myw, timestamp)

#plot a test point around a north western sea mount
x=-161.902
y=24.556
myn=24.61
mys=24.51
myw=-161.93
mye=-161.83
timestamp="NW"
testplot(x, y, myn, mys, mye, myw, timestamp)



x=-161.8
y=25

myn=25
mys=17.1
myw=-161.9
mye=-153
timestamp="OB"
testplot(x, y, myn, mys, mye, myw, timestamp)


