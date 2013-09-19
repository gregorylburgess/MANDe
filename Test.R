source('src/Utility.R')
source('src/Bathy.R')


cells=50
myn = 21.923
myw= -158.869
mys=19.974
mye=-156.912

# Verification point
#myn = 21.923
#myw= -158.869
#mys=s+.001
#mye=-154.96


# Data for the 50m dataset
if(cells==50){
	dpc=0.0005
	n= 24  
	s= 18
	e= -154
	w= -161
	noffset=0
	woffset=0
	filename="b.png"
}

# Data for the 1km dataset
if(cells==1){
	# These are the reported limit values for the himbsyn 1km bathy file
	dpc=0.01
	n= 24.5
	s= 17.5
	e= -154.5
	w= -161.5
	# These are the emprical corrections we found
	noffset=.55
	woffset=.45
	filename="a.png"
}



myn =myn +noffset
mys = mys + noffset
myw = myw+woffset
mye = mye+woffset
errors={}
if(myn>n) {
	errors=c(errors, sprintf("Northern value (%g) must be between %g and %g.\n", myn, s, n-noffset))
}
if(mys<s) {
	errors=c(errors, sprintf("Southern value (%g) must be between %g and %g.\n", mys, s, n-noffset))
}
if(mye>=e) {
	errors=c(errors, sprintf("Eastern value (%g) must be between %g and %g.\n", mye, w, e-woffset))
}
if(myw<=w) {
	errors=c(errors, sprintf("Western value (%g) must be between %g and %g.\n", myw, w, e-woffset))
}
if(myn<=mys) {
	errors=c(errors, "Northern value must be greater than Southern Value.\n")
}
if(mye<=myw) {
	errors=c(errors, "Western value must be greater than Eastern Value.\n")
}

if( length(errors) > 0 ){
	cat(errors)
}else {
	startX=ceiling((myw-w)/dpc)
	XDist=ceiling((mye-myw)/dpc)
	startY=ceiling((mys-s)/dpc)
	YDist=ceiling((myn-mys)/dpc)
	
	print(startX)
	print(XDist)
	print(startY)
	print(YDist)
	if(cells ==1){
		bGrid = getBathy("src/himbsyn.bathytopo.1km.v19.grd/himbsyn.bathytopo.1km.v19.grd", "netcdf", startX, startY, 
			XDist, YDist, 'z', FALSE)
	}
	if(cells ==50) {
		bGrid = getBathy("src/himbsyn.bathy.v19.grd/himbsyn.bathy.v19.grd", "netcdf", startX, startY, 
				XDist, YDist, 'z', FALSE)
	}
	options(width=10000)
	capture.output(print(bGrid), file=filename)
	bGrid = list("bGrid"=bGrid, "cellRatio"=cells)
	if(!("x" %in% names(bGrid))) {
		bGrid$x = (1:dim(bGrid$bGrid)[1])*cells
	}
	if(!("y" %in% names(bGrid))) {
		bGrid$y = (1:dim(bGrid$bGrid)[2])*cells
	}
	result={}
	result$bGrid=bGrid
	png(filename)
	plotGrid(result,type='bGrid',"x","y",plot.bathy=FALSE, plot.sensors=FALSE)
	dev.off()
	print("Done")
}


