## Script for converting from asc format to R
## An asc file containing bathymetry data from Palmyra is located in this zip file: ftp://ftp.soest.hawaii.edu/pibhmc/website/data/pria/bathymetry/Pal_IKONOS.zip
## Place the .asc file in the data folder
## 07.06.2013
rm(list=ls())
filename <- '../data/pal_dbmb.asc'

asd <- scan(file=filename,what='character',nmax=2,sep='')
ncols <- as.numeric(asd[2])
asd <- scan(file=filename,what='character',skip=1,nmax=2,sep='')
nrows <- as.numeric(asd[2])
asd <- scan(file=filename,what='character',skip=2,nmax=2,sep='')
xll <- as.numeric(asd[2])
asd <- scan(file=filename,what='character',skip=3,nmax=2,sep='')
yll <- as.numeric(asd[2])
asd <- scan(file=filename,what='character',skip=4,nmax=2,sep='')
dx <- as.numeric(asd[2])
asd <- scan(file=filename,what='character',skip=5,nmax=2,sep='')
nodata <- as.numeric(asd[2])

dat <- scan(file=filename,skip=6)

dat[dat==nodata] <- NA
bath <- matrix(dat,nrows,ncols,byrow=TRUE)
bath <- bath[nrows:1,]
bath <- t(bath)

## From meta data file
##Spatial_Domain:
## Bounding_Coordinates:
## West_Bounding_Coordinate: -162.177002
## East_Bounding_Coordinate: -161.995141
## North_Bounding_Coordinate: 5.903546
## South_Bounding_Coordinate: 5.855663
lat <- seq(5.855663,5.903546,length=nrows)
lon <- seq(-162.177002, -161.995141,length=ncols)

indlon <- lon > -162.125 & lon < -162.118
indlat <- lat > 5.873 & lat < 5.876

bathstu <- bath[indlon,indlat]
min(bathstu[!is.na(bathstu)])
max(bathstu[!is.na(bathstu)])

bathstu[which(is.na(bathstu))] <- bathstu[which(is.na(bathstu))+1] ## Remove a NA
lonstu <- lon[indlon]
latstu <- lat[indlat]

filled.contour(lonstu,latstu,bathstu)

bathlat <- lat
bathlon <- lon

## Save to compressed R-file
##save(bathlat,bathlon,bath,file='palmyrabath.RData')

##load(file='palmyrabath.RData')
##filled.contour(bathlon,bathlat,bath)
