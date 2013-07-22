\name{var.add.ncdf}
\alias{var.add.ncdf}
\title{Add New netCDF Variable to Existing File}
\description{
 Special purpose routine for adding a new variable to a netCDF file
 that already exists on disk.
}
\usage{
 var.add.ncdf( nc, v, verbose=FALSE, indefine=FALSE )
}
\arguments{
 \item{nc}{The already-existing netCDF file we want to add a new variable to. This must
 be a value of class "ncdf" returned by a call to open.ncdf(...,write=TRUE).}
 \item{v}{The variable to be added to the file.  This must be 
 a value of class "var.ncdf" returned by a call to var.def.ncdf().}
 \item{verbose}{If true, prints diagnostic messages.}
 \item{indefine}{If true, the file is assumed to already be in define mode.}
}
\value{
 A handle to the netCDF file that describes the newly modified file.
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 There are two cases in which you might want to add a variable to a netCDF
 file.  The first, and most common way, is when you are creating a brand
 new netCDF file.
 Usually when you create a netCDF file, you specify what variables you want
 the file to contain.
 This is the method
 most users will use to make netCDF files.  To do this, you do NOT use this
 routine; instead, pass a list of the variables you wish to have created in the output
 file to routine create.ncdf().

 The second, less common, case is when you already have an existing netCDF
 file on disk and wish to add a new variable to it.  In that case, use
 this routine.  First define the variable you want to add to the existing
 file using routine var.def.ncdf(); then add it to the already-existing
 and opened (for writing) netCDF file using this routine.  (This routine
 automatically creates any additional dimensions that are needed in the
 output file to handle the new variable.)

 NOTE that the return value of this routine should replace the old netCDF
 file handle that you were using.  This newly returned value reflects the
 modifications to the file that were accomplished by calling this routine.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf]{dim.def.ncdf}}, \code{\link[ncdf]{create.ncdf}}, 
 \code{\link[ncdf]{var.def.ncdf}}.
}
\examples{
#===========================================================================
# PART 1.  MAKE A TEST NETCDF FILE THAT WE WILL ADD A VARIABLE TO IN PART 2.
#===========================================================================

#----------------
# Make dimensions
#----------------
xvals <- 1:360
yvals <- -90:90

nx <- length(xvals)
ny <- length(yvals)

xdim <- dim.def.ncdf( 'Lon', 'degreesE', xvals )
ydim <- dim.def.ncdf( 'Lat', 'degreesE', yvals )
tdim <- dim.def.ncdf( 'Time', 'days since 1900-01-01', 0, unlim=TRUE )

#---------
# Make var
#---------
mv <- 1.e30     # missing value
var_temp <- var.def.ncdf( 'Temperature', 'K', list(xdim,ydim,tdim), mv )

#---------------------
# Make new output file
#---------------------
output_fname <- 'test_real3d.nc'
ncid_new <- create.ncdf( output_fname, list(var_temp))

#-------------------------------
# Put some test data in the file
#-------------------------------
data_temp <- array(0.,dim=c(nx,ny,1))
for( j in 1:ny )
for( i in 1:nx )
        data_temp[i,j,1] <- sin(i/10)*sin(j/10)

put.var.ncdf( ncid_new, var_temp, data_temp, start=c(1,1,1), count=c(nx,ny,1))

#--------------------------
# Close our new output file
#--------------------------
close.ncdf( ncid_new )

#===========================================================================
# PART 2.  ADD A NEW VARIABLE TO THE FILE
#===========================================================================

#---------------------------------------------------
# Open the existing file we're going to add a var to
#---------------------------------------------------
ncid_old <- open.ncdf( output_fname, write=TRUE )

#------------------------------------------------------------------
# Make a NEW variable to put into the file.  Have this new variable
# use the same dimensions already in the file
#------------------------------------------------------------------
xdim2 <- ncid_old$dim[['Lon']]
ydim2 <- ncid_old$dim[['Lat']]
tdim2 <- ncid_old$dim[['Time']]
mv2 <- 1.e30
var_q <- var.def.ncdf( 'Humidity', 'g/kg', list(xdim2,ydim2,tdim2), mv2 )

ncid_old <- var.add.ncdf( ncid_old, var_q )	# NOTE this returns a modified netcdf file handle 

#-------------------------------------------------------------------
# Make a DIFFERENT new var that will be added to the file.  This var
# uses a dim that does NOT already exist in the file.
#-------------------------------------------------------------------
zdim   <- dim.def.ncdf( 'Level', 'hPa', seq(1000,100,by=-100))
var_cf <- var.def.ncdf( 'CloudFraction', 'percent', list(xdim2,ydim2,zdim,tdim2), mv2 )

ncid_old <- var.add.ncdf( ncid_old, var_cf )

print(ncid_old)

close.ncdf( ncid_old )
}
\keyword{utilities}
