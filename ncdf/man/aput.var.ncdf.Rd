\name{put.var.ncdf}
\alias{put.var.ncdf}
\title{Write data to a netCDF file}
\description{
 Writes data to an existing netCDF file.  The variable to be written to
 must already exist on disk. 
}
\usage{
 put.var.ncdf( nc, varid, vals, start=NA, count=NA, verbose=FALSE ) 
}
\arguments{
 \item{nc}{An object of class \code{ncdf} (as returned by either function
 \code{open.ncdf()}
 or function \code{create.ncdf()}), indicating what file to write to.}
 \item{varid}{What variable to write the data to.  Can be a string with the name
 of the variable, an object of class \code{var.ncdf}, or the "id" field of
 a \code{var.ncdf} object.}
 \item{vals}{The values to be written.}
 \item{start}{A vector of indices indicating where to start writing the passed values
 (starting at 1).  
 The length of this vector must equal the number of dimensions the variable has.
 Order is X-Y-Z-T (i.e., the time dimension is last).  If not specified, writing 
 starts at the beginning of the file (1,1,1,...).}
 \item{count}{A vector of integers indicating the count of values to write
 along each dimension (order is X-Y-Z-T).  
 The length of this vector must equal the number of dimensions the variable has.
 If not specified and the variable does
 NOT have an unlimited dimension, the entire variable is written.  If the variable
 has an unlimited dimension, this must be specified.  As a special case, the value
 "-1" indicates that all entries along that dimension should be written.}
 \item{verbose}{If true, prints information while executing.}
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 This routine writes data values to a variable in a netCDF file.  The file
 should have either been created with \code{\link[ncdf]{create.ncdf}}, or 
 opened with \code{\link[ncdf]{open.ncdf}} called with parameter \code{write=TRUE}..

 Note that the type of the values written to the file is determined when the
 variable is created; in particular, it does not matter what type you pass
 to this function to be written.  In other words, if the variable was created
 with type 'integer', passing double precision values to this routine will 
 still result in integer values being written to disk.

 Values of "NA" are supported; they are converted to the netCDF variable's
 missing value attribute before being written.  See \code{\link[ncdf]{set.missval.ncdf}}
 for more information.

 Data in a netCDF file is conceived as being a multi-dimensional array.
 The number and length of dimensions is determined when the variable is
 created.  The 'start' and 'count' indices that this routine takes indicate
 where the writing starts along each dimension, and the count of values
 along each dimension to write.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf]{dim.def.ncdf}}, \code{\link[ncdf]{create.ncdf}}, 
 \code{\link[ncdf]{get.var.ncdf}}.
}
\examples{
# Make a few dimensions we can use
nx <- 3
ny <- 4
nt <- 5
xvals <- (1:nx)*100.
dimX <- dim.def.ncdf( "X", "meters", xvals )
dimY <- dim.def.ncdf( "Y", "meters", (1:ny)*100. )
dimT <- dim.def.ncdf( "Time", "seconds", (1:nt)/100., unlim=TRUE )

# Make varables of various dimensionality, for illustration purposes
mv <- 1.e30		# missing value to use
var1d <- var.def.ncdf( "var1d", "units", dimX, mv )
var2d <- var.def.ncdf( "var2d", "units", list(dimX,dimY), mv )
var3d <- var.def.ncdf( "var3d", "units", list(dimX,dimY,dimT), mv )

# Create the test file
nc <- create.ncdf( "writevals.nc", list(var1d,var2d,var3d) )

# Write some data to the file
data1d <- runif(nx)
put.var.ncdf( nc, var1d, data1d )	# no start or count: write all values
put.var.ncdf( nc, var1d, 27.5, start=3, count=1 ) # Write a value to the third slot

data2d <- runif(nx*ny)
put.var.ncdf( nc, var2d, data2d )	# no start or count: write all values
# Write a 1-d slice to the 2d var
put.var.ncdf( nc, var2d, data1d, start=c(1,2), count=c(nx,1) )	
# Note how "-1" in the count means "the whole dimension length", 
# which equals nx in this case
put.var.ncdf( nc, var2d, data1d, start=c(1,3), count=c(-1,1) )	

# The 3-d variable has an unlimited dimension.  We will loop over the timesteps,
# writing one 2-d slice per timestep.
for( i in 1:nt) 
	put.var.ncdf( nc, var3d, data2d, start=c(1,1,i), count=c(-1,-1,1) )

close.ncdf(nc)

#----------------------------------------------------------------------
# Illustrate creating a character type variable
#----------------------------------------------------------------------
cnames   <- c("red", "orange", "green", "yellow", "puce", "colorwithverylongname" )
nstrings <- length(cnames)

#--------------------------------------------------------------
# Make dimensions. Setting "dimnchar" to have a length of 12
# means that the maximum color name
# length can be 12.  Longer names will be truncated to this.
#--------------------------------------------------------------
dimnchar   <- dim.def.ncdf("nchar",   "", 1:12 )
dimcolorno <- dim.def.ncdf("colorno", "", 1:nstrings )

varcolors  <- var.def.ncdf("colors", "", list(dimnchar, dimcolorno), 
			NA, prec="char" )

ncid <- create.ncdf( "colornames.nc", list(varcolors) )

put.var.ncdf( ncid, "colors", cnames )

close.ncdf( ncid )
}
\keyword{utilities}
