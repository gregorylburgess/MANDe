\name{sync.ncdf}
\alias{sync.ncdf}
\title{Synchronize (flush to disk) a netCDF File}
\description{
 Flushes any pending operations on a netCDF file to disk.
}
\usage{
 sync.ncdf( nc )
}
\arguments{
 \item{nc}{An object of class \code{ncdf} (as returned by either 
 function \code{open.ncdf()}
 or function \code{create.ncdf()}, indicating what file to read from.}
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 Data in a netCDF file might be cached in memory, for better performance.
 An example of when this might be bad is if a long-running job writes one
 timestep of the output file at a time; if the job crashes near the end, 
 the results of many timesteps might be lost.  In such an event, the user
 can manually force any cached data to be written to disk using this call.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\examples{
# The time you would use the sync.ncdf function is when you have an unlimited
# dimension and are writing to the file timestep-by-timestep. Make a netCDF file 
# that has an unlimited dimension for illustration.
nx <- 5
ny <- 8
dimx <- dim.def.ncdf( "X", "meters", 1:nx )
dimy <- dim.def.ncdf( "Y", "meters", 1:ny )
dimt <- dim.def.ncdf( "Time", "days since 1900-01-01", 0, unlim=TRUE )

vartemp <- var.def.ncdf( "Temperature", "degC", list(dimx,dimy,dimt), 1.e30 )
nc  <- create.ncdf( "temperature.nc", vartemp )

nt <- 10  # Imagine this is actually some very large number of timesteps
for( i in 1:nt ) {
	# Long, slow computation to get the data ... for illustration, we just
	# use the following:
	data <- runif(nx*ny)

	# Write the data to this timestep
	put.var.ncdf( nc, vartemp, data, start=c(1,1,i), count=c(nx,ny,1) )

	# Write the time value for this timestep as well
	timeval <- i*10
	put.var.ncdf( nc, dimt, timeval, start=i, count=1 )

	# Flush this timesteps data to the file so we dont lose it
	# if there is a crash or other problem
	sync.ncdf( nc )
	}

# Always remember to close the file when done!!
close.ncdf(nc)
}
\keyword{utilities}
