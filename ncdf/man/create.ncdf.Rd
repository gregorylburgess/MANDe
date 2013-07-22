\name{create.ncdf}
\alias{create.ncdf}
\title{Create a netCDF File}
\description{
 Creates a new netCDF file, given the variables the new file is to contain. 
}
\usage{
 create.ncdf( filename, vars, verbose=FALSE )
}
\arguments{
 \item{filename}{Name of the  netCDF file to be created.}
 \item{vars}{Either an object of class \code{var.ncdf} describing the variable to be created, or a vector of such objects to be created.}
 \item{verbose}{If TRUE, then information is printed while the file is being created.}
}
\value{
 An object of class \code{ncdf}, which has the fields described in \code{\link[ncdf]{open.ncdf}}.
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 This routine creates a new netCDF file on disk.  It must be given the variables in 
 the file that will be created.
 Keep in mind that the new file may not actually be written to disk until 
 \code{\link[ncdf]{close.ncdf}} is called.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf]{dim.def.ncdf}}, \code{\link[ncdf]{var.def.ncdf}}. 
}
\examples{
# Define an integer dimension 
dimState <- dim.def.ncdf( "StateNo", "count", 1:50 )

# Make an integer variable.  Note that an integer variable can have
# a double precision dimension, or vice versa; there is no fixed
# relationship between the precision of the dimension and that of the
# associated variable.  We just make an integer variable here for
# illustration purposes.
varPop <- var.def.ncdf("Pop", "count", dimState, -1, 
		longname="Population", prec="integer")

# Create a netCDF file with this variable
ncnew <- create.ncdf( "states_population.nc", varPop )

# Write some values to this variable on disk.
popAlabama <- 4447100
put.var.ncdf( ncnew, varPop, popAlabama, start=1, count=1 )

close.ncdf(ncnew)
}
\keyword{utilities}
