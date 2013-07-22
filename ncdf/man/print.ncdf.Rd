\name{print.ncdf}
\alias{print.ncdf}
\title{Print Information About a netCDF File}
\description{
 Prints information about a netCDF file, including the variables and dimensions
 it contains.
}
\usage{
 \method{print}{ncdf}( x, ... )
}
\arguments{
 \item{x}{An object of class "ncdf".}
 \item{...}{Extra arguments are passed to the generic print function.}
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 NetCDF files contain variables, which themselves have dimensions.  
 This routine prints out useful information about a netCDF file's
 variables and dimensions.  It is overloaded on the regular print
 function, so if "nc" is an object of class "ncdf", then just
 calling \code{print(nc)} will suffice.  Objects of class
 "ncdf" are returned from \code{open.ncdf}.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf]{var.def.ncdf}}
}
\examples{
# Open a netCDF file, print information about it
nc <- open.ncdf( "salinity.nc" )
print(nc)
}
\keyword{utilities}
