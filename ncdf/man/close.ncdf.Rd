\name{close.ncdf}
\alias{close.ncdf}
\title{Close a netCDF File}
\description{
 Closes an open netCDF file, which flushes any unwritten data to disk.
}
\usage{
 close.ncdf( con, ... )
}
\arguments{
 \item{con}{An object of class \code{ncdf} (as returned by either function
 \code{open.ncdf()} or function \code{create.ncdf()}, indicating what file to read from.}
 \item{...}{Other arguments passed to or from other methods.}
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 Data in a netCDF file might be cached in memory, for better performance.
 This data is written out when the file is closed.  Therefore, always
 remember to close the file when done with it.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf]{sync.ncdf}}.
}
\examples{
\dontrun{nc <- open.ncdf("salinity.nc")}
\dontrun{data <- get.var.ncdf( nc )	# Read the "only" var in the file}
\dontrun{close.ncdf(nc)}
}
\keyword{utilities}
