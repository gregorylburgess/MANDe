\name{redef.ncdf}
\alias{redef.ncdf}
\title{Puts a netCDF file back into define mode}
\description{
 Puts a netCDF that is not currently in define mode back into define mode.
}
\usage{
 redef.ncdf( nc )
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
 NetCDF files can be in "define mode", at which time dimensions and variables can be
 defined, or new attributes added to a file, or in "data mode", at which time data
 can be read from the file.
 This call puts a file that is currently in data mode back into define mode.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf]{enddef.ncdf}}.
}
\examples{
# This function is for advanced useage only, and will never
# be needed by the typical users R code.
}
\note{
 The typical user will never need this call, nor will ever have to worry about
 "define mode" or "data mode".  THIS CALL IS PROVIDED FOR ADVANCED USERS ONLY!
 If the user goes through this package's standard functional interface,
 the file will always automatically be set to whatever mode it needs to be in
 without the user having to do anything.
 In particular, the call to write an attribute (\code{\link[ncdf]{att.put.ncdf}})
 handles this automatically.
 An example of the kind of situation where you would need this call is if you
 are adding a new variable to an already-existing netCDF file.
 This case is not really handled by this package.
}
\keyword{utilities}
