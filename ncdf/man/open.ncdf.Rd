\name{open.ncdf}
\alias{open.ncdf}
\title{Open a netCDF File}
\description{
 Opens an existing netCDF file for reading (or, optionally, writing).
}
\usage{
 open.ncdf( con, write=FALSE, readunlim=TRUE, verbose=FALSE, ... )
}
\arguments{
 \item{con}{Name of the existing netCDF file to be opened.}
 \item{write}{If FALSE (default), then the file is opened read-only.  If TRUE, then writing to the file is allowed.}
 \item{readunlim}{When invoked, this function reads in the values of all dimensions from 
 the associated variables.  This can be slow for a large file with a long unlimited 
 dimension.  If set to FALSE, the values for the unlimited dimension are not 
 automatically read in (they can be read in later, manually, using \code{get.var.ncdf()}).}
 \item{...}{Arguments passed to or from other methods.}
 \item{verbose}{If TRUE, then messages are printed out during execution of this function.}
}
\value{
 An object of class \code{ncdf} that has the fields described below.
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 This routine opens an existing netCDF file for reading (or, if write=TRUE, for writing).
 To create a new netCDF file, use \code{create.ncdf()} instead.

 In addition to simply opening the file, information about the file and its contents
 is read in and stored in the returned object, which is of class \code{ncdf}.
 This class has the following user-accessible fields, all of which are read-only: 1) filename,
 which is a character string holding the name of the file; 2) ndims, which is an
 integer holding the number of dimensions in the file; 3) nvars, which is an integer
 holding the number of the variables in the file that are NOT coordinate variables
 (aka dimensional variables); 4) natts, which is an integer holding the number of
 global attributes; 5) unlimdimid, which is an integer holding the dimension id of
 the unlimited dimension, or -1 if there is none; 6) dim, which is a list of
 objects of class \code{dim.ncdf}; 7) var, which is a list of objects of class
 \code{var.ncdf}; 8) writable, which is TRUE or FALSE, depending on whether the file
 was opened with write=TRUE or write=FALSE.
 
 The concept behind the R interface to a netCDF file is that the \code{ncdf} object
 returned by this function, as well as the list of \code{dim.ncdf} objects contained
 in the ncdf object's "dim" list and the \code{var.ncdf} objects contained in the
 ncdf object's "var" list, completely describe the netCDF file. I.e., they hold
 the entire contents of the file's metadata.  Therefore, there are no R interfaces
 to the explicit netCDF query functions, such as "nc\_inq\_nvars" or "nc\_inq\_natts".
 The upshot is, look in the ncdf object or its children to get information about
 the netCDF file.  (Note: the \code{dim.ncdf} object is described in the help
 file for \code{\link[ncdf]{dim.def.ncdf}}; the \code{var.ncdf} object is described
 in the help file for \code{\link[ncdf]{var.def.ncdf}}).

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

# Add source info metadata to file
att.put.ncdf( ncnew, 0, "source", "Census 2000 from census bureau web site")

close.ncdf(ncnew)

# Now open the file and read its data
ncold <- open.ncdf("states_population.nc")
data <- get.var.ncdf(ncold)
print("here is the data in the file:")
print(data)
close.ncdf( ncold )
}
\keyword{utilities}
