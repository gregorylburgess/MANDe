\name{dim.def.ncdf}
\alias{dim.def.ncdf}
\title{Define a netCDF Dimension}
\description{
 Defines a netCDF dimension.  This dimension initially only exists in memory.
 It is added to a netCDF variable using \code{var.def.ncdf()}, and written to disk
 using \code{create.ncdf()}.
}
\usage{
 dim.def.ncdf( name, units, vals, unlim=FALSE, create_dimvar=TRUE )
}
\arguments{
 \item{name}{Name of the dimension to be created (character string).}
 \item{units}{The dimension's units (character string; empty character string to suppress creation of the coordinate variable (dimvar)).}
 \item{vals}{The dimension's values (vector of numeric type). If integers are passed, the associated dimensional variable will be integer type; otherwise, it will be double precision.}
 \item{unlim}{If TRUE, this dimension is unlimited.}
 \item{create_dimvar}{If TRUE, a dimensional variable (aka coordinate variable) will be created for this dimension.}
}
\value{
 An object of class \code{dim.ncdf} that can later be passed to 
 \code{var.def.ncdf()}.
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 This routine creates a netCDF dimension in memory.  The dimension can then
 be passed to the routine \code{var.def.ncdf()} when creating a variable.

 Note that this interface to the netCDF library by default includes that more than the
 minimum required by the netCDF standard.  I.e., the netCDF standard allows
 dimensions with no units or values.  This call encourages creating dimensions
 that have units and values,
 as it is useful to ensure that all dimensions have units and
 values, and considerably easier to include them in this call than it is
 to add them later.  The units and values are implemented through "dimensional
 variables," which are variables with the same name as the dimension.  By default, these
 dimensional variables are created automatically -- there is no need for the
 user to create them explicitly.  Dimensional variables are standard practice
 in netCDF files.  To suppress the creation of the dimensional variable for
 the dimension, set passed parameter create\_dimvar to FALSE.  As a
 check, if create\_dimvar is FALSE, you must ALSO pass an empty string ('') as the unit, and the values
 must be simple integers from 1 to the length of the dimension (e.g., 1:10 to
 make a dimension of length 10).
 This empahsizes that without a dimensional variable, a netCDF file cannot
 store a dimension's units or values.

 The dimensional variable is usually created as a double precision floating
 point.  The other possibility is to pass integer values (using \code{as.integer},
 for example), in which case the dimensional variable with be integer.

 The return value of this function is an object of class \code{dim.ncdf}, which
 describes the newly created dimension.
 The \code{dim.ncdf} object is used for more than just creating a new 
 dimension, however.
 When opening an existing file, function \code{open.ncdf()} returns a 
 \code{ncdf} class object, which itself has a list of \code{dim.ncdf} objects
 that describe all the dimensions in that existing file.

 The \code{dim.ncdf} object has the following fields, which are all read only:
 1) name, which is a character string containing the name of the dimension;
 2) units, which is a character string containing the units for the dimension,
 if there are any (technically speaking, this is the "units" attribute of the
 associated coordinate variable); 3) vals, which is a vector containing the
 dimension's values (i.e., the values of the associated coordinate variable,
 or, if there is none, an integer sequence from 1 to the length of the dimension);
 3) len, which is the length of this dimension; 4) unlim, which is a boolean
 indicating whether or not this is an unlimited dimension.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf]{var.def.ncdf}}, \code{\link[ncdf]{create.ncdf}}
}
\examples{
# Define some straightforward dimensions
x <- dim.def.ncdf( "Lon", "degreesE", 0.5:359.5)
y <- dim.def.ncdf( "Lat", "degreesN", as.double(-89:89))
t <- dim.def.ncdf( "Time", "days since 1900-01-01", 1:10, unlim=TRUE)

# Make a variable with those dimensions.  Note order: time is LAST
salinity <- var.def.ncdf("Salinity",    "ppt",  list(x,y,t), 1.e30 )

# Create a netCDF file with this variable
ncnew <- create.ncdf( "salinity.nc", salinity )

close.ncdf(ncnew)

# Now, illustrate some manipulations of the dim.ncdf object.
filename <- "salinity.nc"
nc <- open.ncdf( filename )
print(paste("File",filename,"contains",nc$ndims,"dimensions"))
for( i in 1:nc$ndims ) {
	print(paste("Here is information about dimension number",i,":"))
	d <- nc$dim[[i]]
	print(paste("    Name  :",d$name))
	print(paste("    Units :",d$units))
	print(paste("    Length:",d$len))
	print("    Values:")
	print(d$vals)
	print(paste("    Unlimited:",d$unlim))
	}
}
\note{It is good practice, but not necessary, to pass the dimension's values 
 to this routine when the dimension is created.  It is also possible to 
 write them later with a call to 'put.var.ncdf', using as the dimension
 name as the 'varid' in the call.  This is useful when creating large variables
 with long unlimited dimensions; it can take a long time to write out the
 unlimited dimension's values.  In this case, it can be more efficient to
 step through the file, writing one timestep at a time, and write that 
 timestep's dimensional value at the same time.
}
\keyword{utilities}
