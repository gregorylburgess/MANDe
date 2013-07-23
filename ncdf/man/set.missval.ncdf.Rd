\name{set.missval.ncdf}
\alias{set.missval.ncdf}
\title{Set the Missing Value Attribute For a netCDF Variable}
\description{
 Sets the missing\_value attribute for a netCDF variable.
}
\usage{
 set.missval.ncdf( nc, varid, missval ) 
}
\arguments{
 \item{nc}{An object of class \code{ncdf}, as returned by \code{open.ncdf} or \code{create.ncdf}.}
 \item{varid}{The name, \code{var.ncdf} object, or variable ID whose missing value will be set.}
 \item{missval}{The missing value to set.}
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 Missing values are special values in netCDF files whose value is to be taken
 as indicating the data is "missing".
 This is a convention, and is indicated by the netCDF variable having an
 attribute named "missing\_value" that holds this number.
 This function sets the "missing\_value" attribute for a variable.

 R uses a similar concept to indicate missing values, the "NA" value.
 When the ncdf library reads in data set from a pre-existing file, all data values
 that equal that variable's missing value attribute appear to the R code
 as being "NA" values.
 When the R code writes values to a netCDF variable, any "NA" values are
 set to that variable's missing value before being written out.
 This makes the mapping between netCDF's "missing\_value" attribute and
 R's "NA" values transparent to the user.

 For this to work, though, the user still has to specify a missing value
 for a variable.  Usually this is specified when the variable is created,
 as a required argument to \code{\link[ncdf]{var.def.ncdf}}.
 However, sometimes it is useful to add (or change) a missing value for variable
 that already exists in a disk file.  This function enables that.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf]{var.def.ncdf}}.
}
\examples{
# Make an example netCDF file with a given missing value.  We will
# then change the missing value in the file using set.missval.ncdf.

origMissVal <- -1.
dimX   <- dim.def.ncdf( "X", "meters", 1:7 )
varAlt <- var.def.ncdf( "Altitude", "km", dimX, origMissVal )
ncnew  <- create.ncdf( "transect.nc", varAlt )
data <- c(10.,2.,NA,1.,7.,NA,8.)
put.var.ncdf( ncnew, varAlt, data )
close.ncdf(ncnew)

# At this point, the actual data values in the netCDF
# file will be: 10 2 -1 1 7 -1 8
# because the "NA" values were filled with the missing
# value, -1.  Also, the missing_value attribute of variable
# "varAlt" will be equal to -1.

# Now change the missing value to something else.  Remember
# we have to open the file as writable to be able to change
# the missing value on disk!

newMissVal <- 999.9
nc <- open.ncdf( "transect.nc", write=TRUE )
varname <- "Altitude"
data <- get.var.ncdf( nc, varname )  # data now has: 10., 2., NA, 1., 7., NA, 8.
print(data)
set.missval.ncdf( nc, varname, newMissVal )
put.var.ncdf( nc, varname, data )
close.ncdf(nc)

# Now, the actual data values in the netCDF file will be:
# 10 2 999.9 1 7 999.9 8
# and the variables "missing_value" attributre will be 999.9

# **NOTE** that we had to explicitly read in the data and write
# it out again in order for the on-disk missing values in the
# data array to change!  The on-disk missing_value attribute for
# the variable is set automatically by this function, but it is 
# up to you whether or not you want to read in all the files 
# data and change the values to the new missing value.
}
\keyword{utilities}
