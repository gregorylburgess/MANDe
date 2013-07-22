\name{att.get.ncdf}
\alias{att.get.ncdf}
\title{Get attribute from netCDF file}
\description{
 Reads an attribute from a netCDF file.
}
\usage{
 att.get.ncdf( nc, varid, attname )
}
\arguments{
 \item{nc}{An object of class \code{ncdf} (as returned from \code{open.ncdf}),
 indicating what file to read from.}
 \item{varid}{The variable whose attribute is to be read.  Can be a
 character string with the variable's name, an object of class \code{var.ncdf},
 or an id contained in the "id" field of a \code{var.ncdf} object. As a
 special case, if varid==0, then it is assumed that we are reading a global
 attribute rather than a particular variable's attribute.}
 \item{attname}{Name of the attribute to read.}
}
\value{
 A list with two attributes, "hasatt" and "value".
 "hasatt" is TRUE if the named attribute was found, and FALSE otherwise.
 "value" is the (possibly vector) value of the attribute.
 If the on-disk type of the attribute is short or integer, then an integer
 value is returned.  If the on-disk type is float or double, than a 
 double value is returned.  If the on-disk type is character, than a
 character string is returned.
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 This function gets an attribute from a netCDF variable (or a global attribute
 from a netCDF file, if the passed argument "varid" is zero).
 Multiple attributes are returned in a vector.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf]{att.put.ncdf}}.
}
\examples{
# Make a simple netCDF file
filename <- "atttest_types.nc"
dim <- dim.def.ncdf( "X", "inches", 1:12 )
var <- var.def.ncdf( "Data", "unitless", dim, -1 ) 
ncnew <- create.ncdf( filename, var )

# Define some attributes of various types
attvaldbl <- 3.1415926536
att.put.ncdf( ncnew, var, "testatt_dbl", attvaldbl, prec="double" )
attvalsingle <- c(1.0,4.0,9.0,16.0)
att.put.ncdf( ncnew, var, "testatt_single", attvalsingle )
# varid=0 means it is a global attribute
att.put.ncdf( ncnew, 0, "globalatt_int", 32000, prec="int" ) 
att.put.ncdf( ncnew, 0, "globalatt_short", 7, prec="short" )
att.put.ncdf( ncnew, 0, "description", 
	"this is a test file with attributes of various types")
close.ncdf(ncnew)

# Now illustrate the use of the att.get.ncdf function by reading them back in
doitfor <- function( nc, var, attname ) {
	av <- att.get.ncdf( nc, var, attname )
	if( av$hasatt ) {
		print(paste("File",nc$filename,", var",var,"DOES have attribute",
				attname))
		print(paste("Storage mode:",storage.mode(av$value)))
		print("Attribute value:")
		print(av$value)
	} else {
		print(paste("File",nc$filename,", var",var,"does NOT have",
			"attribute", attname))
		}
}
	
nc <- open.ncdf( filename )
var <- "Data"
doitfor( nc, var, "testatt_dbl" )
doitfor( nc, var, "testatt_single" )
doitfor( nc, var, "testatt_wacko" )
doitfor( nc, 0,   "globalatt_int" )
doitfor( nc, 0,   "globalatt_short" )
doitfor( nc, 0,   "description" )
}
\keyword{utilities}
