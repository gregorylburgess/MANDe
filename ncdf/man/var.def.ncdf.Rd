\name{var.def.ncdf}
\alias{var.def.ncdf}
\title{Define a netCDF Variable}
\description{
 Defines a netCDF variable.  This variable initially only exists in memory.
 It is written to disk using \code{create.ncdf()}.
}
\usage{
 var.def.ncdf( name, units, dim, missval, longname=name, prec="single")
}
\arguments{
 \item{name}{Name of the variable to be created (character string).}
 \item{units}{The variable's units (character string).}
 \item{dim}{The variable's dimension(s) (one or a list of "dim.netcdf" class objects).}
 \item{missval}{The variable's missing value.}
 \item{longname}{Optional longer name for the variable, which is assigned to the variables "long\_name" attribute. For example, a variable named "TS" might have the longname "Surface Temperature"}
 \item{prec}{Precision of the created variable.  Valid options: 'short' 'integer' 'single' 'double' 'char'.}
}
\value{
 An object of class \code{var.ncdf} that can later be passed to 
 \code{create.ncdf()}.
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 This routine creates a netCDF variable in memory.  The variable can then
 be passed to the routine \code{create.ncdf} when writing a file to disk.

 Note that this interface to the netCDF library includes that more than the
 minimum required by the netCDF standard.  I.e., the netCDF standard allows
 variables with no units or missing values.  This call requires units and 
 a missing value,
 as it is useful to ensure that all variables have units and
 missing values, and considerably easier to include them in this call than it is
 to add them later.  The units and missing value are implemented through attributes
 to the variable, named "units" and "missing\_value", respectively.
 This is standard practice in netCDF files.

 After a variable is defined with this call, and created on disk using
 \code{create.ncdf}, then data values for the variable can be written
 to disk using \code{put.var.ncdf}.

 This function returns a \code{var.ncdf} object, which describes the newly-created
 variable.  However, the \code{var.ncdf} object is used for more than just
 creating new variables.  The function \code{open.ncdf} returns a \code{ncdf}
 object that itself contains a list of \code{var.ncdf} objects that describe
 the variables in an existing, on-disk netCDF file.  (Note that coordinate
 variables are NOT included in this list.  Attributes of the coordinate variables
 are kept in the \code{dim.ncdf} object instead.)

 The \code{var.ncdf} object has the following fields, which are all read-only:
 1) name, which is a character string containing the name of the variable;
 2) units, which is a character string containing the contents of the
 variable's "units" attribute; 3) missval, which contains the contents of the
 variable's "missing\_value" attribute; 4) longname, which is the 
 contents of the variable's "long\_name" attribute, or defaults to the name
 of the variable if there is no "long\_name" attribute; 5) ndims, which is the
 number of dimensions this variable has; 6) dim, which is a list of objects of
 class "dim.ncdf" (see \code{\link[ncdf]{dim.def.ncdf}}), and describe this
 variable's dimensions; 7) unlim, which is TRUE if this variable has an unlimited
 dimension and FALSE otherwise; 8) varsize, which is a convenience array
 that gives the shape of the variable (in XYZT ordering).

 Note that the missval attribute does not need to be used much in R, because
 R's special value NA is fully supported.  I.e., when data is read in from an
 existing file, any values equal to the "missing" value are set to NA.  When
 data is written out, any NAs are set equal to the missing value.  If not
 explicitly set by the user, a default value of 1.e30 is used for the missing
 value.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf]{dim.def.ncdf}}, \code{\link[ncdf]{create.ncdf}}, 
 \code{\link[ncdf]{put.var.ncdf}}.
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

# Now illustrate some manipulations of the var.ncdf object
filename <- "states_population.nc"
nc <- open.ncdf(filename)
print(paste("File",nc$filename,"contains",nc$nvars,"variables"))
for( i in 1:nc$nvars ) {
	v <- nc$var[[i]]
	print(paste("Here is information on variable number",i))
	print(paste("   Name: ",v$name))
	print(paste("   Units:",v$units))
	print(paste("   Missing value:",v$missval))
	print(paste("   # dimensions :",v$ndims))
	print(paste("   Variable size:",v$varsize))
	}

# Illustrate creating variables of various types.  You will find
# that the type of the missing_value attribute automatically follows
# the type of the variable.
dimt <- dim.def.ncdf( "Time", "days", 1:3 ) 
missval <- -1
varShort <- var.def.ncdf( "varShort", "meters", dimt, missval, prec="short")
varInt   <- var.def.ncdf( "varInt",   "meters", dimt, missval, prec="integer")
varFloat <- var.def.ncdf( "varFloat", "meters", dimt, missval, prec="single")
varDouble<- var.def.ncdf( "varDouble","meters", dimt, missval, prec="double")
nctypes <- create.ncdf("vartypes.nc", list(varShort,varInt,varFloat,varDouble) )
close.ncdf(nctypes)
}
\keyword{utilities}
