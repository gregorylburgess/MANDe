# This provides an interface to netCDF functions.  As a general
# comment, awkwardness arises from two things:
#	1. R starts counting at 1, and netCDF counting
#	   starts at 0.
#	2. R array subscripts go in Fortran order (XYZT),
#	   while netCDF subscripts go in C order (TZYX).
# We take care of these problems EXCLUSIVELY in the R code,
# NOT in the C interface code!!!   This means that it is
# always the responsibility of the R program to take these
# differences into account.  From the point of view of an
# R program that calls any of these functions, they are 
# strictly R compliant (Fortran order, counting starts at 1).
#
# David W. Pierce
# Climate Research Division
# Scripps Institution of Oceanography
# dpierce@ucsd.edu
# 9-April-2001
#
#-----------------------------------------------------------------
#
# Here are the relevant objects:
#
# class: ncdf is a list with the following fields:
#	filename: name of the file, or "IN-MEMORY"
#	id	: netcdf file id
#	ndims	: integer # of dims in the file
#	nvars	: integer # of vars in the file that are **NOT** dimvars
#	natts	: integer # of global attributes
#	unlimdimid : integer dimension id of unlimited dimension, or -1 if none
#	dim	: a list of dim.ncdf objects
#	var	: a list of var.ncdf objects
#	varid2Rindex : for internal use only; maps a numeric varid stright
#		from the netcdf file to which element in the list of vars this var is.
#	writable: TRUE or FALSE
#
# class: dim.ncdf (returned by dim.def.ncdf, which creates a NEW 
#		netCDF dimension in memory, and part of the list of dims
#		in a ncdf object. NOTE that this is NOT what is returned
#		by dim.inq.ncdf, which is the low-level netCDF dim, not
#		the user-level R version of a netCDF dim.
#     *	name	: character dim name
#	units	: character units in udunits format.  
#	vals	: a vector of dimension values
#     *	len	: size of this dimension
#	id	: dimension ID in the netcdf file, if applicable
#	dimvarid: IFF this dim corresponds to a netCDF dim in an existing file,
#		  then this field will hold the dimvarid, or -1 if no dimvar.
#     *	unlim	: boolean, T or F to indicated unlimited or not
#	create_dimvar : usually TRUE; if FALSE, then no dimvar will be created,
#		  AND the units string must be empty, AND the values must
#		  be simple integers from 1 to the length of the dim.
# *=Indicates element is filled out by the low-level routine 'dim.inq.ncdf'.
#   Other elements are filled out in 'open.ncdf'.
#
# class: var.ncdf (returned by var.def.ncdf, which creates a NEW
#		netCDF variable in memory, and part of the list of vars
#		in a ncdf object.  NOTE that this is NOT what is
#		returned by var.inq.ncdf, which is the low-level netCDF
#		var, NOT the user-level R version of a netCDF var.
#	name	: character var name
#	units	: character units in udunits format
#	missval	: the 'missing_value' attribute, or defaults to default.missval.ncdf(). NOTE: 
#		  if the var has no missing value then this is 'NA'.  This does not mean
#		  that the missing value is NA -- the missing value cannot be NA.  It means
#		  the variable has no missing value.  For instance, character variables have
#		  no default missing value.
#	longname: the 'long_name' attribute, or defaults to name
#	id	: the varid of this variable, IFF it is in a file already
#	ndims	: number of dims this variable has
#	dim	: a list of type dim.ncdf, which is this variable's dims
#	unlim	: boolean, T if this var has an unlimited dim, F otherwise
#	varsize : a convenience array that gives the (X,Y,Z,T) size of
#		  the variable.
#	prec    : The precision of the ON-DISK representation of the
#		  variable.  Can be "byte", "short", "float", "double", 
#	    	  "integer", or "char".
#
# class: vals.ncdf (returned by get.var.ncdf)
#	xvals, yvals, zvals, tvals: the dimensional values, as appropriate
#	vals	: the data values
#
#=================================================================
version.ncdf <- function() {
	
	return("1.6")

}

#=================================================================
# Utility to return a streng of length 'n'; this is used for storage
#
blankstring.ncdf <- function( n ) {

	s10  <- '          '
	if( n <= 10 )
		return(s10)

	s100 <- paste(s10,s10,s10,s10,s10,s10,s10,s10,s10,s10)
	if( n <= 100 )
		return( s100 )

	s1000 <- paste(s100,s100,s100,s100,s100,s100,s100,s100,s100,s100)
	if( n <= 1000 )
		return( s1000 )

	n1000 <- as.integer(n/1000) + 1
	stor <- ''
	for( i in 1:n1000 )
		stor <- paste(stor,s1000)

	return(stor)
}

#=================================================================
# This is where the default missing value for vars is set.
#
default.missval.ncdf <- function()  
{
	return(1.e30)
}

#=================================================================
# Returns -1 if element 'el' is not in the list, and 
# if 'el' IS in the list, returns which number it is.
#
in.list.name.ncdf <- function( el, list ) {
	l <- length(list)
	if( l == 0 ) {
		return(-1)
		}
	for( i in 1:l ) {
		if( el == list[[i]]$name ) {
			return(i)
			}
		}
	return(-1)
}

#===============================================================
# This is the public interface for making a netCDF dimension
# class.  It makes it in memory, not in a file!  Returns a
# "dim.ncdf" object.  The dimvar will have the same precision
# as the passed values.  Therefore, it will always have double
# precision unless the vals are passed like this: as.integer(vals),
# in which case they will be integer.
#
# Example useage:
#	lon <- ncdf.def.dim("Lon", "degreesE", 0:359)
#	lat <- ncdf.def.dim("Lat", "degreesN", -90:90)
#	time <- ncdf.def.dim("time", "days since 1900-01-01", 0, unlim=T)
#
dim.def.ncdf <- function( name, units, vals, unlim=FALSE, create_dimvar=TRUE ) {

	if( ! is.character(name) ) {
		stop("Passed a dim name that is NOT a string of characters!")
		}

	len <- length(vals)
	if( ! create_dimvar ) {
		if( (units != '') || (storage.mode(vals) != "integer" ) || (vals[1] != 1) || (vals[len] != len))
			stop(paste("Error trying to create dimension named",name,": create_dimvar was specified",
				"to be FALSE, which indicates that NO dimensional variable is to be created;",
				"in this case, the unit string MUST be empty ('') and the dimension values MUST",
				"be simple integers from 1 to the length of the dimension (e.g., 1:len)"))
		}
	dim <- list()
	dim$name   	  <- name
	dim$units  	  <- units
	dim$vals   	  <- vals
	dim$len    	  <- len
	dim$id     	  <- -1
	dim$unlim  	  <- unlim
	dim$dimvarid 	  <- -1	# Only exists for on-disk dims
	dim$create_dimvar <- create_dimvar
	attr(dim,"class") <- "dim.ncdf"
	return(dim)
}

#===============================================================
# Produce a more useful listing of the netcdf file than
# just a dump of the object.
#
print.ncdf <- function( x, ... ) {

	nc <- x
	print(paste("file",nc$filename,"has",nc$ndims,"dimensions:"))
	if( nc$ndims > 0 )
		for( i in 1:nc$ndims ) {
			print(paste(nc$dim[[i]]$name, "  Size:", nc$dim[[i]]$len))
			}
	print("------------------------")

	print(paste("file",nc$filename,"has",nc$nvars,"variables:"))
	if( nc$nvars < 1 )
		return;
	for( i in 1:nc$nvars ) {
		nd <- nc$var[[i]]$ndims
		dimstring <- '['
		if( nd > 0 ) {
			for( j in 1:nd ) {
				dimstring <- paste(dimstring,nc$var[[i]]$dim[[j]]$name,sep='')
				if( j < nd )
					dimstring <- paste(dimstring,',',sep='')
				}
			}
		dimstring <- paste(dimstring,'] ',sep='')
		print(paste(nc$var[[i]]$prec, ' ', nc$var[[i]]$name,dimstring,
			' Longname:',nc$var[[i]]$longname,
			' Missval:',nc$var[[i]]$missval,sep=''))
		}
	}

#===============================================================
# This is the public interface for making a netCDF variable
# class.  It makes it in memory, not in a file!  Returns
# a "var.ncdf" object.
#
# Example useage, where "lon" and "lat" are objects of class "dim.ncdf":
#
#	ncvar <- ncdf.def.var( "temp", "degC", list(lon,lat), 1.e30, "Temperature" )
# or
#	ncvar <- ncdf.def.var( "temp", "degC", lon, 1.e30, "Temperature" )
#
# NOTE that the passed dimentions (in this example, "lon" and "lat"
# should have been made by "ncdf.def.dim", and so are of class
# "dim.ncdf".   Argument "prec" will indicate what precision the variable
# is made on the disk.  Allowed values are "short", "integer", "single", 
# "double", "byte", and "char".
# 'missval' is the value to be assigned to the 'missing_value' attribute
# for the variable.  It should be a number representable in the precision-type
# of the variable.  So, for example, for single precision, a floating point
# number with magnitude less than 1.e36 should be used.  As a special value,
# if missval is set to NA, then *no* missing value will be created for the 
# variable.
#
var.def.ncdf <- function( name, units, dim, missval, longname=name, prec="single" ) {

	if( ! is.character(name) ) {
		stop("Passed a var name that is NOT a string of characters!")
		}

	if( storage.mode(missval) == "character" )
		prec <- 'char'

	var <- list()
	var$name     <- name
	var$units    <- units
	var$missval  <- missval
	var$longname <- longname
	var$id       <- -1
	attr(var,"class") <- "var.ncdf"

	if( prec == "float" )
		prec <- 'single'

	if( (prec != "short")   && (prec != "single") && (prec != "double") && 
	    (prec != "integer") && (prec != "char")   && (prec != "byte"))
		stop(paste("var.def.ncdf: error: unknown precision specified:",prec,". Known values: short single double integer char byte"))
	var$prec <- prec

	#-----------------------------------------------------
	# Have to figure out if 'dim' is a dim.ncdf object or 
	# a LIST of dim.ncdf objects.
	#-----------------------------------------------------
	if( is.character(class(dim)) && (class(dim) == "dim.ncdf") )
		dim <- list(dim)
	var$dim  <- dim
	var$ndims <- length(var$dim)
	for( i in 1:var$ndims ) {
		if( class(var$dim[[i]]) != "dim.ncdf" ) {
			print(var)
			stop("Error, passed variable has a dim that is NOT of class dim.ncdf!")
			}
		}

	#-------------------------------------------------
	# A variable is unlimited if any of its dimensions
	# is unlimited
	#-------------------------------------------------
	varunlimited <- FALSE
	if( var$ndims != 0 ) {
		for( i in 1:var$ndims ) {
			if( var$dim[[i]]$unlim ) 
				varunlimited <- TRUE
			}
		}
	var$unlim <- varunlimited

	return(var)
}

#===============================================================
# This is the public interface for opening an already-existing
# netCDF file.  If you want to create a NEW netCDF file, use
# create.ncdf instead.  This returns an object of class 'ncdf',
# used to access netCDF functions.  If you want to modify an 
# already-existing netCDF file, use this function and set write=T.
# By default, this will read in the values of the unlimited dimension
# along with the values of all the other dimensions.  Since this
# can be time comsuming, there is the option of setting 'readunlim' 
# to 'FALSE', which will avoid this behavor.  Then you must read
# in the values yourself with a 'get.var.ncdf' call, instead of
# accessing the ncid$dim[[DIMNAME]]$vals array.
#
open.ncdf <- function( con, write=FALSE, readunlim=TRUE, verbose=FALSE, ... ) {

	if( verbose )
		print(paste("open.ncdf: entering, version=",version.ncdf()))

	rv <- list()

	if( write )
		rv$cmode <- 1
	else
		rv$cmode <- 0

	rv$id    <- -1
	rv$error <- -1
	rv <- .C("R_nc_open",
		as.character(con),
		as.integer(rv$cmode),
		id=as.integer(rv$id),
		error=as.integer(rv$error),
		PACKAGE="ncdf")
	if( rv$error != 0 ) 
		stop(paste("Error in open.ncdf trying to open file",con))
	if( verbose )
		print(paste("open.ncdf: back from call to R_nc_open, ncid=",rv$id))

	#-------------------------------------------------
	# Now we make our elaborate ncdf class object
	#-------------------------------------------------
	nc       <- list()
	nc$id    <- rv$id
	attr(nc,"class") <- "ncdf"

	#---------------------------------------
	# Get general information about the file
	#---------------------------------------
	if( verbose )
		print("open.ncdf: getting general info about the file")
	rv            <- list()
	rv$ndims      <- -1
	rv$nvars      <- -1
	rv$natts      <- -1
	rv$unlimdimid <- -1
	rv$error      <- -1
	rv <- .C("R_nc_inq",
		as.integer(nc$id),
		ndims=as.integer(rv$ndims),
		nvars=as.integer(rv$nvars),
		natts=as.integer(rv$natts),
		unlimdimid=as.integer(rv$unlimdimid),
		error=as.integer(rv$error),
		PACKAGE="ncdf")
	if( rv$error != 0 ) 
		stop(paste("R_nc_inq returned error on file",con,"!"))
	if( verbose )
		print(paste("open.ncdf: back from call to R_nc_inq for id=",nc$id,"   ndims=",rv$ndims,
			"   nvars=",rv$nvars,"  natts=",rv$natts,"  umlimdimid=",rv$unlimdimid))
	nc$ndims        <- rv$ndims
	nc$natts        <- rv$natts
	nc$unlimdimid   <- rv$unlimdimid + 1	# Change from C to R convention
	nc$filename     <- con
	nc$varid2Rindex <- rep(0,rv$nvars)	
	nc$writable     <- write

	#-----------------------------------------------
	# Get all the dimensions that this file has.
	# Get their values as well (for caching).  
	#-----------------------------------------------
	nc$dim   <- list()
	dimnames <- character()
	for( i in 1:nc$ndims ) {	
		if( verbose )
			print(paste("open.ncdf: getting dim info for dim number",i))
		#-----------------------------------------------------------
		# As a general note, the function dim.inq.ncdf does NOT
		# return a full-fledged "dim.ncdf" object. It returns 
		# only a subset of fields that directly correspond to 
		# a low-level, netCDF dimension in a file.  We now fill
		# out the rest of the fields to make it into a real dim.ncdf 
		# object.
		#-----------------------------------------------------------
		d          <- dim.inq.ncdf(nc,i)
		d$id	   <- i
		d$dimvarid <- varid.inq.ncdf(nc,d$name)
		if( verbose )
			print(paste(".....dim name is",d$name,"  len=",d$len,"     dimvarid=",d$dimvarid))
		if( d$dimvarid == -1 ) {	# No dimvar for this dim
			d$vals  <- 1:d$len
			d$units <- ""
			d$create_dimvar <- FALSE	# in case this dim is passed to create.ncdf()
			}
		else {	
			# This dim has a dimvar -- get its properties
			if( verbose )
				print(paste("open.ncdf: getting dimvar info for dim ",d$name))
			attv <- att.get.ncdf( nc, d$dimvarid, "units" )
			if( attv$hasatt )
				d$units <- attv$value
			else
				d$units <- ""
			if( d$unlim && (! readunlim)) # Is unlimited, don't read vals, too slow
				d$vals <- rep(NA,d$len)
			else			# Otherwise, read vals
				d$vals <- get.var.ncdf( nc, forcevarid=d$dimvarid, verbose=verbose )
			d$create_dimvar <- TRUE		# in case this dim is passed to create.ncdf()
			if( verbose )
				{
				print("------------------------------")
				print("Here is new dim:")
				print(paste("name=",d$name,"  len=",d$len,"   unlim=",d$unlim,"   id=",d$id,"   dimvarid=",d$dimvarid,"   units=",d$units))
				print("------------------------------")
				}
			}
		attr(d,"class") <- "dim.ncdf"	# Is a complete dim.ncdf object now
		nc$dim[[i]] <- d
		dimnames[i] <- d$name
		if( verbose )
			print(paste(".......open.ncdf: done processing dim ",d$name))
		}
	attr(nc$dim,"names") <- dimnames
	if(verbose) {
		print("open.ncdf: setting dim$<names> to:")
		print(dimnames)
		}
	
	#-------------------------------------------
	# Get all the vars that this file has.  Note
	# that dimvars are NOT included in the count
	# of vars!!
	#-------------------------------------------
	if( verbose )
		print(paste("open.ncdf: getting var info.  # vars (INCLUDING dimvars)=",rv$nvars))
	nc$nvars <- 0
	nc$var   <- list()
	varnames <- character()
	for( i in 1:rv$nvars ) {
		name <- varname.inq.ncdf( nc, i )
		if( dimid.inq.ncdf( nc, name ) == -1 ) {	# Only process if NOT a dimvar
			if( verbose )
				print(paste("open.ncdf: will process varid=",i,"  name=",name))
			#--------------------------------------
			# No dim with same name as this var, so
			# this var must NOT be a dimvar.
			#--------------------------------------
			v <- var.inq.ncdf( nc, i )
			attr(v,"class") <- "var.ncdf"
			nc$nvars <- nc$nvars + 1
			#--------------------
			# Get this var's dims
			#--------------------
			v$dims   <- list()
			varunlim <- FALSE
			if( v$ndims > 0 ) {
				for( j in 1:v$ndims ) {
					v$dim[[j]] <- nc$dim[[ v$dimids[j] ]]
					if( v$dim[[j]]$unlim )
						varunlim <- TRUE
					v$varsize <- append(v$varsize, v$dim[[j]]$len)
					}
				}
			v$unlim <- varunlim

			#----------------------------------------
			# Get this var's missing value, or set to
			# a default value if it does not have one
			#----------------------------------------
			mv <- att.get.ncdf( nc, i, "missing_value" )
			if( mv$hasatt )
				v$missval <- mv$value
			else if( (v$prec=="float") || (v$prec=="double"))
				v$missval <- default.missval.ncdf()
			else
				v$missval <- NA

			#-------------------------------------------
			# Get add_offset and scale_factor attributes 
			#-------------------------------------------
			ao <- att.get.ncdf( nc, i, "add_offset" )
			if( ao$hasatt ) {
				v$hasAddOffset <- TRUE
				v$addOffset    <- ao$value
				}
			else
				v$hasAddOffset <- FALSE
			sf <- att.get.ncdf( nc, i, "scale_factor" )
			if( sf$hasatt ) {
				v$hasScaleFact <- TRUE
				v$scaleFact    <- sf$value
				}
			else
				v$hasScaleFact <- FALSE
			
			nc$var[[nc$nvars]] <- v
			nc$varid2Rindex[i] <- nc$nvars 
			varnames <- append(varnames,v$name)
			if( verbose ) {
				print("-----------------------")
				print("Here is new var:")
				print(paste("name=",v$name,"  id=",v$id,"   ndims=",v$ndims,"   prec=",v$prec))
				print("size=")
				print(v$size)
				print("dimids=")
				print(v$dimids)
				}
			}
		}
	attr(nc$var,"names") <- varnames

	if( verbose )
		print(paste("open.ncdf: leaving for ncid=",nc$id))

	return(nc)
}

#===============================================================
# This is used to set a missing value to a desired value.  It
# is used when a file has not been created on disk yet, or when
# an existing disk file has been opened to be writable.
# 'varid' may be a character string with the var's name,
# a var.ncdf class, or an integer.  If it is an integer, it is 
# R-style varid, i.e., the varid as reported by a C function call
# but incremented by one (to adjust for R's counting, which starts
# at 1).  Note that this is the number stored in a var.ncdf class
# "id" element.  So, for example, all these will operate on the
# same variable:
#
# If v is of class "var.ncdf", and:
#	v$name == "myvarname", and
#	v$id   == 7,
#
# then all these calls set the missing value for this var to -1:
#
# set.missval.ncdf( ncid, v, -1 )
# set.missval.ncdf( ncid, 7, -1 )
# set.missval.ncdf( ncid, "myvarname", -1 )
#
set.missval.ncdf <- function( nc, varid, missval ) {

	if( class(nc) != "ncdf" ) 
		stop("set.missval.ncdf: passed nc NOT of class ncdf.file!")

	#------------------------------------------------------
	# Can't do this if the file is on disk and not writable
	#------------------------------------------------------
	if( (nc$filename != "IN-MEMORY") && (! nc$writable))
		stop("set.missval.ncdf: the netcdf file was NOT opened in write mode!")

	varid <- vobjtovarid( nc, varid )
	idx   <- nc$varid2Rindex[varid]
	nc$var[[idx]]$missval <- missval

	if( nc$filename != "IN-MEMORY" )
		att.put.ncdf( nc, varid, "missing_value", missval )
}

#===============================================================
# This is the public interface for creating a netCDF file
# on disk.  It takes arguments of class ncdf.var to put in the file.
# It creates the file on disk and returns an object of class "ncdf"
# that can be used to access that file.
#
# Example usage, where "temp" and "salin" are objects of class ncdf.var:
#
#	nc <- ncdf.create( "test.nc", list(temp,salin))
# or
#	nc <- ncdf.create( "test.nc", salin )
#
create.ncdf <- function( filename, vars, verbose=FALSE ) {

	if( ! is.character(filename))
		stop("input filename must be a character string")
	if( nchar(filename) < 1 )
		stop("input filename must be at least 1 character long")

	#---------------------------------------------------
	# Have to tell if the input vars.orig is a single
	# var or a list of vars.   Do it by examining 
	# vars.orig$class.  If vars.orig is a single var, 
	# then this will equal "var.ncdf"; if vars.orig is a 
	# list of vars, this will be NULL.
	#---------------------------------------------------
	if( is.character(class(vars)) && (class(vars) == "var.ncdf") ) {
		vars <- list(vars)
		if( verbose )
			print("create.ncdf: input was a single var")
		}
	else
		{
		if( verbose )
			print("create.ncdf: input was a list of vars")
		}

	nc <- list()

	#----------------
	# Create the file
	#----------------
	nc$cmode    <- 0
	nc$error    <- -1
	nc$id       <- -1
	if( verbose )
		print(paste("Calling R_nc_create for file ",filename))
	nc<-.C("R_nc_create",
		filename,
		as.integer(nc$cmode),
		id=as.integer(nc$id),
		error=as.integer(nc$error),
		PACKAGE="ncdf")
	if( nc$error != 0 )
		stop("Error in create.ncdf!")
	if( verbose )
		print(paste("back from R_nc_create for file ",filename))
	nc$nvars  <- 0
	attr(nc,"class")  <- "ncdf"
	nc$filename <- filename
	nc$writable <- TRUE

	if( verbose )
		print("create.ncdf: about to create the dims")
	nc$ndims  <- 0
	nc$dim    <- list()
	nc$var    <- list()

	max_nc_dims <- 20
	nc$varid2Rindex <- rep(0,(length(vars)+10)*max_nc_dims)	# size is max no. of vars we can have (including dimvars)

	#---------------------------------------------------------------
	# Add the vars to the file.  NOTE that this also adds the unique
	# dims (and hence, dimvars) to the file as a side effect!
	# Note also that the returned 'nc' value is updated each time
	# this subroutine is called.
	#---------------------------------------------------------------
	for(ivar in 1:length(vars)) 
		nc <- var.add.ncdf( nc, vars[[ivar]], verbose=verbose, indefine=TRUE )

	#-----------------------------------------------------------
	# Set the names attribute on the $var and $dim lists so that
	# we can access them by name instead of only by position
	#-----------------------------------------------------------
	varnames <- array('',nc$nvars)
	for( ivar in 1:nc$nvars )
		varnames[ivar] <- nc$var[[ivar]]$name
	attr(nc$var,"names") <- varnames
	dimnames <- array('',nc$ndims)
	for( idim in 1:nc$ndims )
		dimnames[idim] <- nc$dim[[idim]]$name
	attr(nc$dim,"names") <- dimnames

	#-----------------
	# Exit define mode
	#-----------------
	enddef.ncdf( nc )

	return(nc)
}

#===============================================================
dim.same.ncdf <- function( d1, d2 ) {

	if( class(d1) != "dim.ncdf" ) 
		stop("error, class of first passed argument is not dim.ncdf!")
	if( class(d2) != "dim.ncdf" ) 
		stop("error, class of first passed argument is not dim.ncdf!")

	if( d1$name != d2$name )
		return(FALSE)

	if( d1$len != d2$len )
		return(FALSE)

	if( d1$unlim != d2$unlim )
		return(FALSE)

	return(TRUE)
}

#===============================================================
# This is a SPECIAL PURPOSE function ONLY to be used when adding
# an already defined variable (accomplished via "var.def.ncdf")
# to an ALREADY EXISTING netcdf file.  Normally, when making
# a new netcdf file from scratch, a list of vars to be created
# would be passed to "create.ncdf"; this is the preferred method
# of putting vars in a file.  However, sometimes it's necessary
# to add a var to an already-existing file; that's what this
# routine is for.
#
var.add.ncdf <- function( nc, v, verbose=FALSE, indefine=FALSE ) {
	
	if( verbose )
		print(paste("var.add.ncdf: entering with indefine=",indefine))

	if( class(nc) != "ncdf" ) 
		stop("var.add.ncdf: passed nc NOT of class ncdf!  The first arg to var.add.ncdf must be the return value from a call to open.ncdf(...,write=TRUE)")
	if( verbose )
		print(paste("var.add.ncdf: ncid of file to add to=",nc$id,
			"   filename=",nc$filename,"    writable=",nc$writable))

	if( class(v) != "var.ncdf" ) 
		stop("var.add.ncdf: passed var NOT of class var.ncdf! The second arg to var.add.ncdf must be the return value from a call to var.def.ncdf")
	if( verbose )
		print(paste("var.add.ncdf: varname to add=",v$name))

	if( ! indefine ) {
		if( verbose )
			print(paste("var.add.ncdf: about to redef ncid=",nc$id))
		redef.ncdf( nc )	# Go back into define mode
		}

	#-----------------------------------------------------
	# Create the dims for this var.  Harder than it sounds 
	# because we must take care not to repeat making a dim 
	# that occurs in more than one variable.
	#---------------------------------------------------
	nd <- v$ndims
	dimvarids <- array(0,nd)
	if( verbose )
		print(paste("var.add.ncdf: creating",nd,"dims for var",v$name))
	for( idim in 1:nd ) {
		d <- v$dim[[idim]]
		if( verbose )
			print(paste("var.add.ncdf: working on dim >",d$name,"< (number",idim,") for var",v$name))

		#-----------------------------------------------
		# See if we've already made a dim with this name
		#-----------------------------------------------
		place <- -1
		if( length(nc$dim) > 0 ) {
			for( ii in 1:length(nc$dim)) {
				if( nc$dim[[ii]]$name == d$name ) {
					#---------------------------------------------------------------
					# Check to make sure this is REALLY the same dim, even though we
					# know it has the same name as an existing dim!
					#---------------------------------------------------------------
					if( ! dim.same.ncdf( nc$dim[[ii]], d )) {
						paste("Error, when trying to add variable named",
							v$name, "to file",nc$filename,"I found this variable has a dim named",d$name)
						paste("However, the file ALREADY has a dim named",nc$dim[[ii]]$name,"with different characteristics than the new dim with the same name!")
						stop("This is not allowed.")
						}
					place <- ii
					break
					}
				}
			}
		if( place == -1 ) {
			#--------------------------------------------
			# This dim has not been seen before -- create
			#--------------------------------------------
			if( verbose )
				print(paste("create.ncdf: creating dim",d$name))
			ids         <- dim.create.ncdf(nc,d,verbose)	# *** NOTE: makes the dimvar, too! ***
			dimid       <- ids[1]
			dimvarid    <- ids[2]

			newel       <- list()
			attr(newel,"class") <- "dim.ncdf"
			newel$name     <- d$name
			newel$units    <- d$units
			newel$vals     <- d$vals
			newel$lene     <- d$len
			newel$id       <- dimid
			newel$unlim    <- d$unlim
			newel$dimvarid <- dimvarid
			nc$ndims    <- nc$ndims + 1
			nc$dim[[nc$ndims]] <- newel
			}
		else
			dimid <- nc$dim[[place]]$id

		dimvarids[idim] <- dimid
		}

	#----------------------------------------------------
	# Reverse the dimvarids, because R uses Fortran-style
	# ordering and we are using the C netCDF interface
	#----------------------------------------------------
	dimids <- dimvarids
	dimids <- dimids[length(dimids):1]
	newvar       <- list()
	newvar$id    <- -1
	newvar$error <- -1

	#-----------------------------------------------------------
	# Select the routine we will be using to create the variable
	# based on its precision
	#-----------------------------------------------------------
	if( verbose )
		print(paste("create.ncdf: creating",v$prec,"precision var",v$name))
	if( (v$prec == "integer") || (v$prec == "int") )
		funcname <- "R_nc_def_var_int"
	else if( v$prec == "short" )
		funcname <- "R_nc_def_var_short"
	else if( (v$prec == "single" ) || (v$prec == "float"))
		funcname <- "R_nc_def_var_float"
	else if( v$prec == "double" )
		funcname <- "R_nc_def_var_double"
	else if( v$prec == "char" )
		funcname <- "R_nc_def_var_char"
	else if( v$prec == "byte" )
		funcname <- "R_nc_def_var_byte"
	else
		stop(paste("internal error in create.ncdf: var has unknown precision:",v$prec,". Known vals: short single double integer char byte"))

	#---------------------------------
	# Now actually create the variable
	#---------------------------------
	newvar<-.C(funcname,
		as.integer(nc$id),
		v$name,
		as.integer(v$ndims),
		as.integer(dimids-1),	# Change from R to C convention
		id=as.integer(newvar$id),
		error=as.integer(newvar$error),
		PACKAGE="ncdf")
	if( verbose )
		print(paste("create.ncdf: C call returned value",newvar$error))
	if( newvar$error != 0 ) 
		stop("Error in create.ncdf, defining var!")
	newvar$id <- newvar$id + 1	# Change from C to R convention
	v$id      <- newvar$id
	nc$nvars  <- nc$nvars + 1
	nc$var[[nc$nvars]] <- v
	nc$varid2Rindex[newvar$id] <- nc$nvars

	#------------------------------------------------------
	# Add the attributes -- units, missing_value, long_name
	#------------------------------------------------------
	if( (! is.null( v$units )) && (! is.na(v$units))) 
		att.put.ncdf( nc, newvar$id, "units", v$units, definemode=TRUE )

	if( is.null( v$missval ) && ((v$prec=="single") || (v$prec=="float") || (v$prec=="double"))) {
		att.put.ncdf( nc, newvar$id, "missing_value", default.missval.ncdf(), definemode=TRUE )
		}
	else
		{
		if( ! is.na(v$missval) ) {
			att.put.ncdf( nc, newvar$id, "missing_value", v$missval, definemode=TRUE )
			}
		}
	if( v$longname != v$name )
		att.put.ncdf( nc, newvar$id, "long_name", v$longname, definemode=TRUE )

	if( ! indefine )
		enddef.ncdf( nc )	# Exit define mode

	return(nc)
}

#===============================================================
# This is the private interface that actually does the 
# netCDF calls.  User code should never go through this.
# To make a ncdf.dim object, use ncdf.def.dim() instead.
# This makes BOTH the dim AND the dimvar (and RETURNS
# dimid AND dimvarid).
#
dim.create.ncdf <- function( nc, d, verbose=FALSE ) {

	if( class(nc) != "ncdf" ) 
		stop("dim.create.ncdf: passed nc NOT of class ncdf.file!")
	if( verbose )
		print(paste("dim.create.ncdf: entering for ncid=",nc$id))

	if( class(d) != "dim.ncdf" ) 
		stop("dim.create.ncdf: passed d NOT of class ncdf.dim!")
	if( verbose )
		print(paste("dim.create.ncdf: entering for dim",d$name))

	#-------------------
	# Make the dimension
	#-------------------
	ncdim <- list()
	ncdim$error <- -1
	ncdim$id    <- -1
	sizetouse   <- d$len
	if( d$unlim )
		sizetouse <- 0
	if( verbose )
		print(paste("dim.create.ncdf: about to call R_nc_def_dim for dim",d$name))
	ncdim<-.C("R_nc_def_dim",
		as.integer(nc$id),
		d$name,
		as.integer(sizetouse),
		id=as.integer(ncdim$id),
		error=as.integer(ncdim$error),
		PACKAGE="ncdf")
	if( ncdim$error != 0 ) 
		stop("Error in dim.create.ncdf!")
	ncdim$id <- ncdim$id + 1	# Change from C to R convention

	#-----------------------------
	# Make the dimvar if requested
	#-----------------------------
	dimvar<-list()
	if( d$create_dimvar ) {
		if( verbose ) print(paste("dim.create.ncdf: making dimvar for dim",d$name))
		dimvar$id    <- -1
		dimvar$error <- -1
		if( storage.mode(d$vals) == "integer" ) {
			if( verbose )
				print(paste("dim.create.ncdf: about to call R_nc_def_var_int for dimvar",d$name))
			dimvar<-.C("R_nc_def_var_int",
				as.integer(nc$id),
				d$name,
				as.integer(c(1)),
				as.integer(ncdim$id-1),	# Change from R to C convention
				id=as.integer(dimvar$id),
				error=as.integer(dimvar$error),
				PACKAGE="ncdf")
			}
		else
			{
			if( verbose )
				print(paste("dim.create.ncdf: about to call R_nc_def_var_double for dimvar",d$name))
			dimvar<-.C("R_nc_def_var_double",
				as.integer(nc$id),
				d$name,
				as.integer(c(1)),
				as.integer(ncdim$id-1), # Change from R to C convention
				id=as.integer(dimvar$id),
				error=as.integer(dimvar$error),
				PACKAGE="ncdf")
			}
		if( dimvar$error != 0 ) 
			stop("Error defining dimvar in routine dim.create.ncdf")
		dimvar$id <- dimvar$id + 1	# Change from C to R convention
		dimvarid  <- dimvar$id

		#---------------------------------
		# Put in the dimvals as specified.
		#---------------------------------
		enddef.ncdf( nc )	# Must exit define mode for this
		rv <- list()
		rv$error <- -1
		start <- 0		# Use C convention
		count <- length(d$vals)
		if( count > 0 ) {
			if( storage.mode(d$vals) == "integer" ) {
				if( verbose )
					print(paste("dim.create.ncdf: about to call R_nc_put_vara_int dimvals for dimvar",d$name))
				rv <- .C("R_nc_put_vara_int",
					as.integer(nc$id),
					as.integer(dimvar$id-1),	# Change from R to C convention
					as.integer(start),
					as.integer(count),
					as.integer(d$vals),
					error=as.integer(rv$error),
					PACKAGE="ncdf")
				}
			else if( storage.mode(d$vals) == "double" ) {
				if( verbose )
					print(paste("dim.create.ncdf: about to call R_nc_put_vara_double dimvals for dimvar",d$name))
				rv <- .C("R_nc_put_vara_double",
					as.integer(nc$id),
					as.integer(dimvar$id-1),	# Change from R to C convention
					as.integer(start),
					as.integer(count),
					as.double(d$vals),
					error=as.integer(rv$error),
					PACKAGE="ncdf")
				}
			else
				stop(paste("dim.create.ncdf: unknown storage mode:",storage.mode(d$vals),"for dim",d$name))
			if( rv$error != 0 )
				stop("Error in dim.create.ncdf, while writing dimvar values!")
			}
		redef.ncdf( nc )	# Go back into define mode

		#----------------------------------------------------
		# Set the dimension's (dimvar's, actually) attributes
		#----------------------------------------------------
		att.put.ncdf( nc, dimvar$id, "units", d$units, definemode=TRUE )

		}	# end of "if(create_dimvar)"
	else
		{
		if( verbose ) print(paste("dim.create.ncdf: NOT making dimvar for dim",d$name))
		#----------------------------------------------------------
		# if we were NOT asked to create the dimvar (via an empty
		# units string) than make sure NO dim values were specified
		# except simple integers from 1 to len!
		#----------------------------------------------------------
		if( (storage.mode( d$vals ) != "integer" ) || (d$vals[1] != 1) || (d$vals[d$len] != d$len))
			stop(paste("Error trying to create dimension named",d$name,": the passed units string",
				"was empty, which indicates that NO dimensional variable is to be created.",
				"In this case, the dimension values MUST be simple integers from 1 to the length",
				"of the dimension (e.g., 1:len)"))
		dimvar$id = -1
		}

	if( verbose )
		print(paste("dim.create.ncd: exiting for ncid=",nc$id,"  dim=",d$name))

	#----------------------------------------------
	# Return the dimvar ID of the newly created dim
	#----------------------------------------------
	return(c(ncdim$id,dimvar$id))
}

#===============================================================
# This returns a list; first element of list (named "hasatt")
# is TRUE if the variable had an attribute with name "attname", and
# is FALSE otherwise.  Second element of the list (named "value")
# holds the value IFF the variable had an attribute of that name.
# "varid" can be an integer (R-type varid with counting starting
# at 1), object of class var.ncdf, or character string with a var's name.
# If the attribute type is short or integer, an integer value is
# returned.  If the attribute type is float or double, a double
# value is returned.  If the attribute type is text, a character
# value is returned.
#
att.get.ncdf <- function( nc, varid, attname ) {

	verbose <- FALSE

	if( verbose )
		print(paste("att.get.ncdf: entering with varid=",varid," attname=",attname))
	
	if( class(nc) != "ncdf" ) 
		stop("Error: the first argument to att.get.ncdf is not of class ncdf!")

	varid <- vobjtovarid( nc, varid, allowdimvar=TRUE )

	retval <- list()

	#----------------------------------------------------
	# Find out if the attribute exists for this variable, 
	# and, if so, what type and length it is.
	#----------------------------------------------------
	rv0 <- list()
	rv0$error  <- -1
	rv0$attlen <- -1
	rv0$type   <- -1
	rv0 <- .C("R_nc_inq_att",
		as.integer(nc$id),
		as.integer(varid-1),	# Convert from R to C convention
		as.character(attname),
		type=as.integer(rv0$type), # 1=short 2=int 3=single 4=double 5=text  6=byte
		attlen=as.integer(rv0$attlen),
		error=as.integer(rv0$error),
		PACKAGE="ncdf")
	if( rv0$error != 0 ) {
		#---------------------------------------------------------
		# This variable did NOT have an attribute named 'attname',
		# or it is of a type not handled.
		#---------------------------------------------------------
		retval$hasatt <- FALSE
		retval$value  <- 0
		return(retval)
		}

	retval$hasatt <- TRUE

	rv <- list()
	rv$error     <- -1

	if( (rv0$type == 1) || (rv0$type == 2) || (rv0$type == 6)) {
		#--------------------
		# Short, Int, or Byte
		#--------------------
		rv$attribute <- rep(as.integer(0),rv0$attlen)
		rv <- .C("R_nc_get_att_int",
			as.integer(nc$id),
			as.integer(varid-1),	# go from R to C convention
			as.character(attname),
			attribute=as.integer(rv$attribute),
			error=as.integer(rv$error),
			PACKAGE="ncdf")
		}
	else if( (rv0$type == 3) || (rv0$type == 4)) {
		#-----------------
		# Single or Double
		#-----------------
		rv$attribute <- rep(0.0,rv0$attlen)
		rv <- .C("R_nc_get_att_double",
			as.integer(nc$id),
			as.integer(varid-1),	# go from R to C convention
			as.character(attname),
			attribute=as.double(rv$attribute),
			error=as.integer(rv$error),
			PACKAGE="ncdf")
		}
	else if( rv0$type == 5 ) {
		#------------------------------------------------------
		# This is a string NC_MAX_LEN long, to provide storage.
		# I'm not sure if this is needed or not....
		#------------------------------------------------------
		rv$attribute <- blankstring.ncdf( rv0$attlen )
		rv <- .C("R_nc_get_att_text",
			as.integer(nc$id),
			as.integer(varid-1),	# go from R to C convention
			as.character(attname),
			attribute=as.character(rv$attribute),
			error=as.integer(rv$error),
			PACKAGE="ncdf")
		}
	else
		stop("error, unhandled attribute type!")

	if( rv$error != 0 ) {
		#---------------------------------------------
		# ? Got some strange error -- return as if the
		# attribute did not exist
		#---------------------------------------------
		retval$hasatt <- FALSE
		retval$value  <- 0
		return(retval)
		}

	retval$value <- rv$attribute
	return(retval)
}

#===============================================================
# Put an attribute into a netCDF file.  This puts the file into
# define mode, then takes it back out of define mode when done (unless
# definemode=TRUE, in which case the file is ASSUMED to be in define
# mode already, and is left in define mode as well).
#
# NOTE this does not work with IN-MEMORY ncdf files.
# (Note: you could extend this to work with in-memory files by
# having a list of variable (or file) attributes in the appropriate
# R classes, but I haven't bothered to do this since it does not
# come up very much for me.)
#
# "varid" can be an integer (R-type varid with counting starting
# at 1), object of class var.ncdf, or character string with a var's name.
# If varid==0, then a global attribute is set instead of a variable's
# attribute.
#
# Precision: ordinarily the precision (type) of the attribute will always
# follow the precision (type) of the var that this is an attribute of.
# However, you can explicitly override this by setting "prec" to the
# desired precision, which can be one of: short single double text( or character) integer.
# In the event of a global attribute, which of course has no associated
# variable, the storage mode of the passed attval will be used to
# determine the precision of the attribute, UNLESS "prec" is set.
# If "prec" is set, it always determines the created attribute type.
#
att.put.ncdf <- function( nc, varid, attname, attval, prec=NA, 
	verbose=FALSE, definemode=FALSE ) {

	if( verbose )
		print(paste("Making attribute",attname,"with value",attval,"for ncid=",nc$id))

	rvdim    <- vobjtodimname( nc, varid, verbose )
	isdimvar <- rvdim$isdim

	#-------------------------------------------------------
	# Can't do this if the file is in memory or not writable
	#-------------------------------------------------------
	if( (nc$filename == "IN-MEMORY") || (! nc$writable))
		stop("att.put.ncdf: the netcdf file has not been written to disk yet, or was not opened in write mode!")

	if( isdimvar ) {
		varid <- rvdim$name
		if( verbose ) print(paste("Passed obj is a dimension .. .using name=",varid))
		}

	varid <- vobjtovarid( nc, varid, allowdimvar=TRUE )

	if( varid == 0 )
		global <- TRUE
	else
		global <- FALSE
	if( verbose ) print(paste("varid to use:",varid))
	if( verbose ) print(paste("global:",global))

	#----------------------------------------------------------
	# Note there are TWO types here.  One is the storage mode
	# of the passed attval.  The netCDF routine to call is based
	# on this stoarge mode. The second type is the type of 
	# attribute to create.  This is passed as a parameter to
	# the netcDF routine.
	#----------------------------------------------------------

	#---------------------------------------------------------
	# Get the netCDF function to call ... this always depends
	# exclusively on the storage mode of the attval
	#---------------------------------------------------------
	if( storage.mode(attval) == "integer" ) 
		funcname <- "R_nc_put_att_int"
	else if( storage.mode(attval) == "double" )
		funcname <- "R_nc_put_att_double"
	else if( storage.mode(attval) == "character")
		funcname <- "R_nc_put_att_text"
	else
		stop(paste("att.put.ncdf: error, passed an attribute with a storage mode not handled.  Att name:",attname,"Att value:",attval,"Storage mode passed:",storage.mode(attval),".  Handled types: integer double character"))
	if( verbose ) print(paste("using function",funcname))

	#-------------------------------------------------------------
	# Get the type of attribute to create.  This follows the var's
	# type, in general, but can be manually overridden.
	#-------------------------------------------------------------
	atttypeShort <- 1  # These MUST match the values in the C code
	atttypeInt   <- 2 
	atttypeFloat <- 3
	atttypeDbl   <- 4
	atttypeText  <- 5
	atttypeByte  <- 6
	typetocreate <- -1
	if( (length(prec)==1) && is.na(prec) ) {
		if( global ) {
			prec <- storage.mode(attval)
			}
		else if( storage.mode(attval) == "character" ) {
			prec <- "character"
			}
		else
			{
			if( isdimvar ) {
				varidtouse <- nc$dim[[rvdim$name]]$dimvar
				if( varidtouse == -1 ) 
					prec <- storage.mode(attval) 
				else
					prec <- nc$var[[varidtouse]]$prec
				}
			else
				{
				vobj <- nc$var[[ nc$varid2Rindex[varid] ]]
				if( verbose )
					print(paste("getting precision from vobs for var",vobj$name))
				prec <- vobj$prec
				}
			}
		}

	if( verbose ) print(paste("prec to ccreate:",prec))
	if( (prec == "single") || (prec == "float"))
		typetocreate <- atttypeFloat
	else if( prec == "short" )
		typetocreate <- atttypeShort
	else if( prec == "byte" )
		typetocreate <- atttypeByte
	else if( prec == "double" )
		typetocreate <- atttypeDbl
	else if( (prec == "integer" ) || (prec == "int"))
		typetocreate <- atttypeInt
	else if( (prec == "text") || (prec == "character") || (prec == "char"))
		typetocreate <- atttypeText
	else
		stop(paste("Error in att.put.ncdf: unknown prec type specified:",prec,". Known values: short integer single double character"))

	if( ! definemode )
		redef.ncdf(nc)

	rv     <- list()
	rv$error <- -1
	rv <- .C(funcname,
		as.integer(nc$id),
		as.integer(varid-1),	# Change from R to C convention
		as.character(attname), 
		as.integer(typetocreate),
		as.integer(length(attval)),
		attval,
		error=as.integer(rv$error),
		PACKAGE="ncdf")
	if( rv$error != 0 ) {
		print(paste("Error in att.put.ncdf, while writing attribute",
			attname,"with value",attval))
		stop(paste("Error return from C call",funcname,"for attribute",attname))
		}

	if( ! definemode )
		enddef.ncdf(nc)
}

#===============================================================
# Internal use only.
#
varname.inq.ncdf <- function( nc, varid ) {

	if( class(nc) != "ncdf" ) 
		stop("Error: the first argument to varname.inq.ncdf (nc) is not of class ncdf!")

	#------------------------------------------------------
	# This is a string NC_MAX_LEN long, to provide storage.
	# I'm not sure if this is needed or not....
	#------------------------------------------------------
	str.nc.max.name <- "12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678"
	rv <- list()
	rv$varname <- str.nc.max.name
	rv$error   <- -1
	rv <- .C("R_nc_inq_varname",
		as.integer(nc$id),
		as.integer(varid-1),	# Change from R to C convention
		varname=as.character(rv$varname),
		error=as.integer(rv$error),
		PACKAGE="ncdf")
	if( rv$error != 0 ) {
		stop(paste("using ncid ",nc$id," varid ",varid))
		}
	return( rv$varname )
}

#=============================================================
# This takes as input a varid that may be the real R-type variable
# ID (an integer), a variable's name (a character string),
# or an object of class 'var.ncdf'.  It returns an integer
# variable id.  If the input varid is a NA, then we look to
# see if there is only one var in the file, and return the 
# varid of that if true, and generate an error if false.
# If "allowdimvar" is true, and if the passed string does not
# match a variable's name, then we ALSO check if it matches
# a dimvar's name.
#
vobjtovarid <- function( nc, varid, verbose=FALSE, allowdimvar=TRUE) {

	if( verbose )
		print(paste("vobjtovarid: entering with varid=",varid,sep=""))

	#---------------------------------------------------------------------------
	# If varid is NA, then return the only var in the file (if there IS only one
	# var in the file).  If there is more than one var in the file, return the
	# one with the most dimensions, IF that highest-dimensionality var has more
	# dimensions than any other var in the file.  Otherwise, generate an error.
	#---------------------------------------------------------------------------
	if( (length(varid)==1) && is.na(varid)) {
		if( nc$nvars == 1 ) {
			varToUse   <- 1
			}
		else
			{
			# Choose the most complicated var, if there is one, otherwise
			# halt with an error
			varToUse   <- -1
			ndimsItHas <- -1
			for( ii in 1:nc$nvars ) {
				if( nc$var[[ii]]$ndims > ndimsItHas ) {
					varToUse <- ii
					ndimsItHas <- nc$var[[ii]]$ndims
					}
				}
			for( ii in 1:nc$nvars ) {
				if( (ii != varToUse) && (nc$var[[ii]]$ndims == ndimsItHas)) {
					stop(paste("File",nc$filename,"has more than one variable, so you must explicitly specify which one you want"))
					}
				}
			}
		if( verbose )
			print(paste("vobjtovarid: returning with only var in file; id=",nc$var[[varToUse]]$id))
		return( nc$var[[varToUse]]$id )
		}

	origvarid <- varid
	if( ! is.numeric(varid) ) {
		if( is.character(varid)) {	# we were given a variable's name
			origvarid <- varid
			#--------------------------------------------
			# See if any vars in this file have this name
			#--------------------------------------------
			varidInteger <- -1
			for( kk in 1:nc$nvars ) {
				if( origvarid == nc$var[[kk]]$name ) 
					varidInteger <- nc$var[[kk]]$id
				}

			if( varidInteger != -1 ) {
				if(verbose)
					print(paste("Variable named",origvarid,"found in file with varid=",varidInteger))
				varid <- varidInteger
				}
			else
				{
				if(verbose)
					print(paste("Variable named",origvarid,"NOT found in file; looking for a dimvar with this name"))
				#---------------------------------------------------------------
				# A var with this name was NOT found in the file.  But, it could
				# be the name of a dimvar in the file.  Check to see if we are
				# allowed to return dimvars in this case.
				#---------------------------------------------------------------
				if( allowdimvar ) {
					for( i in 1:nc$ndims )
						if( origvarid == nc$dim[[i]]$name ) {
							#---------------------
							# Yes, it IS a dimvar!
							#---------------------
							varid <- nc$dim[[i]]$dimvarid 
							if( verbose )
								print(paste("vobjtovarid: returning with DIMvarid deduced from name; varid=",varid))
							return(varid)
							}
					}
				print("vobjtovarid: error: I could not find the var whose name was passed in the file!")
				print(paste("var name:",origvarid))
				print(paste("file name:",nc$filename))
				if( ! allowdimvar ) 
					print("Note: I was NOT allowed to check to see if this was a dimvar name")
				print(nc)
				stop("Variable not found")
				}
			if( verbose )
				print(paste("vobjtovarid: returning with varid deduced from name; varid=",varid))
			}
		else if( class(varid) == "var.ncdf" ) {
			if(verbose)
				print(paste("vobjtovarid: passed a var.ncdf class, name=",varid$name))
			varid <- nc$var[[varid$name]]$id # Note we do NOT use varid$id in case var is from different file (but names are same)
			varidOK <- ((varid>=0) && (varid<=100000))
			if( is.na(varidOK) || (!varidOK)) {
				print("vobjtovarid: I was passed a var.ncdf object, BUT this object does NOT refer to any valid var in the netcdf file!")
				print(paste("This happened for netCDF filename:",nc$filename))
				print("Here are the vars in the netCDF file:")
				for( ii in 1:nc$nvars )
					print(paste(ii,": ",nc$var[[ii]]$name, sep='' ))
				print(paste("The passed varid object (which does NOT exist in that file) is:"))
				print(origvarid)
				print(paste("Hint: make SURE the variable was not only defined with a call to def.var.ncdf, but also included in the list passed to create.var.ncdf"))
				stop("stopping")
				}
			
			if(verbose)
				print(paste("vobjtovarid: returning varid=",varid,"  (OK=",varidOK,")"))
			}
		else if( allowdimvar && (class(varid) == "dim.ncdf") ) {
			if(verbose)
				print(paste("vobjtovarid: passed a dim.ncdf class, name=",varid$name))
			varid <- nc$var[[varid$name]]$id # Note we do NOT use varid$id in case var is from different file (but names are same)
			varidOK <- ((varid>=0) && (varid<=100000))
			if( is.na(varidOK) || (!varidOK))
				stop("vobjtovarid: I was passed a dim.ncdf object, BUT this object does NOT refer to any valid dimvar in the netcdf file!")
			if(verbose)
				print(paste("vobjtovarid: returning varid=",varid))
			}
		else
			{
			print("Error in vobjtovarid: second argument (varid) must be a integer varid, character string variable name, or var.ncdf object")
			print("here is what I was passed:")
			print(varid)
			print("--------------------------")
			stop("second argument (varid) must be a integer varid, character string variable name, or var.ncdf object")
			}
		}

	return(varid)
}

#===============================================================
# Internal use only.  Turns an integer varid in R-format (1-based
# counting) into a var.ncdf object.
#
var.inq.ncdf <- function( nc, varid ) {

	if( class(nc) != "ncdf" ) 
		stop("Error: the first argument to var.inq.ncdf (nc) is not of class ncdf!")

	#------------------------------------------------------
	# This is a string NC_MAX_LEN long, to provide storage.
	# I'm not sure if this is needed or not....
	#------------------------------------------------------
	str.nc.max.name <- "12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678"

	rv <- list()
	rv$name    <- str.nc.max.name
	rv$varlen  <- -1
	rv$error   <- -1
	rv$type    <- -1
	rv$ndims   <- -1
	rv$natts   <- -1
	rv$precint <- -1 # INTEGER (not character) form of precision. 1=SHORT, 2=INT, 3=FLOAT, 4=DOUBLE, 5=CHAR 6=BYTE.  Must match C code values!!
	rv$dimids  <- integer(varndims.ncdf( nc, varid ))
	rv <- .C("R_nc_inq_var",
		as.integer(nc$id),
		as.integer(varid-1),	# Change from R to C convention
		name=as.character(rv$name),
		type=as.integer(rv$type),
		ndims=as.integer(rv$ndims),
		dimids=as.integer(rv$dimids),
		natts=as.integer(rv$natts),
		precint=as.integer(rv$precint),
		error=as.integer(rv$error),
		PACKAGE="ncdf")
	if( rv$error != 0 )
		stop("call to C function R_nc_inq_var failed")
	
	var       <- list()
	var$id    <- varid
	var$name  <- rv$name
	var$ndims <- rv$ndims
	var$natts <- rv$natts
	var$size  <- varsize.ncdf( nc, varid )
	if( rv$precint == 1 )
		var$prec <- "short"
	else if( rv$precint == 2 )
		var$prec <- "int"
	else if( rv$precint == 3 )
		var$prec <- "float"
	else if( rv$precint == 4 )
		var$prec <- "double"
	else if( rv$precint == 5 )
		var$prec <- "char"
	else if( rv$precint == 6 )
		var$prec <- "byte"
	else
		stop(paste("Error, unrecognized type code of variable returned from C call:",rv$precint,". I currently know about the following types: byte short int float double char."))

	#---------------------------------------
	# Convert dimids from C to R conventions
	#---------------------------------------
	var$dimids <- rv$dimids + 1
	var$dimids <- var$dimids[ var$ndims:1 ]

	#--------------------------
	# Get this var's attributes
	#--------------------------
	attu <- att.get.ncdf( nc, varid, "units" )
	if( attu$hasatt )
		var$units <- attu$value
	else
		var$units <- ""
	attu <- att.get.ncdf( nc, varid, "long_name" )
	if( attu$hasatt ) 
		var$longname <- attu$value
	else
		var$longname <- var$name

	return(var)
}

#===============================================================
# NOTE that this does NOT return a full-fledged "dim.ncdf" object,
# it is the lower-level interface and only returns the portions
# of a dim.ncdf object that directly correspond to a real netCDF
# dimension.  The dim.ncdf object, by contrast, also has information
# from the dimvar.
#
# Internal use only.
#
dim.inq.ncdf <- function( nc, dimid ) {

	if( class(nc) != "ncdf" ) 
		stop("Error: the first argument to dim.inq.ncdf (nc) is not of class ncdf!")

	#------------------------------------------------------
	# This is a string NC_MAX_LEN long, to provide storage.
	# I'm not sure if this is needed or not....
	#------------------------------------------------------
	str.nc.max.name <- "12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678"

	rv <- list()
	rv$dimname <- str.nc.max.name
	rv$dimlen  <- -1
	rv$error   <- -1
	rv <- .C("R_nc_inq_dim",
		as.integer(nc$id),
		as.integer(dimid-1),	# Change from R to C convention
		dimname=as.character(rv$dimname),
		dimlen=as.integer(rv$dimlen),
		error=as.integer(rv$error),
		PACKAGE="ncdf")
	if( rv$error != 0 ) {
		stop(paste("using ncid ",nc$id," dimid ",dimid))
		}
	d <- list()
	d$name  <- rv$dimname
	d$len   <- rv$dimlen

	#---------------------------------
	# Is this the unlimited dimension?
	#---------------------------------
	if( dimid == unlimdim.ncdf(nc))
		d$unlim <- TRUE
	else
		d$unlim <- FALSE

	return(d)
}

#===============================================================
# This returns the dimid of the unlimited dimension, or -1 if 
# there is NO unlimited dimension in the file.  The dimid is
# returned in R conventional 1-based counting.
#
unlimdim.ncdf <- function( nc ) {
	if( class(nc) != "ncdf")
		stop("unlimdim.ncdf passed something that is NOT class ncdf!")

	rv            <- list()
	rv$unlimdimid <- -1
	rv$error      <- -1
	rv <- .C("R_nc_inq_unlimdim",
		as.integer(nc$id),
		unlimdimid=as.integer(rv$unlimdimid),
		error=as.integer(rv$error),
		PACKAGE="ncdf")
	if( rv$error != 0 ) 
		stop("unlimdim.ncdf returned error")

	if( rv$unlimdimid != -1 )
		rv$unlimdimid <- rv$unlimdimid + 1	# Switch from C to R convention
	return(rv$unlimdimid)
}

#=====================================================================================
# This turns a passed 'var' object into the dimension name, IFF it refers to
# a dimvar rather than a regular var.  
# Return value:
#	a list with logical 'isdim', and character string 'name'
#
vobjtodimname <- function( nc, varid, verbose=FALSE ) {

	if( verbose ) 
		print(paste("entering vobjtodimname with varid=",varid))

	if( (length(varid)==1) && is.na(varid)) {
		if( verbose ) 
			print(paste("vobjtodimname: is NA"))
		retval <- list()
		retval$isdim <- FALSE
		return(retval)
		}

	if( ((length(varid)>1)||! is.na(varid)) && (! is.numeric(varid))     && 
	    (! is.character(varid)) && (! is.null(class(varid))) &&
	    ( class(varid) == "dim.ncdf") ) {
	    	if( verbose ) {
			print("vobjtodimname: passed varid is a dim.ncdf object!")
			print(varid)
			}
		retval <- list()
		retval$isdim <- TRUE
		retval$name  <- varid$name  # NOTE: this is actually the dim name since varid is a dim
		return(retval)
		}

	if( is.character(varid) ) {
		if( verbose ) 
			print(paste("vobjtodimname: is a character type varid.  This file has",nc$ndims,"dims"))
		#----------------------------------------------------------------------
		# See if this string is a dim's name (note carefully: NOT a var's name)
		#----------------------------------------------------------------------
		for( i in 1:nc$ndims )
			if( varid == nc$dim[[i]]$name ) {
				#---------------------
				# Yes, it IS a dimvar!
				#---------------------
				if( verbose )
					print("put.var.ncdf: passed varid is the name of a dimvar!")
				retval <- list()
				retval$isdim <- TRUE
				retval$name  <- varid
				return(retval)
				}
		}

	if( verbose ) 
		print(paste("vobjtodimname: no cases found, returning FALSE"))

	retval <- list()
	retval$isdim <- FALSE
	return(retval)
}
	
#===============================================================
# Writes a vector of values to a netCDF file.  'start' and 'count'
# are given in R convention, i.e., starting at 1, and
# in XYZT order.  If they are omitted, the full array is written.
# 'varid' can be the variable's name, an object of class var.ncdf,
# or the integer varid.  If varid is NA, then the "only" var in 
# the file is assumed to be the selected one.  Values that are NA's 
# in the input data are converted to that variable's 'missing_value' 
# attribute before being written out to the file.
#
put.var.ncdf <- function( nc, varid=NA, vals=NA, start=NA, count=NA, verbose=FALSE ) {

	if( class(nc) != "ncdf" )
		stop("first argument is not of class ncdf!")

	if( (length(vals)==1) && is.na(vals) ) 
		stop("requires a vals argument to be set, and not NA (to set a single NA, use c(NA))")
		
	if( verbose ) {
		print(paste("put.var.ncdf: entering with filename",
			nc$filename," and varid:"))
		print(varid)
		}

	if( ! nc$writable ) 
		stop(paste("trying to write to file",nc$filename,"but it was not opened with write=TRUE"))

	#-------------------------------------------------------------
	# First check to see if we are ACTUALLY putting dimvar values.  
	#-------------------------------------------------------------
	if( verbose )
		print('Checking to see if passed varid is ACTUALLY a dimension')
	rvdim    <- vobjtodimname( nc, varid, verbose )
	isdimvar <- rvdim$isdim
	if( verbose ) {
		if( isdimvar ) 
			print('...YES, passed obj WAS a dimension')
		else
			print('...NO, passed obj was NOT a dimension')
		}

	if( isdimvar ) {
		if(verbose) print("put.var.ncdf: putting values to a DIMVAR") 
		varid <- nc$dim[[rvdim$name]]$dimvarid
		if( is.null(varid) || (varid == -1)) {
			print("put.var.ncdf: error -- trying to write values to a dimvar, BUT the file does not have the")
			print("specified dimvar in it! (It probably has the dimension, but not the dimension *variable*)" );
			print(paste("Passed nc filename:",nc$filename))
			print(paste("Passed varid:",varid))
			stop("Dimvar does not exist")
			}
		}
	else
		{
		if(verbose) print("put.var.ncdf: about to call vobjtovarid")
		varid <- vobjtovarid( nc, varid, verbose=verbose )
		if(verbose) print(paste("put.var.ncdf: vobjtovarid returned: >",
					varid,"<"))
		}

	if( verbose )
		print(paste("put.var.ncdf: ending up using varid=",varid))

	varsize <- varsize.ncdf ( nc, varid )
	ndims   <- varndims.ncdf( nc, varid )

	#--------------------------------------------------------
	# Fix up start and count to use (in R convention for now)
	#--------------------------------------------------------
	if( (length(start)==1) && is.na(start) )
		start <- rep(1,ndims)	# Note: use R convention for now
	else
		{
		if( length(start) != ndims ) 
			stop(paste("'start' should specify",ndims,
				"dims but actually specifies",length(start)))
		}
	if( (length(count)==1) && is.na(count)) 
		count <- varsize - start + 1	
	else
		{
		if( length(count) != ndims ) 
			stop(paste("'count' should specify",ndims,
				"dims but actually specifies",length(count)))
		count <- ifelse( (count == -1), varsize-start+1, count)
		}

	#------------------------------
	# Switch from R to C convention
	#------------------------------
	c.start <- start[ ndims:1 ] - 1
	c.count <- count[ ndims:1 ]

	#--------------------------------------------
	# Change NA's to the variable's missing value
	#--------------------------------------------
	if( verbose )
		print("about to change NAs to variables missing value")
	if( isdimvar )
		mv <- default.missval.ncdf()
	else
		mv <- varid.to.missing.value( nc, varid )
	vals <- ifelse( is.na(vals), mv, vals)

	#---------------------------------
	# Get the correct type of variable
	#---------------------------------
	precint <- vartype.ncdf( nc, varid ) # 1=short, 2=int, 3=float, 4=double, 5=char, 6=byte
	if( verbose )
		print(paste("Putting var of type",precint," (1=short, 2=int, 3=float, 4=double, 5=char, 6=byte)"))

	#----------------------------------------------------------
	# Sanity check to make sure we have at least as many values 
	# in the data array as we are writing.  Chars are a special
	# case because typically they are defined with an extra
	# "nchar" dim that is not included in the passed array.
	#----------------------------------------------------------
	n2write <- prod(count)
	if( (precint != 5) && (length(vals) != n2write)) {
		if( length(vals) > n2write ) 
			print(paste("put.var.ncdf: warning: you asked to write",n2write,
				"values, but the passed data array has",length(vals),
				"entries!"))
		else
			stop(paste("put.var.ncdf: error: you asked to write",n2write,
				"values, but the passed data array only has",length(vals),
				"entries!"))
		}

	rv <- list()
	rv$error <- -1

	if( verbose ) {
		print("put.var.ncdf: calling C routines with C-style count=")
		print(c.count)
		print("and C-style start=")
		print(c.start)
		}
	if( (precint == 1) || (precint == 2) || (precint == 6)) {
		#--------------------
		# Short, Int, or Byte
		#--------------------
		rv <- .C("R_nc_put_vara_int", 
			as.integer(nc$id),
			as.integer(varid-1),	# Switch from R to C convention
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			data=as.integer(vals),
			error=as.integer(rv$error),
			PACKAGE="ncdf")
		if( rv$error != 0 ) 
			stop("C function R_nc_put_var_int returned error")
		}

	else if( (precint == 3) || (precint == 4)) {
		#----------------
		# Float or double
		#----------------
		rv <- .C("R_nc_put_vara_double", 
			as.integer(nc$id),
			as.integer(varid-1),	# Switch from R to C convention
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			data=as.double(vals),
			error=as.integer(rv$error),
			PACKAGE="ncdf")
		if( rv$error != 0 ) 
			stop("C function R_nc_put_var_double returned error")
		}

	else if( precint == 5 ) {
		#----------
		# Character
		#----------
		rv <- .C("R_nc_put_vara_text", 
			as.integer(nc$id),
			as.integer(varid-1),	# Switch from R to C convention
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			data=as.character(vals),
			error=as.integer(rv$error),
			PACKAGE="ncdf")
		if( rv$error != 0 ) 
			stop("C function R_nc_put_var_double returned error")
		}

	else
		stop(paste("Internal error in put.var.ncdf: unhandled variable type=",precint,". Types I know: 1=short 2=int 3=float 4=double 5=char"))
}

#===============================================================
# Returns data values from a netCDF file.  'start' and 'count'
# are given in R convention, i.e., starting at 1, and
# in XYZT order.  If they are omitted, the full array is read in.
# 'varid' can be the variable's name, an object of class var.ncdf,
# or the integer varid.  'varid' can also be omitted entirely,
# in which case the only variable in the file is identified and
# that one read in (if no such variable can be identified, an
# error is generated).  The 'count' array can have -1's, which 
# indicate that all the values in that dim are to be read (subject
# to the start array).  Missing values in the source file (i.e.,
# values that match that variable's 'missing_value' attribute)
# are set to NA's.
# Argument 'signedbyte' can be TRUE for bytes to be interpreted as 
# signed, or FALSE to be unsigned.
#
get.var.ncdf <- function( nc, varid=NA, start=NA, count=NA, verbose=FALSE, signedbyte=TRUE, forcevarid=NA ) {

	if( verbose ) {
		if( !is.na(forcevarid)) {
			print("get.var.ncdf: entering with forcevarid set to:")
			print(forcevarid)
			}
		else
			{
			print("get.var.ncdf: entering. Here is varid:")
			print(varid)
			}
		}

	have_start = (length(start)>1) || ((length(start)==1) && (!is.na(start)))
	have_count = (length(count)>1) || ((length(count)==1) && (!is.na(count)))

	if( class(nc) != "ncdf" )
		stop("first argument (nc) is not of class ncdf!")

	if( signedbyte )
		byte_style = 1	# 1=signed
	else
		byte_style = 2	# 2=unsigned

	if( is.na(forcevarid) ) {
		#-------------------------------------------------------------
		# First check to see if we are ACTUALLY getting dimvar values.  
		#-------------------------------------------------------------
		if( verbose ) 
			print("checking to see if passed varid is actually a dimvar")
		rvdim    <- vobjtodimname( nc, varid, verbose )
		isdimvar <- rvdim$isdim
		if( isdimvar )
			dimidtouse <- rvdim$name

		if( verbose )
			print(paste("get.var.ncdf: isdimvar:",isdimvar))

		if( isdimvar ) {
			if( verbose )
				print(paste("get.var.ncdf: dimname:",dimidtouse))
			varid <- nc$dim[[dimidtouse]]$dimvarid
			#--------------------------------------------------------
			# Here we return default integers for dims with no dimvar
			#--------------------------------------------------------
			if( varid == -1 ) {
				if( ! have_start )
					start <- 1
				if( ! have_count )
					count <- nc$dim[[dimidtouse]]$len
				if( count == 1 )
					return( start )
				else
					return( start:(start+count-1) )
				}
			}
		else
			{
			varid <- vobjtovarid  ( nc, varid, verbose=verbose )
			}
		}
	else
		{
		varid <- forcevarid
		isdimvar <- TRUE
		}

	if( verbose ) 
		print(paste("get.var.ncdf: ending up using varid=",varid))
	varsize <- varsize.ncdf ( nc, varid )
	ndims   <- varndims.ncdf( nc, varid )
	if( verbose ) {
		print(paste("ndims:",ndims))
		print("get.var.ncdf: varsize:")
		print(varsize)
		}

	#------------------------------
	# Fix up start and count to use
	#------------------------------
	if( ndims == 0 ) {
		start <- 1
		count <- 1
		}
	else
		{
		if( ! have_start )
			start <- rep(1,ndims)	# Note: use R convention for now
		if( ! have_count )
			count <- varsize - start + 1	
		else
			{
			#------------------
			# Take care of -1's
			#------------------
			count <- ifelse( (count == -1), varsize-start+1, count)
			}
		}
	if( verbose ) {
		print("get.var.ncdf: start:")
		print(start)
		print("get.var.ncdf: count:")
		print(count)
		}

	if( ndims > 0 ) {
		if( length(start) != ndims ) 
			stop(paste("Error: variable has",ndims,"dims, but start has",length(start),"entries.  They must match!"))
		if( length(count) != ndims ) 
			stop(paste("Error: variable has",ndims,"dims, but count has",length(count),"entries.  They must match!"))
		}

	#----------------------------------------
	# Need to know how much space to allocate
	#----------------------------------------
	totvarsize <- prod(count)
	if( verbose )
		print(paste("get.var.ncdf: totvarsize:",totvarsize))
	
	#------------------------------
	# Switch from R to C convention
	#------------------------------
	c.start <- start[ ndims:1 ] - 1
	c.count <- count[ ndims:1 ]

	rv <- list()
	rv$error <- -1

	#---------------------------------
	# Get the correct type of variable
	#---------------------------------
	precint <- vartype.ncdf( nc, varid ) # 1=short, 2=int, 3=float, 4=double, 5=char, 6=byte
	if( verbose )
		print(paste("Getting var of type",precint," (1=short, 2=int, 3=float, 4=double, 5=char, 6=byte)"))
	if( (precint == 1) || (precint == 2) || (precint == 6)) {
		#--------------------
		# Short, Int, or Byte
		#--------------------
		rv$data  <- integer(totvarsize)
		rv <- .C("R_nc_get_vara_int", 
			as.integer(nc$id),
			as.integer(varid-1),	# Switch from R to C convention
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			as.integer(byte_style), # 1=signed, 2=unsigned
			data=as.integer(rv$data),
			error=as.integer(rv$error),
			PACKAGE="ncdf",
			DUP=FALSE)
		if( rv$error != 0 ) 
			stop("C function R_nc_get_var_int returned error")
		}
	else if( (precint == 3) || (precint == 4)) {
		#----------------
		# Float or double
		#----------------
		rv$data  <- double(totvarsize)
		rv <- .C("R_nc_get_vara_double", 
			as.integer(nc$id),
			as.integer(varid-1),	# Switch from R to C convention
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			data=as.double(rv$data),
			error=as.integer(rv$error),
			PACKAGE="ncdf",
			DUP=FALSE)
		if( rv$error != 0 ) 
			stop("C function R_nc_get_vara_double returned error")
		}
	else if( precint == 5 ) {
		strndims <- ndims - 1
		strlen   <- count[1] + 1
		strdim   <- 1
		if( strndims >= 1 ) {
			strdim <- count[2:ndims]
			nstr   <- prod(strdim)
			}
		else
			nstr <- 1
		if(verbose)
			print(paste("ndims:",ndims,"strndims:",strndims,"strlen:",strlen,"nstr:",nstr))

		#----------------------------------------------
		# Make a character string of the specified size
		#----------------------------------------------
		stor     <- blankstring.ncdf( totvarsize )
		stordata <- blankstring.ncdf(strlen)
		if( verbose )
			print(paste("length of stor string:",nchar(stor)))
		rv$tempstore <- stor
		rv$data      <- array(stordata, dim=strdim)

		rv <- .C("R_nc_get_vara_text", 
			as.integer(nc$id),
			as.integer(varid-1),	# Switch from R to C convention
			as.integer(c.start),	# Already switched to C convention...
			as.integer(c.count),	# Already switched to C convention...
			tempstore=as.character(rv$tempstore),
			data=as.character(rv$data),
			error=as.integer(rv$error),
			PACKAGE="ncdf")
		if( rv$error != 0 ) 
			stop("C function R_nc_get_var_text returned error")

		dim(rv$data) <- strdim
		}
	else
		{
		stop(paste("Trying to get variable of an unhandled type code: ",precint))
		}
	if( verbose )
		print(paste("get.var.ncdf: C call returned",rv$error))

	#--------------------------------------------------------
	# Set our dims...but collapse degenerate dimensions first
	#--------------------------------------------------------
	if( ndims > 0 ) {
		count.nodegen <- vector()
		foundone <- 0
		for( i in 1:ndims )
			if( count[i] > 1 ) {
				count.nodegen <- append(count.nodegen, count[i])
				foundone <- 1
				}
		if( foundone == 0 ) 
			dim(rv$data) <- (1)
		else
			{
			if( verbose )
				print(paste("count.nodegen:",count.nodegen,"   Length of data:",length(rv$data)))
			if( precint != 5 )
				dim(rv$data) <- count.nodegen
			}
		if( verbose ) {
			print("get.var.ncdf: final dims of returned array:")
			print(dim(rv$data))
			}
		}

	if( verbose ) {
		print(paste("varid:",varid))
		print(paste("nc$varid2Rindex:",nc$varid2Rindex))
		print(paste("nc$varid2Rindex[varid]:",nc$varid2Rindex[varid]))
		}

	#----------------------------------------------------------
	# Change missing values to "NA"s.  Note that 'varid2Rindex'
	# is NOT filled out for dimvars, so skip this if a dimvar
	#----------------------------------------------------------
	if( (!isdimvar) && (precint != 5)) {
		if( verbose ) 
			print("get.var.ncdf: setting missing values to NA")
		if( (precint==1) || (precint==2)) {
			#---------------------
			# Short, Int, and Byte
			#---------------------
			mv  <- nc$var[[ nc$varid2Rindex[varid] ]]$missval
			if( ! is.na(mv) ) {
				if( verbose )
					print(paste("missval:",mv))
				rv$data[rv$data==mv] <- NA
				}
			}
		else if( (precint==3) || (precint==4)) {
			#-----------------
			# Float and Double
			#-----------------
			mv  <- nc$var[[ nc$varid2Rindex[varid] ]]$missval
			if( ! is.na(mv) ) {
				tol <- abs(mv*1.e-5)
				if( verbose )
					print(paste("missval:",mv,"  tol:",tol))
				rv$data[abs(rv$data-mv)<tol] <- NA
				}
			}
		}

	#--------------------------------------
	# Implement add_offset and scale_factor
	#--------------------------------------
	if( ! isdimvar ) {
		if( verbose ) 
			print(paste("get.var.ncdf: implementing add_offset (",
				nc$var[[ nc$varid2Rindex[varid] ]]$hasAddOffset,
				") and scale_factor (",
				nc$var[[ nc$varid2Rindex[varid] ]]$hasScaleFact, ")" ))
		if( nc$var[[ nc$varid2Rindex[varid] ]]$hasAddOffset &&
		    nc$var[[ nc$varid2Rindex[varid] ]]$hasScaleFact ) {
		    	if(verbose)
				print(paste("var has BOTH add_offset (",
					nc$var[[ nc$varid2Rindex[varid] ]]$addOffset,
					") and scale_fact (",
					 nc$var[[ nc$varid2Rindex[varid] ]]$scaleFact, ")" ))
			rv$data <- rv$data * nc$var[[ nc$varid2Rindex[varid] ]]$scaleFact +
				nc$var[[ nc$varid2Rindex[varid] ]]$addOffset
			}

		else if( nc$var[[ nc$varid2Rindex[varid] ]]$hasAddOffset ) {
		    	if(verbose)
				print(paste("var has add_offset (only):",
					nc$var[[ nc$varid2Rindex[varid] ]]$addOffset))
			rv$data <- rv$data + nc$var[[ nc$varid2Rindex[varid] ]]$addOffset
			}

		else if( nc$var[[ nc$varid2Rindex[varid] ]]$hasScaleFact ) {
		    	if(verbose)
				print(paste("var has scale_factor (only):",
					nc$var[[ nc$varid2Rindex[varid] ]]$scaleFact))
			rv$data <- rv$data * nc$var[[ nc$varid2Rindex[varid] ]]$scaleFact
			}
		else
			{
		    	if(verbose)
				print("var has NEITHER add_offset nor scale_factor")
			}
		}

	return(rv$data)
}

#===============================================================
# Returns a vector of the size of the variable, in
# R order (XYZT).
#
# Internal use only.  Use v$varsize (where v is an object of
# class "var.ncdf") if you want this info.
#
varsize.ncdf <- function( nc, varid ) {
	
	if( class(nc) != "ncdf" )
		stop("the first passed argument (nc) does not have class ncdf!")

	ndims <- varndims.ncdf( nc, varid )
	if( ndims == 0 )
		return(vector())

	rv         <- list()
	rv$error   <- -1
	rv$varsize <- integer(ndims)
	rv <- .C("R_nc_varsize",
		as.integer(nc$id),
		as.integer(varid-1),	# Switch from R to C convention
		varsize=as.integer(rv$varsize),
		error=as.integer(rv$error),
		PACKAGE="ncdf")
	if( rv$error != 0 ) 
		stop("error returned from C routine R_nc_varsize")

	#-------------------------------------
	# Switch order from C-style to R-style
	#-------------------------------------
	rv$varsize <- rv$varsize[ ndims:1 ]

	return(rv$varsize)
}

#===============================================================
# Internal use only.   Input: integer varid.  Output: one of the
# integer R type codes (1=short, 2=int, 3=float, 4=double,
# 5=char, 6=byte).
#
vartype.ncdf <- function( nc, varid ) {

	rv         <- list()
	rv$error   <- -1
	rv$precint <- -1

	rv <- .C("R_nc_inq_vartype", 
		as.integer(nc$id),
		as.integer(varid-1),	# Switch from R to C convention
		precint=as.integer(rv$precint),
		error=as.integer(rv$error),
		PACKAGE="ncdf")
	if( rv$error != 0 ) 
		stop("error returned from C call")
	return( rv$precint )
}

#===============================================================
# Internal use only.  Use v.ndims if you want this info.
#
varndims.ncdf <- function( nc, varid ) {
	rv <- list()
	rv$error <- -1
	rv$ndims <- -1
	rv <- .C("R_nc_inq_varndims", 
		as.integer(nc$id),
		as.integer(varid-1),	# Switch from R to C convention
		ndims=as.integer(rv$ndims),
		error=as.integer(rv$error),
		PACKAGE="ncdf")
	if( rv$error != 0 ) 
		stop("error returned from C call")
	return( rv$ndims )
}

#===============================================================
sync.ncdf <- function( nc ) {
	.C("R_nc_sync", as.integer(nc$id), PACKAGE="ncdf")
}

#===============================================================
redef.ncdf <- function( nc ) {
	.C("R_nc_redef", as.integer(nc$id), PACKAGE="ncdf")
}

#===============================================================
# Returns -1 if the dim is NOT found in the file, and the
# dimid of the dim otherwise.
#
dimid.inq.ncdf <- function( nc, dimname ) {
	rv       <- list()
	rv$dimid <- -1
	rv <- .C("R_nc_inq_dimid", 
		as.integer(nc$id),
		as.character(dimname),
		dimid=as.integer(rv$dimid),
		PACKAGE="ncdf")
	if( rv$dimid != -1 )
		rv$dimid <- rv$dimid + 1	# Switch from C to R convention
	return(rv$dimid)
}

#===============================================================
# Returns -1 if the var is NOT found in the file, and the
# varid of the var otherwise.
#
varid.inq.ncdf <- function( nc, varname ) {
	rv       <- list()
	rv$varid <- -1
	rv <- .C("R_nc_inq_varid", 
		as.integer(nc$id),
		as.character(varname),
		varid=as.integer(rv$varid),
		PACKAGE="ncdf")
	if( rv$varid != -1 )
		rv$varid <- rv$varid + 1	# Switch from C to R convention
	return(rv$varid)
}

#===============================================================
enddef.ncdf <- function( nc ) {
	.C("R_nc_enddef", as.integer(nc$id), PACKAGE="ncdf")
	sync.ncdf( nc )
}

#===============================================================
varid.to.missing.value <- function( nc, varid ) {
	idx <- nc$varid2Rindex[varid]
	return( nc$var[[idx]]$missval )
}

#===============================================================
close.ncdf <- function( con, ... ) {
	if(class(con) != "ncdf")
		stop("Error, close.ncdf passed something NOT of class ncdf!")
	.C("R_nc_close", as.integer(con$id), PACKAGE="ncdf")
}

#dyn.load("ncdf.so")

#x <- dim.def.ncdf( "Lon", "degreesE", 0.5:359.5)
#y <- dim.def.ncdf( "Lat", "degreesN", as.double(-89:89))
#t <- dim.def.ncdf( "Time", "days since 1900-01-01", 1:10, unlim=TRUE)
#temp  <- var.def.ncdf( "Temperature", "degC", list(x,y,t), 1.e30 )
#salin <- var.def.ncdf( "Salinity",    "ppt",  list(x,y,t), 1.e30 )
#ncnew <- create.ncdf( "test2.nc", list(temp,salin), verbose=TRUE)

# Note: to put the value of an unlimited dim, use:
#put.var.ncdf( ncid, DIMNAME, value, start=(), count=() )

#nc <- open.ncdf("t.atm.nc")
#var <- nc$var[["SNOWH"]]
#vals <- get.var.ncdf(nc,var)
#ncnew <- create.ncdf("test9.nc",var)
#put.var.ncdf(ncnew,var,vals)
#close.ncdf(ncnew)

# Test making different precision dims and vars
#dimint <- dim.def.ncdf("IntegerDim", "count", as.integer(1:10))
#dimdbl <- dim.def.ncdf("DblDim", "count", c(1,2,3) )
#intvar <- var.def.ncdf( "IntVar", "m", list(dimint, dimdbl), as.integer(999), prec="integer")
#floatvar1 <- var.def.ncdf( "FloatVar1", "m", list(dimint, dimdbl), 1.e30)
#floatvar2 <- var.def.ncdf( "FloatVar2", "m", list(dimint, dimdbl), 1.e30)
#dblvar <- var.def.ncdf( "DblVar", "m", list(dimint, dimdbl), 1.e30, prec="double")
#ncidnew <- create.ncdf("testprec.nc", list(intvar, floatvar1, floatvar2, dblvar), verbose=TRUE)

# Test putting attributes of various types
#att.put.ncdf( ncidnew, intvar, "intattdefault", 3. )
#att.put.ncdf( ncidnew, floatvar1, "floatattdefault", 3. )
#att.put.ncdf( ncidnew, dblvar, "dblattdefault", 3. )
#
#att.put.ncdf( ncidnew, dblvar, "intattforce", 3, prec="integer" )
#att.put.ncdf( ncidnew, dblvar, "floatattforce", 3., prec="single" )
#att.put.ncdf( ncidnew, intvar, "dblattforce", 3., prec="double" )
#
#att.put.ncdf( ncidnew, intvar, "textxtt", "This is a tesst attribute") 
#att.put.ncdf( ncidnew, 0, "title", "This is a test TITLE global attribute")
#
#close.ncdf(ncidnew)
