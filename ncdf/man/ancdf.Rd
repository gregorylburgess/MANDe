\name{ncdf}
\alias{ncdf}
\title{Read, write, and create netCDF files (v1.6)}
\description{
 Read from or write to existing netCDF format files, or create new ones.
}
\references{
 http://www.unidata.ucar.edu/packages/netcdf/
}
\details{
 The netCDF data file format from Unidata is a platform-independent, binary
 file that also contains metadata describing the contents and format of the data
 in the file.  NetCDF files contain one or more variables, which are 
 structured as regular N-dimensional arrays.  They also contain dimensions, which 
 describe the extent of the variables' arrays.  Data can be read from
 or written to variables in arbitrary hyperslabs.  The R package 'ncdf' allows
 reading from, writing to, and creation of netCDF files.  Note that the
 netCDF library must already be installed on your machine for this R interface
 to the library to work.

 If you are absolutely new to netCDF files, they can be a little overwhelming,
 so here is a brief sketch of what documentation you need to read next.  
 
 If you want to READ data 
 from an already-existing netCDF file, first call \code{\link[ncdf]{open.ncdf}} 
 to open the file, then call \code{\link[ncdf]{get.var.ncdf}} to read the data
 from the file.  
 
 If you want to WRITE data to a new netCDF file, first call
 \code{\link[ncdf]{dim.def.ncdf}} to define the dimensions that your data
 exists along (for example, perhaps latitude and longitude), then call
 \code{\link[ncdf]{var.def.ncdf}} to define a variable in the netCDF
 file that will hold your data, then call \code{\link[ncdf]{create.ncdf}}
 to create the netCDF file, then call \code{\link[ncdf]{put.var.ncdf}}
 to write your data to the newly created netCDF file.

 This is version 1.6 of the ncdf library.
}
\seealso{ 
 \code{\link[ncdf]{att.put.ncdf}}, \code{\link[ncdf]{att.get.ncdf}},
 \code{\link[ncdf]{close.ncdf}}, \code{\link[ncdf]{create.ncdf}},
 \code{\link[ncdf]{dim.def.ncdf}},
 \code{\link[ncdf]{get.var.ncdf}},
 \code{\link[ncdf]{put.var.ncdf}},
 \code{\link[ncdf]{open.ncdf}},
 \code{\link[ncdf]{print.ncdf}},
 \code{\link[ncdf]{set.missval.ncdf}},
 \code{\link[ncdf]{sync.ncdf}},
 \code{\link[ncdf]{var.def.ncdf}}.
 \code{\link[ncdf]{redef.ncdf}}.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\keyword{utilities}
