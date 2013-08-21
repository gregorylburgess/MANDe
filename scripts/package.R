#' @name change
#' @title Reads in the input file, and replaces all matching pattern with value.
#' @param input Path to the file to read.
#' @param output Path to the file to output.
#' @param pattern The pattern to search for and replace.
#' @param value The value to replace matched patterns with.
#' @return N/A
change <- function(input, output, pattern, value, debug=TRUE) {
	conn <- file(input, "r")
	options(warn=-1)
	out =  readLines(conn, 1)
	while(length(line <- readLines(conn, 1)) > 0) {
		out = paste(out, line, sep="\n")
	}
	close(conn)
	conn <- file(output, "w+")
	out = gsub(pattern, value, out)
	cat(out,file=output)
	close(conn)
	if(debug) {
		cat(sprintf("Changed all references to \'%s\' to \'%s\' in \'%s\'\n", pattern, value, output))
	}
}

library(roxygen2)
# name of the package to create (Fill in whatever you want)
packageName = "acoustic"
# Names of source files.  These must exist in 'src/' folder of the working directory.
files = c("Description.R","Bathy.R", 
		"FishModel.R", 
		"Utility.R", 
		"ShapeFunctions.R", 
		"Main.R")

# Delete any old packages of the same name
unlink(packageName, TRUE)

# A bug in package.skeleton() requires us to import these...
library(methods)
library(utils)

code_files = {}
for (file in files) {
	code_files = c(code_files, paste("src/", file, sep=""))
}

# Sets up the default package directories and files
package.skeleton(packageName, 
		code_files=code_files)

# Changes all file paths in source() calls to their new values.
for (file in files) {
	name = paste(packageName, "/R/", file, sep="")
	change(name, name, "src/", "")
}

# Call Roxygen to make the .Rd files
print('--- Roxygenize')
roxygenize(packageName, copy=FALSE)

# Changes all file paths in source() calls to their new values.
for (file in files) {
	name = paste(packageName, "/R/", file, sep="")
	change(name, name, paste("source\\(\'", sep=""),paste("source\\(","\'R/", sep=""))
}

# Generate HTML files from Rd Files
rdFiles = list.files(paste(packageName, "/man", sep=""))
rdFiles = gsub(".Rd", "", rdFiles)

# write a list of files to file for the help.html file to read.
cat(rdFiles, file="pages/help_pages/index.txt", sep="\n")

outPath = paste(packageName, "/man/", sep="")
for (file in rdFiles) {
	command = paste("R CMD Rdconv -t html -o ", paste(outPath, file, ".html", sep=""),
					" ", packageName, "/man/", sep="")
	print(paste(command, file, ".Rd", sep=""))
	system(command=paste(command, file, ".Rd", sep=""))
	# Copy files to the pages/help_pages folder
	file.copy(paste(outPath, file, ".html", sep=""), paste("pages/help_pages/", file, ".html", sep=""),
			overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
}
# file.remove("pages/help_pages/acoustic-package.html")
# Delete the examples{} section of the rd file
# gsub won't let me match an open curly brace...
path = paste(packageName, "/man/", packageName, "-package.Rd", sep="")
conn <- file(path, "r")
options(warn=-1)
out =  readLines(conn, 1)
lines = 0
while(length(line <- readLines(conn, 1)) > 0 && lines < 38) {
	lines = lines + 1
	out = paste(out, line, sep="\n")
}
close(conn)
conn <- file(path, "w+")
cat(out,file=path)
close(conn)

