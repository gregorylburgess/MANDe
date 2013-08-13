change <- function(input, output, pattern, value) {
	conn <- file(input, "r")
	out =  readLines(conn, 1)
	while(length(line <- readLines(conn, 1)) > 0) {
		out = paste(out, line, sep="\n")
	}
	close(conn)
	conn <- file(output, "w+")
	out = gsub(pattern, value, out)
	cat(out,file=output)
	close(conn)
	return()
}

change("acoustic/R/Bathy.R", "acoustic/R/Bathy.R", "src/", "")
change("acoustic/R/FishModel.R", "acoustic/R/FishModel.R", "src/", "")
change("acoustic/R/Utility.R", "acoustic/R/Utility.R", "src/", "")
change("acoustic/R/ShapeFunctions.R", "acoustic/R/ShapeFunctions.R", "src/", "")
change("acoustic/R/Main.R", "acoustic/R/Main.R", "src/", "")