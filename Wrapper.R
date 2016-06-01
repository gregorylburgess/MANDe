source("src/Main.R")
library("rjson")

#' Parses JSON objects recieved from the client.
#' @param params A JSON string without the outter curly braces.
#' @return an R dictionary containing the key/value pairs given.
parseJSON <- function(params) {
	print(params)
	params = as.character(paste("{",params,"}", sep=""))
	print(params)
	parser = newJSONParser()
	parser$addData(params)
	return(parser$getObject())
}

args = commandArgs(trailingOnly = TRUE)
#params = parseJSON(args[1])
print(params)
#acousticRun(params)

