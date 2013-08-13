#' Defines handlers for Rook URIs.
#' 
#' 
source("Acoustic/R/Main.R")
library("rjson")

#' The main function that calls the webapp with parameters
#' @param env The Rook environment object.
#' @return none.
query <- function(env) {
	req = Rook::Request$new(env)
	res = Rook::Response$new()
	params = parseJSON(req$params())
	results=run(params)
	res$write(toJSON(results, method="C"))
	res$finish()
}

#' Parses JSON objects recieved from the client.
#' @param params A JSON string without the outter curly braces.
#' @return an R dictionary containing the key/value pairs given.
parseJSON <- function(params) {
	params = as.character(paste("{",params,"}", sep=""))
	parser = newJSONParser()
	parser$addData(params)
	return(parser$getObject())
}
