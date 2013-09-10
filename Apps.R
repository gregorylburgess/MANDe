#' Defines handlers for Rook URIs.
#' 
#' 
source("src/Main.R")
library("rjson")

#' The main function that calls the webapp with parameters
#' @param env The Rook environment object.
#' @return none.
query <- function(env) {
	req = Rook::Request$new(env)
	res = Rook::Response$new()
	params = parseJSON(req$params())
	# Try and run an asynchrynous request if multicore is present
	if(require('multicore')) {
		library('multicore')
		parallel(acousticRun(params))
		res$write("processing")
	}
	# Run a serial request otheriwse.
	else {
		acousticRun(params)
		res$write("finished")
	}
	res$finish()
}

#' Parses JSON objects recieved from the client.
#' @param params A JSON string without the outter curly braces.
#' @return an R dictionary containing the key/value pairs given.
parseJSON <- function(params) {
	params = paste("{",params,"}", sep="")
	parser = newJSONParser()
	parser$addData(params)
	return(parser$getObject())
}
