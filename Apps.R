#' Defines handlers for Rook URIs.
#' 
#' 
source("src/Main.R")
#' @include src/Main.R
library("rjson")

#' The main function that calls the webapp with parameters
#' @param env The Rook environment object.
#' @return none.
query <- function(env) {
	req = Rook::Request$new(env)
	res = Rook::Response$new()
	
	params = parseJSON(req$params())
	data = params$data
	
	if(params$part==1) {
		queries[toString(params$timestamp)] <<- data
	}
	else {
		queries[toString(params$timestamp)] <<- paste(queries[toString(params$timestamp)], data, sep="") 
	}

	if(params$part == params$complete){
		parameters = parseJSON(gsub("\'", "\"", queries[toString(params$timestamp)]))
		# Try and run an asynchrynous request if multicore is present
		if(require('multicore')) {
			library('multicore')
			jobs[toString(params$timestamp)] <<- parallel(acousticRun(parameters))
			res$write("processing")
		}
		# Run a serial request otheriwse.
		else {
			acousticRun(parameters)
			res$write("finished")
		}
	}
	res$finish()
}

#' Returns the status of a particular request
#' @param env The Rook environment object.
#' @return none.
getStatus <- function(env) {
	req = Rook::Request$new(env)
	res = Rook::Response$new()
	params = parseJSON(req$params())
	for(job in jobs) {
		print(collect(job, wait, 1))
	}
	jobStatus = checkStatus(toString(params$timestamp))
	res$write(jobStatus)
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
