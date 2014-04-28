#' Defines handlers and helper methods for Rook URIs.

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
	
	#if this is the first package of data, create a new index for it
	if(params$part==1) {
		acousticQueries[toString(params$timestamp)] <<- data
	}
	#otherwise append the data to an existing index
	else {
		acousticQueries[toString(params$timestamp)] <<- paste(acousticQueries[toString(params$timestamp)], data, sep="") 
	}

	#Once we have all the packets, run the request.
	if(params$part == params$complete) {
		parameters = parseJSON(gsub("\'", "\"", acousticQueries[toString(params$timestamp)]))
		# Try and run an asynchrynous request if multicore is present
		if(require('multicore')) {
			library('multicore')
			acousticJobs[toString(params$timestamp)] <<- parallel(execute(res, parameters))
		}
		# Run a serial request otheriwse.
		else {
			execute(res, parameters)
		}
		
	}
	
	#Finish the request
	res$finish()
}

#' Calls acousticRun and writes a json file.  Required for a parallel call.
#' @param res The Rook response object.
#' @param parameters A dictionary of parameters to pass to acousticRun().
execute = function (res, parameters) {
	res$write("processing...")
	result = acousticRun(parameters)
	writeJSON(result$filenames$jsonFile, result)
	res$write("finished!")
}

#' Writes a JSON file with the provided data.
#' @param data The data to write.
#' @param filename The output filepath.
writeJSON = function (filename, data) {
	# Write results to a json file
	jsonFile = data$filenames$jsonFile
	file.create(jsonFile)
	cat(toJSON(data), file=jsonFile, append=FALSE)
}


#' Returns the status of a particular request
#' @param env The Rook environment object.
#' @return none.
getStatus <- function(env) {
	req = Rook::Request$new(env)
	res = Rook::Response$new()
	params = parseJSON(req$params())
	for(job in jobs) {
		print(collect(job, wait=FALSE, timeout=1))
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
