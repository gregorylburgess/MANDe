source("src/Main.R")
library("rjson")

## The main function that calls the webapp with parameters
query <- function(env) {
	req = Rook::Request$new(env)
	res = Rook::Response$new()
	params = parseJSON(req$params())
	results=run(params)
	res$write(toJSON(results, method="C"))
	res$finish()
}

## Parses JSON objects recieved from the client.
parseJSON <- function(params) {
	params = as.character(paste("{",params,"}", sep=""))
	parser = newJSONParser()
	parser$addData(params)
	return(parser$getObject())
}
