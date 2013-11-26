library('Rook')
source('Apps.R')

myPort <- 80
myInterface <- "0.0.0.0"
status <- -1
queries <<- {}
status <- .Call(tools:::startHTTPD, myInterface, myPort)
if (status == 0) {
	print("Trying to Start Server...")
	unlockBinding("httpdPort", environment(tools:::startDynamicHelp))
	assign("httpdPort", myPort, environment(tools:::startDynamicHelp))
	
	rook <- Rhttpd$new()
	rook$listenAddr <- myInterface
	rook$listenPort <- myPort
	
	# Change this line to your own application. You can add more than one
	# application if you like
	rook$add(name="base", Redirect$new("/../static/pages/index.html"))
	rook$add(query, name="query")
	rook$add(getStatus, name="status")
	## Define static content that should be exposed via url
	rook$add(name="static", 
			Builder$new( 
					Static$new(
							urls = c('/css','/img','/js', '/pages', '/txt', '/zip'),
							root = '.'
					),
					Static$new(urls='/plots',root=tempdir()), 
					Brewery$new(url='/brew',root='.'),
					App$new(function(env) {
								req <- Request$new(env)
								res <- Response$new()
								res$redirect(req$to_url('/pages/404.html'))
								res$finish()
							}
					)
			)
	)
	print("Server is Running!")
	while(TRUE) Sys.sleep(1)
	#message("Press Return To Stop Server")
	#invisible(readLines("stdin", n=1))
	#rook$stop()
	return(rook)
}

# If we get here then the web server didn't start up properly
warning(paste("Oops! Couldn't start Rook app, status = ", status, sep=""))