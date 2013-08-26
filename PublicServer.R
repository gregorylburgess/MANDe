library('Rook')
source('Apps.R')

myPort <- 80
myInterface <- "0.0.0.0"
status <- -1

# R 2.15.1 uses .Internal, but the next release of R will use a .Call.
# Either way it starts the web server.
if (as.integer(R.version[["svn rev"]]) > 59600) {
	status <- .Call(tools:::startHTTPD, myInterface, myPort)
} else {
	status <- .Internal(startHTTPD(myInterface, myPort))
}

if (status == 0) {
	unlockBinding("httpdPort", environment(tools:::startDynamicHelp))
	assign("httpdPort", myPort, environment(tools:::startDynamicHelp))
	
	rook <- Rhttpd$new()
	rook$listenAddr <- myInterface
	rook$listenPort <- myPort
	
	# Change this line to your own application. You can add more than one
	# application if you like
	rook$add(name="base", Redirect$new("/../static/pages/index.html"))
	rook$add(query, name="query")
	
	## Define static content that should be exposed via url
	rook$add(name="static", 
			Builder$new( 
					Static$new(
							urls = c('/css','/img','/js', '/pages', '/txt'),
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
	while (TRUE) Sys.sleep(24 * 60 * 60)
}

# If we get here then the web server didn't start up properly
warning("Oops! Couldn't start Rook app")