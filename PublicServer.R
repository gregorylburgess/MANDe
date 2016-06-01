library('Rook')
source('Apps.R')

myPort <- 80
myInterface <- "0.0.0.0"
#Rook Status
status <- -1

# Stores partial parameter sections
acousticQueries <<- {}
# Status dictionary for parallel jobs
acousticJobs <<-{}

status <- .Call(tools:::startHTTPD, myInterface, myPort)
if (status == 0) {
	print("Trying to Start Server...")
	unlockBinding("httpdPort", environment(tools:::startDynamicHelp))
	assign("httpdPort", myPort, environment(tools:::startDynamicHelp))
	
	rook <- Rhttpd$new()
	rook$listenAddr <- myInterface
	rook$listenPort <- myPort
	# Home page
	rook$add(name="base", app=Redirect$new("/../static/pages/index.html"))
	# Submit a job
	rook$add(app=query, name="query")
	# Check on the status of a job
	rook$add(app=getStatus, name="status")
	## Expose static directories via URL.
	rook$add(name="static", 
		app=Builder$new( 
			Static$new(
				urls = c('/css','/img','/js', '/pages', '/txt', '/zip'),
				root = '.'
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
