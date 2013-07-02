require('Rook')
source('Apps.R')

rook <- Rhttpd$new()
rook$add(index,'index')
rook$add(helloWorld,'helloWorld')
rook$add(summary,'summary')
rook$add(rookTestApp,'rookTestApp')
# rook$start('127.0.0.1', 8000)
rook$start(listen="0.0.0.0", port=as.numeric(Sys.getenv("PORT")))

# check if we should stop
#while(TRUE) {
#	y<-scan(n=1, what=character())
	# s stops the server
#	if('s' %in% y || 'q' %in% y) {
#		print("Exiting...")
#		rook$stop()
#		break
#	}
#}