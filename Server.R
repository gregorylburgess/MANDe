require('Rook')
source('Apps.R')
rook <- Rhttpd$new()

helloWorld <- function(env){
	req <- Rook::Request$new(env)
	res <- Rook::Response$new()
	res$write("HELLO WORLD")
	res$finish()
}
rook$add(index,'index')
rook$add(helloWorld,'helloWorld')
rook$add(summary,'summary')
rook$add(rookTestApp,'rookTestApp')
#rook$start(listen="0.0.0.0", port=as.numeric(Sys.getenv("PORT")))
rook$start(listen="127.0.0.1", port=8000)
Redirect$new("custom/index")

while(T) {
	Sys.sleep(10000)
}