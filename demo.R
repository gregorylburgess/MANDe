require('Rook')

rook <- Rhttpd$new()
helloWorld <- function(env){
	req <- Rook::Request$new(env)
	res <- Rook::Response$new()
	res$write("HELLO WORLD")
	res$finish()
}

rook$add(helloWorld,'helloWorld')
rook$start(listen="0.0.0.0", port=as.numeric(Sys.getenv("PORT")))
while(T) {
	Sys.sleep(10000)
}