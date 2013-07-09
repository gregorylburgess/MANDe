require('Rook')
source('Apps.R')
rook <- Rhttpd$new()

Heroku = TRUE

rook$add(name="base", Redirect$new("/../static/pages/index.html"))
rook$add(index,'index')
rook$add(helloWorld,'helloWorld')
rook$add(summary,'summary')
rook$add(rookTestApp,'rookTestApp')
rook$add(name="static", 
		Builder$new( 
				Static$new(
							urls = c('/css','/img','/js', '/pages'),
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


if (Heroku) {
	rook$start(listen="0.0.0.0", port=as.numeric(Sys.getenv("PORT")))
	while(T) {
		Sys.sleep(10000)
	}
} else {
	rook$start(listen="127.0.0.1", port=8000)

	while(TRUE) {
	  y<-scan(n=1, what=character())
	  # s stops the server
	  if('s' %in% y ) {
		print("Exiting...")
	    rook$stop()
		break
	  }
	  if('r' %in% y) {
		rook$stop()
		rook$start(listen="127.0.0.1", port=8000)
	  }
	}
}