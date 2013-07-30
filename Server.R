require('Rook')
source('Apps.R')

rook <- Rhttpd$new()
## Indicates whether this is a Heroku deployment (=TRUE) or a local server (=FALSE).
Heroku = FALSE

## Add the handler functions
rook$add(name="base", Redirect$new("/../static/pages/index.html"))
rook$add(query, name="query")

## Define static content that should be exposed via url
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
## Sets the default Heroku port and url bindings, and then starts an infinite loop.
if (Heroku) {
	rook$start(listen="0.0.0.0", port=as.numeric(Sys.getenv("PORT")))
	while(T) {
		Sys.sleep(10000)
	}
## Starts a local server.
} else {
	rook$start(listen="127.0.0.1", port=8000)

}