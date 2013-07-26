require('Rook')
source('Apps.R')

rook <- Rhttpd$new()
Heroku = FALSE

rook$add(name="base", Redirect$new("/../static/pages/index.html"))
rook$add(query, name="query")
rook$add(testquery, name="testquery")
rook$add(index,name='index')
rook$add(helloWorld,name='helloWorld')
rook$add(summary,name='summary')
rook$add(rookTestApp,name='rookTestApp')
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

}