

query <- function(env) {
	req <- Rook::Request$new(env)
	res <- Rook::Response$new()
	res$write(req$params())
	res$finish()
}

index <- function(env){
	req <- Rook::Request$new(env)
	res <- Rook::Response$new()
	res$write("<h1>Index</h1></br>")
	res$write("<a href='/custom/base/'>Bootstrap</a></br>")
	res$write("<a href='/custom/helloWorld/'>HelloWorld</a></br>")
	res$write("<a href='/custom/summary/'>Forms</a></br>")
	res$write("<a href='/custom/rookTestApp/'>Post Data</a>")
	res$finish()
}

helloWorld <- function(env){
	req <- Rook::Request$new(env)
	res <- Rook::Response$new()
	res$write("HELLO WORLD")
	res$finish()
}

summary <- function(env) {
	req <- Rook::Request$new(env)
	res <- Rook::Response$new()
	res$write('Choose a CSV file:\n')
	res$write('<form method="POST" enctype="multipart/form-data">\n')
	res$write('<input type="file" name="data">\n')
	res$write('<input type="submit" name="Upload">\n</form>\n<br>')
	
	if (!is.null(req$POST())){
		data <- req$POST()[['data']]
		res$write("<h3>Summary of Data</h3>");
		res$write("<pre>")
		res$write(paste(capture.output(summary(read.csv(data$tempfile,stringsAsFactors=FALSE)),file=NULL),collapse='\n'))
		res$write("</pre>")
		res$write("<h3>First few lines (head())</h3>");
		res$write("<pre>")
		res$write(paste(capture.output(head(read.csv(data$tempfile,stringsAsFactors=FALSE)),file=NULL),collapse='\n'))
		res$write("</pre>")
	}
	res$finish()
}

rookTestApp <- function(env){
	req <- Request$new(env)
	res <- Response$new()
	envstr <- paste(capture.output(str(as.list(env)),file=NULL),collapse='\n')
	poststr <- paste(capture.output(str(req$POST()),file=NULL),collapse='\n')
	getstr <- paste(capture.output(str(req$GET()),file=NULL),collapse='\n')
	randomString <- function() paste(letters[floor(runif(10,0,26))],collapse='')
	randomNumber <- function() runif(1,0,26)
	res$write( c(
					'<HTML><head><style type="text/css">\n',
					'table { border: 1px solid #8897be; border-spacing: 0px; font-size: 10pt; }',
					'td { border-bottom:1px solid #d9d9d9; border-left:1px solid #d9d9d9; border-spacing: 0px; padding: 3px 8px; }',
					'td.l { font-weight: bold; width: 10%; }\n',
					'tr.e { background-color: #eeeeee; border-spacing: 0px; }\n',
					'tr.o { background-color: #ffffff; border-spacing: 0px; }\n',
					'</style></head><BODY>',
					'<img alt="rook logo" src="http://wiki.rapache.net/static/rook.png">',
					'<H1>Welcome to Rook</H1>\n',
					sprintf('<form enctype="multipart/form-data" method=POST action="%s/%s?called=%s">',env$SCRIPT_NAME,randomString(),randomNumber()),
					'Enter a string: <input type=text name=name value=""><br>\n',
					'Enter another string: <input type=text name=name2 value=""><br>\n',
					'Upload a file: <input type=file name=fileUpload><br>\n',
					'Upload another file: <input type=file name=anotherFile><br>\n',
					'<input type=submit name=Submit><br><br>',
					'Environment:<br><pre>',envstr,'</pre><br>',
					'Get:<br><pre>',getstr,'</pre><br>',
					'Post:<br><pre>',poststr, '</pre><br><br>'
			))
	res$finish()
}
