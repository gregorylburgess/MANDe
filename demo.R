require('Rook')
rook <- Rhttpd$new()
rook$start(listen="0.0.0.0", port=as.numeric(Sys.getenv("PORT")))

paste("Hello World")

