require('Rook')

library(Rook)
library(rjson)

rook = Rhttpd$new()
rook$add(
		name ="summarize",
		app = function(env) {
			req = Rook::Request$new(env)
			numbers = as.numeric(unlist(strsplit(req$params()$numbers, ",")))
			results = list()
			results$mean = mean(numbers)
			results$sd = sd(numbers)
			
			res = Rook::Response$new()
			res$write(toJSON(results))
			res$finish()
		}
)

rook$browse("summarize")