library('Rook')
source('Apps.R')
rook=Rhttpd$new()
rook$add(query,"query")
rook$add(getStatus,"status")
status<<-{}
queries<<-{}
jobs<<-{}
