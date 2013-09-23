library("rjson")
source("src/Main.R")

#' Parses JSON objects recieved from the client.
#' @param params A JSON string without the outter curly braces.
#' @return an R dictionary containing the key/value pairs given.
parseJSON <- function(params) {
	params = paste("{",params,"}", sep="")
	parser = newJSONParser()
	parser$addData(params)
	return(parser$getObject())
}


testString = "\"startX\":\"308\",\"bias\":\"2\",\"startY\":\"452\",\"XDist\":\"55\",\"YDist\":\"44\",\"inputFile\":\"src/himbsyn.bathytopo.1km.v19.grd/himbsyn.bathytopo.1km.v19.grd\",\"inputFileType\":\"netcdf\",\"seriesName\":\"z\",\"cellSize\":\"1000\",\"numSensors\":\"10\",\"sensorElevation\":\"1\",\"shapeFcn\":\"shape.gauss\",\"peak\":\".98\",\"detectionRange\":\"1000\",\"north\":\"21.923\",\"south\":\"21.479\",\"west\":\"-158.869\",\"east\":\"-158.316\",\"resolution\":\"src/himbsyn.bathytopo.1km.v19.grd/himbsyn.bathytopo.1km.v19.grd\",\"suppressionRangeFactor\":\"2\",\"suppressionFcn\":\"suppression.scale\",\"maxsuppressionValue\":\"1\",\"minsuppressionValue\":\".5\",\"timestamp\":\"1379916128897\""


params = parseJSON(testString)
# Try and run an asynchrynous request if multicore is present
acousticRun(params)