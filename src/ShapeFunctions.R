## This file contains shape functions (also known as detection functions)
## Input: distance from receiver, shape parameters
## Output: detection probability
## Todo: Possibly add support for Gompertz, linear, student's t, sigmoidal, and other functions

#' @name shape.gauss
#' @title  Detection function with the shape of a half Gaussian distribution peaking
#' at distance zero, and declining with increasing distances.
#' 
#' @param dist The distance from the tag to the reciever. 
#' @param params A dictionary containing the keys 'sd' and 'peak', where sd defines a standard deviation, and peak defiens the peak value of the curve (maximum value).
#' @return The probability of detecting a tag at a given dist with the given parameters.
shape.gauss = function(dist,params){
	sd = params$sd
	peak = params$peak
    return (peak*dnorm(dist/sd,0,1)/dnorm(0,0,1))
}
