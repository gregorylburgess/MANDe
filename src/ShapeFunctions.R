## This file contains shape functions (also known as detection functions)
## Input: distance from receiver, shape parameters
## Output: detection probability
## Todo: Possibly add Gompertz, linear, and other functions

shape.t <- function(dist,par){
    sd <- par$sd
    peak <- par$peak
    ## Shape based on t-distribution
    return(peak*dt(dist/sd,df=1)/dt(0,df=1))
}

shape.gauss <- function(dist,par){
	sd <- par$sd
	peak <- par$peak
    ## Shape based on normal distribution
    return (peak*dnorm(dist/sd,0,1)/dnorm(0,0,1))
}

shape.sigmoidal <- function(dist,par){
    pmax <- par[1]
    D50 <- par[2]
    D95 <- par[3]
    ## Detection probability function from How and de Lestang 2012, eqn 3, note error in eqn 3, it should +exp not -exp
    return(pmax/(1 + exp(log(19)*(dist-D50)/(D95-D50))))
}