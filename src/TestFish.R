source('src/FishModel.R')
source('src/Utility.R')
options(digits=14)
tolerance = 1e-10
TestFish.fish = function() {

	sol = list(
			s1 = matrix(.04,
						nrow=5,
						ncol=5),
		
			s2 = matrix(c(0.023246839878294, 0.033823952439922, 0.038327559383904, 0.033823952439922, 
						  0.023246839878294, 0.033823952439922, 0.049213560408541, 0.055766269846849, 
						  0.049213560408541, 0.033823952439922, 0.038327559383904, 0.055766269846849, 
						  0.063191462410265, 0.055766269846849, 0.038327559383904, 0.033823952439922, 
						  0.049213560408541, 0.055766269846849, 0.049213560408541, 0.033823952439922, 
						  0.023246839878294, 0.033823952439922, 0.038327559383904, 0.033823952439922, 
						  0.023246839878294),
						nrow=5,
						ncol=5),
		  
		  	s3 = matrix(c(0.000000000000000, 0.000000000000000, 0.000000000000000, 0.077542768547716, 
						  0.053294313470615, 0.077542768547716, 0.112824062502653, 0.127846411893597, 
						  0.112824062502653, 0.077542768547716, 0.087867468226364, 0.127846411893597, 
						  0.144868963867373, 0.000000000000000, 0.000000000000000, 0, 0, 0, 0, 0, 0, 
						  0, 0, 0, 0),
				  		nrow=5,
				  		ncol=5),
		  
		   s4 = matrix(c(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 
						 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
				  		nrow=5,
				  		ncol=5)
		 )
	bGrid =  {}
	fGrid = {}
	bGrid$bGrid = matrix(-1:-25,nrow=5,ncol=5)
	
	params = list(fishmodel="rw", numSensors=2, bias=1)
	params$cellSize=50
	bGrid$x = (1:dim(bGrid$bGrid)[1])*params$cellSize 
	bGrid$y = (1:dim(bGrid$bGrid)[2])*params$cellSize

	#rw
	fGrid$f1 = fish(params,bGrid)
	
	#ou
	params$fishmodel="ou"
	params$mux=.5
	params$muy=.5
	params$oucor=0
	params$ousdx=100
	params$ousdy=100
	params = checkParams(params)
	fGrid$f2 = fish(params,bGrid)
	
	
	#ou with vertical habitat restriction
	params$mindepth = -3
	params$maxdepth = -14
	fGrid$f3 = fish(params,bGrid)
	
	#rw with vertical habitat restriction
	params$fishmodel ="rw"
	fGrid$f4 = fish(params,bGrid)
	
	for(i in 1:length(fGrid)) {
		#check that we got the right data
		if(!isTRUE(all.equal(fGrid[[i]], sol[[i]], tolerance=tolerance))){
			print("Expected:")
			print(sol[[i]])
			print("Recieved:")
			print(fGrid[[i]])
			print(split(fGrid[[i]], rep(1:ncol(fGrid[[i]]), each = nrow(fGrid[[i]]))))
			print("Diff:")
			print(fGrid[[i]]-sol[[i]])
			stop(sprintf("[TestFish.fish %i]: FAIL",i))
		}
	}
	print("[fish]: Pass")
}


#' Runs a battery of tests.
#' 
#' @return N/A
#' @export
runFishTests = function () {
	print("--- Testing Fish Functions ---")
	TestFish.fish()
	print("Success! All Fish Tests Passed.")
}


#Executes the tests, and stops if an error occurs
tryCatch({
			runFishTests()
		}, warning = function(e) {
			print(e[1])
		}, error = function(e) {
			print(e[1])
		}, finally = function(e){})
