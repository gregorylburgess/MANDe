args <- commandArgs(trailingOnly = TRUE)
dependencies = c("multicore",
				 "mvtnorm", 
				 "raster",
				 "rgdal",
				 "rjson", 
                 "Rook", 
                 "sp",
				 "roxygen2")
if (length(args) > 0 && tolower(args[1]) == "uninstall") {
    print("Uninstalling R packages")
    uninstall.packages(dependencies)
}else {
    print("Installing R packages.")
    install.packages(dependencies, repos='http://cran.cnr.Berkeley.edu')
}
print("Done")