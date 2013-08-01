args <- commandArgs(trailingOnly = TRUE)
dependencies = c("rjson", 
                  "Rook", 
                   "mvtnorm", 
                   "ncdf", 
                   "sp", 
                   "rgdal", 
                   "raster")
if (length(args) > 0 && tolower(args[1]) == "uninstall") {
    print("Uninstalling R packages")
    uninstall.packages(dependencies)
}else {
    print("Installing R packages.")
    install.packages(dependencies)
}
print("Done")