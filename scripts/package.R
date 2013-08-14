# Delete any old packages
unlink("acoustic", TRUE)

# A bug in package.skeleton() makes me import these...
library(methods)
library(utils)
# Frame the package
package.skeleton("acoustic", code_files=c("src/Bathy.R", "src/FishModel.R", "src/Utility.R", "src/ShapeFunctions.R", "src/Main.R"))

source("scripts/change.R")

library(roxygen2)
roxygenize("acoustic", copy=FALSE, roclets = c("collate", "namespace", "rd"))
