source('Bathy.R')
ncdf = list("ncdf")
arcgis = list("rgdal", "raster", "sp")
prereqs = list(ncdf=ncdf, arcgis=arcgis)

## Check that required packages installed
for (prereq in names(prereqs)) {
	for (requirement in prereqs[prereq]) {
		if(!(require(requirement, character.only=TRUE))) {
			print(sprintf("Package %s not installed!  Unable to parse %s files!", requirement, prereq))
		}
		else {
			print(sprintf("Package %s was successfully loaded.",requirement))
		}
	}
}