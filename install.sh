sudo apt-get update
sudo apt-get install r-base-dev
sudo apt-get install git
sudo apt-get install libnetcdf-dev
git clone https://github.com/gregorylburgess/acoustic-deploy.git
cd acoustic-deploy
cd ncdf
export NETCDF_INCLUDE=/usr/include/
export NETCDF_LIB=/user/lib/
Rscript install.R