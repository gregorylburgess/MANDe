sudo apt-get update
sudo apt-get install r-base-dev
sudo apt-get install git
sudo apt-get install libnetcdf-dev
git clone https://github.com/gregorylburgess/acoustic-deploy.git
cd acoustic-deploy
sudo R CMD INSTALL --configure-args="-with-netcdf-include=/usr/include/ -with-netcdf-lib=/usr/lib/" ncdf_1.6.6.tar.gz
Rscript install.R