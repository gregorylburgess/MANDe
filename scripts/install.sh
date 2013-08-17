echo "deb http://cran.rstudio.com/bin/linux/ubuntu precise/" | sudo tee -a /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo add-apt-repository ppa:marutter/rdev
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install r-base-dev
sudo apt-get install git
sudo apt-get install libgdal-dev
sudo apt-get install libproj-dev
sudo apt-get install libnetcdf-dev
sudo apt-get update
cd acoustic-deploy
sudo R CMD INSTALL --configure-args="-with-netcdf-include=/usr/include/ -with-netcdf-lib=/usr/lib/" ncdf_1.6.6.tar.gz
sudo Rscript scripts/install.R