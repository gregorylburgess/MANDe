#sudo apt-get -y update
#sudo apt-get -y upgrade
sudo apt-get -y install r-base-dev
sudo apt-get -y install git
sudo apt-get -y install libgdal-dev
sudo apt-get -y install libproj-dev
sudo apt-get -y install libnetcdf-dev
sudo apt-get -y install texlive-base-bin
sudo apt-get -y install texinfo
sudo apt-get -y install texlive-latex-extra
sudo apt-get -y install texlive-fonts-extra
sudo apt-get -y install texlive-math-extra
sudo texhash
sudo apt-get -y update
sudo R CMD INSTALL --configure-args="-with-netcdf-include=/usr/include/ -with-netcdf-lib=/usr/lib/" ncdf_1.6.6.tar.gz
sudo Rscript scripts/install.R
mkdir img
mkdir txt
