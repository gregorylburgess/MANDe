mkdir pages
mkdir pages/help_pages
sudo Rscript scripts/package.R
sudo R CMD check acoustic
sudo R CMD build acoustic
