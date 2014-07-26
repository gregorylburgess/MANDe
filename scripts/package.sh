mkdir pages
mkdir pages/help_pages
sudo Rscript scripts/package.R
R CMD Rd2pdf acoustic
sudo mv acoustic.pdf acoustic/
sudo R CMD check acoustic
sudo R CMD build acoustic
