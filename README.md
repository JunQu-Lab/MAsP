# MAsP, R Shiny-based web application

[![DOI](https://zenodo.org/badge/499283690.svg)](https://zenodo.org/badge/latestdoi/499283690)

Micro-scaffold Assisted Spatial Proteomics (MASP) is a novel method for spatial proteomic studies. MASP precisely procures spatial information by precise tissue spatial compartmentalization using a 3D-printed micro-scaffold, and accurately/sensitively quantifies thousands of proteins in a whole-tissue slice.

![](https://github.com/JunQu-Lab/UHR-IonStar/blob/master/manual_flowchart-1.png)

MAsP is a R Shiny-based interactive application designed for processing, visualization and analysis of quantitative proteomics data generated by MASP.

## MAsP Installation (less than 5 minutes)
R (Version 4.0.5 or above) is required for Windows 10 or MacOS.

Users can install MAsP with the following line in R:
```
install.packages("devtools")
library(devtools)
devtools::install_github("JunQu-Lab/MAsP")
```
## Troubleshooting
- R will occasionally ask users whether they want to update any dependent packages that have a new version. We recommend that you update them all.
- Rtools is also required to install R packages.
- During the installation of dependent packages, sometimes R would inquire, "Do you want to install from sources the packages which need compliation?" If users select "Yes" and receive an error, selecting "No" is equally acceptable.

## Run MAsP
If no error pops up, the MAsP web app could be started with the following codes:
```
library(MAsP)
MAsP::MAsPShiny()
```

## Manual
User can download the manual either at the mainpage of MAsP app or at the github directory [MAsP/inst/shiny/MAsP](https://github.com/JunQu-Lab/MAsP/tree/master/inst/shiny/MAsP).

## Test data
The files provided in this [Google Drive link](https://drive.google.com/drive/folders/16oierixQPBpj_b1WYT689bUHbGNDt0OU?usp=sharing) can be uploaded in "Data upload" section in the App, please refer to the manual for details.

The files are:
- "MASP_data.csv", the abundance ratios of the 5019 proteins quantified by MASP with spatial locations;
- "Locations.csv", the spatial locations of the micro-specimans;
- "Brain_cover.png", a brain cover image.


## Related Articles
[Whole-tissue Mapping of >5000 proteins by Micro-scaffold Assisted Spatial Proteomics(MASP), preprint.](https://www.researchsquare.com/article/rs-1786070/v1)

[Shen, Xiaomeng, et al. "IonStar enables high-precision, low-missing-data proteomics quantification in large biological cohorts." Proceedings of the National Academy of Sciences 115.21 (2018): E4767-E4776.](https://www.pnas.org/content/115/21/E4767.short)

[Wang, Xue, et al. "Ultra-High-Resolution IonStar Strategy Enhancing Accuracy and Precision of MS1-Based Proteomics and an Extensive Comparison with State-of-the-Art SWATH-MS in Large-Cohort Quantification." Analytical chemistry 93.11 (2021): 4884-4893.](https://pubs.acs.org/doi/abs/10.1021/acs.analchem.0c05002?casa_token=12l8WRigfZ0AAAAA:0qwzMnfjpE2stVCpMYKICmvqwofN15Q6ItzDZ7ATFY3m3aFI6oSzB1z20CJGzzwASyaegR5POgS8xA)

[Fonville, Judith M., et al. "Robust data processing and normalization strategy for MALDI mass spectrometric imaging." Analytical chemistry 84.3 (2012): 1310-1319.](https://pubs.acs.org/doi/full/10.1021/ac201767g)
