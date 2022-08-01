# Environmental Context Dependency in Species Interactions
This repository accompanies the manuscript "Environmental Context Dependency in Species Interactions".

A full description and step-by-step tutorial is provided in `edm_analysis.html`, with all of the underlying code provided in the RMarkdown file `edm_analysis.Rmd`. To download and run a fully reproducible version of the analysis, this repository can be forked or cloned. The file `edm_analysis_code_only.R` is simply an extraction of the code from `edm_analysis.Rmd`, but viewers are encouraged to use the `html` and `Rmd` to most easily follow the steps of the analysis.

Raw data and all associated data preparation scripts are in the `data` folder. The main analysis (`edm_analysis`) begins by calling scripts that perform these data processing steps.

Please contact the author with any questions.

![](ccm_network_w_icons.png)


# Package versions and system requirements

The required packages are loaded at the beginning of each script. Once packages are installed, the entire analysis should take less than five minutes to run (tested using Windows 7 x64 and R version `3.5.1`). In the analysis in the manuscript, the following package versions were used for analysis and data visualization:

```
rEDM_1.11.0
here_1.0.1
ncdf4_1.19
knitr_1.38
kableExtra_1.3.4
ggsci_2.9
igraph_1.3.0
RANN_2.6.1
gridExtra_2.3
plot3D_1.4
quantreg_5.88
fields_13.3
tidyverse_1.3.1
ggplot2_3.3.5
purrr_0.3.4
tibble_3.1.6
tidyr_1.2.0
stringr_1.4.0
dplyr_1.0.8
readr_2.1.2
forcats_0.5.1
lubridate_1.8.1
```

## Installation of packages

If there is a need to install any of the above packages, they can be installed directly in R by calling `install.packages("PKG NAME")`, e.g. `install.packages('tidyverse')`. Individual packages are small and should not take more than a few minutes to install on a regular desktop computer. After installation of packages, the scripts mentioned above should be able to be run immediately.
Note: Installing the `tidyverse` package automatically installs many of the other packages (all those under it in the list above).