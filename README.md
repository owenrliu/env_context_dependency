# Environmental Context Dependency in Species Interactions
This repository accompanies the manuscript "Environmental Context Dependency in Species Interactions".

A full description and step-by-step tutorial is provided in `edm_analysis.html`, with all of the underlying code provided in the RMarkdown file `edm_analysis.Rmd`. To download and run a fully reproducible version of the analysis, this repository can be forked or cloned. The file `edm_analysis_code_only.R` is simply an extraction of the code from `edm_analysis.Rmd`, but viewers are encouraged to use the `html` and `Rmd` to most easily follow the steps of the analysis.

Raw data and all associated data preparation scripts are in the `data` folder. The main analysis (`edm_analysis`) begins by calling scripts that perform these data processing steps.

Please contact the author with any questions.

![](ccm_network_w_icons.png)


# Package versions and system requirements

The required packages are loaded at the beginning of each script. Once packages are installed, the entire analysis should take less than five minutes to run (tested using Windows 7 x64 and R version `3.5.1`). In the analysis in the manuscript, the following package versions were used for analysis and data visualization:

```
rEDM_0.7.4
here_0.1
ncdf4_1.16
knitr_1.21
kableExtra_0.9.0
ggsci_2.9
igraph_1.2.2
RANN_2.6.1
gridExtra_2.3
plot3D_1.1.1
quantreg_5.38
fields_9.6
tidyverse_1.2.1
  ggplot2_3.1.0
  purrr_0.2.5
  tibble_2.0.0
  tidyr_0.8.2
  stringr_1.3.1
  dplyr_0.7.8
  readr_1.3.1
  forcats_0.3.0
  lubridate_1.7.4
```

## Installation of packages

If there is a need to install any of the above packages, they can be installed directly in R by calling `install.packages("PKG NAME")`, e.g. `install.packages('tidyverse')`. Individual packages are small and should not take more than a few minutes to install on a regular desktop computer. After installation of packages, the scripts mentioned above should be able to be run immediately.
Note: Installing the `tidyverse` package automatically installs many of the other packages (all those under it in the list above).