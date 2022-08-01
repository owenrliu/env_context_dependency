### Join benthic swath and physical variables ###

## The data preparation scripts are nested:
## Each data set is processed separately, then the biological and physical data are joined
## and finally, this script joins the biological and physical data.

## This script produces the final dataset ready for the main EDM analysis

## Normalized time series from benthic swath data and physical variables
require(dplyr)
require(here)

## Run the scripts that processed normalized time series
source(here::here('data','data preparation scripts','physical_data_preparation.R'))
source(here::here('data','data preparation scripts','benthic_data_preparation.R'))


## Now we simply join the two datasets by their common monitoring period identifier
## Physical variables will be duplicated to attach to each of the benthic swath replicates

westend <- westend.norm %>% left_join(phys.dat.norm,by="period")

# Linearly interpolate missing values
fill_vec1<- function(arrvec){
  x <- 1:63
  z <- approx(x,arrvec,xout=x,method="linear")$y
  z
}

# interpolate 

westend <- westend %>% 
  group_by(site) %>% 
  arrange(period) %>% 
  mutate(across(lam:sst,fill_vec1)) %>% 
  ungroup() %>% arrange(site,period)

# Remove excess data from the workspace
rm(list=setdiff(ls(), 'westend'))

### FINISH ###