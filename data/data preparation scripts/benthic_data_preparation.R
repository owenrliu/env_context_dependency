# Benthic swath data preparation
# Data of benthic densities of key species of interest for EDM
# Michael C. Kenner, James A. Estes, M. Tim Tinker, James L. Bodkin, Robert K. Cowen, 
# Christopher Harrold, Brian B. Hatfield, Mark Novak, Andrew Rassweiler, and Daniel C. Reed. 2013. 
# A multi-decade time series of kelp forest community structure at San Nicolas Island, California (USA). 
# Ecology 94:2654. http://dx.doi.org/10.1890/13-0561.1

# This function takes the raw benthic density data in their raw form, edits the species names, 
# and fills in missing values to create complete timeseries.
library(tidyverse)
library(here)
here <- here::here

#### Species Key ####
# Kenner et al. published a species key that we use to match species code to species name
# All species, all surveys filtered to just benthic density survey
spp.key <- read_csv(here('data','raw','species_name_key.csv'),col_types = cols()) %>% 
  filter(DataSet=="Benthic density")

# better shorthand names
spp.key$name <- c("pat","red","purp","wavy","astro","derm","halc","halr","limp","paras","pis","pyc","cys","lam","ymac","youn","eis","pter","mac")


#species of interest
study_spp <- c("mac","ymac","lam","pter","red","purp")

#study sites of interest (permanent transects at the "West End" of San Nicolas Island)
westend.sites <- c("2_10L", '2_22L', "2_32R", "2_39L", "2_45L", "3_10R" ,"3_22R", "3_32L", "3_39R", "3_45L")

##### Load, filter, and normalize the data ####
westend.norm <- read_csv('data/raw/Benthic density raw data.csv',col_types="iiccid") %>%
  
  #fix date column
  mutate(Date=as.POSIXct(Date,format="%m/%d/%y")) %>%

# Replace species code with species name and rename some variables
  rename(station=Station,period=Period,date=Date,swath=Swath,dens=Density) %>%
  left_join(select(spp.key,SpeciesCode,name),by="SpeciesCode") %>%
  
  # remove unneeded columns
  select(-SpeciesCode,-date) %>%
  rename(spp=name) %>%
  arrange(spp,period,swath) %>%
  
  # join station and swath columns into one unique 'site' identifier
  unite(site,station,swath) %>% 
  
  # filter out just species of interest, then normalize the time series from each site
  # to zero mean, unit variance
  filter(spp %in% study_spp, site %in% westend.sites) %>%
  group_by(spp,site) %>%
  # remove linear trends by using residuals of simple linear model
  mutate(dens_res= residuals(lm(dens~period))) %>% 
  
  # normalize to zero mean, unit variance
  mutate(norm=(dens_res-mean(dens_res,na.rm=T))/sd(dens_res,na.rm=T)) %>%

  ungroup() %>%
  distinct() %>%
  
  # Cast from long to wide-form data (periods by species) for use in later analyses
  select(-dens_res,-dens) %>%
  spread(key=spp,value=norm) %>%
  
  # Finally, fill in empty monitoring periods with NA values
  group_by(site) %>% complete(period=full_seq(x=period,1)) %>% ungroup()

