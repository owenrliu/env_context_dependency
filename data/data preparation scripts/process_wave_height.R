# Wave height data from Tom Bell
# updated modeled data from the California Coastal Data Information Program CDIP MOPv1.1

require(tidyverse)
require(lubridate)
require(here)
library(ncdf4)

hs <- read_csv(here('data','raw','swh.csv'),col_types = "iiiidd") %>%
  # Site 2 is the West End of SNI
  select(Year,Month,Day,Mean_2,Max_2) %>%
  unite("date",Year,Month,Day,sep="-",remove=F)%>%
  mutate(date=as.Date(date)) %>%
  rename(maxHs=Max_2)%>%

  rename(year=Year,month=Month,day=Day)%>%
  ungroup()
hs$maxHs[is.na(hs$maxHs)]<-NA

### USGS GFDL for gap filling ###

# http://www.sccoos.org/data/waves/
# USGS Coastal and Marine Geology Program
# http://cmgwindwave.usgsportals.net/
# options chosen in map data: 
# Model Grid: ENP (Eastern North Pacific)
# Model output: US Geophysical Fluid Dynamics Laboratory GFDL-ESM2M
# Format: netCDF
# Location: iooos/USGS station n46219
# Latitude: 33.221 Longitude -119.882

modwaves <- nc_open(here('data','raw','wave_significant_height_enp_gfdlhistorical.nc'))
##[Not run] view structure of NCDF file
# print(modwaves)

times<-ncvar_get(modwaves,"time")

# NCDF is pretty simple (only one longitude/latitude, so the variable wave height is one-dimensional)
heights <- ncvar_get(modwaves,"sea_surface_wave_significant_height")

# Maximum wave height defined as 
# average height of the one third highest waves in the record over the time period

# Convert to data frame
wavesdf <- data_frame(time=times,Hs=heights,date=as_date(as.POSIXct(time,origin="1970-01-01 00:00:00"))) %>%
  mutate(year=year(date),month=month(date),day=day(date)) %>%
  
  # Calculate the max Hs for each day
  group_by(year,month,day) %>%
  summarise(maxHs_usgs=max(Hs)) %>%
  ungroup()

# Remove unneeded data from workspace
rm(modwaves,times,heights)

# join the 2 datasets #
wavesdf <- hs %>%
  full_join(wavesdf,by=c("year","month","day")) %>%

  ungroup()

# coalesce data (fill in Bell data with USGS data)
# Monthly mean of daily Max Hs
wavesdf <- wavesdf %>%
  group_by(year,month) %>%
  mutate(maxHs_comb=coalesce(maxHs,maxHs_usgs)) %>%
  group_by(year,month) %>%
  summarise(maxHs=mean(maxHs_comb,na.rm=T)) %>%
  ungroup()

rm(hs)
### FINISH ###