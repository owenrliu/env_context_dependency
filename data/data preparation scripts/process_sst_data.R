### Extract appropriate SNI SST data

# Extracts the data from the .nc file, and casts to a long form dataframe

require(ncdf4)
require(fields)
require(lubridate)
require(dplyr)
require(here)
here <- here::here

extr_ncdf_SNI <- function(fname) {
  nc = nc_open(fname)
  ncdates = nc$dim$time$vals
  ncdates = as.Date(ncdates,origin = '1800-1-1')
  nlon = length(nc$dim$lon$vals)
  lons= nc$dim$lon$vals
  nlat = length(nc$dim$lat$vals)
  lats=nc$dim$lat$vals
  ndates = length(ncdates)
  
  sstout = array(data = NA, dim = c(nlon,nlat,ndates))
  # Extract the data from the NetCDF file
  sstout[,,] = ncvar_get(nc, varid = 'sst',
                         count = c(nlon,nlat,ndates))
  # The output array sstout will be arranged with longitudes in rows, 
  # increasing in an easterly direction as you move down a row (larger 
  # longitude values), and latitudes in columns, increasing in latitude (more 
  # northerly) as you move across columns. The 3rd dimension represents 
  # different dates.
  
  # If there are missing data in the NetCDF, they should appear as 32767. 
  # Replace that value with NA if it occurs anywhere.
  sstout = ifelse(sstout == 32767, NA, sstout)
  
  # Rearrange the output matrix or array so that latitudes run from north to
  # south down the rows, and longitudes run from west to east across columns.
  dims = dim(sstout) # get size of array
  # Make new output array to hold rearranged data. The dimension names will
  # match the newly rearranged latitude and longitude values
  sst2 = array(data = NA, dim = c(dims[2],dims[1],dims[3]),
               dimnames = list(Lat = rev(lats),
                               Long = lons,
                               Date = as.character(ncdates)))
  # Step through each page of array and rearrange lat/lon values
  for (i in 1:dims[3]){
    # Extract one day's worth of lat/lon pairs
    temp = as.matrix(sstout[,,i])
    temp = t(temp) # transpose lon/lat to lat/lon
    temp = temp[nrow(temp):1,] # reverse row order to reverse latitudes
    sst2[,,i] = temp # write data to sst2 array
  }	
  
  require(reshape2)
  sst3 <- melt(sst2,value.name="sst")
  
  return(sst3)
}

## Pull out all data for this 1x1deg grid around SNI, 1981-2011
sst.dat.all <- data.frame(Lat=numeric(),Long=numeric(),Date=character(),sst=numeric())
for(i in 1981:2014) {
  fname <- paste0(here('data','raw','ncdf sst'),'/sst.day.mean.',i,'.nc')
  temp <- extr_ncdf_SNI(fname)
  sst.dat.all <- rbind(sst.dat.all,temp)
}


sst.dat.all <- sst.dat.all %>% mutate(month=month(Date),year=year(Date),day=day(Date))

## Final output data is monthly mean SST based on the above daily SSTs

noaa.sst <- sst.dat.all %>% group_by(year, month, Lat, Long) %>%
  summarise(mean.monthly.sst = mean(sst,na.rm=T)) %>% 
  ungroup() %>% 
  # pick the closest pixel of SST data to San Nicolas Island
  filter(Lat==33.125,Long==240.625)

## Remove all unneeded data from workspace ##
rm(fname,i,extr_ncdf_SNI,temp,sst.dat.all)