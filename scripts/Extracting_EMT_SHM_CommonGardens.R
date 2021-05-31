# EMT, SHM and MCMT in the three common gardens (Pierroton, Asturias and Fundao)

# EMT = extreme minium temperature over the study period (1901-1950)
# SHM = summer heat moiture index
# MCMT = mean coldest month temperature

library(raster)
library(tidyverse)


df <- readRDS(file="data/ClonapinData/PhenoDataNovember2019_AnnualTraits.rds") %>% 
  dplyr::select(site,longitude_site,latitude_site) %>% 
  distinct()

xysite <- SpatialPoints(df[,c("longitude_site","latitude_site")], 
                        proj4string=CRS("+proj=longlat +datum=WGS84")) %>% 
          spTransform(crs(rast)) # reproject in ETRS89 - LAEA 

for(x in c("EMT","SHM","MCMT")){
  rast <- raster(paste0("data/Climate/Raw/FromMaurizio/Annualndices/",x,".tif"))
  df[,x] <- raster::extract(rast,xysite)
}

