# EMT and SHM in the three common gardens (Pierroton, Asturias and Fundao)

library(raster)
library(tidyverse)


df <- readRDS(file="data/ClonapinData/PhenoDataNovember2019_AnnualTraits.rds") %>% 
  dplyr::select(site,longitude_site,latitude_site) %>% 
  distinct()

xysite <- SpatialPoints(df[,c("longitude_site","latitude_site")], 
                        proj4string=CRS("+proj=longlat +datum=WGS84")) %>% 
          spTransform(crs(rast)) # reproject in ETRS89 - LAEA 

for(x in c("EMT","SHM")){
  rast <- raster(paste0("data/Climate/Raw/FromMaurizio/Annualndices/",x,".tif"))
  df[,x] <- raster::extract(rast,xysite)
}

