# ---------------------------------------------------------------------
# Table with the values of the climatic variables in the common gardens
# ---------------------------------------------------------------------

library(raster)
library(tidyverse)
library(kableExtra)

# Table with three row, one for each common garden:
df <- readRDS(file="data/ClonapinData/PhenoDataNovember2019_AnnualTraits.rds") %>% 
  dplyr::select(site,longitude_site,latitude_site) %>% 
  distinct()

# Coordinates of the common gardens in in ETRS89 - LAEA 
xysiteLAEA <- SpatialPoints(df[,c("longitude_site","latitude_site")], 
                        proj4string=CRS("+proj=longlat +datum=WGS84")) %>% 
          spTransform(crs(rast))

# Annual climatic variables:
for(x in c("EMT","MCMT","MWMT","PRCsum","SHM","SPR","TD")){
  rast <- raster(paste0("data/Climate/Raw/FromMaurizio/Annualndices/",x,".tif"))
  df[,x] <- raster::extract(rast,xysiteLAEA)
}

# Monthly climatic variables:
vars <- list.files(path="data/Climate/Raw/FromMaurizio/MonthlyIndices",pattern=".tif")
vars <- vars[-grep("PET",vars)] # remove monthly PET
vars <- str_sub(vars,1,-5) # extract variable names

for(x in vars){
  rast <- raster(paste0("data/Climate/Raw/FromMaurizio/MonthlyIndices/",x,".tif"))
  df[,x] <- raster::extract(rast,xysiteLAEA)
}

rast <- raster("/data/GenomicOffsetPinpin/Topography/TRI/TifsLAEA/TRI_LAEA.tif")
df[,"TRI"] <- raster::extract(rast,xysiteLAEA)


# Soil variables

# Rasters in WGS84, so we extract the coordinates of the common gardens in WGS84:
xysiteWGS84 <- SpatialPoints(df[,c("longitude_site","latitude_site")], 
                        proj4string=CRS("+proj=longlat +datum=WGS84"))

# Load the soil rasters in a stack
soil.rast <- list(raster("../../Pinpin_Clonapin/data/raw-data/soil/DERIVED_LAYERS/STU_EU_DEPTH_ROOTS_WGS84.tif"),
                  raster("../../Pinpin_Clonapin/data/raw-data/soil/DERIVED_LAYERS/STU_EU_T_CLAY_WGS84.tif"),
                  raster("../../Pinpin_Clonapin/data/raw-data/soil/DERIVED_LAYERS/STU_EU_T_SAND_WGS84.tif"),
                  raster("../../Pinpin_Clonapin/data/raw-data/soil/DERIVED_LAYERS/STU_EU_T_SILT_WGS84.tif"))
names(soil.rast) <- c("DepthRoots","Clay","Sand","Silt")

# Extract the soil variables
for(x in c("DepthRoots","Clay","Sand","Silt")){
  df[,x] <- raster::extract(soil.rast[[x]],xysiteWGS84)
}

options(knitr.table.format = "latex")
df %>%  
  gather(var, value, -site) %>% 
  spread(site, value) %>% 
  dplyr::rename_with(str_to_title, everything()) %>% 
  kable(digits=2,escape = FALSE)
