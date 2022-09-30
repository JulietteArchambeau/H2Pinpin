# Juliette Archambeau
# September 2022

# In this script, we generate the input files for the FEEM analysis
# See the tutorial https://github.com/NovembreLab/feems/blob/main/docsrc/notebooks/getting-started.ipynb

library(genio)
library(compare)
library(tidyverse)
library(sf)




# Genomic data #### 
###################

# The genomic data as to be in PLINK format

geno_names <- read.delim2("../../GenomicOffset/GenomicOffsetPinPin/data/ClonapinBlups523IndPiMASSJuly2019.txt", 
                          row.names=1) # file with the genotype names (clone names)
geno <- read.csv("../../GenomicOffset/GenomicOffsetPinPin/data/5165snps523genotypesNA.txt", 
                 header=FALSE, row.names=1) # # genomic matrix (5,165 SNPs x 523 genotypes)
geno <- geno[,3:dim(geno)[[2]]] # remove the first two columns with allele info (A,T, G or C)
colnames(geno) <- rownames(geno_names) # attrribute the genotype name for each column of the genomic matrix

head(geno[,1:10])


# Imputing missing data based on the most common allele of the gene pool
# ======================================================================

prop <- readRDS(file="data/ClonapinData/PhenoDataNovember2019_AnnualTraits.rds") %>% 
  dplyr::select(clon,prov,paste0(rep("Q",6),1:6),max.Q) %>%  
  unique() %>% 
  drop_na()

compare(prop$clon[order(prop$clon)],colnames(geno)) # check that we have the same clones in the population structure data and in the genomic matrix

mainGP <- prop %>% dplyr::select(clon,max.Q) # extract the main gene pool for each clone.

geno.imp <- geno %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("clon") %>%  
  left_join(mainGP,by="clon") %>% 
  dplyr::select(clon,max.Q,everything())

geno.imp[1:10,1:10]

for(i in unique(geno.imp$max.Q)){
  subset <- geno.imp[geno.imp$max.Q==i,]
  subset <- apply(subset[,3:ncol(subset)], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
  geno.imp[geno.imp$max.Q==i,3:ncol(geno.imp)] <- subset
}

# Generating the PLINK format files (bim fam bed)
# ===============================================
geno.imp <- geno.imp %>% 
  column_to_rownames(var="clon") %>% 
  dplyr::select(-max.Q) %>% 
  t() %>% 
  as.matrix()

geno_plink <- write_plink("data/FEEMs/GenotypeMatrixPlink",geno.imp)


# Clone coordinates ####
########################

data <- readRDS(file="data/ClonapinData/PhenoDataNovember2019_AnnualTraits.rds")
clon.coord <- unique(data[,c("clon","longitude_prov","latitude_prov")])
clon.coord <- clon.coord[order(clon.coord$clon),]
write.table(clon.coord[,2:3], file = "data/FEEMs/clon_coord.coord", sep = "\t",row.names=FALSE,col.names = FALSE)




# Central latitude and longitude
################################


outer <- read.table("~/Documents/H2Pinpin/H2Pinpin/data/FEEMs/OuterPolygonMaritimePinePopulations.outer", 
                    quote="\"", comment.char="")
apply(outer,2,function(x) min(x) + (max(x) - min(x))/2)



### Create a grid
#################

library(dggridR)

# https://cran.microsoft.com/snapshot/2017-08-01/web/packages/dggridR/vignettes/dggridR.html

#Generate a dggs specifying an intercell spacing of ~25 miles
dggs      <- dgconstruct(spacing=25, metric=FALSE, resround='nearest',projection="ISEA",aperture=4,topology="TRIANGLE")
dgearthgrid(dggs, savegrid = "data/FEEMs/MyOwnGrid_XXresolution.shp") # can be very long!

dgearthgrid(dgconstruct(res=10, projection = "ISEA", aperture = 4, topology = "TRIANGLE"),
                        savegrid = "data/FEEMs/MyOwnGrid_triangle_res8.shp")
