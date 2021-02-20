############################################################################################"
##################  Script from Juliette Archambeau - July 2019      #######################"
##################         CLEANING PHENOTYPIC DATABASE              #######################"
############################################################################################"


# Libraries
library(readr)
library(tidyverse)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

#                  Import phenotypic data                           ####
# ====================================================================="



# Asturias
dataA <- read_csv("data/ClonapinData/CLONAPIN_data20191119_asturias.csv")

dataA <- dataA %>% dplyr::select(region =Region,    # mediterranean/atlantic
                               site,                # trial site
                               metap,               # metapopulation
                               prov,                # provenance 
                               clon_original,       # clone name in trial site raw data
                               block_original,      # block number in trial site
                               clon,                # unified and between sites comparable clone name
                               block,               # unified and between sites comparable bloc name
                               clones522,           # selection of 522 clones done by I. Rodriguez-Quilon
                               contains("_sur"),    # survival
                               contains("_ht"),     # height
                               -contains("YEMA"))   # remove Asturias phenoloty estimation by the height of whorl 2012



# Bordeaux
dataB <- read_csv("data/ClonapinData/CLONAPIN_data20191119_bordeaux.csv")

dataB <- dataB %>% dplyr::select(region =Region,    # mediterranean/atlantic
                               site,                # trial site
                               metap,               # metapopulation
                               prov,                # provenance 
                               clon_original,       # clone name in trial site raw data
                               block_original,      # block number in trial site
                               clon=clone,          # unified and between sites comparable clone name
                               block,               # unified and between sites comparable bloc name
                               clones522,           # selection of 522 clones done by I. Rodriguez-Quilon
                               contains("_sur"),    # survival
                               contains("_ht"),     # height
                               contains("BDX_s1"),  # Bordeaux phenology degrees-day stage 1
                               contains("BDX_s3"),  # Bordeaux phenology degrees-day stage 3
                               contains("BDX_s4"),  # Bordeaux phenology degrees-day stage 4
                               contains("BDX_lbb"), # Bordeaux overall duration of the whole growth cycle
                               -contains("YEMA")) %>% 
                  dplyr::mutate(BDX_dbb2013 = BDX_s42013 - BDX_s12013)


# Portugal
dataP <- read_csv("data/ClonapinData/CLONAPIN_data20191119_portugal.csv")

dataP <- dataP %>% dplyr::select(region =Region,      # mediterranean/atlantic
                                 site,                # trial site
                                 metap,               # metapopulation
                                 prov,                # provenance 
                                 clon_original,       # clone name in trial site raw data
                                 block_original,      # block number in trial site
                                 clon,                # unified and between sites comparable clone name
                                 block,               # unified and between sites comparable bloc name
                                 clones522,           # selection of 522 clones done by I. Rodriguez-Quilon
                                 contains("_sur"),    # survival
                                 contains("_ht"),     # height
                                 contains("d13c"),    # Portugal isotopic composition of carbon 13
                                 contains("POR_SLA"), # Portugal Specific Leaf Area
                                 -contains("YEMA"))   # remove height in of the tree in spring 2012 (for phenology estimation)

data <- bind_rows(dataB,dataP,dataA)



# Checking the phenology calculation (it's ok)

# # 2014
# df <- dataB
# df$diff <- df$BDX_s42014 - df$BDX_s12014
# df$diffabc <- ifelse(df$diff > df$BDX_lbb2014, 'A',ifelse(df$diff < df$BDX_lbb2014, 'B', 'C'))
# unique(df$diffabc)
# df[df$diffabc=="A"& !is.na(df$diffabc),c("BDX_s12014","BDX_s42014","BDX_lbb2014","diff","diffabc")]
# 
# # 2015
# df$diff <- df$BDX_s42015 - df$BDX_s12015
# df$diffabc <- ifelse(df$diff > df$BDX_lbb2015, 'A',ifelse(df$diff < df$BDX_lbb2015, 'B', 'C'))
# unique(df$diffabc)
# df[df$diffabc=="A"& !is.na(df$diffabc),c("BDX_s12015","BDX_s42015","BDX_lbb2015","diff","diffabc")]
# 
# # 2017
# df$diff <- df$BDX_S4_2017 - df$BDX_S1_2017
# df$diffabc <- ifelse(df$diff > df$BDX_lbb_2017, 'A',ifelse(df$diff < df$BDX_lbb_2017, 'B', 'C'))
# unique(df$diffabc)
# df[df$diffabc=="A"& !is.na(df$diffabc),c("BDX_S1_2017","BDX_S4_2017","BDX_lbb_2017","diff","diffabc")]

#*************************************************************************************************************************


#            Format the data                            ####
# ========================================================="


# A column to define a different IDs per tree (ID = clon_block)
data <- data %>% 
  mutate(tree = NA) %>% 
  unite(tree,clon,block, sep="_",remove=FALSE) 


# !!! UNITS DIFFERENCES BETWEEN SITES !!!
# In the file "CLONAPIN_data20191119.csv", height measurements are:
# in cm for bordeaux and asturias
# in mm for caceres, madrid and portugal
# We are going to use mm for all sites:
col.cm <- c("AST_htdic11_cm",
            "AST_htnov12_cm",
            "AST_htmar14_cm",
            "BDX_htnov13_cm",
            "BDX_htnov14_cm",
            "BDX_htnov15_cm",
            "BDX_htnov18_cm")
data[,col.cm] <- data[,col.cm] * 10


# remane columns
data <- data %>% rename(AST_survdec11 = AST_survdic11,
                        AST_htdec11 = AST_htdic11_cm, 
                        AST_htnov12 = AST_htnov12_cm, 
                        AST_htmar14 = AST_htmar14_cm,
                        BDX_htnov13 = BDX_htnov13_cm,
                        BDX_htnov14 = BDX_htnov14_cm,
                        BDX_htnov15 = BDX_htnov15_cm,
                        BDX_htnov18 = BDX_htnov18_cm,
                        POR_htjan12 = POR_htjan12_mm,
                        POR_htmay12 = POR_htmay12_mm,
                        POR_htoct12 = POR_htoct12_mm,
                        POR_htmay13 = POR_htmay13_mm,
                        BDX_s32017 = BDX_S3_2017,
                        BDX_s42017 = BDX_S4_2017,
                        BDX_s12017 = BDX_S1_2017,
                        BDX_dbb2014 = BDX_lbb2014,
                        BDX_dbb2015 = BDX_lbb2015,
                        BDX_dbb2017 = BDX_lbb_2017)


# !! OUTLIERS !!
# ggplot(data, aes(x= BDX_htnov13)) + geom_histogram(binwidth = 20) +  theme_bw()
# ggplot(data, aes(x= BDX_htnov14)) + geom_histogram(binwidth = 20) +  theme_bw()
# ggplot(data, aes(x= BDX_htnov15)) + geom_histogram(binwidth = 20) +  theme_bw()
# ggplot(data, aes(x= BDX_htnov18)) + geom_histogram(binwidth = 20) +  theme_bw()
# ggplot(data, aes(x= BDX_s32013)) + geom_histogram(binwidth = 10) +  theme_bw()
# ggplot(data, aes(x= BDX_s32014)) + geom_histogram(binwidth = 10) +  theme_bw()
# ggplot(data, aes(x= BDX_s32015)) + geom_histogram(binwidth = 10) +  theme_bw()
# ggplot(data, aes(x= BDX_s32017)) + geom_histogram(binwidth = 10) +  theme_bw()
# ggplot(data, aes(x= POR_htjan12)) + geom_histogram(binwidth = 10) +  theme_bw()
# ggplot(data, aes(x= POR_htmay12)) + geom_histogram(binwidth = 10) +  theme_bw()
# ggplot(data, aes(x= POR_htoct12)) + geom_histogram(binwidth = 10) +  theme_bw()
# ggplot(data, aes(x= POR_htmay13)) + geom_histogram(binwidth = 10) +  theme_bw()
# ggplot(data, aes(x= d13C)) + geom_histogram(binwidth = 0.05) +  theme_bw()
# ggplot(data, aes(x= AST_htdec11)) + geom_histogram(binwidth = 10) +  theme_bw()
# ggplot(data, aes(x= AST_htnov12)) + geom_histogram(binwidth = 10) +  theme_bw()
# ggplot(data, aes(x= AST_htmar14)) + geom_histogram(binwidth = 10) +  theme_bw() # two outliers here


data[data$tree=="VAL16_7","AST_htmar14"] <- data[data$tree=="VAL16_7","AST_htmar14"] / 10
data[data$tree=="VER12_3","AST_htmar14"] <- data[data$tree=="VER12_3","AST_htmar14"] / 10



#*************************************************************************************************************************


#                    Relative contribution of each gene pool                            #####
# =========================================================================================="



geno <- read.delim("data/ClonapinData/ClonapinBlups523IndPiMASSJuly2019.txt")
geno <- geno[,c("X","prov","Q1","Q2","Q3","Q4","Q5")]
geno <- geno %>% rename(clon=X)

# Checking NAs
sapply(geno, function(x) sum(is.na(x)))

# Adding Q6
geno$Q6 <- 1 - (geno$Q1 + geno$Q2 + geno$Q3 + geno$Q4 + geno$Q5)

# Adding the cluster the most impt for each clone
for (i in 1:length(geno$clon)){
  geno[i,"max.Qvalue"] <- max(geno[i,3:8])
  geno[i,"max.Q"] <- names(geno[i,3:8])[which.max(apply(geno[i,3:8],MARGIN=2,max))]
}
geno$max.Q <- as.factor(geno$max.Q)


# 18 clones have no structure geno
setdiff(unique(data %>% pull(clon)),unique(geno[,"clon"]))
setdiff(unique(geno[,"clon"]),unique(data %>% pull(clon)))
length(setdiff(unique(data %>% pull(clon)),unique(geno[,"clon"])))


### Let's look at these clones that has no structure data:

# !! CAS12 is probably CAR12 !!
data[data$clon_original=="Carb_F12_P1",]
data$prov[data$clon=="CAS12"] <- "CAR"
data$tree[data$clon=="CAS12"] <- "CAR12_7"
data$clon[data$clon=="CAS12"] <- "CAR12"

data[data$clon=="SAL6",] # is only in Bordeaux
data[data$prov=="SAC"&data$site=="portugal",] # the original clone SSanC_F2_P5 is alone 
data[data$clon=="PIE27",] # a lot of this clone in Asturias and Portugal, but no geno data
data[data$clon=="PIA8",]  # a lot of this clone in Asturias and Portugal, but no geno data
data[data$clon=="OLO17",]  
data[data$clon=="OLBsn",]  # only in Bordeaux

# !! CAD14 is probably CAR14 !!
data[data$clon_original=="CAR_14",] 
data$prov[data$clon=="CAD14"] <- "CAR"
data$tree[data$clon=="CAD14"] <- "CAR14_13"
data$clon[data$clon=="CAD14"] <- "CAR14"

# Email of Santi the 16/02/2021: he thinks MAD of Bordeaux and MAD1 of Portugal adn Asturias are the same
data[data$clon=="MAD"|data$clon=="MAD1",] %>% tbl_df %>% print(n = Inf)
data$tree[data$clon=="MAD"] <- paste0(stri_sub(g,1,3),1,stri_sub(g,4,-1))
data$clon[data$clon=="MAD"] <- "MAD1"


data[data$clon=="BAY18",] # a lot of this clone in Asturias, Bordeaux and Portugal


data <- merge(data,geno,by=c("clon","prov"),all.x=TRUE)

# Removing clones for which there is no population structure data
sapply(data, function(x) sum(is.na(x)))
data <- data[!(is.na(data$Q1)),]
sapply(data, function(x) sum(is.na(x)))



#*************************************************************************************************************************


#                       Resurrected trees                           ####
# ====================================================================="




## REMOVING THE "RESURRECTED TREES"
######## Some trees are noted as dead at a given measurement, but as alive in some of the following measurements. 
######## For these trees, we replace 0 (which means that the tree is dead) with 1 (the tree is alive) because we 
######## hypothesize that the trees were considered dead and might not be seen in the experiment (many tall grasses hiding the trees), 
######## but were actually alive (and noted as alive in the following measures). See Marina's email of 07/01/2019.

# BORDEAUX -> no trees concerned
df <- subset(data, site=="bordeaux"&(BDX_surv12==0|is.na(BDX_surv12)|
                                     BDX_surv13==0|is.na(BDX_surv13)|
                                     BDX_surv14==0|is.na(BDX_surv15))&BDX_surv15==1)
df <- subset(data, site=="bordeaux"&(BDX_surv12==0|is.na(BDX_surv12)|
                                     BDX_surv13==0|is.na(BDX_surv13))&BDX_surv14==1)
df <- subset(data, site=="bordeaux"&(BDX_surv12==0|is.na(BDX_surv12))&BDX_surv13==1)

# ASTURIAS
df <- subset(data, site=="asturias"&(AST_survdec11==0|is.na(AST_survdec11)|
                                     AST_survnov12==0|is.na(AST_survnov12))&AST_survmar14==1)
data <- anti_join(data,df)
df$AST_survdec11<- 1
df$AST_survnov12 <- 1
data <- rbind(data,df)

df <- subset(data, site=="asturias"&(AST_survdec11==0|is.na(AST_survdec11))&AST_survnov12==1)
# data <- anti_join(data,df)
# df$AST_survdec11<- 1
# data <- rbind(data,df)



# PORTUGAL
df <- subset(data, site=="portugal"&(POR_survjan12==0|is.na(POR_survjan12)|
                                     POR_survmay12==0|is.na(POR_survmay12)|
                                     POR_survoct12==0|is.na(POR_survoct12))&POR_survmay13==1)
# data <- anti_join(data,df)
# df$POR_survjan12 <- 1
# df$POR_survmay12<- 1
# df$POR_survoct12 <- 1
# data <- rbind(data,df)

df <- subset(data, site=="portugal"&(POR_survjan12==0|is.na(POR_survjan12)|
                                       POR_survmay12==0|is.na(POR_survmay12))&POR_survoct12==1) # no trees concerned

df <- subset(data, site=="portugal"&(POR_survjan12==0|is.na(POR_survjan12))&POR_survmay12==1)
# data <- anti_join(data,df)
# df$POR_survjan12 <- 1
# data <- rbind(data,df)


# Trees were planted:
# 02/2011 in Asturias
# 04/2011 in Caceres
# 11/2011 in Bordeaux
# 11/2010 in Madrid
# 02/2011 in Portugal



#*************************************************************************************************************************

dataSave <- data 

#                        Dataframe with one column per annual trait             ####
#=================================================================================="

# Three rows in which there is only NAs
dataNA <- data %>% select(d13C,contains("BDX"),contains("POR"),contains("AST"))
which(rowSums(is.na(dataNA))==ncol(dataNA))
dataNA <- dataNA[which(rowSums(is.na(dataNA))==ncol(dataNA)),]
dataNA
data <- anti_join(data,dataNA)


### Remove prov "SID" that is only in Bordeaux 
data <- data[!(data$prov=="SID"),]

### Remove the provenance ROD as it had not been genotyped
sapply(data, function(x) sum(is.na(x)))
data <- data[!(is.na(data$Q1)),]
sapply(data, function(x) sum(is.na(x)))


### ADDING LATITUDE AND LONGITUDE
# Import & clean provenance coordinates
prov.coord <- read_csv("data/ClonapinData/coordinates_provenances.csv")
prov.coord <- prov.coord[,c("CODE","ALTITUDE","LATITUDE","LONGITUDE")] %>%
  rename(prov=CODE,altitude_prov=ALTITUDE,latitude_prov=LATITUDE, longitude_prov=LONGITUDE) %>%
  slice(1:35) # remove the last line with only NAs (which remained from the original table)

# Merge provenance coordinates & phenotypic data
data <- merge(data,prov.coord,by="prov",all.x=T)
data <- as_tibble(data)

# Import & clean site coordinates
site.coord <- read_csv("data/ClonapinData/coordinates_sites.csv")

# Merge site coordinates & phenotypic data
data <- merge(data,site.coord,by="site")
data <- as_tibble(data)



### Save ###################"
saveRDS(data,file="data/ClonapinData/PhenoDataNovember2019_AnnualTraits.rds")
############################"

#*************************************************************************************************************************

data <- dataSave

#                            Zombie trees                           ####
# ====================================================================="


### REMOVING ZOMBIE TREES
## Here, we delete measurements where trees were noted as "dead" 
# when they were already noted as dead in previous measurements. 

### BORDEAUX
df <- subset(data, site=="bordeaux"&BDX_surv12==0&BDX_surv13==0&BDX_surv14==0&BDX_surv15==0&BDX_surv18==0)
data <- anti_join(data,df)
df$BDX_surv13 <- NA
df$BDX_surv14 <- NA
df$BDX_surv15 <- NA
data <- rbind(data,df)

df <- subset(data, site=="bordeaux"&BDX_surv13==0&BDX_surv14==0&BDX_surv15==0&BDX_surv18==0)
data <- anti_join(data,df)
df$BDX_surv14 <- NA
df$BDX_surv15 <- NA
data <- rbind(data,df)

df <- subset(data, site=="bordeaux"&BDX_surv14==0&BDX_surv15==0&BDX_surv18==0)
data <- anti_join(data,df)
df$BDX_surv15 <- NA
df$BDX_surv18 <- NA
data <- rbind(data,df)

df <- subset(data, site=="bordeaux"&BDX_surv15==0&BDX_surv18==0) # EMPTY
# data <- anti_join(data,df)
# df$BDX_surv18 <- NA
# data <- rbind(data,df)

### ASTURIAS
df <- subset(data, site=="asturias"&AST_survdec11==0&AST_survnov12==0&AST_survmar14==0)
data <- anti_join(data,df)
df$AST_survnov12<- NA
df$AST_survmar14 <- NA
data <- rbind(data,df)

df <- subset(data, site=="asturias"&AST_survnov12==0&AST_survmar14==0)
data <- anti_join(data,df)
df$AST_survmar14 <- NA
data <- rbind(data,df)

### PORTUGAL
df <- subset(data, site=="portugal"&POR_survjan12==0&POR_survmay12==0&POR_survoct12==0&POR_survmay13==0)
data <- anti_join(data,df)
df$POR_survmay12<- NA
df$POR_survoct12 <- NA
df$POR_survmay13 <- NA
data <- rbind(data,df)

df <- subset(data, site=="portugal"&POR_survmay12==0&POR_survoct12==0&POR_survmay13==0)
data <- anti_join(data,df)
df$POR_survoct12 <- NA
df$POR_survmay13 <- NA
data <- rbind(data,df)

df <- subset(data, site=="portugal"&POR_survoct12==0&POR_survmay13==0)
data <- anti_join(data,df)
df$POR_survmay13 <- NA
data <- rbind(data,df)





#                        Dataframe with one column per trait                    ####
#=================================================================================="



###############################  > ASTURIAS  ##################################################
####### december 2011
AST_dec11 <-  data[!(is.na(data$AST_survdec11)&is.na(data$AST_htdec11)),]
# 216 trees for which there was a value for survival but not for height. 
AST_dec11$survival <- AST_dec11$AST_survdec11
AST_dec11$height <- AST_dec11$AST_htdec11
AST_dec11$age <- 10
AST_dec11 <- AST_dec11 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### november 2012
AST_nov12 <-  data[!(is.na(data$AST_survnov12)&is.na(data$AST_htnov12)),]
AST_nov12$survival <- AST_nov12$AST_survnov12
AST_nov12$height <- AST_nov12$AST_htnov12
AST_nov12$age <- 21
AST_nov12 <- AST_nov12 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### march 2014
AST_mar14 <-  data[!(is.na(data$AST_survmar14)&is.na(data$AST_htmar14)),]
AST_mar14$survival <- AST_mar14$AST_survmar14
AST_mar14$height <- AST_mar14$AST_htmar14
AST_mar14$age <- 37
AST_mar14 <- AST_mar14 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)



###############################  > BORDEAUX  ##################################################
####### November 2012: survival
BDX_nov12 <-  data[!(is.na(data$BDX_surv12)),]
BDX_nov12$survival <- BDX_nov12$BDX_surv12
BDX_nov12$height <- NA
BDX_nov12$age <- 13
BDX_nov12 <- BDX_nov12 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### November 2013: survival and height
BDX_nov13 <-  data[!(is.na(data$BDX_surv13)&is.na(data$BDX_htnov13)),]
BDX_nov13$survival <- BDX_nov13$BDX_surv13
BDX_nov13$height <- BDX_nov13$BDX_htnov13
BDX_nov13$age <- 25
BDX_nov13 <- BDX_nov13 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)


####### Phenology 2013: duration of bud burst and bud burst
BDX_pheno13 <-  data[!(is.na(data$BDX_s32013)&is.na(data$BDX_dbb2013)),]
BDX_pheno13$bb <- BDX_pheno13$BDX_s32013
BDX_pheno13$dbb <- BDX_pheno13$BDX_dbb2013
BDX_pheno13$age <- 15 # age in January 2013
BDX_pheno13 <- BDX_pheno13 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, bb, dbb)


####### November 2014: survival and height
BDX_nov14 <-  data[!(is.na(data$BDX_surv14)&is.na(data$BDX_htnov14)),]
BDX_nov14$survival <- BDX_nov14$BDX_surv14
BDX_nov14$height <- BDX_nov14$BDX_htnov14
BDX_nov14$age <- 37
BDX_nov14 <- BDX_nov14 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### Phenology 2014: duration of bud burst and bud burst
BDX_pheno14 <-  data[!(is.na(data$BDX_s32014)&is.na(data$BDX_dbb2014)),]
BDX_pheno14$bb <- BDX_pheno14$BDX_s32014
BDX_pheno14$dbb <- BDX_pheno14$BDX_dbb2014
BDX_pheno14$age <- 27 # age in January 2014
BDX_pheno14 <- BDX_pheno14 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                    clones522, age, bb, dbb)

####### November 2015: survival and height
BDX_nov15 <-  data[!(is.na(data$BDX_surv15)&is.na(data$BDX_htnov15)),]
BDX_nov15$survival <- BDX_nov15$BDX_surv15
BDX_nov15$height <- BDX_nov15$BDX_htnov15
BDX_nov15$age <- 49
BDX_nov15 <- BDX_nov15 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### Phenology 2015: duration of bud burst and bud burst
BDX_pheno15 <-  data[!(is.na(data$BDX_s32015)&is.na(data$BDX_dbb2015)),]
BDX_pheno15$bb <- BDX_pheno15$BDX_s32015
BDX_pheno15$dbb <- BDX_pheno15$BDX_dbb2015
BDX_pheno15$age <- 29 # age in January 2015
BDX_pheno15 <- BDX_pheno15 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                      clones522, age, bb, dbb)

####### Phenology 2017: duration of bud burst and bud burst
BDX_pheno17 <-  data[!(is.na(data$BDX_s32017)&is.na(data$BDX_dbb2017)),]
BDX_pheno17$bb <- BDX_pheno17$BDX_s32017
BDX_pheno17$dbb <- BDX_pheno17$BDX_dbb2017
BDX_pheno17$age <- 53 # age in January 2017
BDX_pheno17 <- BDX_pheno17 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                      clones522, age, bb, dbb)

####### November 2018: survival and height 
BDX_nov18 <-  data[!(is.na(data$BDX_surv18)&is.na(data$BDX_htnov18)),]
BDX_nov18$survival <- BDX_nov18$BDX_surv18
BDX_nov18$height <- BDX_nov18$BDX_htnov18
BDX_nov18$age <- 85
BDX_nov18 <- BDX_nov18 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)


###############################  > CACERES  ##################################################"
####### december 2011
# CAC_dec11 <-  data[!(is.na(data$CAC_survdec11)&is.na(data$CAC_htdec11)),]
# CAC_dec11$survival <- CAC_dec11$CAC_survdec11
# CAC_dec11$height <- CAC_dec11$CAC_htdec11
# CAC_dec11$age <- 8
# CAC_dec11 <- CAC_dec11 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
#                                   clones522, age, survival, height)
# 
# 
# ###############################  > MADRID  ##################################################"
# ####### december 2011
# MAD_dec11 <-  data[!(is.na(data$MAD_survdec11)&is.na(data$MAD_htdec11)),]
# MAD_dec11$survival <- MAD_dec11$MAD_survdec11
# MAD_dec11$height <- MAD_dec11$MAD_htdec11
# MAD_dec11$age <- 13
# MAD_dec11 <- MAD_dec11 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
#                                   clones522, age, survival, height)


############################### > PORTUGAL  ##################################################
####### january 2012
POR_jan12 <-  data[!(is.na(data$POR_survjan12)&is.na(data$POR_htjan12)),]
POR_jan12$survival <- POR_jan12$POR_survjan12
POR_jan12$height <- POR_jan12$POR_htjan12
POR_jan12$age <- 11
POR_jan12 <- POR_jan12 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### may 2012
POR_may12 <-  data[!(is.na(data$POR_survmay12)&is.na(data$POR_htmay12)),]
POR_may12$survival <- POR_may12$POR_survmay12
POR_may12$height <- POR_may12$POR_htmay12
POR_may12$age <- 15
POR_may12 <- POR_may12 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### october 2012
POR_oct12 <-  data[!(is.na(data$POR_survoct12)&is.na(data$POR_htoct12)),]
POR_oct12$survival <- POR_oct12$POR_survoct12
POR_oct12$height <- POR_oct12$POR_htoct12
POR_oct12$age <- 20
POR_oct12 <- POR_oct12 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### may 2013
POR_may13 <-  data[!(is.na(data$POR_survmay13)&is.na(data$POR_htmay13)),]
POR_may13$survival <- POR_may13$POR_survmay13
POR_may13$height <- POR_may13$POR_htmay13
POR_may13$age <- 27
POR_may13 <- POR_may13 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)


####### d13C
POR_d13C <-  data[!(is.na(data$d13C)),]
POR_d13C$age <- NA
POR_d13C <- POR_d13C %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                      clones522, age, d13C)

####### SLA
POR_SLA <-  data[!(is.na(data$POR_SLA)),]
POR_SLA$SLA <- POR_SLA$POR_SLA
POR_SLA$age <- NA
POR_SLA <- POR_SLA %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                clones522, age, SLA)


##########################################"
##########################################"

### Merge all
total <- bind_rows(AST_dec11,AST_nov12,AST_mar14,
               BDX_nov12,BDX_nov13,BDX_nov14,BDX_nov15,BDX_nov18,
               BDX_pheno13,BDX_pheno14,BDX_pheno15,BDX_pheno17,
               #CAC_dec11,MAD_dec11, 
               POR_jan12,POR_may12,POR_oct12,POR_may13,
               POR_d13C,POR_SLA)

### quick check
table(total$age,total$site)

### A column to define a different IDs per obs (several osb per tree)
total <- total %>% mutate(obs = tree, num=NA) %>% group_by(tree) %>%  mutate(num = 1:n()) %>% 
  unite(obs, obs, num, sep="_") %>% select(1:5,7,10,6,9,8,obs,11:ncol(total))
length(unique(total$obs))

######################################"
## MISSING VALUES
## 8 trees have only NAs -> have not been included during the merging
names.total <- unique(total$tree) 
names.data <- unique(data$tree)
diff <- setdiff(names.data,names.total)

missing.values <- as.data.frame(matrix(NA,length(diff),ncol(data),dimnames = list(c(diff),c(colnames(data)))))
for (i in diff){
  missing.values[i,] <- data[data$tree==i,]
}
######################################"


### Remove prov "SID" that is only in Bordeaux 
total <- total[!(total$prov=="SID"),]


### ADDING LATITUDE AND LONGITUDE
# Import & clean provenance coordinates
prov.coord <- read_csv("data/ClonapinData/coordinates_provenances.csv")
prov.coord <- prov.coord[,c("CODE","ALTITUDE","LATITUDE","LONGITUDE")] %>%
  rename(prov=CODE,altitude_prov=ALTITUDE,latitude_prov=LATITUDE, longitude_prov=LONGITUDE) %>%
  slice(1:35) # remove the last line with only NAs (which remained from the original table)

# Merge provenance coordinates & phenotypic data
total <- merge(total,prov.coord,by="prov",all.x=T)
total <- as_tibble(total)

# Import & clean site coordinates
site.coord <- read_csv("data/ClonapinData/coordinates_sites.csv")

# Merge site coordinates & phenotypic data
total <- merge(total,site.coord,by="site")
total <- as_tibble(total)


### Save ###################"
saveRDS(total,file="data/ClonapinData/PhenoDataNovember2019.rds")
############################"



