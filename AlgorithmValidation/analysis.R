library(dplyr)
library(googledrive)
library(MASS)
library(mvtnorm)
library(plotly)
library(pROC)
library(stringr)
library(tidyr)

# Load helper functions
source('AlgorithmValidation/functions.R')

# Get our master data from Google Drive
IW_folder <- '16XRyxU7PXpvGI4eYro9HNJLUBC1kg-z5'
MAD_folder <- '1B7TmgwQ96OabvtzstCZ53FttEPN2rRzG'

IWfiles <- drive_ls(as_id(IW_folder))
IWid <- as_id(filter(IWfiles, grepl("IW.rds", name)))
MADfiles <- drive_ls(as_id(MAD_folder))
MADid <- as_id(filter(MADfiles, grepl("MAD.rds", name)))

drive_download(file = IWid, path = 'AlgorithmValidation/data/IW.rds', overwrite = TRUE)
drive_download(file = MADid, path = 'AlgorithmValidation/data/MAD.rds', overwrite = TRUE)

#randomForest
forest.out <- randomForest(as.factor(Change) ~ rcvmax_z+cv_z+nbr_z+ndsi_z+ndvi_z,
                           data = IW,
                           importance = TRUE, ntree = 400, mtry = 5, keep.forest = FALSE)


#Create ROC curves for generic IW and MAD
generic <- lda_analysis("", "Bare|None")
plot_ROC_crv(generic$IW, generic$MAD)

shrub <- lda_analysis("Shrub", "")
plot_ROC_crv(shrub$IW, shrub$MAD)
shrub_bare <- lda_analysis("Shrub", "Bare|None")
shrub_solar <- lda_analysis("Shrub", "Solar|None")

forest <- lda_analysis("Forest", "")
plot_ROC_crv(forest$IW, forest$MAD)

desert <- lda_analysis("Desert", "")
plot_ROC_crv(desert$IW, desert$MAD)
desert_bare <- lda_analysis("Desert","Bare|None")
desert_solar <- lda_analysis("Desert", "Solar|None")
desert_residential <- lda_analysis("Desert", "Residential|None")

grassland <- lda_analysis("Grassland", "")
plot_ROC_crv(grassland$IW, grassland$MAD)
grassland_bare <- lda_analysis("Grassland", "Bare|None")
grassland_solar <- lda_analysis("Grassland", "Solar|None")

grassland_a <- lda_analysis("Grassland", "[^Bare]")
plot_ROC_crv(grassland_a$IW, grassland_a$MAD)


wetland <- lda_analysis("Wetland", "")
plot_ROC_crv(wetland$IW, wetland$MAD)
