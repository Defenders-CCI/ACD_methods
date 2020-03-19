library(dplyr)
library(googledrive)
library(MASS)
library(mvtnorm)
library(plotly)
library(pROC)
library(stringr)
library(tidyr)

##Authorize GDRive
drive_auth()

## READ FROM GOOGLE DRIVE
IW_folder <- '16XRyxU7PXpvGI4eYro9HNJLUBC1kg-z5'
MAD_folder <- '1B7TmgwQ96OabvtzstCZ53FttEPN2rRzG'

IWfiles <- drive_ls(as_id(IW_folder))
IWfiles <- filter(IWfiles, grepl('.csv', name))
IWids <- as_id(IWfiles)

MADfiles <- drive_ls(as_id(MAD_folder))
MADfiles <- filter(MADfiles, grepl('.csv', name))
MADids <- as_id(MADfiles)

IW <- readRDS(file = 'AlgorithmValidation/data/IW.rds')
MAD <- readRDS(file = 'AlgorithmValidation/data/MAD.rds')
files <- readRDS(file = 'AlgorithmValidation/data/files.rds')

# Create/update a
files <- bind_rows(MADfiles[!MADfiles$id %in% files$id],
                   IWfiles[!IWfiles$id %in% files$id])%>%
  mutate(loc = str_split_fixed(name, "_", 4)[,1],
         algorithm = str_split_fixed(name, "_", 4)[,2])%>%
  select(name, id, loc, algorithm)%>%
  bind_rows(files)

#Process IW files to standardize variables and combine
IW <- data.frame()
for(i in IWids[!IWids %in% files$id]){
  temp <- drive_download(file = as_id(i), path = 'temp.csv', overwrite = TRUE)
  f <- read.csv(file = 'temp.csv', sep = ",", header = TRUE, stringsAsFactors = FALSE)
  colnames(f) <- trimws(
    tolower(
      gsub("\\.", "", colnames(f))
    )
  )
  f$Loc <- strsplit(temp$name, "_")[[1]][1]
  f$Train <- rbinom(nrow(f), 1, 0.5)
  f$habitat <- str_to_title(f$habitat)
  f$disturbance <- str_to_title(f$disturbance)
  if(!'change' %in% colnames(f)){
    f$change <- str_split_fixed(f$systemindex, pattern = "_", n = 3)[,1]
    f$change[f$change == 2] <- 0
  }
  f$change <- as.character(f$change)
  f <- select(f, Loc, habitat, disturbance, cv_z, nbr_z, ndsi_z, ndvi_z, ndwi_z, rcvmax_z, change, Train)
  IW <- bind_rows(f, IW)
}

#Process MAD files to standardize variables and combine
MAD <- data.frame()
for(i in MADids[!MADids %in% files$id]){
  temp <- drive_download(file = as_id(MADids[i]), path = 'temp.csv', overwrite = TRUE)
  f <- read.csv(file = 'temp.csv', sep = ",", header = TRUE, stringsAsFactors = FALSE)
  colnames(f) <- trimws(
    tolower(
      gsub("\\.", "", colnames(f))
    )
  )
  f$Loc <- strsplit(temp$name, "_")[[1]][1]
  f$Train <- rbinom(nrow(f), 1, 0.5)
  f$habitat <- str_to_title(f$habitat)
  f$disturbance <- str_to_title(f$disturbance)
  if(!'change' %in% colnames(f)){
    f$change <- str_split_fixed(f$systemindex, pattern = "_", n = 3)[,1]
    f$change[f$change == 2] <- 0
  }
  f$change <- as.character(f$change)
  f <- select(f, Loc, habitat, disturbance, v1, v2, v3, v4, v5, v6, chi, change, Train)
  MAD <- bind_rows(f, MAD)
}

# Save data to rds files locally and on GDrive
saveRDS(files, file = 'AlgorithmValidation/data/files.rds')
saveRDS(MAD, file = "AlgorithmValidation/data/MAD.rds")
saveRDS(IW, file = "AlgorithmValidation/data/IW.rds")

drive_upload('AlgorithmValidation/data/IW.rds', path = as_id(IW_folder), name = 'IW.rds')
drive_upload('AlgorithmValidation/data/MAD.rds', path = as_id(MAD_folder), name = 'MAD.rds')

#Inspect the data
grid <- expand.grid(habitat = unique(IW$habitat), disturbance = unique(IW$disturbance))

plot_IW_zs <- function(i){
  a <- grid$habitat[i]
  b <- grid$disturbance[i]
  df <- gather(IW, "Metric", "Z", 2:6)%>%
    filter(habitat == a, disturbance == b)
  plot_ly(type = "box")%>%
    add_trace(data = df[df$Change == 1,], y = ~Z, x = ~Metric, name = "Yes")%>%
    add_trace(data = df[df$Change == 2,], y = ~Z, x = ~Metric, name = "No")%>%
    layout(xaxis = list(title = a),
           yaxis = list(title = b))
}

plot_MAD_zs <- function(i){
  a <- grid$habitat[i]
  b <- grid$disturbance[i]
  df <- filter(MAD, habitat == a, disturbance == b)
  plot_ly(type = "box")%>%
    add_trace(data = df[df$Change == 1,], y = ~chi,
              name = "Yes", legendgroup = ~Change)%>%
    add_trace(data = df[df$Change == 2,], y = ~chi,
              name = "No", legendgroup = ~Change)%>%
    layout(xaxis = list(title = a),
           yaxis = list(title = b))
}

s1 <- subplot(lapply(1:nrow(grid), plot_MAD_zs),
              nrows = 3, shareY = TRUE, titleY = TRUE, shareX =TRUE,titleX= TRUE)


