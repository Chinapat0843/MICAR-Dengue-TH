options(dplyr.summarise.inform = FALSE)
options("max.print"=5000)
library(INLA)

########################
## Some necessary packages (This only installs packages not yet installed)
########################
## Package names
packages <- c("spdep", "spData", "rgdal","fastDummies","reshape2","spatialreg","dplyr","sp","RColorBrewer","geodata","ggplot2","readxl","tidyverse")
## Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)## Packages loading
invisible(lapply(packages, library, character.only = TRUE))



# read data dengue
data <- read_csv("C:/Users/Public/senior/DataDengue_2011_5factor.csv")



J <- length(unique(data$type))
S <- length(unique(data$area))
T <- length(unique(data$time))

datafilter <- data[(nrow(data) - J*S)+1 :nrow(data),]
datafilter <- na.omit(datafilter)
datafilter$OBS <- NA
datafilter$time <- datafilter$time +1
data_combine <- rbind(data, datafilter)
#data_origin <- data_combine[1:(nrow(data_combine) - J*S),]
#J <- length(unique(data$type))
#S <- length(unique(data$area))
#T <- length(unique(data$time))


####################################
## 2) Fit spatio-temporal M-model ##
####################################

#load initial M-model 
source("C:/Users/Public/senior/model/Mmodel_icar.R")  ## Intrinsic multivariate CAR latent effect
#source("./INLA_Mmodel_lcar.R")  
#source("./INLA_Mmodel_pcar.R")  

## Fit a spatial-temporal multivariate Poisson mixed model to areal count data, 

## through the use of M-models Pipeline
source("C:/Users/Public/senior/model/MCAR_INLA_ST.R")
source("C:/Users/Public/senior/model/MCAR_INLA_ST_Prediction.R")
#source("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/R/Model_data_inla.R")
#load("results/MODELS_inla_icar_5y5f.RData")


############################################
## 2.1) Multivariate intrinsic CAR models ##
############################################

#Load geospatial Thailand map data
thai_map <- read_sf("C:/Users/Public/senior/THA_adm/THA_adm1.shp")
colnames(thai_map)[5] <- "area" #change column name in Thailand map data
#proj4string(thai_map) <- CRS("+proj=longlat +ellps=clrk66")


# run Icar
#MODELS.inla.icar <- MCAR_INLA_ST(carto=thai_map, data=data, ID.area="area", ID.year="time", ID.disease="type",
                       # O="OBS", E="EXP", spatial="intrinsic", temporal="rw1", interaction="TypeI",
                        #strategy="simplified.laplace")

MODELS.inla.icar <- MCAR_INLA_ST(carto=thai_map, data=data_combine,  ID.area="area", ID.year="time", ID.disease="type",
                                 O="OBS", E="EXP", spatial="intrinsic", temporal="rw1", interaction="none",
                                 strategy="simplified.laplace")

MODELS.inla.icar.prediction.2 <- MCAR_INLA_ST_Prediction(carto=thai_map, data=data_combine,  ID.area="area", ID.year="time", ID.disease="type",
                                                       O="OBS", E="EXP", spatial="intrinsic", temporal="rw2", interaction="TypeI",
                                                       strategy="simplified.laplace")



MODELS.inla.lcar.prediction.2 <- MCAR_INLA_ST_Prediction(carto=thai_map, data=data_combine,  ID.area="area", ID.year="time", ID.disease="type",
                                                       O="OBS", E="EXP", spatial="Leroux", temporal="rw2", interaction="TypeI",
                                                       strategy="simplified.laplace")



MODELS.inla.pcar.prediction <- MCAR_INLA_ST_Prediction(carto=thai_map, data=data_combine,  ID.area="area", ID.year="time", ID.disease="type",
                                                       O="OBS", E="EXP", spatial="proper", temporal="rw1", interaction="TypeI",
                                                       strategy="simplified.laplace")






# run Lcar
#Lcar.t1.1 <- MCAR_INLA_ST(carto=thai_map, data=data, ID.area="NAME_1", ID.year="time", ID.disease="type",
                          #O="OBS", E="EXP", spatial="Leroux", temporal="rw1", interaction="TypeI",
                          #strategy="simplified.laplace")

# run Pcar
#pcar.t1.1 <- MCAR_INLA_ST(carto=thai_map, data=data, ID.area="NAME_1", ID.year="time", ID.disease="type",
                          #O="OBS", E="EXP", spatial="proper", temporal="rw1", interaction="TypeI",
                          #strategy="simplified.laplace"


#MODELS.inla.n <- list(icar_TypeI_RW1=icar.t1.1, icar_TypeI_RW2 = icar.t1.2, lcar_TypeI_RW1=lcar.t1.1, lcar_TypeI_RW2=lcar.t1.2, pcar_TypeI_RW1=pcar.t1.1)


# Save model results into local host
#MODELS.inla.icar <- list(TypeI_RW1=icar.t1.1.ff)
save(MODELS.inla.pcar.prediction, file=paste0("./results/", gsub("\\.", "_", "MODELS.inla.pcar.prediction.RW1.2011"), ".RData"))
save(MODELS.inla.pcar.prediction.2, file=paste0("./results/", gsub("\\.", "_", "MODELS.inla.pcar.prediction.RW2.2011"), ".RData"))
.
load("./results/MODELS_inla_icar_prediction_2016.RData")

Model <- MODELS.inla.icar.prediction

summary_fitted_values <-  Model$summary.fitted.values[(nrow(data_combine) - J*S)+1 :nrow(data_combine),][,1:5]
summary_fitted_values <- na.omit(summary_fitted_values)


summary_fitted_values <- round(summary_fitted_values, digits = 0)
summary_fitted_values

type_list <- vector()
for (j in 1:J) {
  for (s in 1:S) {
    type_list <- c(type_list, rep(j))
  }
}

type_S <- vector()
for (j in 1:J) {
  for (s in 1:S) {
    type_S <- c(type_S, rep(s))
  }
}

summary_fitted_values$ID <- type_S
summary_fitted_values$type <- type_list
num_col <- length((summary_fitted_values))
summary_fitted_values <- summary_fitted_values[,2:num_col]
thai_map$ID <- seq(1:77)
name_area <- thai_map[,c('area','ID')]

data_predict <- merge(x = summary_fitted_values , y = name_area , by.x = "ID" , by.y = "ID")

data_predict <- data_predict[,c('ID','area','mean','type')]

