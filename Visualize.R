################################################################################
## Title: Bayesian inference in multivariate spatio-temporal areal models     ##
##        using INLA: analysis of gender-based violence in small areas        ##
##                                                                            ##
## Authors: Vicente, G. - Goicoa, T. -  Ugarte, M.D.                          ##
##                                                                            ##
## https://doi.org/10.1007/s00477-020-01808-x                                 ##
##                                                                            ##
################################################################################
##                          Spatio-temporal M-models                          ##
################################################################################
rm(list=ls())

library(grid)
library(INLA)
library(RColorBrewer)
library(sf)
library(tmap)
library(tmaptools)

options(dplyr.summarise.inform = FALSE)
options("max.print"=5000)


########################
## Some necessary packages (This only installs packages not yet installed)
########################
## Package names
packages <- c("spdep", "spData", "rgdal","fastDummies","reshape2","spatialreg","INLA","dplyr","sp","RColorBrewer","geodata","ggplot2","readxl","tidyverse")
## Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

## Packages loading
invisible(lapply(packages, library, character.only = TRUE))



data <- read.csv("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/data_dengue/data_ntime.csv")

load("results/MODELS_inla_icar_test.RData")
#load("results/MODELS_inla_lcar.RData")
#load("results/MODELS_inla_pcar.RData")
J <- length(unique(data$type))
S <- length(unique(data$area))
T <- length(unique(data$time))

t.from <- min(data$time)
t.to <- max(data$time)

##################################################################################
################                    Figures                       ################
##################################################################################

## Folder to save the figures ##
if(!file.exists("figures")) {dir.create("figures")}

Model_cluster <- MODELS.inla.icar_B

thai_map <- readOGR("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/shp_thai/gadm41_THA_1.shp")

## Maps of posterior exceedence probabilities ##
probs <- matrix(1-Model_cluster$summary.fitted.values$`1 cdf`,S*T,J,byrow=F)

paleta <- brewer.pal(8,"OrRd")[-1]
values <- c(0,0.95,1)

carto.DF <- thai_map
carto.DHF <- thai_map
carto.DSS <- thai_map
for(i in seq(T)){
  carto.DF$var <- matrix(probs[,1],S,T,byrow=F)[,i]
  carto.DHF$var <- matrix(probs[,2],S,T,byrow=F)[,i]
  carto.DSS$var <- matrix(probs[,3],S,T,byrow=F)[,i]
  
  names(carto.DF)[ncol(carto.DF)] <- paste("Year",i,sep=".")
  names(carto.DHF)[ncol(carto.DHF)] <- paste("Year",i,sep=".")
  names(carto.DSS)[ncol(carto.DSS)] <- paste("Year",i,sep=".")
}

Map.prob1 <- tm_shape(carto.DF) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="", main.title.position="left", panel.label.size=1,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=2, ncol=12)

Map.prob2 <- tm_shape(carto.DHF) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="", main.title.position="left", panel.label.size=1.5,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=10, ncol=24)

Map.prob3 <- tm_shape(carto.DSS) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="", main.title.position="left", panel.label.size=1.5,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=10, ncol=12)


## File: figure_6.pdf
tmap_save(tmap_arrange(Map.prob1, Map.prob2,Map.prob3), file="figures/prob_DF_all.jpg")
tmap_save( Map.prob1, width = 2000 ,height = 1000, file="figures/prob_DF.jpg")
tmap_save( Map.prob2, width = 1000 ,height = 800, file="figures/prob_DHF.jpg")
tmap_save( Map.prob3, width = 1000 ,height = 800, file="figures/prob_DSS.jpg")


Model = icar.t1.1.ff$association_model

predicted.p.value <- c()
n <- length(data[,1])
for(i in (1:n)) {
  predicted.p.value[i] <- inla.pmarginal(q=data$OBS[i],
                                         marginal=icar.t1.1.f$marginals.fitted.values[[i]])
}

plot(data$OBS[5082:length(data)],Model[["summary.fitted.values"]][["mean"]][5082:length(data)]*data$EXP[5082:length(data)],
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")

plot(data$OBS[5082:length(data)],Model[["summary.fitted.values"]][["mean"]][5082:length(data)]*data$EXP[5082:length(data)],
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")


plot(cor(data$OBS[5082:length(data)],Model[["summary.fitted.values"]][["mean"]][5082:length(data)]*data$EXP[5082:length(data)]))

mse_calc <- function(x,y){
  result = 0
  for (i in 1:length(x)){
    result = result + (x[i] - y[i])^2
  }
  result = result/(length(x))
  return(result)
}

data$OBS[5082:length(data)]
icar.t1.1.f$summary.fitted.values$mean[5082:length(data)]

print(mse_calc(data$OBS[5082:length(data)],Model$summary.fitted.values$mean[5082:length(data)] * data$EXP[5082:length(data)]))
cor(data$OBS,Model$summary.fitted.values$mean * data$EXP)



Model <- icar.t1.1.ff$association_model


t.from <- min(data$time)
t.to <- max(data$time)
## Selected districts ##
id.area<- c(3)

data <- data[
  order( data[,"area"], data[,"type"] ),]

start <- 1+(60*(as.numeric(id.area)-1))
stop <- 60+(60*(as.numeric(id.area)-1))

RR.DF <- Model$summary.fitted.values[(Model$.args$data$ID.area %in% id.area) & Model$.args$data$ID.disease==1,] 
RR.DHF <- Model$summary.fitted.values[(Model$.args$data$ID.area %in% id.area) & Model$.args$data$ID.disease==2,]
RR.DSS <- Model$summary.fitted.values[(Model$.args$data$ID.area %in% id.area) & Model$.args$data$ID.disease==3,]

RR.DF <- Model$summary.fitted.values[(Model$.args$data$ID.area %in% id.area) & Model$.args$data$ID.disease==1,] * data$EXP[start:stop]
RR.DHF <- Model$summary.fitted.values[(Model$.args$data$ID.area %in% id.area) & Model$.args$data$ID.disease==2,] * data$EXP[60+start:60+stop]
RR.DSS <- Model$summary.fitted.values[(Model$.args$data$ID.area %in% id.area) & Model$.args$data$ID.disease==3,] * data$EXP[120+start:120 +stop ]


OBS.DF <-  data$OBS[start:stop]
OBS.DHF <-  data$OBS[60+start:60+stop]
OBS.DSS <-  data$OBS[120+start:120 +stop ]

OBS_line <- c(OBS.DF,OBS.DHF,OBS.DSS)

aux <- lapply(list(DF=RR.DF,DHF=RR.DHF,DSS=RR.DSS), function(x){
  mean <- data.frame(mean=matrix(x[,"mean"],ncol=T,byrow=F))
  q1 <- data.frame(mean=matrix(x[,"0.025quant"],ncol=T,byrow=F))
  q2 <- data.frame(mean=matrix(x[,"0.975quant"],ncol=T,byrow=F))
  colnames(mean) <- colnames(q1) <- colnames(q1) <- seq(1,T)
  
  list(mean=mean, q1=q1, q2=q2)
})



inf <- min(unlist(aux))-0.05
top <- max(unlist(aux))+1.5
selected_colors <- c(rgb(154,192,205,alpha=150, maxColorValue=255),
                     rgb(69,139,116,alpha=150, maxColorValue=255),
                     rgb(10,50,20,alpha=150, maxColorValue=255))



graphics.off()

pdf("figures/figure_8.3.pdf", height=12, width=12, onefile=FALSE)
par(mfrow=c(1,1))

x <- seq(1,T)
main.title <- unique(paste(Model$.args$data$Area," (ID area ", Model$.args$data$ID.area,")",sep="")[Model$.args$data$ID.area %in% id.area])

for(i in seq(length(id.area))){
  
  plot(range(x), c(inf, top), type="n", xlab="Months", ylab="Number of cases",
       xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.main=2, main=main.title[i])
  
  axis(1, at=seq(1,T,2), labels=seq(t.from,t.to,2), las=0, cex.axis=1)
  
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(aux$DF$q1[i,], tail(aux$DF$q2[i,1]), rev(aux$DF$q2[i,]), aux$DF$q1[i,1])
  polygon(X.Vec, Y.Vec, col=selected_colors[1], border = NA)
  lines(1:T, aux$DF$mean[i,], pch = 18, col = "blue", type = "b", lty = 2)
  lines(1:T, OBS.DF, frame = FALSE, pch = 19, 
        col = "red", xlab = "x",type = "b", ylab = "y")
  
  
  #X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  #Y.Vec <- c(aux$DHF$q1[i,], tail(aux$DHF$q2[i,1]), rev(aux$DHF$q2[i,]), aux$DHF$q1[i,1])
  #polygon(X.Vec, Y.Vec, col=selected_colors[2], border = NA)
  #lines(1:T, aux$DHF$mean[i,], pch = 18, col = "blue", type = "b", lty = 2)
  #lines(1:T, OBS.DHF, ype = "b", frame = FALSE, pch = 19, 
  # col = "red", xlab = "x",type = "b", ylab = "y")
  
  
  #X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  #Y.Vec <- c(aux$DSS$q1[i,], tail(aux$DSS$q2[i,1]), rev(aux$DSS$q2[i,]), aux$DSS$q1[i,1])
  #polygon(X.Vec, Y.Vec, col=selected_colors[3], border = NA)
  #lines(1:T, aux$DSS$mean[i,], lwd=1.5)
  
  abline(h=1,lty=2)
  legend("topleft", legend=c("Observaion of case", "Prediction"),
         col=c("red", "blue"), lty = 2:4, cex=1.6)
  #legend("topleft", inset=.02, c("DFF","DHF","DSS"), fill=selected_colors, horiz=FALSE, cex=2.2, box.lty=0)
}

dev.off() 


map <- ggplot() + geom_polygon(data = thai_map, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

map + theme_void()

thai_map <- readOGR("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/shp_thai/gadm41_THA_1.shp")
thai_map@data$ID_area <- seq(78,154)

summary_random <- Model[["summary.random"]][[3]][78:154,]


#carto <- sf::st_as_sf(thai_map)
#carto <- carto[order(sf::st_set_geometry(carto, NULL)[,"NAME_1" ]),]


#data <- merge(x = carto[,c("NAME_1","ID_area")] , y = summary_random , by.x = "ID_area" , by.y = "ID")
thai_map@data <- merge(x = thai_map@data , y = summary_random , by.x = "ID_area" , by.y = "ID")

thai_map@data <- thai_map@data %>% mutate(Significant =
                                            case_when( 
                                              `0.025quant` / `0.975quant` < 0 ~ "no",
                                              `0.025quant` / `0.975quant` > 0 ~ "yes")
)


shp_df <- broom::tidy(thai_map, region = "Significant")
head(shp_df)



thai_map@data <- thai_map@data %>% mutate(id = row.names(.))
shp_df <- broom::tidy(thai_map, region = "id")
shp_df <- shp_df %>% left_join(thai_map@data, by = c("id"="id"))
map <- ggplot() + geom_polygon(data = shp_df, aes(x = long, y = lat, group = group, fill = Significant), colour = "black") + theme_void()

col <- c("NAME_1","mean","sd","0.025quant","0.5quant","0.975quant","Significant","type")
shp_info <- thai_map@data[,names(thai_map@data) %in% col]

map 
#data$geometry <- NULL
#data[,ID.year] <- paste(sprintf("%02d", as.numeric(as.character(data[,ID.year]))))
#data[,ID.disease] <- paste(sprintf("%02d", as.numeric(as.character(data[,ID.disease]))))
#data <- data[order(data[,ID.disease],data[,ID.year],data[,ID.area]),
