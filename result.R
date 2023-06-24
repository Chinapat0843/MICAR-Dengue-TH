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

#data <- read.csv("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/data_dengue/Final_Data_Dengue_10y.csv")
#J <- length(unique(data$Severity))
#S <- length(unique(data$NAME_1))
#T <- length(unique(data$year))


########################################
## 1) Load data and cartography files ##
########################################
#load("./dataMmodel.RData")

#head(carto_UP)
#str(data_UP)
data <- read_excel("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/data_dengue/2016-4f.xlsx")
#source("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/R/Data_e.R")
#J <- length(data_UP)
#S <- length(unique(data_UP[[1]]$NAME_1))
#T <- length(unique(data_UP[[1]]$year))

#data <- do.call(rbind,data_UP)
#data$Severity <- rep(1:J,each=S*T)


######################################
## 2) Load previously fitted models ##
######################################

## INLA ##
load("./results/MODELS_inla_icar_prediction_2016.RData")
#load("results/MODELS_inla_lcar.RData")
#load("results/MODELS_inla_pcar.RData")


#t.from <- min(data_UP[[1]]$year)
#t.to <- max(data_UP[[1]]$year)

t.from <- min(data$time)
t.to <- max(data$time)
## WinBUGS ##
#load("resul/results_winbugs_bym_fe.RData")
#load("resul/results_winbugs_bym_re.RData")
#load("resul/results_winbugs_icar_fe.RData")
#load("resul/results_winbugs_icar_re.RData")
#load("resul/results_winbugs_pcar_fe.RData")
#load("resul/results_winbugs_pcar_re.RData")
#load("resul/results_winbugs_lcar_fe.RData")
#load("resul/results_winbugs_lcar_re.RData")


##################################################################################
################                    Figures                       ################
##################################################################################

## Folder to save the figures ##
if(!file.exists("figures")) {dir.create("figures")}




####################################
##  Figure 2. Evolution of the crude rates (per 100000 women) of rapes and
##            dowry deaths in Uttar Pradesh in the period 2001-2014
####################################
crude_rates <- data.frame(DF=1e+5*aggregate(data_UP$DF$OBS, by=list(data_UP$DF$year), sum)$x/aggregate(data_UP$DF$pop, by=list(data_UP$DF$year), sum)$x,
                          DHF=1e+5*aggregate(data_UP$DHF$OBS, by=list(data_UP$DHF$year), sum)$x/aggregate(data_UP$DHF$pop, by=list(data_UP$DHF$year), sum)$x,
                          DSS=1e+5*aggregate(data_UP$DSS$OBS, by=list(data_UP$DSS$year), sum)$x/aggregate(data_UP$DSS$pop, by=list(data_UP$DSS$year), sum)$x)

inf <- round(min(crude_rates))
top <- round(max(crude_rates))

selected_colors <- c("chocolate1", "tomato4","chocolate2")


## File: figure_2.pdf
pdf("figures/figure_2.pdf", height=5, width=7.5, onefile=FALSE)
plot(t.from:t.to, crude_rates$DHF, type="l", xlab ="",ylab ="", 
     ylim=c(inf, top), col=selected_colors[1], cex.axis=1, lwd=4)
lines(t.from:t.to, crude_rates$DF, col=selected_colors[2],lwd=4)
lines(t.from:t.to, crude_rates$DSS, col=selected_colors[3],lwd=4)
legend("topleft",  c("DF","DHF","DSS"), ncol=1, pch=c("-", "-"), 
       col=selected_colors, bty="n",lwd=c(4,4), cex=1.2)
dev.off()


Model <- MODELS.inla.icar$TypeI

carto <- thai_map

## Posterior means of district-specific spatial risks ##
spatial <- matrix(unlist(lapply(Model$marginals.random$idx, function(x) inla.emarginal(exp,x))),S,J,byrow=F)

inf <- min(spatial)
top <- max(spatial)
values <- c(round(seq(inf,1,length.out=5),2), round(seq(1,top,length.out=5),2)[-1])

paleta <- brewer.pal(8,"YlOrRd")
carto$spatial.DF <- spatial[,1]
carto$spatial.DHF <- spatial[,2]
carto$spatial.DSS <- spatial[,3]

Map.spatial <- tm_shape(carto) + 
  tm_polygons(col=c("spatial.DF","spatial.DHF","spatial.DSS"), palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") + 
  tm_layout(panel.labels=c("DF","DHF","DSS"), legend.position=c("right","top","left"))


## Posterior exceedence probabilities ##
probs <- matrix(unlist(lapply(Model$marginals.random$idx, function(x){1-inla.pmarginal(0, x)})),S,J,byrow=F)

paleta <- brewer.pal(5,"PuBu")
values <- c(0,0.1,0.2,0.8,0.9,1)

carto$prob.DF <- probs[,1]
carto$prob.DHF <- probs[,2]
carto$prob.DSS <- probs[,3]

Map.prob <- tm_shape(carto) + 
  tm_polygons(col=c("prob.DF","prob.DHF","prob.DSS"), palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left",
              labels=c("[0-0.1)","[0.1-0.2)","[0.2-0.8)","[0.8-0.9)","[0.9-1]")) + 
  tm_layout(panel.labels=c("DF","DHF","DSS"), legend.position=c("right","top","left"))


## File: figure_4.pdf
tmap_save(tmap_arrange(Map.spatial, Map.prob, ncol=2), file="figures/figure_4.pdf")


####################################
##  Figure 6. Map of estimated incidence risks for rape (top) and posterior probabilities
##            that the relative risk is greater than 1 (P(Ritj>1|O)) (bottom) in Uttar 
##            Pradesh
##  Figure 7. Map of estimated incidence risks for dowry deaths (top) and posterior 
##            probabilities that the relative risk isgreater than one (P(Ritj > 1|O)) 
##            (bottom) in Uttar Pradesh
####################################
Model <-  MODELS.inla.icar$TypeI_RW1

## Maps of posterior mean estimates of relative risks ##
#risks <- matrix(Model$summary.fitted.values$mean,S*T,J,byrow=F)
SMR <- data$SMR

inf <- min(SMR)
top <- max(SMR)
values <- c(round(seq(inf,1,length.out=5),2), round(seq(1,top,length.out=5),2)[-1])

paleta <- brewer.pal(8,"YlOrRd")

carto.DF <- thai_map
carto.DHF <- thai_map
carto.DSS <- thai_map
for(i in seq(T)){
  carto.DF$var <- matrix(data$SMR,S,T,byrow=F)[,i]
  carto.DHF$var <- matrix(data$SMR[,2],S,T,byrow=F)[,i] 
  carto.DSS$var <- matrix(data$SMR[,3],S,T,byrow=F)[,i]
  
  names(carto.DF)[ncol(carto.DF)] <- paste("Year",i,sep=".")
  names(carto.DHF)[ncol(carto.DHF)] <- paste("Year",i,sep=".")
  names(carto.DSS)[ncol(carto.DSS)] <- paste("Year",i,sep=".")
}

Map.risk1 <- tm_shape(carto.DF) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="DF", main.title.position="center", panel.label.size=1,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=1, ncol=12)

Map.risk2 <- tm_shape(carto.DHF) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="DHF", main.title.position="center", panel.label.size=1,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=4, ncol=3)

Map.risk3 <- tm_shape(carto.DSS) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="DSS", main.title.position="center", panel.label.size=1,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=4, ncol=3)


thai_map <- readOGR("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/shp_thai/gadm41_THA_1.shp")
load("results/MODELS_inla_icar.RData")
Model <-  MODELS.inla.icar_B
## Maps of posterior exceedence probabilities ##
probs <- matrix(1-Model$summary.fitted.values$`1 cdf`,S*T,J,byrow=F)

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
  tm_layout(main.title="", main.title.position="left", panel.label.size=1.5,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=1, ncol=12)

Map.prob2 <- tm_shape(carto.DHF) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="", main.title.position="left", panel.label.size=1.5,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=1, ncol=12)

Map.prob3 <- tm_shape(carto.DSS) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="", main.title.position="left", panel.label.size=1.5,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=1, ncol=12)


## File: figure_6.pdf
tmap_save(tmap_arrange(Map.prob1, Map.prob2,Map.prob3), file="figures/prob_DF_all.jpg")
tmap_save( Map.prob1, width = 2400 ,height = 800, file="figures/prob_DF.jpg")
tmap_save( Map.prob2, width = 2400 ,height = 800, file="figures/prob_DHF.jpg")
tmap_save( Map.prob3, width = 2400 ,height = 800, file="figures/prob_DSS.jpg")
#tmap_save(tmap_arrange(Map.risk2, Map.prob2, nrow=2), file="figures/figure_7.pdf")

## File: figure_7.pdf
#tmap_save(tmap_arrange(Map.risk2, Map.prob2, nrow=2), file="figures/figure_7.pdf")
#tmap_save( Map.prob2, width = 2400 ,height = 3000, file="figures/figure_DHF_cluster.pdf")

#tmap_save(tmap_arrange(Map.risk3, Map.prob3, nrow=2), file="figures/figure_7_2.pdf")
#tmap_save( Map.prob3, width = 2400 ,height = 3000, file="figures/figure_DSS_cluster.pdf")
#tmap_save( Map.risk3, width = 2400 ,height = 3000, file="figures/figure_DSS_maprisk.pdf")
####################################
##  Figure 5. Temporal pattern of incidence risks (posterior means of exp(gamma_tj))
##            for rape and dowry deaths in Uttar Pradesh
####################################
Model <- MODELS.inla.icar$TypeI

## Posterior means of year-specific temporal risks ##
temporal <- matrix(unlist(lapply(Model$marginals.random$idy, function(x) inla.emarginal(exp,x))),T,J,byrow=F)

aux <- lapply(Model$marginals.random$idy, function(x) inla.tmarginal(exp,x))
q1 <- matrix(unlist(lapply(aux, function(x) inla.qmarginal(0.025,x))),T,J,byrow=F) 
q2 <- matrix(unlist(lapply(aux, function(x) inla.qmarginal(0.975,x))),T,J,byrow=F) 

inf <- min(q1)-0.05
top <- max(q2)+0.05
selected_colors <- c(rgb(154,192,205,alpha=150, maxColorValue=255),
                     rgb(69,139,116,alpha=150, maxColorValue=255),
                     rgb(200,100,50,alpha=150, maxColorValue=255))


## File: figure_5.pdf
graphics.off()

pdf("figures/figure_5.pdf", height=5, width=7.5, onefile=FALSE)
plot(range(x), c(inf, top), type="n", xlab="Year", ylab="",
     xaxt="n", cex.lab=1, cex.axis=1, cex.main=1, main=NULL)
title(ylab=expression(exp(gamma[t])), line=2.5, cex.lab=1.1)
axis(1, at=seq(1,T), labels=seq(t.from,t.to), las=0, cex.axis=1)

x <- seq(1,T)
for(i in seq(J)){
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(q1[,i], tail(q2[,i], 1), rev(q2[,i]), q1[1,i])
  polygon(X.Vec, Y.Vec, col = selected_colors[i], border = NA)
  lines(temporal[,i])
  abline(h=1,lty=2)
}
legend("topleft", inset=.02, legend=c("DF","DHF","DSS"), fill=selected_colors, horiz=FALSE, cex=1,box.lty=0)
dev.off()

####################################
##  Table 1.  Descriptive statistics. Minimun, first quartile (q1), mean, third
##            quartile (q3), maximun, standard desviation (sd), and coefficient of
##            variation of the number of rapes and dowry deaths in the districts 
##            of Uttar Pradesh per year
####################################

table.1 <- matrix(NA, nrow=T, ncol=15)
colnames(table.1) <- c("Year", "min", "q1", "mean", "q3", "max","sd", "|cv|",
                       "min", "q1", "mean", "q3", "max","sd", "|cv|")


k <- 1
for(i in seq(t.from,t.to)){
  table.1[k,]<- c(i, summary(data$OBS[data$year==i & data$Severity==1])[-3],
                  sd(data$OBS[data$year==i & data$Severity==1]), NA,
                  
                  summary(data$OBS[data$year==i & data$Severity==2])[-3],
                  sd(data$OBS[data$year==i & data$Severity==2]), NA)
                  
  k <- k+1
}
table.1[,8] <- abs(table.1[7]/ table.1[,4])
table.1[,15]<- abs(table.1[,14]/ table.1[,11])
table.1 <- round(table.1, 1)
print(table.1)


table.2 <- matrix(NA, nrow=T, ncol=10)
colnames(table.2) <- c("Year", "Total Case","Mean Per Month", "SD", "Total case" , "Mean Per Month" , "SD","Total case" , "Mean Per Month" , "SD")


k <- 1
for(i in seq(t.from,t.to)){
  table.2[k,]<- c(i, sum(data$OBS[data$year==i & data$Severity==1])[-3],
                  sum(data$OBS[data$year==i & data$Severity==1])/12,
                  sd(data$OBS[data$year==i & data$Severity==1]), 
                  sum(data$OBS[data$year==i & data$Severity==2])[-3],
                  sum(data$OBS[data$year==i & data$Severity==2])/12,
                  sd(data$OBS[data$year==i & data$Severity==2]),
                  sum(data$OBS[data$year==i & data$Severity==3])[-3],
                  sum(data$OBS[data$year==i & data$Severity==3])/12,
                  sd(data$OBS[data$year==i & data$Severity==3]))
                  
  
  k <- k+1
}

table.2 <- round(table.2)
print(table.2)

write.csv(table.2,"C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/data_dengue/EDA_table.csv", row.names = FALSE)



## latex
latex_table<-xtable::xtable(table.2,
                            caption="Descriptive statistics. Minimun, first quartile ($q_1$), mean, third quartile ($q_3$), maximun, standard desviation (sd), and coefficient of variation of the number of rapes and dowry deaths in the districts of Uttar Pradesh per year.",
                            label="t_Severity", digits=c(1),
                            display=c("d", "d", "d","d","d","d","d","d","d", "d","d"))
xtable::print.xtable(latex_table, include.rownames = FALSE, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))
##################################################################################
################                    Tables                        ################
##################################################################################

####################################
##  Table 3.  Model selection criteria, DIC, WAIC and LS for INLA models
####################################
Model <-  Model_association
predicted.p.value <- c()
n <- length(data[,1])
for(i in (1:n)) {
  predicted.p.value[i] <- inla.pmarginal(q=data$OBS[i],
                                         marginal=Model$marginals.fitted.values[[i]])
}

plot(data$OBS,Model$summary.fitted.values$mean * data$EXP,
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")


mse_calc <- function(x,y){
  result = 0
  for (i in 1:length(x)){
    result = result + (x[i] - y[i])^2
  }
  result = result/(length(x))
  return(result)
}

bias_calc <- function(x,y){
  result = 0
  for (i in 1:length(x)){
    result = result + (x[i] - y[i])
  }
  result = result/(length(x))
  return(result)
}

train_mse_1 <- mse_calc(Model$summary.fitted.values$mean,data$OBS)
#bias_mse_1 <- bias_calc(Model$summary.fitted.values$mean,data$OBS)

plot((Model$summary.fitted.values$mean * data$EXP),data$OBS)




table3.iCAR <- do.call(rbind,lapply(Model, function(x) data.frame(DIC=Model$dic$dic, WAIC=Model$waic$waic, LS=-sum(log(Model$cpo$cpo)) ,MSE = mse_calc(Model$summary.fitted.values$mean[5082:length(data)]* data$EXP[5082:length(data)],data$OBS[5082:length(data)]), Error = bias_calc(Model$summary.fitted.values$mean[5082:length(data)]* data$EXP[5082:length(data)],data$OBS[5082:length(data)]), Correlation = cor(Model$summary.fitted.values$mean[5082:length(data)]* data$EXP[5082:length(data)],data$OBS[5082:length(data)]))))
rownames(table3.iCAR) <- paste("iCAR",rownames(table3.iCAR),sep=".")

table3.LCAR <- do.call(rbind,lapply(Model, function(x) data.frame(DIC=x$dic$dic, WAIC=x$waic$waic, LS=-sum(log(x$cpo$cpo)) ,MSE = mse_calc(x$summary.fitted.values$mean* data$EXP,data$OBS),Error = bias_calc(x$summary.fitted.values$mean* data$EXP,data$OBS), Correlation = cor(x$summary.fitted.values$mean* data$EXP,data$OBS))))
rownames(table3.LCAR) <- paste("LCAR",rownames(table3.LCAR),sep=".")

table3.pCAR <- do.call(rbind,lapply(Model, function(x) data.frame(DIC=x$dic$dic, WAIC=x$waic$waic, LS=-sum(log(x$cpo$cpo)) ,MSE = mse_calc(x$summary.fitted.values$mean* data$EXP,data$OBS),Error = bias_calc(x$summary.fitted.values$mean* data$EXP,data$OBS),Correlation = cor(x$summary.fitted.values$mean* data$EXP,data$OBS))))
rownames(table3.pCAR) <- paste("pCAR",rownames(table3.pCAR),sep=".")

table.3 <- rbind(table3.iCAR,table3.LCAR,table3.pCAR)
print(table.3)

write.csv(table3.iCAR,"C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/data_dengue/Modelcompar_table_Pcar_pred.csv", row.names = FALSE)

## latex
x_table<-xtable::xtable(table3.iCAR, 
                        caption="Model selection criteria, DIC, WAIC and LS, for different models.", 
                        label="t_dic", digits=c(3))

xtable::print.xtable(x_table, include.rownames = T, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))


plot(MODELS.inla$summary.fitted.values$mean,data$OBS)

cor_model

