
#rm(list=ls())

library(INLA)
library(dplyr)
library(sp)
library(rgdal)
library(RColorBrewer)
library(geodata)
library(spdep)
library(leaflet)
library(readxl)
## install 'webshot' package
library(devtools)
#install_github("wch/webshot")

## load packages
#library(leaflet)
library(htmlwidgets)
library(webshot)
library(mapview)
# data1 <- read.csv("D:/Desktop/Year 4/Sr project MU PSCM/Model/data_2011.csv")
data <- read_excel("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/data_dengue/2011-4f.xlsx")


# source("D:/Desktop/Year 4/Sr project MU PSCM/Model/main/main.R") 

#load("./results/MODELS_inla_icar_5y5f.RData")
# load("C:/Users/tussa/OneDrive/Documents/MSApp/results/MODELS_icar_t1.RData")

Model_cluster <- MODELS.inla.icar_B$cluster_model

J <- length(unique(data$type))
S <- length(unique(data$area))
T <- length(unique(data$time))

t.from <- min(data$time)
t.to <- max(data$time)

# thai_map <- readOGR("D:/Desktop/Year 4/Sr project MU PSCM/Model/gadm40_THA_1.shp")
map <- readOGR("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/shp_thai/gadm41_THA_1.shp")
proj4string(map) <- CRS("+proj=longlat +ellps=clrk66")
map <- spTransform(map,CRS("+proj=longlat +datum=WGS84 +no_defs"))
#colnames(thai_map@data)[4] <- "area"


## Maps of posterior exceedence probabilities ##
prob <- matrix(1-Model_cluster$summary.fitted.values$`1 cdf`,S*T,J,byrow=F)

prob = c(prob[,1],prob[,2],prob[,3])

data <- data[order(data$type),]
data_p <- cbind(data, prob)


num_month = 8
num_year = 1
start = 1 + ((231 * num_month*num_year)-231) 

stop = 77 +  ((231 * num_month*num_year)-231)
# Add data to map
#datafiltered <- data_p[77 + ((231 * num_month*num_year)-231) : 154 +  ((231 * num_month*num_year)-231),]
data_p2 <- data_p[order(data_p$time),]
row.names(data_p2) <- NULL
datafiltered <- data_p2[start:stop,]

ordercounties <- match(map@data$NAME_1,datafiltered$area)
map@data <- datafiltered[ordercounties, ]

# map@data <- merge(x = thai_map@data , y = datafiltered , by.x = "NAME_1" , by.y = "area")
# map@data <- merge(y = thai_map@data , x = datafiltered , by.x = "area" , by.y = "area")
# map@data <- map@data[order(map@data$time),]

for(i in seq(T)){
  map$exceedance_prob <- matrix(map@data$prob,S,T,byrow=F)[,i]
}

print(unique(map$exceedance_prob))




# Create a new variable 'var1_group' that groups the values of 'var1' into two ranges
map$exceedance_prob_group <- cut(map$exceedance_prob, c(0, 0.95, 1), labels = c("0 - 0.95", "0.95 - 1"))
l <- leaflet(map) %>% addTiles()
# Create a new color palette for the two groups
pal <- colorFactor(palette = c("#BCBCBC", "#FF2D00"), domain = map$exceedance_prob_group)
labels <- sprintf("<strong> %s </strong> <br/> Exceedance Probability: %s ", map$area, map$exceedance_prob) %>%
  lapply(htmltools::HTML)
# Use 'var1_group' instead of 'var1' in the fillColor argument
l %>%
  addPolygons(
    color = "white", weight = 1,
    fillColor = ~ pal(exceedance_prob_group), fillOpacity = 0.8,
    highlightOptions = highlightOptions(weight = 4),
    label = labels,
    labelOptions = labelOptions(
      style = list(
        "font-weight" = "normal",
        padding="3px 8px"
      ),
      textsize = "15px", direction = "auto"
    )
  ) %>%
  # Use 'var1_group' instead of 'var1' in the values argument
  addLegend(
    pal = pal, values = ~exceedance_prob_group, opacity = 0.7,
    title = "Exceedance Probability", position = "bottomright"
  )

install.packages("mapview", repo="http://cran.r-project.org", dep=T)

#saveWidget(l, "temp.html", selfcontained = FALSE)
#webshot("temp.html", file=paste0("./figures/cluster_report/df.1.1.png")
      # )#
#m = mapview(l)
mapshot(l, file = paste0(getwd(), "./figures/cluster_report/df.1.1.png"),
        remove_controls = c("homeButton", "layersControl"))
