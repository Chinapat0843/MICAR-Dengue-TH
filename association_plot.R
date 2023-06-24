
#rm(list=ls())

library(INLA)
library(dplyr)
library(sp)
library(rgdal)
library(RColorBrewer)
library(geodata)
library(spdep)
library(leaflet)
library(tidyverse)
library(maptools)
library(ggplot2)
library(readxl)
# data1 <- read.csv("D:/Desktop/Year 4/Sr project MU PSCM/Model/data_2011.csv")
data <- read_excel("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/data_dengue/2011-2012-4f.xlsx")
# source("D:/Desktop/Year 4/Sr project MU PSCM/Model/main/main.R") 

๒load("./results/MODELS_inla_icar_2011_2012.RData")
Model_association <- MODELS.inla.icar_B$association_model

map <- readOGR("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/shp_thai/gadm41_THA_1.shp")
proj4string(map) <- CRS("+proj=longlat +ellps=clrk66")
map <- spTransform(map,CRS("+proj=longlat +datum=WGS84 +no_defs"))



J <- length(unique(data$type))
S <- length(unique(data$area))

summary_random <- Model_association[["summary.random"]]
summary_random <- summary_random[!(names(summary_random) %in% c("idx", "idy"))]

# initialize empty vector to store type_list values
type_list <- vector()

for (j in 1:J) {
  for (s in 1:S) {
    type_list <- c(type_list, rep(j))
  }
}


# loop over sublists in summary_random
for (i in 1:length(summary_random)) {
  # add type_list as a new column to the i-th sublist
  summary_random[[i]]$type <- type_list
}



map@data$ID_area <- seq(1,77)
summary_random <- summary_random[[1]][1:77,]


# datafiltered <- summary_random
# ordercounties <- match(map@data$ID_area,datafiltered$ID)
# map@data <- datafiltered[ordercounties, ]

map@data <- merge(x = map@data , y = summary_random , by.x = "ID_area" , by.y = "ID")

map@data <- map@data %>% mutate(Significant =
                                            case_when(`0.025quant` / `0.975quant` > 0 ~ "yes", 
                                                      `0.025quant` / `0.975quant` < 0 ~ "no")
)


# shp_df <- broom::tidy(map, region = "Significant")
# 
# map@data <- map@data %>% mutate(id = row.names(.))
# shp_df <- broom::tidy(map, region = "id")
# shp_df <- shp_df %>% left_join(map@data, by = c("id"="id"))



col <- c("NAME_1","mean","sd","0.025quant","0.5quant","0.975quant","Significant","type")
shp_info <- map@data[,names(map@data) %in% col]
colnames(shp_info)[1] <- "area"


map <- readOGR("C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/shp_thai/gadm41_THA_1.shp")
proj4string(map) <- CRS("+proj=longlat +ellps=clrk66")
map <- spTransform(map,CRS("+proj=longlat +datum=WGS84 +no_defs"))

# datafiltered <- data[which(data$level == 1), ]
datafiltered <- shp_info
ordercounties <- match(map@data$NAME_1, datafiltered$area)
map@data <- datafiltered[ordercounties, ]

# Create leaflet
l <- leaflet(map) %>% addTiles()

# Convert Significant column to factor
map$Significant <- factor(map$Significant, levels = c("no", "yes"))

pal <- colorFactor(palette = "YlOrRd", domain = map$Significant)
labels <- sprintf("<strong> %s </strong> <br/> Significant: %s ",
                  map$area, map$Significant
) %>%
  lapply(htmltools::HTML)

l %>%
  addPolygons(
    color = "grey", weight = 1,
    fillColor = ~ pal(Significant), fillOpacity = 0.7,
    highlightOptions = highlightOptions(weight = 4),
    label = labels,
    labelOptions = labelOptions(
      style = list(
        "font-weight" = "normal",
        padding = "3px 8px"
      ),
      textsize = "15px", direction = "auto"
    )
  ) %>%
  addLegend(
    pal = pal, values = ~Significant, opacity = 0.7,
    title = "Significant", position = "bottomright"
  )


if (is.null(names(data)))
  xd <- character(0)

xd<- names(data) 
xd2 <- xd[!(names(data) %in% c("area", "time", "OBS", "EXP", "pop", "type"))]

if (is.null(xd))
  xd <- character(0)

xd2<- c("-", xd)

write.csv(shp_info,"C:/Users/skybl/OneDrive/Documents/Dengue_spatiotemporal/data_dengue/association_table.csv", row.names = FALSE)
