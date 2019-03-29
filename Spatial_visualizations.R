# loading the required packages and data
library(ggplot2)
library(ggmap)
library(maps)
lake_level_eff <- read_csv("Documents/GitHub/Hg_Fish/lake_level_eff.csv")
lake_level_eff

# getting the map
mapmn <- get_map(location = c(left = -98, bottom = 43, right = -89, top = 50))
# plotting the map with some points on it
ggmap(mapmn) +
  geom_point(data = lake_level_eff, aes(x = Longitude, y =Latitude, color = Random_Lake_Effect, alpha = .5), size = 1) +
   guides(fill=FALSE, alpha=FALSE, size=FALSE) +scale_colour_gradient2(low = "white", high = "black", midpoint = -1.7)


