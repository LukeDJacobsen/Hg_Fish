################################################################################
#Data import
################################################################################

library(tidyverse)
library(readr)
fish_clean <- read_csv("Documents/Hg_Fish/fish_clean.csv")

################################################################################
#Data formatting
################################################################################

#Getting rid of composite samples because they will bias
#variance estimates

no_comp <- fish_clean %>% filter(composite == "N")

#make a species variable that's fewer levels
library(forcats)
table(no_comp$SPEC)
no_comp$SPEC_clumped <- fct_lump(no_comp$SPEC, n = 50)

table(no_comp$SPEC_clumped)
library(skimr)
skim(no_comp)

################################################################################
#Model
################################################################################
library(lme4)
names(no_comp)
no_comp$Year_F <- as.factor(no_comp$Year)
lme1 <- lmer(LN_HGPPM ~  + SPEC_clumped + SPEC_clumped:LGTHIN + 
       (1|Year_F) + (1|DOWID), data = no_comp)
summary(lme1)
lake_ranef <- ranef(lme1)

################################################################################
#Organize data
################################################################################
#goal: df with column for renef estimate, dowid, lat, and long

#making df with mean lat and long by lake.
fish_clean$lat <- ifelse(fish_clean$LAT == 0, NA, fish_clean$LAT)
fish_clean$long <- ifelse(fish_clean$LONG == 0, NA, fish_clean$LONG)


lake_extra <- fish_clean %>% filter(composite == "N") %>% select(DOWID, lat, long) %>%
  group_by(DOWID) %>% summarise(mean_LAT = mean(lat, na.rm =T), mean_LONG = mean(long, na.rm = T))

#get rid of duplicated rows
lake <- lake_extra[!duplicated(lake_extra),]

#awkward data munging stuff
a <- as.data.frame(lake_ranef[[1]])
names(a) <- c("Random_Lake_Effect")
DOWID <- rownames(a)
a$DOWID <- DOWID
a <- as_tibble(a)

lake_level_eff <- inner_join(a, lake, by = "DOWID")

#appears to be a data issue with lake 03003000. Changed to match what it should be
lake_level_eff$mean_LONG <- ifelse(lake_level_eff$DOWID == "03003000", 952649, lake_level_eff$mean_LONG)

################################################################################
#Spatial covariance modeling and visualization
################################################################################
ggplot(lake_level_eff, aes(mean_LONG, mean_LAT)) + geom_point(aes(color = Random_Lake_Effect)) 
lake_level_eff %>% filter(mean_LONG < 870000)


#now what does plot look like
ggplot(lake_level_eff, aes(-mean_LONG, mean_LAT)) + geom_point(aes(color = Random_Lake_Effect)) +
  scale_colour_gradient2(low = "dark green", high = "dark red")


#models
names(lake_level_eff)

mod_loc2 <- lm(Random_Lake_Effect ~ mean_LAT + mean_LONG + 
                 I(mean_LAT^2) + I(mean_LONG^2), data = lake_level_eff)
summary(mod_loc2)


mod_loc3 <- lm(Random_Lake_Effect ~ mean_LAT + mean_LONG + 
                 I(mean_LAT^2) + I(mean_LONG^2) +   I(mean_LAT^3) + I(mean_LONG^3), data = lake_level_eff)
summary(mod_loc3)

mod_loc4 <- lm(Random_Lake_Effect ~ mean_LAT + mean_LONG + 
                 I(mean_LAT^2) + I(mean_LONG^2) +   I(mean_LAT^3) + I(mean_LONG^3) +
                 I(mean_LAT^4) + I(mean_LONG^4), data = lake_level_eff)
summary(mod_loc4)


############################

#distance between lake i and j matrix
distance <- matrix(nrow = nrow(lake_level_eff), ncol = nrow(lake_level_eff))
for (i in 1:nrow(lake_level_eff)){
  for (j in 1:nrow(lake_level_eff)){
    distance[i,j] <- sqrt((lake_level_eff$mean_LAT[i] - lake_level_eff$mean_LAT[j])^2 +
                            (lake_level_eff$mean_LONG[i] - lake_level_eff$mean_LONG[j])^2)
  }
}

#change in random effect estimate between lake i and j
change_eff <- matrix(nrow = nrow(lake_level_eff), ncol = nrow(lake_level_eff))
for (i in 1:nrow(lake_level_eff)){
  for (j in 1:nrow(lake_level_eff)){
    change_eff[i,j] <- abs(lake_level_eff$Random_Lake_Effect[i] - lake_level_eff$Random_Lake_Effect[j])
  }
}
plot(distance[1201,], change_eff[1201,])
distance <- ifelse(distance == 0, NA, distance)
dist_change <- tibble(as.vector(distance), as.vector(change_eff))
names(dist_change) <- c("distance", "change_effect")
dist_effect <- ggplot(dist_change, aes(distance, change_effect)) + geom_smooth()
dist_effect + coord_cartesian(xlim = c(0,50000))
